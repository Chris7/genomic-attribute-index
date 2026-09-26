use std::{
    collections::HashMap,
    fmt,
    sync::{Arc, Mutex},
    time::{Duration, Instant},
};

#[cfg(unix)]
use pprof::{ProfilerGuardBuilder, Report};
use tracing::{Id, Subscriber, span::Attributes};
use tracing_subscriber::{
    layer::{Context, Layer, SubscriberExt},
    registry::{LookupSpan, Registry},
};

#[derive(Clone, Default)]
pub(crate) struct Profiler {
    stats: Arc<Mutex<HashMap<String, Stat>>>,
}

#[derive(Default)]
struct Stat {
    calls: u64,
    total: Duration,
}

struct SpanState {
    stack: String,
    starts: Vec<Instant>,
}

impl Profiler {
    pub(crate) fn run<T>(self, f: impl FnOnce() -> T) -> T {
        let subscriber = Registry::default().with(ProfilerLayer {
            profiler: self.clone(),
        });
        tracing::subscriber::with_default(subscriber, || {
            let result = f();
            self.print_report();
            result
        })
    }

    fn record(&self, stack: String, duration: Duration) {
        let mut stats = self
            .stats
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let entry = stats.entry(stack).or_default();
        entry.calls = entry.calls.saturating_add(1);
        entry.total += duration;
    }

    fn print_report(&self) {
        let stats = self
            .stats
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let mut rows = stats
            .iter()
            .map(|(stack, stat)| (stack_frames(stack), stat.calls, stat.total))
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| right.2.cmp(&left.2).then_with(|| left.0.cmp(&right.0)));
        let stack_width = rows
            .iter()
            .map(|(frames, _, _)| stack_label_width(frames))
            .max()
            .unwrap_or(0)
            .max("Stack".len());

        eprintln!("Profile results");
        eprintln!(
            "{:<stack_width$} {:>10} {:>14} {:>14}",
            "Stack", "Calls", "Total (ms)", "Avg (us)"
        );
        eprintln!(
            "{:<stack_width$} {:>10} {:>14} {:>14}",
            "-----", "-----", "----------", "--------"
        );
        for (frames, calls, total) in rows {
            let average_micros = if calls == 0 {
                0.0
            } else {
                total.as_secs_f64() * 1_000_000.0 / calls as f64
            };
            print_stacked_row(
                &frames,
                stack_width,
                format_args!(
                    "{:>10} {:>14.3} {:>14.3}",
                    calls,
                    total.as_secs_f64() * 1000.0,
                    average_micros
                ),
            );
        }
    }
}

struct ProfilerLayer {
    profiler: Profiler,
}

impl<S> Layer<S> for ProfilerLayer
where
    S: Subscriber + for<'span> LookupSpan<'span>,
{
    fn on_new_span(&self, attributes: &Attributes<'_>, id: &Id, context: Context<'_, S>) {
        if let Some(span) = context.span(id) {
            let metadata = attributes.metadata();
            let frame = format!("{}::{}", metadata.target(), metadata.name());
            let stack = match span.parent() {
                Some(parent) => parent
                    .extensions()
                    .get::<SpanState>()
                    .map(|state| format!("{} -> {}", state.stack, frame))
                    .unwrap_or(frame),
                None => frame,
            };
            span.extensions_mut().insert(SpanState {
                stack,
                starts: Vec::new(),
            });
        }
    }

    fn on_enter(&self, id: &Id, context: Context<'_, S>) {
        if let Some(span) = context.span(id)
            && let Some(state) = span.extensions_mut().get_mut::<SpanState>()
        {
            state.starts.push(Instant::now());
        }
    }

    fn on_exit(&self, id: &Id, context: Context<'_, S>) {
        if let Some(span) = context.span(id) {
            let snapshot = span
                .extensions_mut()
                .get_mut::<SpanState>()
                .and_then(|state| {
                    state
                        .starts
                        .pop()
                        .map(|started| (state.stack.clone(), started.elapsed()))
                });
            if let Some((stack, duration)) = snapshot {
                self.profiler.record(stack, duration);
            }
        }
    }
}

#[cfg(unix)]
#[derive(Clone, Default)]
pub(crate) struct SamplingProfiler;

#[cfg(unix)]
impl SamplingProfiler {
    pub(crate) fn run<T>(self, f: impl FnOnce() -> T) -> Result<T, String> {
        let guard = ProfilerGuardBuilder::default()
            .frequency(100)
            .blocklist(&["libc", "libgcc", "pthread", "vdso"])
            .build()
            .map_err(|error| format!("could not start sample profiler: {error}"))?;

        let result = f();
        let report = guard
            .report()
            .build()
            .map_err(|error| format!("could not collect sample profile: {error}"))?;
        self.print_report(&report);

        Ok(result)
    }

    fn print_report(&self, report: &Report) {
        let mut stats: HashMap<String, u64> = HashMap::new();
        let mut total_samples = 0_u64;
        let sample_frequency = report.timing.frequency.max(1) as f64;

        for (frames, count) in &report.data {
            if *count <= 0 {
                continue;
            }

            let sample_count = *count as u64;
            total_samples = total_samples.saturating_add(sample_count);
            if let Some(stack) = sample_stack(frames) {
                let entry = stats.entry(stack).or_default();
                *entry = entry.saturating_add(sample_count);
            }
        }

        let mut rows = stats
            .into_iter()
            .map(|(stack, samples)| (stack_frames(&stack), samples))
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        let stack_width = rows
            .iter()
            .map(|(frames, _)| stack_label_width(frames))
            .max()
            .unwrap_or(0)
            .max("Stack".len());

        eprintln!("Sampling profile results");
        eprintln!(
            "{:<stack_width$} {:>10} {:>14} {:>10}",
            "Stack", "Samples", "Time (ms)", "Pct"
        );
        eprintln!(
            "{:<stack_width$} {:>10} {:>14} {:>10}",
            "-----", "-------", "---------", "---"
        );
        for (frames, sample_count) in rows {
            let estimated_ms = sample_count as f64 / sample_frequency * 1000.0;
            let percentage = if total_samples == 0 {
                0.0
            } else {
                sample_count as f64 * 100.0 / total_samples as f64
            };
            print_stacked_row(
                &frames,
                stack_width,
                format_args!(
                    "{:>10} {:>14.3} {:>9.2}%",
                    sample_count, estimated_ms, percentage
                ),
            );
        }
    }
}

#[cfg(unix)]
fn sample_stack(frames: &pprof::Frames) -> Option<String> {
    let mut stack = Vec::new();
    for frame in frames.frames.iter().rev() {
        if let Some(name) = frame.iter().find_map(owned_frame_name) {
            stack.push(name);
        }
    }

    if stack.is_empty() {
        return None;
    }

    trim_sample_prefix(&mut stack);
    (!stack.is_empty()).then(|| stack.join(" -> "))
}

#[cfg(unix)]
fn owned_frame_name(frame: &pprof::Symbol) -> Option<String> {
    let name = frame.name();
    if name.starts_with("gai::") || name.starts_with("r#gai::") {
        Some(name.strip_prefix("r#").unwrap_or(&name).to_string())
    } else {
        None
    }
}

#[cfg(unix)]
fn trim_sample_prefix(frames: &mut Vec<String>) {
    if let Some(index) = frames.iter().position(|frame| frame == "gai::run") {
        frames.drain(..index);
    }
}

fn stack_frames(stack: &str) -> Vec<String> {
    stack.split(" -> ").map(ToString::to_string).collect()
}

fn stack_label_width(frames: &[String]) -> usize {
    match frames.split_last() {
        Some((leaf, parents)) => parents.len() + leaf.len(),
        None => "<empty>".len(),
    }
}

fn print_stacked_row(frames: &[String], stack_width: usize, tail: fmt::Arguments<'_>) {
    match frames.split_last() {
        Some((leaf, parents)) => {
            for (depth, frame) in parents.iter().enumerate() {
                eprintln!("{}{}", " ".repeat(depth), frame);
            }
            let label = format!("{}{}", " ".repeat(parents.len()), leaf);
            eprintln!("{label:<stack_width$} {tail}");
        }
        None => eprintln!("{:<stack_width$} {tail}", "<empty>"),
    }
}
