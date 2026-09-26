use crate::*;

#[path = "index/bed.rs"]
pub(crate) mod bed;
#[path = "index/gff.rs"]
mod gff;

#[derive(Clone, Debug)]
pub(crate) struct ExtractedRecord {
    pub(crate) span: SpanKey,
    pub(crate) terms: Vec<String>,
}

fn read_index_records<'a, R>(
    reader: &'a mut R,
    format: SortFormat,
    reference_ids: &'a HashMap<String, u32>,
    configured: &'a HashSet<String>,
    case_sensitive: bool,
) -> impl Iterator<Item = Result<ExtractedRecord>> + 'a
where
    R: BufRead,
{
    let mut line = Vec::new();
    let mut line_number = 0_usize;
    let mut stopped = false;

    std::iter::from_fn(move || {
        if stopped {
            return None;
        }

        loop {
            line.clear();
            match reader.read_until(b'\n', &mut line) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(error) => return Some(Err(Error::Io(error))),
            }
            line_number += 1;

            let mut raw = line.as_slice();
            if let Some(stripped) = raw.strip_suffix(b"\n") {
                raw = stripped;
            }
            if let Some(stripped) = raw.strip_suffix(b"\r") {
                raw = stripped;
            }

            match format {
                SortFormat::Gff => {
                    if raw == b"##FASTA" {
                        stopped = true;
                        return None;
                    }
                    if raw.first() == Some(&b'#') {
                        continue;
                    }
                    return Some(gff::extract_gff_line(
                        raw,
                        line_number,
                        reference_ids,
                        configured,
                        case_sensitive,
                    ));
                }
                SortFormat::Bed => {
                    if raw.first() == Some(&b'#')
                        || raw.starts_with(b"track ")
                        || raw.starts_with(b"browser ")
                    {
                        continue;
                    }
                    return Some(bed::extract_bed_record(
                        raw,
                        line_number,
                        reference_ids,
                        case_sensitive,
                    ));
                }
            }
        }
    })
}

fn scan_index_reader<R, F>(
    mut reader: R,
    format: SortFormat,
    reference_ids: &HashMap<String, u32>,
    configured: &HashSet<String>,
    case_sensitive: bool,
    process: &mut F,
) -> Result<([u8; 32], u64)>
where
    R: HashingInput,
    F: FnMut(ExtractedRecord) -> Result<()>,
{
    {
        let records = read_index_records(
            &mut reader,
            format,
            reference_ids,
            configured,
            case_sensitive,
        );
        for record in records {
            process(record?)?;
        }
    }
    reader.drain_and_finish().map_err(Error::from)
}

struct ScanIndexContext<'a> {
    format: SortFormat,
    bgzf_threads: usize,
    bytes_read: Arc<AtomicU64>,
    reference_ids: &'a HashMap<String, u32>,
    configured: &'a HashSet<String>,
    case_sensitive: bool,
}

fn scan_index_path<F>(
    path: &Path,
    context: ScanIndexContext<'_>,
    mut process: F,
) -> Result<([u8; 32], u64)>
where
    F: FnMut(ExtractedRecord) -> Result<()>,
{
    let mut source = File::open(path)?;
    let mut magic = [0_u8; 2];
    let magic_length = source.read(&mut magic)?;
    source.rewind()?;
    let hashing = HashingReader::with_counter(source, context.bytes_read);
    if magic_length == magic.len() && magic == [0x1f, 0x8b] {
        if context.bgzf_threads > 1 {
            let workers = NonZeroUsize::new(context.bgzf_threads)
                .ok_or_else(|| Error::InvalidInput("BGZF worker count must be positive".into()))?;
            scan_index_reader(
                bgzf::io::MultithreadedReader::with_worker_count(workers, hashing),
                context.format,
                context.reference_ids,
                context.configured,
                context.case_sensitive,
                &mut process,
            )
        } else {
            scan_index_reader(
                bgzf::io::Reader::new(hashing),
                context.format,
                context.reference_ids,
                context.configured,
                context.case_sensitive,
                &mut process,
            )
        }
    } else {
        scan_index_reader(
            BufReader::new(hashing),
            context.format,
            context.reference_ids,
            context.configured,
            context.case_sensitive,
            &mut process,
        )
    }
}

/// Builds a deterministic GAI beside a BGZF or plain GFF3 source.
pub fn build_name_index(
    gff_path: impl AsRef<Path>,
    coordinate_index_path: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    options: &NameIndexOptions,
) -> Result<IndexStats> {
    build_name_index_with_options(
        gff_path,
        coordinate_index_path,
        destination,
        options,
        &BuildOptions::default(),
    )
}

/// Builds a deterministic GAI with explicit resource and progress controls.
pub fn build_name_index_with_options(
    gff_path: impl AsRef<Path>,
    coordinate_index_path: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    options: &NameIndexOptions,
    build_options: &BuildOptions,
) -> Result<IndexStats> {
    build_name_index_with_options_and_span_block_size(
        gff_path,
        coordinate_index_path,
        destination,
        options,
        DEFAULT_SPANS_PER_BLOCK,
        build_options,
    )
}

/// Builds a GAI with an explicit reference-specific span-block row target.
/// The default [`build_name_index`] value is 4,096; this variant is provided
/// for reproducible compression experiments and deployment tuning.
pub fn build_name_index_with_span_block_size(
    gff_path: impl AsRef<Path>,
    coordinate_index_path: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    options: &NameIndexOptions,
    spans_per_block: usize,
) -> Result<IndexStats> {
    build_name_index_with_options_and_span_block_size(
        gff_path,
        coordinate_index_path,
        destination,
        options,
        spans_per_block,
        &BuildOptions::default(),
    )
}

/// Builds a GAI with explicit span-block and resource controls.
pub fn build_name_index_with_options_and_span_block_size(
    gff_path: impl AsRef<Path>,
    coordinate_index_path: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    options: &NameIndexOptions,
    spans_per_block: usize,
    build_options: &BuildOptions,
) -> Result<IndexStats> {
    if !(1..=1_000_000).contains(&spans_per_block) {
        return Err(Error::InvalidInput(
            "span block size must be between 1 and 1,000,000".into(),
        ));
    }
    if build_options.memory_budget_bytes == 0 {
        return Err(Error::InvalidInput("memory budget must be positive".into()));
    }
    if build_options.compression_threads == 0 {
        return Err(Error::InvalidInput(
            "compression worker count must be positive".into(),
        ));
    }
    if build_options.bgzf_threads == 0 {
        return Err(Error::InvalidInput(
            "BGZF worker count must be positive".into(),
        ));
    }
    let started = Instant::now();
    let source_path = gff_path.as_ref();
    let source_format = SortFormat::from_path(source_path)?;
    let coordinate_index_path = coordinate_index_path.as_ref();
    let destination = destination.as_ref();
    let name_options = match source_format {
        SortFormat::Gff => {
            NameIndexOptions::new(options.attributes.clone(), options.case_sensitive)?
        }
        SortFormat::Bed => NameIndexOptions::bed(options.case_sensitive),
    };
    let destination_parent = destination.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(destination_parent)?;
    let (coordinate_index, coordinate_index_fingerprint) =
        read_coordinate_index_with_fingerprint(coordinate_index_path)?;
    let dictionary = coordinate_index.dictionary(source_format)?;
    let reference_ids = dictionary
        .names
        .iter()
        .enumerate()
        .map(|(index, name)| {
            Ok((
                name.clone(),
                u32::try_from(index)
                    .map_err(|_| Error::InvalidInput("too many reference sequences".into()))?,
            ))
        })
        .collect::<Result<HashMap<_, _>>>()?;

    let configured: HashSet<String> = name_options.attributes.iter().cloned().collect();
    let mut collector = SpillCollector::new(
        build_options.memory_budget_bytes,
        destination_parent,
        &destination
            .file_name()
            .map(|name| name.to_string_lossy())
            .unwrap_or_else(|| "gai".into()),
    );
    let bytes_read = Arc::new(AtomicU64::new(0));
    let mut records_processed = 0_u64;
    let mut records_indexed = 0_u64;
    let mut span_observations = 0_u64;
    let mut pair_observations = 0_u64;
    let mut spill_duration = Duration::ZERO;
    let mut scan_spill_duration = Duration::ZERO;
    let scan_started = Instant::now();
    let mut process_record = |extracted: ExtractedRecord| -> Result<()> {
        records_processed = records_processed
            .checked_add(1)
            .ok_or(Error::InvalidCoordinate)?;
        if extracted.terms.is_empty() {
            if records_processed.is_multiple_of(PROGRESS_RECORD_INTERVAL) {
                report_progress(
                    build_options,
                    BuildPhase::Scan,
                    records_processed,
                    &bytes_read,
                    started,
                );
            }
            return Ok(());
        }
        records_indexed = records_indexed
            .checked_add(1)
            .ok_or(Error::InvalidCoordinate)?;
        span_observations = span_observations
            .checked_add(1)
            .ok_or(Error::InvalidCoordinate)?;
        pair_observations = pair_observations
            .checked_add(extracted.terms.len() as u64)
            .ok_or(Error::InvalidCoordinate)?;
        let should_spill = collector.add(extracted.span, &extracted.terms)?;
        if should_spill {
            let spill_started = Instant::now();
            collector.spill()?;
            let elapsed = spill_started.elapsed();
            spill_duration = spill_duration.saturating_add(elapsed);
            scan_spill_duration = scan_spill_duration.saturating_add(elapsed);
            report_progress(
                build_options,
                BuildPhase::Spill,
                records_processed,
                &bytes_read,
                started,
            );
        } else if records_processed.is_multiple_of(PROGRESS_RECORD_INTERVAL) {
            report_progress(
                build_options,
                BuildPhase::Scan,
                records_processed,
                &bytes_read,
                started,
            );
        }
        Ok(())
    };
    let (source_fingerprint, source_bytes) = scan_index_path(
        source_path,
        ScanIndexContext {
            format: source_format,
            bgzf_threads: build_options.bgzf_threads,
            bytes_read: Arc::clone(&bytes_read),
            reference_ids: &reference_ids,
            configured: &configured,
            case_sensitive: name_options.case_sensitive,
        },
        &mut process_record,
    )?;
    let scan_duration = scan_started.elapsed();
    let final_spill_started = Instant::now();
    collector.finish()?;
    spill_duration = spill_duration.saturating_add(final_spill_started.elapsed());
    report_progress(
        build_options,
        BuildPhase::Scan,
        records_processed,
        &bytes_read,
        started,
    );

    let merge_started = Instant::now();
    let coordinate_spans = merge_unique_spans(collector.run_paths())?;
    let merge_duration = merge_started.elapsed();
    report_progress(
        build_options,
        BuildPhase::Merge,
        records_processed,
        &bytes_read,
        started,
    );

    let postings_started = Instant::now();
    let posting_encoding = encode_postings_from_runs(
        collector.run_paths(),
        &coordinate_spans,
        build_options.compression_threads,
    )?;
    let (terms, posting_directory, postings_data, postings_before, unique_postings, term_count) =
        posting_encoding;
    let postings_duration = postings_started.elapsed();
    report_progress(
        build_options,
        BuildPhase::EncodePostings,
        records_processed,
        &bytes_read,
        started,
    );

    let spans_started = Instant::now();
    let (span_directory, starts_data, lengths_data, span_encoding_stats) =
        encode_span_blocks_with_threads(
            &coordinate_spans,
            spans_per_block,
            build_options.compression_threads,
        )?;
    let spans_duration = spans_started.elapsed();
    let span_directory_bytes = encode_span_directory(&span_directory)?;
    report_progress(
        build_options,
        BuildPhase::EncodeSpans,
        records_processed,
        &bytes_read,
        started,
    );

    let serialize_started = Instant::now();
    let serialized = serialize_index(
        &name_options.attributes,
        name_options.case_sensitive,
        source_fingerprint,
        coordinate_index_fingerprint,
        dictionary.fingerprint,
        term_count,
        coordinate_spans.len() as u64,
        unique_postings,
        posting_directory.len() as u64,
        span_directory.len() as u64,
        u32::try_from(dictionary.names.len())
            .map_err(|_| Error::InvalidInput("too many reference sequences".into()))?,
        terms,
        &posting_directory,
        postings_data,
        span_directory_bytes,
        starts_data,
        lengths_data,
        spans_per_block,
    )?;
    let index_bytes = serialized.len() as u64;
    atomic_write(destination, &serialized)?;
    let serialize_duration = serialize_started.elapsed();
    report_progress(
        build_options,
        BuildPhase::Serialize,
        records_processed,
        &bytes_read,
        started,
    );
    let timings = BuildTimings {
        scan: scan_duration.saturating_sub(scan_spill_duration),
        spill: spill_duration,
        merge: merge_duration,
        encode_postings: postings_duration,
        encode_spans: spans_duration,
        serialize: serialize_duration,
        total: started.elapsed(),
    };
    report_progress(
        build_options,
        BuildPhase::Complete,
        records_processed,
        &bytes_read,
        started,
    );
    let postings_after = posting_directory
        .iter()
        .map(|entry| u64::from(entry.compressed_length))
        .sum();
    let stats = IndexStats {
        records_processed,
        records_indexed,
        distinct_terms: term_count,
        unique_spans: coordinate_spans.len() as u64,
        postings: unique_postings,
        duplicate_postings_removed: pair_observations.saturating_sub(unique_postings),
        duplicate_spans_removed: span_observations.saturating_sub(coordinate_spans.len() as u64),
        index_bytes,
        postings_bytes_before_compression: postings_before,
        postings_bytes_after_compression: postings_after,
        span_bytes_fixed_width: span_encoding_stats.fixed_width_bytes,
        span_bytes_structural: span_encoding_stats
            .starts_structural_bytes
            .checked_add(span_encoding_stats.lengths_structural_bytes)
            .ok_or(Error::InvalidCoordinate)?,
        span_bytes_after_compression: span_encoding_stats
            .starts_compressed_bytes
            .checked_add(span_encoding_stats.lengths_compressed_bytes)
            .ok_or(Error::InvalidCoordinate)?,
        span_starts_bytes_before_compression: span_encoding_stats.starts_structural_bytes,
        span_starts_bytes_after_compression: span_encoding_stats.starts_compressed_bytes,
        span_lengths_bytes_before_compression: span_encoding_stats.lengths_structural_bytes,
        span_lengths_bytes_after_compression: span_encoding_stats.lengths_compressed_bytes,
        delta_start_blocks: span_encoding_stats.delta_start_blocks,
        length_varint_blocks: span_encoding_stats.length_varint_blocks,
        length_for_blocks: span_encoding_stats.length_for_blocks,
        bytes_per_term: index_bytes as f64 / (term_count.max(1) as f64),
        bytes_per_posting: index_bytes as f64 / (unique_postings.max(1) as f64),
        bytes_per_unique_span: index_bytes as f64 / (coordinate_spans.len().max(1) as f64),
        timings,
        peak_working_set_bytes: collector.peak_working_set_bytes(),
    };
    debug_assert_eq!(source_bytes, bytes_read.load(Ordering::Relaxed));
    Ok(stats)
}
