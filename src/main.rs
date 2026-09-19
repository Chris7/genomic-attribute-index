use std::{path::PathBuf, process::ExitCode};

use clap::{Args, Parser, Subcommand};
use gni::{BuildOptions, IndexedGff, NameIndexOptions, Result, build_name_index_with_options};

#[derive(Debug, Parser)]
#[command(name = "gni", about = "GFF3 indexing and exact attribute-name queries")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Build a deterministic GFF Name Index (GNI).
    #[command(name = "build-index")]
    Build(BuildIndexArgs),
    /// Query configured attribute values through TBI/CSI and print GFF3 records.
    #[command(name = "query-index")]
    Query(QueryIndexArgs),
    /// Display GNI format, normalization, fingerprint, and block metadata.
    #[command(name = "inspect-index")]
    Inspect(InspectArgs),
}

#[derive(Debug, Args)]
struct BuildIndexArgs {
    /// BGZF GFF3 source.
    input: PathBuf,
    /// Repeatable configured GFF3 attribute tag. At least one is required.
    #[arg(long = "attribute", required = true)]
    attributes: Vec<String>,
    /// Explicit TBI or CSI path. If omitted, discover an unambiguous sibling.
    #[arg(long = "coordinate-index")]
    coordinate_index: Option<PathBuf>,
    /// Destination GNI path. Defaults to <input>.gni.
    #[arg(long)]
    output: Option<PathBuf>,
    /// Preserve case after trimming values.
    #[arg(long)]
    case_sensitive: bool,
    /// Approximate bounded scan working-set budget in bytes.
    #[arg(long, default_value_t = 64 * 1024 * 1024)]
    memory_budget: usize,
    /// Number of independent block-compression workers.
    #[arg(long)]
    compression_threads: Option<usize>,
    /// Number of BGZF decompression workers.
    #[arg(long)]
    bgzf_threads: Option<usize>,
}

#[derive(Debug, Args)]
struct QueryIndexArgs {
    /// BGZF GFF3 source.
    input: PathBuf,
    /// Query term before normalization.
    term: String,
    /// Explicit TBI or CSI path. If omitted, discover an unambiguous sibling.
    #[arg(long = "coordinate-index")]
    coordinate_index: Option<PathBuf>,
    /// Explicit GNI path. Defaults to <input>.gni.
    #[arg(long)]
    gni: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct InspectArgs {
    /// GNI path.
    input: PathBuf,
}

fn discover_coordinate_index(input: &std::path::Path) -> Result<PathBuf> {
    let tbi = PathBuf::from(format!("{}.tbi", input.display()));
    let csi = PathBuf::from(format!("{}.csi", input.display()));
    match (tbi.exists(), csi.exists()) {
        (true, false) => Ok(tbi),
        (false, true) => Ok(csi),
        (true, true) => Err(gni::Error::InvalidInput(format!(
            "both {} and {} exist; pass --coordinate-index explicitly",
            tbi.display(),
            csi.display()
        ))),
        (false, false) => Err(gni::Error::InvalidInput(format!(
            "could not discover {}.tbi or {}.csi",
            input.display(),
            input.display()
        ))),
    }
}

fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Build(arguments) => {
            let coordinate_index = arguments
                .coordinate_index
                .unwrap_or(discover_coordinate_index(&arguments.input)?);
            let output = arguments
                .output
                .unwrap_or_else(|| PathBuf::from(format!("{}.gni", arguments.input.display())));
            let options = NameIndexOptions::new(arguments.attributes, arguments.case_sensitive)?;
            let mut build_options = BuildOptions::default()
                .with_memory_budget(arguments.memory_budget)
                .with_progress(|progress| {
                    eprintln!(
                        "gni: {:?}: records={} bytes={} elapsed={:.1}s",
                        progress.phase,
                        progress.records_processed,
                        progress.bytes_read,
                        progress.elapsed.as_secs_f64()
                    );
                });
            if let Some(threads) = arguments.compression_threads {
                build_options = build_options.with_compression_threads(threads);
            }
            if let Some(threads) = arguments.bgzf_threads {
                build_options = build_options.with_bgzf_threads(threads);
            }
            let stats = build_name_index_with_options(
                &arguments.input,
                &coordinate_index,
                &output,
                &options,
                &build_options,
            )?;
            println!("wrote {}", output.display());
            println!("records processed: {}", stats.records_processed);
            println!("records indexed: {}", stats.records_indexed);
            println!("distinct terms: {}", stats.distinct_terms);
            println!("unique spans: {}", stats.unique_spans);
            println!("postings: {}", stats.postings);
            println!("index bytes: {}", stats.index_bytes);
            eprintln!(
                "gni: complete: {:.0} records/s; phases scan={:.3}s spill={:.3}s merge={:.3}s postings={:.3}s spans={:.3}s serialize={:.3}s total={:.3}s; peak working set={} bytes",
                stats.records_processed as f64
                    / stats.timings.total.as_secs_f64().max(f64::MIN_POSITIVE),
                stats.timings.scan.as_secs_f64(),
                stats.timings.spill.as_secs_f64(),
                stats.timings.merge.as_secs_f64(),
                stats.timings.encode_postings.as_secs_f64(),
                stats.timings.encode_spans.as_secs_f64(),
                stats.timings.serialize.as_secs_f64(),
                stats.timings.total.as_secs_f64(),
                stats.peak_working_set_bytes,
            );
        }
        Command::Query(arguments) => {
            let coordinate_index = arguments
                .coordinate_index
                .unwrap_or(discover_coordinate_index(&arguments.input)?);
            let gni = arguments
                .gni
                .unwrap_or_else(|| PathBuf::from(format!("{}.gni", arguments.input.display())));
            let mut indexed = IndexedGff::open(&arguments.input, coordinate_index, gni)?;
            for record in indexed.query_name(&arguments.term)? {
                println!("{}", record.raw_line);
            }
        }
        Command::Inspect(arguments) => {
            let metadata = gni::NameIndexReader::open(arguments.input)?.inspect();
            println!(
                "GNI version: {}.{}",
                metadata.major_version, metadata.minor_version
            );
            println!("coordinate convention: zero-based half-open start + length");
            println!("attributes: {}", metadata.attributes.join(", "));
            println!(
                "normalization: {}",
                if metadata.case_sensitive {
                    "trim Unicode whitespace"
                } else {
                    "trim Unicode whitespace; ASCII lowercase"
                }
            );
            println!("terms: {}", metadata.term_count);
            println!("reference sequences: {}", metadata.reference_count);
            println!("unique spans: {}", metadata.unique_span_count);
            println!("postings: {}", metadata.posting_count);
            println!("postings blocks: {}", metadata.postings_block_count);
            println!("span blocks: {}", metadata.span_block_count);
            println!("span block target: {}", metadata.span_block_size);
            println!(
                "span encodings (delta/FOR starts, varint/FOR lengths): {}/{}, {}/{}",
                metadata.delta_start_blocks,
                metadata.for_start_blocks,
                metadata.varint_length_blocks,
                metadata.for_length_blocks
            );
            println!(
                "attribute section bytes: {}",
                metadata.attribute_section_bytes
            );
            println!("FST bytes: {}", metadata.term_dictionary_bytes);
            println!(
                "postings directory/data bytes: {}/{}",
                metadata.postings_directory_bytes, metadata.postings_data_bytes
            );
            println!(
                "postings compression (uncompressed/compressed): {}/{} bytes ({:.2}x)",
                metadata.postings_uncompressed_bytes,
                metadata.postings_data_bytes,
                compression_ratio(
                    metadata.postings_uncompressed_bytes,
                    metadata.postings_data_bytes
                )
            );
            println!(
                "span directory/data bytes: {}/{}",
                metadata.span_directory_bytes, metadata.span_data_bytes
            );
            println!(
                "span compression (uncompressed/compressed): {}/{} bytes ({:.2}x)",
                metadata.span_uncompressed_bytes,
                metadata.span_data_bytes,
                compression_ratio(metadata.span_uncompressed_bytes, metadata.span_data_bytes)
            );
            println!(
                "zstd blocks (postings/spans): {}/{}",
                metadata.compressed_postings_blocks, metadata.compressed_span_blocks
            );
            println!("file bytes: {}", metadata.file_size);
            println!("GFF SHA-256: {}", hex(&metadata.gff_fingerprint));
            println!(
                "coordinate-index SHA-256: {}",
                hex(&metadata.coordinate_index_fingerprint)
            );
            println!(
                "reference-dictionary SHA-256: {}",
                hex(&metadata.reference_dictionary_fingerprint)
            );
        }
    }
    Ok(())
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn compression_ratio(uncompressed: u64, compressed: u64) -> f64 {
    if compressed == 0 {
        0.0
    } else {
        uncompressed as f64 / compressed as f64
    }
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("gni: {error}");
            ExitCode::from(1)
        }
    }
}
