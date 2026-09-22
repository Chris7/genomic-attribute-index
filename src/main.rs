use std::{path::PathBuf, process::ExitCode};

use clap::{Args, Parser, Subcommand, ValueEnum};
use gai::{
    BuildOptions, IndexedGff, MatchMode, NameIndexOptions, Result, build_name_index_with_options,
    sort_file,
};

#[derive(Debug, Parser)]
#[command(name = "gai", about = "GFF3 indexing and attribute-name queries")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Sort GFF/GFF3 or BED records and write the result to stdout.
    #[command(name = "sort")]
    Sort(SortArgs),
    /// Build a deterministic Genomic Attribute Index (GAI) for configured GFF3 attributes.
    #[command(name = "build-index")]
    Build(BuildIndexArgs),
    /// Query configured attribute values through TBI/CSI and print GFF3 records.
    #[command(name = "query-index")]
    Query(QueryIndexArgs),
    /// Display GAI format, normalization, fingerprint, and block metadata.
    #[command(name = "inspect-index")]
    Inspect(InspectArgs),
}

#[derive(Debug, Args)]
struct SortArgs {
    /// Input .gff, .gff3, or .bed path; sorted output is written to stdout.
    input: PathBuf,
    /// Sort on disk; use for files that may fill memory
    #[arg(long, alias = "ds")]
    disk_sort: bool,
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
    /// Destination GAI path. Defaults to <input>.gai.
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
    /// Explicit GAI path. Defaults to <input>.gai.
    #[arg(long)]
    gai: Option<PathBuf>,
    /// Match complete values exactly or stream values beginning with the query.
    #[arg(long = "match", value_enum, default_value_t = QueryMatch::Exact)]
    match_mode: QueryMatch,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum QueryMatch {
    Exact,
    Prefix,
}

impl From<QueryMatch> for MatchMode {
    fn from(value: QueryMatch) -> Self {
        match value {
            QueryMatch::Exact => Self::Exact,
            QueryMatch::Prefix => Self::Prefix,
        }
    }
}

#[derive(Debug, Args)]
struct InspectArgs {
    /// GAI path.
    input: PathBuf,
}

fn discover_coordinate_index(input: &std::path::Path) -> Result<PathBuf> {
    let tbi = PathBuf::from(format!("{}.tbi", input.display()));
    let csi = PathBuf::from(format!("{}.csi", input.display()));
    match (tbi.exists(), csi.exists()) {
        (true, false) => Ok(tbi),
        (false, true) => Ok(csi),
        (true, true) => Err(gai::Error::InvalidInput(format!(
            "both {} and {} exist; pass --coordinate-index explicitly",
            tbi.display(),
            csi.display()
        ))),
        (false, false) => Err(gai::Error::InvalidInput(format!(
            "could not discover {}.tbi or {}.csi",
            input.display(),
            input.display()
        ))),
    }
}

fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Sort(arguments) => {
            sort_file(arguments.input, arguments.disk_sort, std::io::stdout())?;
        }
        Command::Build(arguments) => {
            let coordinate_index = arguments
                .coordinate_index
                .unwrap_or(discover_coordinate_index(&arguments.input)?);
            let output = arguments
                .output
                .unwrap_or_else(|| PathBuf::from(format!("{}.gai", arguments.input.display())));
            let options = NameIndexOptions::new(arguments.attributes, arguments.case_sensitive)?;
            let mut build_options = BuildOptions::default()
                .with_memory_budget(arguments.memory_budget)
                .with_progress(|progress| {
                    eprintln!(
                        "gai: {:?}: records={} bytes={} elapsed={:.1}s",
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
                "gai: complete: {:.0} records/s; phases scan={:.3}s spill={:.3}s merge={:.3}s postings={:.3}s spans={:.3}s serialize={:.3}s total={:.3}s; peak working set={} bytes",
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
            let gai = arguments
                .gai
                .unwrap_or_else(|| PathBuf::from(format!("{}.gai", arguments.input.display())));
            let mut indexed = IndexedGff::open(&arguments.input, coordinate_index, gai)?;
            for record in
                indexed.query_name_with_mode(&arguments.term, arguments.match_mode.into())?
            {
                println!("{}", record.raw_line);
            }
        }
        Command::Inspect(arguments) => {
            let metadata = gai::NameIndexReader::open(arguments.input)?.inspect();
            println!(
                "GAI version: {}.{}",
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
                "span encodings (delta-varint starts, varint/FOR lengths): {}, {}/{}",
                metadata.delta_start_blocks,
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
                "span directory/starts/lengths bytes: {}/{}/{}",
                metadata.span_directory_bytes,
                metadata.starts_data_bytes,
                metadata.lengths_data_bytes
            );
            println!(
                "starts compression (uncompressed/compressed): {}/{} bytes ({:.2}x)",
                metadata.starts_uncompressed_bytes,
                metadata.starts_data_bytes,
                compression_ratio(
                    metadata.starts_uncompressed_bytes,
                    metadata.starts_data_bytes
                )
            );
            println!(
                "lengths compression (uncompressed/compressed): {}/{} bytes ({:.2}x)",
                metadata.lengths_uncompressed_bytes,
                metadata.lengths_data_bytes,
                compression_ratio(
                    metadata.lengths_uncompressed_bytes,
                    metadata.lengths_data_bytes
                )
            );
            println!(
                "zstd blocks (postings/starts/lengths): {}/{}/{}",
                metadata.compressed_postings_blocks,
                metadata.compressed_start_blocks,
                metadata.compressed_length_blocks
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
            eprintln!("gai: {error}");
            ExitCode::from(1)
        }
    }
}
