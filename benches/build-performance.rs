//! Reproducible end-to-end builder benchmark.
//!
//! This is intentionally an ignored custom benchmark. It compares the
//! single-pass no-spill working-set baseline with the bounded spill/merge
//! configuration
//! on a deterministic 100,000-record BGZF GFF3/TBI fixture.

use std::{fs::File, io::Write, time::Instant};

use gai::{BuildOptions, NameIndexOptions, build_name_index_with_options};
use noodles::{
    bgzf,
    core::Position,
    csi::{self, binning_index::index::reference_sequence::bin::Chunk},
    tabix,
};
use tempfile::tempdir;

const RECORD_COUNT: usize = 100_000;

fn main() {
    if cfg!(debug_assertions) {
        return;
    }

    let directory = tempdir().expect("should create benchmark directory");
    let source_path = directory.path().join("build-performance.gff3.gz");
    let coordinate_index_path = directory.path().join("build-performance.gff3.gz.tbi");
    let mut writer = File::create(&source_path)
        .map(bgzf::io::Writer::new)
        .expect("should create BGZF source");
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
    let mut uncompressed_bytes = "##gff-version 3\n".len() as u64;
    writeln!(writer, "##gff-version 3").expect("should write GFF header");

    let fixture_started = Instant::now();
    for index in 0..RECORD_COUNT {
        let reference_id = index / 25_000;
        let row = index % 25_000;
        let reference = format!("chr{}", reference_id + 1);
        let start = row as u64 * 7 + 1;
        let end = start + 4;
        let line = format!(
            "{reference}\tbenchmark\tgene\t{}\t{}\t.\t+\t.\tName=term{:05};Alias=shared{:03}",
            start,
            end,
            index % 10_000,
            index % 128
        );
        uncompressed_bytes += line.len() as u64 + 1;
        let start_position = writer.virtual_position();
        writeln!(writer, "{line}").expect("should write GFF record");
        let end_position = writer.virtual_position();
        indexer
            .add_record(
                &reference,
                Position::try_from(start as usize).expect("valid start"),
                Position::try_from(end as usize).expect("valid end"),
                Chunk::new(start_position, end_position),
            )
            .expect("should add TBI record");
    }
    writer.finish().expect("should finish BGZF source");
    let index = indexer.build();
    let mut index_writer = File::create(&coordinate_index_path)
        .map(tabix::io::Writer::new)
        .expect("should create TBI");
    index_writer.write_index(&index).expect("should write TBI");
    drop(index_writer);

    let options = NameIndexOptions::new(["Name", "Alias"], false)
        .expect("should configure benchmark attributes");
    let baseline_path = directory.path().join("baseline.gai");
    let bounded_path = directory.path().join("bounded.gai");
    let baseline_options = BuildOptions::default()
        .with_memory_budget(1 << 30)
        .with_compression_threads(4)
        .with_bgzf_threads(4);
    let bounded_options = BuildOptions::default()
        .with_memory_budget(1 << 20)
        .with_compression_threads(4)
        .with_bgzf_threads(4);

    let baseline_started = Instant::now();
    let baseline = build_name_index_with_options(
        &source_path,
        &coordinate_index_path,
        &baseline_path,
        &options,
        &baseline_options,
    )
    .expect("should build baseline GAI");
    let baseline_elapsed = baseline_started.elapsed();

    let bounded_started = Instant::now();
    let bounded = build_name_index_with_options(
        &source_path,
        &coordinate_index_path,
        &bounded_path,
        &options,
        &bounded_options,
    )
    .expect("should build bounded GAI");
    let bounded_elapsed = bounded_started.elapsed();

    println!("fixture records: {RECORD_COUNT}");
    println!("fixture uncompressed GFF bytes: {uncompressed_bytes}");
    println!(
        "fixture compressed GFF bytes: {}",
        source_path.metadata().unwrap().len()
    );
    println!("fixture generation: {:?}", fixture_started.elapsed());
    println!(
        "single-pass no-spill baseline: elapsed {:?}, records/sec {:.1}, peak working-set proxy {}, spill {:?}, scan {:?}, merge {:?}, postings {:?}, spans {:?}, serialize {:?}",
        baseline_elapsed,
        baseline.records_processed as f64 / baseline_elapsed.as_secs_f64(),
        baseline.peak_working_set_bytes,
        baseline.timings.spill,
        baseline.timings.scan,
        baseline.timings.merge,
        baseline.timings.encode_postings,
        baseline.timings.encode_spans,
        baseline.timings.serialize,
    );
    println!(
        "bounded spill (1 MiB): elapsed {:?}, records/sec {:.1}, peak working-set proxy {}, spill {:?}, scan {:?}, merge {:?}, postings {:?}, spans {:?}, serialize {:?}",
        bounded_elapsed,
        bounded.records_processed as f64 / bounded_elapsed.as_secs_f64(),
        bounded.peak_working_set_bytes,
        bounded.timings.spill,
        bounded.timings.scan,
        bounded.timings.merge,
        bounded.timings.encode_postings,
        bounded.timings.encode_spans,
        bounded.timings.serialize,
    );
    println!(
        "bounded vs baseline: throughput {:.2}x, peak working-set {:.2}x, byte-identical {}",
        baseline_elapsed.as_secs_f64() / bounded_elapsed.as_secs_f64(),
        bounded.peak_working_set_bytes as f64 / baseline.peak_working_set_bytes.max(1) as f64,
        std::fs::read(&baseline_path).expect("read baseline")
            == std::fs::read(&bounded_path).expect("read bounded"),
    );
}
