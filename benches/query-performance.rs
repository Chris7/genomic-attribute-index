//! Real-fixture indexed query benchmark.
//!
//! This intentionally uses the repository's GENCODE fixture and never writes
//! to it. Run with `cargo bench --bench query-performance --all-features`.

use std::{path::Path, time::Instant};

use gni::IndexedGff;

fn main() {
    if cfg!(debug_assertions) {
        return;
    }

    let source = Path::new("fixtures/gencode_sorted.gff.gz");
    let coordinate_index = Path::new("fixtures/gencode_sorted.gff.gz.tbi");
    let gni = Path::new("fixtures/gencode_sorted.gff.gz.gni");

    if !(source.exists() && coordinate_index.exists() && gni.exists()) {
        eprintln!(
            "query-performance: skipping; expected fixture files are absent under {}",
            source.display()
        );
        return;
    }

    let open_started = Instant::now();
    let mut indexed = IndexedGff::open(source, coordinate_index, gni)
        .expect("should open the real GENCODE fixture");
    let open_elapsed = open_started.elapsed();

    let query_started = Instant::now();
    let (records, stats) = indexed
        .query_name_with_stats("brca1")
        .expect("BRCA1 query should succeed");
    let query_elapsed = query_started.elapsed();

    println!("term: brca1");
    println!("records: {}", records.len());
    println!("open/fingerprint elapsed: {:?}", open_elapsed);
    println!("query elapsed: {:?}", query_elapsed);
    println!(
        "records/sec: {:.1}",
        records.len() as f64 / query_elapsed.as_secs_f64()
    );
    println!("requested spans: {}", stats.requested_spans);
    println!(
        "distinct span blocks decoded: {}",
        stats.distinct_span_blocks_decoded
    );
    println!("exact interval queries: {}", stats.exact_interval_queries);
    println!("raw chunks: {}", stats.raw_chunks);
    println!("merged chunks: {}", stats.merged_chunks);
    println!(
        "unique candidate records parsed: {}",
        stats.unique_candidate_records
    );
    println!("matching records: {}", stats.matching_records);
    println!("uncompressed bytes read: {}", stats.bytes_read);
}
