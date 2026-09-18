use std::{
    collections::BTreeSet,
    fs::File,
    io::Write,
    time::{Duration, Instant},
};

use gni::{IndexedGff, NameIndexOptions, build_name_index_with_span_block_size};
use noodles::{
    bgzf,
    core::Position,
    csi::{self, binning_index::index::reference_sequence::bin::Chunk},
    tabix,
};
use tempfile::tempdir;

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Span {
    reference: u32,
    start: u64,
    length: u64,
}

struct SyntheticRecord {
    span: Span,
    name_term: String,
    alias_term: String,
}

#[derive(Clone, Copy)]
struct BlockMeasurement {
    block_size: usize,
    structural_bytes: u64,
    compressed_bytes: u64,
    block_count: u64,
    gni_bytes: u64,
    lookup_bytes: u64,
    lookup_latency: Duration,
}

fn varint(mut value: u64, output: &mut Vec<u8>) {
    while value >= 0x80 {
        output.push(value as u8 | 0x80);
        value >>= 7;
    }
    output.push(value as u8);
}

fn for_stream(values: &[u64]) -> Vec<u8> {
    let base = values.iter().copied().min().unwrap_or(0);
    let maximum = values
        .iter()
        .map(|value| value.saturating_sub(base))
        .max()
        .unwrap_or(0);
    let width = if maximum == 0 {
        0
    } else {
        64 - maximum.leading_zeros()
    } as usize;
    let bit_count = values.len() * width;
    let byte_count = bit_count / 8 + usize::from(!bit_count.is_multiple_of(8));
    let mut output = vec![0_u8; byte_count];
    for (index, value) in values.iter().enumerate() {
        let relative = value.saturating_sub(base);
        for bit in 0..width {
            if relative & (1_u64 << bit) != 0 {
                let target = index * width + bit;
                output[target / 8] |= 1 << (target % 8);
            }
        }
    }
    output
}

fn measure_queries(
    source_path: &std::path::Path,
    coordinate_index_path: &std::path::Path,
    gni_path: &std::path::Path,
    query_terms: &[String],
) -> (u64, Duration) {
    let mut reader = IndexedGff::open(source_path, coordinate_index_path, gni_path)
        .expect("should open synthetic indexed GFF");
    let started = Instant::now();
    let mut lookup_bytes = 0_u64;
    for term in query_terms {
        let (_, stats) = reader
            .name_index()
            .lookup_spans_with_stats(term)
            .expect("should resolve a known synthetic term");
        lookup_bytes = lookup_bytes
            .checked_add(stats.postings_bytes_decompressed)
            .and_then(|value| value.checked_add(stats.span_bytes_decompressed))
            .expect("lookup byte count should fit in u64");
        let records = reader
            .query_name(term)
            .expect("should query a known synthetic term");
        assert!(!records.is_empty(), "sampled term must be present: {term}");
    }
    let elapsed = started.elapsed();
    (
        lookup_bytes / query_terms.len().max(1) as u64,
        elapsed / query_terms.len().max(1) as u32,
    )
}

fn main() {
    // Cargo includes custom `harness = false` bench binaries in
    // `cargo test --all-targets`; keep that ordinary test command from
    // executing this intentionally expensive workload. `cargo bench` builds
    // with optimizations and reaches the measurement body below.
    if cfg!(debug_assertions) {
        return;
    }
    let mut seed = 0x4753_4e49_u64;
    let mut records: Vec<SyntheticRecord> = Vec::with_capacity(7_200);
    for index in 0..7_200 {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let reference = if index < 6_000 {
            0
        } else if index < 6_800 {
            1
        } else {
            2
        };
        let start = if index % 2 == 0 {
            (index / 2) as u64 * 37 + (seed % 11)
        } else {
            (seed % 10_000_000) + 1_000_000
        };
        let length = if index % 5 == 0 {
            1_000 + seed % 50_000
        } else {
            10 + seed % 200
        };
        let span = Span {
            reference,
            start,
            length,
        };
        let duplicate = index > 0 && index % 257 == 0;
        let duplicate_term = duplicate.then(|| format!("duplicate-span-{:04}", index - 1));
        let name_term = if let Some(duplicate_term) = &duplicate_term {
            duplicate_term.clone()
        } else if index % 3 == 0 {
            format!("term{:04}", index)
        } else {
            format!("random-{:016x}", seed)
        };
        if let Some(duplicate_term) = duplicate_term {
            // Keep the original and duplicate records under one known term so
            // the benchmark includes one-to-many postings for an identical
            // span, not merely duplicate coordinates in the encoding input.
            records[index - 1].name_term = duplicate_term.clone();
        }
        records.push(SyntheticRecord {
            span: if duplicate {
                records[index - 1].span
            } else {
                span
            },
            name_term,
            alias_term: format!("shared-{:03}", index % 128),
        });
    }
    // Source records are coordinate sorted, including the intentionally
    // identical records. Stable sorting keeps their source order visible to
    // the query measurements.
    records.sort_by_key(|record| record.span);
    let unique_spans = records
        .iter()
        .map(|record| record.span)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();

    let started = Instant::now();
    let mut fixed = Vec::new();
    for span in &unique_spans {
        fixed.extend_from_slice(&span.reference.to_le_bytes());
        fixed.extend_from_slice(&span.start.to_le_bytes());
        fixed.extend_from_slice(&span.length.to_le_bytes());
    }
    // Per-record interleaved tuples retain duplicate observations. This is a
    // useful contrast with the globally deduplicated row representation below.
    let mut interleaved = Vec::new();
    for record in &records {
        let span = record.span;
        varint(span.reference as u64, &mut interleaved);
        varint(span.start, &mut interleaved);
        varint(span.length, &mut interleaved);
    }
    let mut row = Vec::new();
    for span in &unique_spans {
        varint(span.reference as u64, &mut row);
        varint(span.start, &mut row);
        varint(span.length, &mut row);
    }
    let mut columnar = Vec::new();
    let mut previous_reference = None;
    let mut previous_start = 0_u64;
    for span in &unique_spans {
        if previous_reference != Some(span.reference) {
            previous_start = span.start;
            previous_reference = Some(span.reference);
            // The first start for each reference is represented in the
            // reference-local block directory, matching the GNI layout.
            continue;
        }
        varint(span.start - previous_start, &mut columnar);
        previous_start = span.start;
    }
    for span in &unique_spans {
        varint(span.length, &mut columnar);
    }
    let mut frame = Vec::new();
    let mut reference_start = 0;
    while reference_start < unique_spans.len() {
        let reference = unique_spans[reference_start].reference;
        let mut reference_end = reference_start + 1;
        while reference_end < unique_spans.len()
            && unique_spans[reference_end].reference == reference
        {
            reference_end += 1;
        }
        let starts = unique_spans[reference_start..reference_end]
            .iter()
            .map(|span| span.start)
            .collect::<Vec<_>>();
        let lengths = unique_spans[reference_start..reference_end]
            .iter()
            .map(|span| span.length)
            .collect::<Vec<_>>();
        frame.extend_from_slice(&for_stream(&starts));
        frame.extend_from_slice(&for_stream(&lengths));
        reference_start = reference_end;
    }
    let compressed = zstd::bulk::compress(&frame, 3).expect("zstd should compress benchmark data");
    println!(
        "workload: {} records, {} unique spans, seed 0x47534e49",
        records.len(),
        unique_spans.len()
    );
    println!("interleaved fixed-width triples: {}", fixed.len());
    println!(
        "interleaved unsigned varints (records): {}",
        interleaved.len()
    );
    println!(
        "row-oriented unsigned varints (unique spans): {}",
        row.len()
    );
    println!("columnar start delta-varints: {}", columnar.len());
    println!("columnar frame-of-reference: {}", frame.len());
    println!("columnar FOR + zstd blocks: {}", compressed.len());

    let directory = tempdir().expect("should create benchmark directory");
    let source_path = directory.path().join("synthetic.gff3.gz");
    let coordinate_index_path = directory.path().join("synthetic.gff3.gz.tbi");
    let mut writer = File::create(&source_path)
        .map(bgzf::io::Writer::new)
        .expect("should create synthetic BGZF");
    let mut uncompressed_gff_bytes = 0_u64;
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
    let mut query_terms = Vec::new();
    for (index, record) in records.iter().enumerate() {
        let reference = format!("chr{}", record.span.reference + 1);
        let start = record.span.start + 1;
        let end = start + record.span.length - 1;
        let line = format!(
            "{reference}\tsynthetic\tgene\t{start}\t{end}\t.\t+\t.\tName={};Alias={}",
            record.name_term, record.alias_term
        );
        uncompressed_gff_bytes += line.len() as u64 + 1;
        let chunk_start = writer.virtual_position();
        writeln!(writer, "{line}").expect("should write synthetic GFF");
        let chunk_end = writer.virtual_position();
        indexer
            .add_record(
                &reference,
                Position::try_from(start as usize).expect("valid synthetic start"),
                Position::try_from(end as usize).expect("valid synthetic end"),
                Chunk::new(chunk_start, chunk_end),
            )
            .expect("should add synthetic record to TBI");
        if index % 149 == 0 {
            query_terms.push(record.name_term.clone());
        }
        if index % 997 == 0 {
            query_terms.push(record.alias_term.clone());
        }
    }
    // A duplicate-span term is always sampled even if sorting moves its
    // records away from a step boundary.
    if let Some(record) = records
        .iter()
        .find(|record| record.name_term.starts_with("duplicate-span-"))
    {
        query_terms.push(record.name_term.clone());
    }
    query_terms.sort();
    query_terms.dedup();
    writer.finish().expect("should finish synthetic BGZF");
    let index = indexer.build();
    let mut index_writer = File::create(&coordinate_index_path)
        .map(tabix::io::Writer::new)
        .expect("should create synthetic TBI");
    index_writer
        .write_index(&index)
        .expect("should write synthetic TBI");
    drop(index_writer);

    let options = NameIndexOptions::new(["Name", "Alias"], false)
        .expect("should configure synthetic attributes");
    let mut block_measurements = Vec::new();
    let mut selected_stats = None;
    let mut selected_metadata = None;
    for block_size in [1_024_usize, 4_096, 16_384] {
        let block_path = directory.path().join(format!("synthetic-{block_size}.gni"));
        let block_stats = build_name_index_with_span_block_size(
            &source_path,
            &coordinate_index_path,
            &block_path,
            &options,
            block_size,
        )
        .expect("should build synthetic GNI");
        let block_metadata = gni::NameIndexReader::open(&block_path)
            .expect("should open synthetic block-size GNI")
            .inspect();
        let (lookup_bytes, lookup_latency) = measure_queries(
            &source_path,
            &coordinate_index_path,
            &block_path,
            &query_terms,
        );
        block_measurements.push(BlockMeasurement {
            block_size,
            structural_bytes: block_stats.span_bytes_structural,
            compressed_bytes: block_stats.span_bytes_after_compression,
            block_count: block_metadata.span_block_count,
            gni_bytes: block_stats.index_bytes,
            lookup_bytes,
            lookup_latency,
        });
        if block_size == 4_096 {
            selected_stats = Some(block_stats);
            selected_metadata = Some(block_metadata);
        }
    }
    for measurement in block_measurements {
        println!(
            "span blocks {}: structural {}, zstd {}, blocks {}, GNI bytes {}, avg lookup bytes {}, avg lookup latency {:?}",
            measurement.block_size,
            measurement.structural_bytes,
            measurement.compressed_bytes,
            measurement.block_count,
            measurement.gni_bytes,
            measurement.lookup_bytes,
            measurement.lookup_latency
        );
    }
    let stats = selected_stats.expect("selected block size should be measured");
    let metadata = selected_metadata.expect("selected block metadata should be available");
    println!("query terms (all known): {}", query_terms.len());
    println!("uncompressed GFF bytes: {}", uncompressed_gff_bytes);
    println!(
        "BGZF GFF bytes: {}",
        source_path.metadata().expect("source metadata").len()
    );
    println!(
        "TBI bytes: {}",
        coordinate_index_path
            .metadata()
            .expect("TBI metadata")
            .len()
    );
    println!("GNI bytes: {}", stats.index_bytes);
    println!("FST bytes: {}", metadata.term_dictionary_bytes);
    println!("postings bytes: {}", metadata.postings_data_bytes);
    println!("span-table bytes: {}", metadata.span_data_bytes);
    println!(
        "block-directory bytes: {}",
        metadata.postings_directory_bytes + metadata.span_directory_bytes
    );
    println!(
        "bytes per term: {:.2}",
        stats.index_bytes as f64 / stats.distinct_terms.max(1) as f64
    );
    println!(
        "bytes per posting: {:.2}",
        stats.index_bytes as f64 / stats.postings.max(1) as f64
    );
    println!(
        "bytes per unique span: {:.2}",
        stats.index_bytes as f64 / stats.unique_spans.max(1) as f64
    );
    println!("elapsed: {:?}", started.elapsed());
}
