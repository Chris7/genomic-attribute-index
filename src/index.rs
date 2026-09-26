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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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

#[cfg(test)]
pub(crate) mod tests {
    use std::sync::{Arc, Mutex};

    use noodles::{
        bgzf,
        core::Position,
        csi::{
            self,
            binning_index::{
                Indexer,
                index::reference_sequence::{bin::Chunk, index::BinnedIndex},
            },
        },
        tabix,
    };
    use tempfile::tempdir;

    use super::*;
    #[tracing::instrument(level = "trace", skip_all)]
    pub(crate) fn fixture_records() -> [&'static str; 6] {
        [
            "##gff-version 3",
            "##sequence-region chr1 1 1000",
            "chr1\tsrc\tgene\t100\t150\t.\t+\t.\tName=BRCA1;Alias=BRCC1,RNF53;gene_name=BRCA1",
            "chr1\tsrc\tgene\t100\t150\t.\t+\t.\tAlias=BRCA1;ID=without_special_status",
            "chr1\tsrc\tgene\t200\t210\t.\t-\t.\tName=Other%20Gene",
            "chr2\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1",
        ]
    }
    #[tracing::instrument(level = "trace", skip_all)]
    pub(crate) fn write_fixture(directory: &Path) -> (PathBuf, PathBuf) {
        let source_path = directory.join("fixture.gff3.gz");
        let index_path = directory.join("fixture.gff3.gz.tbi");
        let records = fixture_records();
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .expect("should create BGZF source");
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
        for (record_index, line) in records.iter().enumerate() {
            if line.starts_with('#') {
                writeln!(writer, "{line}").expect("should write directive");
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let start = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
            let end = Position::try_from(fields[4].parse::<usize>().unwrap()).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").expect("should write record");
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    fields[0],
                    start,
                    end,
                    Chunk::new(start_position, end_position),
                )
                .expect("should index record");
            assert!(record_index > 0);
        }
        writer.finish().expect("should finish BGZF");
        let index = indexer.build();
        let mut index_writer = File::create(&index_path)
            .map(tabix::io::Writer::new)
            .expect("should create TBI");
        index_writer.write_index(&index).expect("should write TBI");
        (source_path, index_path)
    }
    #[tracing::instrument(level = "trace", skip_all)]
    pub(crate) fn write_tbi_lines(
        directory: &Path,
        stem: &str,
        lines: &[&str],
    ) -> (PathBuf, PathBuf) {
        let source_path = directory.join(format!("{stem}.gff3.gz"));
        let index_path = directory.join(format!("{stem}.gff3.gz.tbi"));
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .expect("should create BGZF source");
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
        for line in lines {
            if line.starts_with('#') {
                writeln!(writer, "{line}").expect("should write directive");
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let start = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
            let end = Position::try_from(fields[4].parse::<usize>().unwrap()).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").expect("should write record");
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    fields[0],
                    start,
                    end,
                    Chunk::new(start_position, end_position),
                )
                .expect("should index record");
        }
        writer.finish().expect("should finish BGZF");
        let index = indexer.build();
        let mut index_writer = File::create(&index_path)
            .map(tabix::io::Writer::new)
            .expect("should create TBI");
        index_writer.write_index(&index).expect("should write TBI");
        drop(index_writer);
        (source_path, index_path)
    }
    #[tracing::instrument(level = "trace", skip_all)]
    pub(crate) fn write_csi_fixture(directory: &Path) -> (PathBuf, PathBuf) {
        let source_path = directory.join("fixture-csi.gff3.gz");
        let index_path = directory.join("fixture-csi.gff3.gz.csi");
        let records = fixture_records();
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .expect("should create BGZF source");
        let mut names = csi::binning_index::index::header::ReferenceSequenceNames::new();
        names.insert("chr1".into());
        names.insert("chr2".into());
        let header = csi::binning_index::index::Header::builder()
            .set_format(Format::Generic(CoordinateSystem::Gff))
            .set_reference_sequence_names(names)
            .build();
        let mut indexer = Indexer::<BinnedIndex>::default().set_header(header);
        for line in records {
            if line.starts_with('#') {
                writeln!(writer, "{line}").expect("should write directive");
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let reference_id = if fields[0] == "chr1" { 0 } else { 1 };
            let start = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
            let end = Position::try_from(fields[4].parse::<usize>().unwrap()).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").expect("should write record");
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    Some((reference_id, start, end, true)),
                    Chunk::new(start_position, end_position),
                )
                .expect("should index record");
        }
        writer.finish().expect("should finish BGZF");
        let index = indexer.build(2);
        let mut index_writer =
            csi::io::Writer::new(File::create(&index_path).expect("should create CSI file"));
        index_writer.write_index(&index).expect("should write CSI");
        (source_path, index_path)
    }
    #[tracing::instrument(level = "trace", skip_all)]
    pub(crate) fn write_bed_fixture(directory: &Path) -> (PathBuf, PathBuf) {
        let source_path = directory.join("fixture.bed.gz");
        let index_path = directory.join("fixture.bed.gz.tbi");
        let records = [
            "chr1\t0\t10\tAlpha",
            "chr1\t20\t25\tBeta",
            "chr2\t5\t7\tAlpha",
            "chr2\t8\t9",
        ];
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .expect("should create BGZF BED source");
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::bed().build());
        for line in records {
            let fields = line.split('\t').collect::<Vec<_>>();
            let bed_start = fields[1].parse::<usize>().unwrap();
            let bed_end = fields[2].parse::<usize>().unwrap();
            let start = Position::try_from(bed_start + 1).unwrap();
            let end = Position::try_from(bed_end).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").expect("should write BED record");
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    fields[0],
                    start,
                    end,
                    Chunk::new(start_position, end_position),
                )
                .expect("should index BED record");
        }
        writer.finish().expect("should finish BGZF BED source");
        let index = indexer.build();
        let mut index_writer = File::create(&index_path)
            .map(tabix::io::Writer::new)
            .expect("should create BED TBI");
        index_writer
            .write_index(&index)
            .expect("should write BED TBI");
        (source_path, index_path)
    }
    #[tracing::instrument(level = "trace", skip_all)]
    pub(crate) fn mutate_section(
        bytes: &mut [u8],
        kind: SectionKind,
        mutate: impl FnOnce(&mut [u8]),
    ) {
        let directory_offset =
            usize::try_from(u64::from_le_bytes(bytes[176..184].try_into().unwrap())).unwrap();
        let section_count = u32::from_le_bytes(bytes[24..28].try_into().unwrap()) as usize;
        for index in 0..section_count {
            let entry = directory_offset + index * DIRECTORY_ENTRY_SIZE;
            let entry_kind = u32::from_le_bytes(bytes[entry..entry + 4].try_into().unwrap());
            if entry_kind != kind as u32 {
                continue;
            }
            let section_offset = usize::try_from(u64::from_le_bytes(
                bytes[entry + 8..entry + 16].try_into().unwrap(),
            ))
            .unwrap();
            let section_length = usize::try_from(u64::from_le_bytes(
                bytes[entry + 16..entry + 24].try_into().unwrap(),
            ))
            .unwrap();
            mutate(&mut bytes[section_offset..section_offset + section_length]);
            let section_checksum =
                checksum(&bytes[section_offset..section_offset + section_length]);
            bytes[entry + 32..entry + 36].copy_from_slice(&section_checksum.to_le_bytes());
            return;
        }
        panic!("missing section {kind:?}");
    }
    #[tracing::instrument(level = "trace", skip_all)]
    pub(crate) fn mutate_section_directory_item_count(
        bytes: &mut [u8],
        kind: SectionKind,
        count: u64,
    ) {
        let directory_offset =
            usize::try_from(u64::from_le_bytes(bytes[176..184].try_into().unwrap())).unwrap();
        let section_count = u32::from_le_bytes(bytes[24..28].try_into().unwrap()) as usize;
        for index in 0..section_count {
            let entry = directory_offset + index * DIRECTORY_ENTRY_SIZE;
            let entry_kind = u32::from_le_bytes(bytes[entry..entry + 4].try_into().unwrap());
            if entry_kind == kind as u32 {
                bytes[entry + 24..entry + 32].copy_from_slice(&count.to_le_bytes());
                return;
            }
        }
        panic!("missing section {kind:?}");
    }
    #[tracing::instrument(level = "trace", skip_all)]
    pub(crate) fn mutate_section_directory_offset(
        bytes: &mut [u8],
        kind: SectionKind,
        offset: u64,
    ) {
        let directory_offset =
            usize::try_from(u64::from_le_bytes(bytes[176..184].try_into().unwrap())).unwrap();
        let section_count = u32::from_le_bytes(bytes[24..28].try_into().unwrap()) as usize;
        for index in 0..section_count {
            let entry = directory_offset + index * DIRECTORY_ENTRY_SIZE;
            let entry_kind = u32::from_le_bytes(bytes[entry..entry + 4].try_into().unwrap());
            if entry_kind == kind as u32 {
                bytes[entry + 8..entry + 16].copy_from_slice(&offset.to_le_bytes());
                return;
            }
        }
        panic!("missing section {kind:?}");
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_span_blocks_are_reference_local_and_deterministic() {
        let mut spans = Vec::new();
        for index in 0..5_000_u64 {
            spans.push(SpanKey {
                reference_id: 0,
                start: index * 10,
                length: 5,
            });
        }
        for index in 0..3_u64 {
            spans.push(SpanKey {
                reference_id: 1,
                start: index * 10,
                length: 7,
            });
        }
        let first = encode_span_blocks_with_threads(&spans, 1_024, 1).unwrap();
        let second = encode_span_blocks_with_threads(&spans, 1_024, 1).unwrap();
        assert_eq!(first.0.len(), 6);
        assert_eq!(first.0, second.0);
        assert_eq!(first.1, second.1);
        assert_eq!(first.2, second.2);
        assert_eq!(
            first
                .0
                .iter()
                .map(|entry| entry.span_count)
                .collect::<Vec<_>>(),
            vec![1024, 1024, 1024, 1024, 904, 3]
        );
        assert!(first.0.windows(2).all(|entries| {
            entries[0].reference_id <= entries[1].reference_id
                && entries[0].first_span_id + u64::from(entries[0].span_count)
                    == entries[1].first_span_id
        }));
        assert_eq!(first.0[4].reference_id, 0);
        assert_eq!(first.0[5].reference_id, 1);
        assert!(encode_span_blocks_with_threads(&spans, 0, 1).is_err());
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_spill_budgets_and_compression_threads_are_byte_deterministic() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let options = NameIndexOptions::new(["Name", "Alias", "gene_name"], false).unwrap();
        let first_destination = directory.path().join("deterministic-first.gai");
        let second_destination = directory.path().join("deterministic-second.gai");
        let first_options = BuildOptions::default()
            .with_memory_budget(1)
            .with_compression_threads(1)
            .with_bgzf_threads(1);
        let second_options = BuildOptions::default()
            .with_memory_budget(1024)
            .with_compression_threads(4)
            .with_bgzf_threads(2);
        let first_stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &first_destination,
            &options,
            &first_options,
        )
        .unwrap();
        let second_stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &second_destination,
            &options,
            &second_options,
        )
        .unwrap();
        assert_eq!(
            fs::read(&first_destination).unwrap(),
            fs::read(&second_destination).unwrap()
        );
        assert_eq!(first_stats.distinct_terms, second_stats.distinct_terms);
        assert_eq!(first_stats.unique_spans, second_stats.unique_spans);
        assert!(first_stats.peak_working_set_bytes >= first_options.memory_budget_bytes as u64);
        let first_phase_total = first_stats.timings.scan
            + first_stats.timings.spill
            + first_stats.timings.merge
            + first_stats.timings.encode_postings
            + first_stats.timings.encode_spans
            + first_stats.timings.serialize;
        let second_phase_total = second_stats.timings.scan
            + second_stats.timings.spill
            + second_stats.timings.merge
            + second_stats.timings.encode_postings
            + second_stats.timings.encode_spans
            + second_stats.timings.serialize;
        assert!(first_phase_total <= first_stats.timings.total);
        assert!(second_phase_total <= second_stats.timings.total);
        assert!(
            fs::read_dir(directory.path())
                .unwrap()
                .filter_map(|entry| entry.ok())
                .all(|entry| !entry.file_name().to_string_lossy().contains("-run-"))
        );
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_streaming_progress_and_fasta_fingerprint() {
        let directory = tempdir().expect("should create temp directory");
        let source = directory.path().join("fasta.gff3.gz");
        let coordinate_index = directory.path().join("fasta.gff3.gz.tbi");
        let destination = directory.path().join("fasta.gai");
        let record = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=fasta-term";

        let write_source = |sequence: &str| {
            let mut writer = File::create(&source)
                .map(bgzf::io::Writer::new)
                .expect("should create FASTA source");
            let mut indexer = tabix::index::Indexer::default();
            indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
            writeln!(writer, "##gff-version 3").unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{record}").unwrap();
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    "chr1",
                    Position::try_from(10).unwrap(),
                    Position::try_from(20).unwrap(),
                    Chunk::new(start_position, end_position),
                )
                .unwrap();
            writeln!(writer, "##FASTA").unwrap();
            writeln!(writer, ">chr1").unwrap();
            writeln!(writer, "{sequence}").unwrap();
            writer.finish().unwrap();
            let index = indexer.build();
            let mut index_writer = File::create(&coordinate_index)
                .map(tabix::io::Writer::new)
                .unwrap();
            index_writer.write_index(&index).unwrap();
        };
        write_source("ACGTACGT");

        let progress = Arc::new(Mutex::new(Vec::<BuildProgress>::new()));
        let progress_for_callback = Arc::clone(&progress);
        let build_options = BuildOptions::default()
            .with_memory_budget(1)
            .with_bgzf_threads(2)
            .with_progress(move |event| {
                progress_for_callback.lock().unwrap().push(event);
            });
        let stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
            &build_options,
        )
        .unwrap();
        assert_eq!(stats.records_processed, 1);
        assert_eq!(
            NameIndexReader::open(&destination)
                .unwrap()
                .lookup_span_ids("fasta-term")
                .unwrap()
                .len(),
            1
        );
        let events = progress.lock().unwrap().clone();
        assert!(events.iter().any(|event| event.phase == BuildPhase::Scan));
        assert!(
            events
                .iter()
                .any(|event| event.phase == BuildPhase::Complete)
        );
        assert!(events.windows(2).all(|window| {
            window[0].records_processed <= window[1].records_processed
                && window[0].bytes_read <= window[1].bytes_read
                && window[0].elapsed <= window[1].elapsed
        }));
        let source_bytes = fs::metadata(&source).unwrap().len();
        assert_eq!(events.last().unwrap().bytes_read, source_bytes);

        write_source("TTTTCCCC");
        assert!(matches!(
            IndexedSource::open(&source, &coordinate_index, &destination),
            Err(Error::Stale(message)) if message.contains("source fingerprint")
        ));
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_scan_progress_is_rate_limited() {
        let directory = tempdir().expect("should create temp directory");
        let mut storage = Vec::with_capacity(20_001);
        storage.push("##gff-version 3".to_string());
        for index in 0..20_000 {
            let position = index + 1;
            storage.push(format!(
                "chr1\tsrc\tgene\t{position}\t{position}\t.\t+\t.\tName=progress-{index}"
            ));
        }
        let lines = storage.iter().map(String::as_str).collect::<Vec<_>>();
        let (source, coordinate_index) = write_tbi_lines(directory.path(), "progress", &lines);
        let destination = directory.path().join("progress.gai");
        let progress = Arc::new(Mutex::new(Vec::<BuildProgress>::new()));
        let progress_for_callback = Arc::clone(&progress);
        let options = BuildOptions::default()
            .with_memory_budget(1 << 30)
            .with_bgzf_threads(1)
            .with_progress(move |event| {
                progress_for_callback.lock().unwrap().push(event);
            });
        let stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
            &options,
        )
        .unwrap();
        assert_eq!(stats.records_processed, 20_000);
        let events = progress.lock().unwrap().clone();
        assert_eq!(
            events
                .iter()
                .filter(|event| event.phase == BuildPhase::Scan)
                .count(),
            1,
            "the final scan event should be the only scan callback below the threshold"
        );
        assert_eq!(
            events
                .iter()
                .filter(|event| event.phase == BuildPhase::Complete)
                .count(),
            1
        );
        assert!(
            events.len() <= 8,
            "unexpected progress event count: {}",
            events.len()
        );
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_spill_run_fan_in_compaction_is_equivalent() {
        let directory = tempdir().expect("should create temp directory");
        let mut storage = Vec::new();
        for index in 0..200 {
            storage.push(format!(
                "chr1\tsrc\tgene\t{}\t{}\t.\t+\t.\tName=run-{index:03}",
                index + 1,
                index + 1
            ));
        }
        let mut lines = vec!["##gff-version 3"];
        lines.extend(storage.iter().map(String::as_str));
        let (source, coordinate_index) = write_tbi_lines(directory.path(), "fan-in", &lines);
        let options = NameIndexOptions::new(["Name"], false).unwrap();
        let compacted = directory.path().join("compacted.gai");
        let baseline = directory.path().join("fan-in-baseline.gai");
        build_name_index_with_options(
            &source,
            &coordinate_index,
            &compacted,
            &options,
            &BuildOptions::default()
                .with_memory_budget(1)
                .with_compression_threads(1)
                .with_bgzf_threads(1),
        )
        .unwrap();
        build_name_index(&source, &coordinate_index, &baseline, &options).unwrap();
        assert_eq!(fs::read(&compacted).unwrap(), fs::read(&baseline).unwrap());
        assert!(
            fs::read_dir(directory.path())
                .unwrap()
                .filter_map(|entry| entry.ok())
                .all(|entry| !entry.file_name().to_string_lossy().contains("-run-"))
        );
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_valid_empty_index_has_no_terms_or_spans() {
        let directory = tempdir().unwrap();
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("empty.gai");
        let options = NameIndexOptions::new(["not_present"], false).unwrap();
        let stats = build_name_index(&source, &coordinate_index, &destination, &options).unwrap();
        assert_eq!(stats.distinct_terms, 0);
        assert_eq!(stats.unique_spans, 0);
        let reader = NameIndexReader::open(&destination).unwrap();
        assert_eq!(reader.metadata().term_count, 0);
        assert!(reader.lookup_span_ids("anything").unwrap().is_empty());
        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination).unwrap();
        let (records, stats) = indexed.query_name_with_stats("anything").unwrap();
        assert!(records.is_empty());
        assert_eq!(stats, QueryStats::default());
    }
}
