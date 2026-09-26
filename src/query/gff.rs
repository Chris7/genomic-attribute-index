use super::QueryReadContext;
use crate::*;
#[tracing::instrument(level = "trace", skip_all)]
fn parsed_record_from_record_buf(
    record: &gff::feature::RecordBuf,
    raw_line: String,
    reference_ids: &HashMap<String, u32>,
) -> Result<ParsedRecord> {
    let reference_sequence_name = String::from_utf8(record.reference_sequence_name().to_vec())
        .map_err(|_| Error::InvalidInput("reference sequence name is not UTF-8".into()))?;
    if !reference_ids.contains_key(&reference_sequence_name) {
        return Err(Error::InvalidInput(format!(
            "GFF reference sequence {reference_sequence_name:?} is absent from coordinate index"
        )));
    }
    let start = record.start().get() as u64;
    let end = record.end().get() as u64;
    gff_to_span(start, end)?;
    let source = String::from_utf8(record.source().to_vec())
        .map_err(|_| Error::InvalidInput("GFF source is not UTF-8".into()))?;
    let ty = String::from_utf8(record.ty().to_vec())
        .map_err(|_| Error::InvalidInput("GFF type is not UTF-8".into()))?;
    let attributes = record
        .attributes()
        .as_ref()
        .iter()
        .map(|(tag, value)| {
            let tag = String::from_utf8(<_ as AsRef<[u8]>>::as_ref(tag).to_vec())
                .map_err(|_| Error::InvalidInput("GFF attribute tag is not UTF-8".into()))?;
            let values = value
                .iter()
                .map(|value| {
                    String::from_utf8(<_ as AsRef<[u8]>>::as_ref(value).to_vec())
                        .map_err(|_| Error::InvalidInput("GFF attribute value is not UTF-8".into()))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok((tag, values))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(ParsedRecord {
        record: GffRecord {
            reference_sequence_name,
            source,
            ty,
            start,
            end,
            score: record
                .score()
                .map(|score| score.to_string())
                .unwrap_or_else(|| ".".to_string()),
            strand: match record.strand() {
                gff::feature::record::Strand::None => ".".to_string(),
                gff::feature::record::Strand::Forward => "+".to_string(),
                gff::feature::record::Strand::Reverse => "-".to_string(),
                gff::feature::record::Strand::Unknown => "?".to_string(),
            },
            phase: record
                .phase()
                .map(|phase| match phase {
                    gff::feature::record::Phase::Zero => "0".to_string(),
                    gff::feature::record::Phase::One => "1".to_string(),
                    gff::feature::record::Phase::Two => "2".to_string(),
                })
                .unwrap_or_else(|| ".".to_string()),
            attributes,
            raw_line,
        },
    })
}
#[tracing::instrument(level = "trace", skip_all)]
fn feature_record_span<R>(record: &R, reference_ids: &HashMap<String, u32>) -> Result<SpanKey>
where
    R: gff::feature::Record + ?Sized,
{
    let reference_sequence_name = String::from_utf8(record.reference_sequence_name().to_vec())
        .map_err(|_| Error::InvalidInput("GFF reference sequence name is not UTF-8".into()))?;
    let reference_id = *reference_ids.get(&reference_sequence_name).ok_or_else(|| {
        Error::InvalidInput(format!(
            "GFF reference sequence {reference_sequence_name:?} is absent from coordinate index"
        ))
    })?;
    let start = record
        .feature_start()
        .map_err(|error| Error::InvalidInput(format!("invalid GFF start: {error}")))?
        .get() as u64;
    let end = record
        .feature_end()
        .map_err(|error| Error::InvalidInput(format!("invalid GFF end: {error}")))?
        .get() as u64;
    let (start, length) = gff_to_span(start, end)?;
    Ok(SpanKey {
        reference_id,
        start,
        length,
    })
}
#[tracing::instrument(level = "trace", skip_all)]
fn feature_record_matches_term<R>(
    record: &R,
    configured_attributes: &HashSet<String>,
    case_sensitive: bool,
    matcher: &CompiledMatch,
) -> Result<bool>
where
    R: gff::feature::Record + ?Sized,
{
    for result in record.attributes().iter() {
        let (tag, value) = result
            .map_err(|error| Error::InvalidInput(format!("invalid GFF attributes: {error}")))?;
        let tag = std::str::from_utf8(tag.as_ref())
            .map_err(|_| Error::InvalidInput("GFF attribute tag is not UTF-8".into()))?;
        if !configured_attributes.contains(tag) {
            continue;
        }
        for value in value.iter() {
            let value = value.map_err(|error| {
                Error::InvalidInput(format!("invalid GFF attribute value: {error}"))
            })?;
            let value = std::str::from_utf8(value.as_ref())
                .map_err(|_| Error::InvalidInput("GFF attribute value is not UTF-8".into()))?;
            let normalized_value = normalize_value(value, case_sensitive);
            if matcher.matches_normalized(&normalized_value) {
                return Ok(true);
            }
        }
    }
    Ok(false)
}
#[tracing::instrument(level = "trace", skip_all)]
fn parsed_record_from_feature_record<R>(
    record: &R,
    raw_line: String,
    reference_ids: &HashMap<String, u32>,
) -> Result<ParsedRecord>
where
    R: gff::feature::Record + ?Sized,
{
    let record = gff::feature::RecordBuf::try_from_feature_record(record)
        .map_err(|error| Error::InvalidInput(format!("invalid GFF record: {error}")))?;
    parsed_record_from_record_buf(&record, raw_line, reference_ids)
}
#[tracing::instrument(level = "trace", skip_all)]
pub(super) fn read_gff_query_chunks(
    source: File,
    chunks: &[Chunk],
    context: &QueryReadContext<'_>,
    stats: &mut QueryStats,
) -> Result<Vec<GffRecord>> {
    let mut reader = gff::io::Reader::new(bgzf::io::Reader::new(source));
    let mut line = gff::Line::default();
    let mut positions = HashSet::new();
    let mut records = Vec::new();
    for chunk in chunks {
        reader.get_mut().seek_to_virtual_position(chunk.start())?;
        loop {
            let source_position = reader.get_ref().virtual_position();
            if source_position >= chunk.end() {
                break;
            }
            let length = reader.read_line(&mut line)?;
            if length == 0 {
                break;
            }
            stats.bytes_read = stats
                .bytes_read
                .checked_add(length as u64)
                .ok_or(Error::InvalidCoordinate)?;
            let Some(result) = line.as_record() else {
                continue;
            };
            if !positions.insert(u64::from(source_position)) {
                continue;
            }
            stats.unique_candidate_records += 1;
            let record = result
                .map_err(|error| Error::InvalidInput(format!("invalid GFF record: {error}")))?;
            let span = feature_record_span(&record, context.reference_ids)?;
            if !context.requested_spans.contains(&span) {
                continue;
            }
            if !feature_record_matches_term(
                &record,
                context.configured_terms,
                context.case_sensitive,
                context.matcher,
            )? {
                continue;
            }
            let raw_line = String::from_utf8(line.as_ref().to_vec())
                .map_err(|_| Error::InvalidInput("GFF record is not UTF-8".into()))?;
            let parsed =
                parsed_record_from_feature_record(&record, raw_line, context.reference_ids)?;
            stats.matching_records += 1;
            records.push(parsed.record);
        }
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use super::*;
    use crate::{
        index::tests::{fixture_records, write_csi_fixture, write_fixture, write_tbi_lines},
        query::{merge_query_chunks, read_query_chunks},
    };

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_name_index_round_trip_and_query() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("fixture.gai");
        let options = NameIndexOptions::new(["Name", "Alias", "gene_name", "Name"], false)
            .expect("should validate attributes");
        let stats = build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build GAI");
        let second_destination = directory.path().join("fixture-second.gai");
        build_name_index(&source, &coordinate_index, &second_destination, &options)
            .expect("should rebuild GAI deterministically");
        assert_eq!(
            fs::read(&destination).unwrap(),
            fs::read(&second_destination).unwrap()
        );
        assert_eq!(stats.records_processed, 4);
        assert_eq!(stats.records_indexed, 4);
        assert_eq!(stats.distinct_terms, 4);
        assert_eq!(stats.unique_spans, 3);
        let reader = NameIndexReader::open(&destination).expect("should open GAI");
        assert_eq!(reader.metadata().minor_version, 0);
        let mapped = NameIndexReader::open_mmap(&destination).expect("should mmap GAI");
        assert_eq!(mapped.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
        assert_eq!(reader.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
        assert_eq!(reader.lookup_span_ids("brcc1").unwrap(), vec![0]);
        assert_eq!(reader.lookup_span_ids("rnf53").unwrap(), vec![0]);
        assert_eq!(reader.lookup_span_ids("other gene").unwrap(), vec![1]);
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("br", MatchMode::Prefix)
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("brcc", MatchMode::Prefix)
                .unwrap(),
            vec![0]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode(" RCA ", MatchMode::Contains)
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("^brc(a1|c1)$", MatchMode::Regex)
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode(r"^brc[a-z][0-9]$", MatchMode::Regex)
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode(r"\S", MatchMode::Regex)
                .unwrap(),
            vec![0, 1, 2]
        );
        assert!(
            reader
                .lookup_span_ids_with_mode("", MatchMode::Contains)
                .unwrap()
                .is_empty()
        );
        assert!(matches!(
            reader.lookup_span_ids_with_mode("missing[", MatchMode::Regex),
            Err(Error::InvalidInput(message)) if message.contains("regex")
        ));
        assert_eq!(reader.resolve_span_id(0).unwrap().start, 99);
        assert_eq!(reader.resolve_span_id(0).unwrap().length, 51);
        assert!(
            reader
                .lookup_span_ids("without_special_status")
                .unwrap()
                .is_empty()
        );

        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open indexed source");
        let records = indexed.query_name(" BRCA1 ").expect("should query name");
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].start, 100);
        assert_eq!(records[2].reference_sequence_name, "chr2");
        let (prefix_records, prefix_stats) = indexed
            .query_name_with_mode_and_stats(" br ", MatchMode::Prefix)
            .expect("should query prefix");
        assert_eq!(prefix_records.len(), 3);
        assert_eq!(prefix_stats.requested_spans, 2);
        assert_eq!(prefix_stats.exact_interval_queries, 2);
        assert_eq!(prefix_stats.matching_records, 3);
        assert_eq!(
            prefix_records
                .iter()
                .map(|record| record.raw_line.as_str())
                .collect::<Vec<_>>(),
            vec![
                fixture_records()[2],
                fixture_records()[3],
                fixture_records()[5]
            ]
        );
        let (contains_records, contains_stats) = indexed
            .query_name_with_mode_and_stats(" RCA ", MatchMode::Contains)
            .expect("should query literal substring");
        assert_eq!(contains_records.len(), 3);
        assert_eq!(contains_stats.requested_spans, 2);
        assert_eq!(contains_stats.exact_interval_queries, 2);
        assert_eq!(contains_stats.matching_records, 3);
        let (regex_records, regex_stats) = indexed
            .query_name_with_mode_and_stats("^BR C[A-Z]$", MatchMode::Regex)
            .expect("should query regex");
        assert!(regex_records.is_empty());
        assert_eq!(regex_stats.requested_spans, 0);
        let (regex_records, regex_stats) = indexed
            .query_name_with_mode_and_stats("^BRC(A1|C1)$", MatchMode::Regex)
            .expect("should query anchored regex alternation");
        assert_eq!(regex_records.len(), 3);
        assert_eq!(regex_stats.requested_spans, 2);
        assert_eq!(regex_stats.matching_records, 3);
        let (escape_records, _) = indexed
            .query_name_with_mode_and_stats(r"\S", MatchMode::Regex)
            .expect("should preserve uppercase regex escape");
        assert_eq!(escape_records.len(), 4);
        assert!(indexed.query_name("missing").unwrap().is_empty());
        assert!(
            indexed
                .query_name_with_mode("missing", MatchMode::Prefix)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_configured_attributes_and_disjoint_query_spans() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("id-only.gai");
        let options = NameIndexOptions::new(["ID"], false).unwrap();
        build_name_index(&source, &coordinate_index, &destination, &options).unwrap();
        let reader = NameIndexReader::open(&destination).unwrap();
        assert!(reader.lookup_span_ids("brca1").unwrap().is_empty());
        assert_eq!(
            reader.lookup_span_ids("without_special_status").unwrap(),
            vec![0]
        );
        let duplicate_options = NameIndexOptions::new(["Name", "Name", "Alias"], false).unwrap();
        assert_eq!(duplicate_options.attributes, vec!["Name", "Alias"]);
        assert!(NameIndexOptions::new([""], false).is_err());
        assert!(NameIndexOptions::new(std::iter::empty::<&str>(), false).is_err());

        // The two records below share a term but not a span.  This exercises
        // exact disjoint postings instead of a bounding interval.
        let source_path = directory.path().join("disjoint.gff3.gz");
        let index_path = directory.path().join("disjoint.gff3.gz.tbi");
        let lines = [
            "##gff-version 3",
            "chr1\tsrc\tgene\t10\t10\t.\t+\t.\tName=shared",
            "chr1\tsrc\tgene\t100\t100\t.\t+\t.\tName=shared",
        ];
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .unwrap();
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
        for line in lines {
            if line.starts_with('#') {
                writeln!(writer, "{line}").unwrap();
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let position = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").unwrap();
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    fields[0],
                    position,
                    position,
                    Chunk::new(start_position, end_position),
                )
                .unwrap();
        }
        writer.finish().unwrap();
        let index = indexer.build();
        let mut index_writer = File::create(&index_path)
            .map(tabix::io::Writer::new)
            .unwrap();
        index_writer.write_index(&index).unwrap();
        drop(index_writer);
        let destination = directory.path().join("disjoint.gai");
        build_name_index(
            &source_path,
            &index_path,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .unwrap();
        let reader = NameIndexReader::open(&destination).unwrap();
        assert_eq!(reader.lookup_span_ids("shared").unwrap(), vec![0, 1]);
        let mut indexed = IndexedSource::open(&source_path, &index_path, &destination).unwrap();
        let (records, stats) = indexed.query_name_with_stats("shared").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(stats.requested_spans, 2);
        assert_eq!(stats.exact_interval_queries, 2);
        assert_eq!(stats.matching_records, 2);
        assert_eq!(stats.unique_candidate_records, 2);
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_query_exact_spans_preserves_identical_records_and_overlap_order() {
        let directory = tempdir().expect("should create temporary directory");
        let first = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=overlap";
        let overlapping = "chr1\tsrc\tgene\t15\t25\t.\t+\t.\tName=overlap";
        let (source, coordinate_index) = write_tbi_lines(
            directory.path(),
            "overlap",
            &["##gff-version 3", first, first, overlapping],
        );
        let destination = directory.path().join("overlap.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build overlap GAI");

        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open overlap GAI");
        let (records, stats) = indexed
            .query_name_with_stats("overlap")
            .expect("should query overlap term");
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].raw_line, first);
        assert_eq!(records[1].raw_line, first);
        assert_eq!(records[2].raw_line, overlapping);
        assert_eq!(records[0].start, 10);
        assert_eq!(records[1].start, 10);
        assert_eq!(records[2].start, 15);
        assert_eq!(stats.requested_spans, 2);
        assert_eq!(stats.exact_interval_queries, 2);
        assert_eq!(stats.distinct_span_blocks_decoded, 1);
        assert_eq!(stats.unique_candidate_records, 3);
        assert_eq!(stats.matching_records, 3);
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_contains_filters_nonmatching_records_at_a_shared_span() {
        let directory = tempdir().expect("should create temporary directory");
        let matching = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=Alpha";
        let nonmatching = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=Beta";
        let (source, coordinate_index) = write_tbi_lines(
            directory.path(),
            "contains-shared-span",
            &["##gff-version 3", matching, nonmatching],
        );
        let destination = directory.path().join("contains-shared-span.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build shared-span GAI");

        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open shared-span GAI");
        let (records, stats) = indexed
            .query_name_with_mode_and_stats("ph", MatchMode::Contains)
            .expect("should query interior substring");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].raw_line, matching);
        assert_eq!(stats.requested_spans, 1);
        assert_eq!(stats.distinct_span_blocks_decoded, 1);
        assert_eq!(stats.exact_interval_queries, 1);
        assert_eq!(stats.unique_candidate_records, 2);
        assert_eq!(stats.matching_records, 1);
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_query_chunk_union_and_virtual_position_deduplication() {
        let merged = merge_query_chunks(vec![
            Chunk::new(
                bgzf::VirtualPosition::from(10),
                bgzf::VirtualPosition::from(20),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(15),
                bgzf::VirtualPosition::from(25),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(30),
                bgzf::VirtualPosition::from(35),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(31),
                bgzf::VirtualPosition::from(32),
            ),
        ]);
        assert_eq!(
            merged,
            vec![
                Chunk::new(
                    bgzf::VirtualPosition::from(10),
                    bgzf::VirtualPosition::from(25)
                ),
                Chunk::new(
                    bgzf::VirtualPosition::from(30),
                    bgzf::VirtualPosition::from(35)
                ),
            ]
        );

        let directory = tempdir().expect("should create temp directory");
        let first = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=dedup";
        let (source, coordinate_index) = write_tbi_lines(
            directory.path(),
            "virtual-position",
            &["##gff-version 3", first, first],
        );
        let (index, _) = read_coordinate_index_with_fingerprint(&coordinate_index).unwrap();
        let spans = [SpanKey {
            reference_id: 0,
            start: 9,
            length: 11,
        }];
        let reference_ids = HashMap::from([(String::from("chr1"), 0_u32)]);
        let configured_attributes = HashSet::from([String::from("Name")]);
        let reference_name = "chr1";
        let query_start = Position::try_from(10).unwrap();
        let query_end = Position::try_from(20).unwrap();
        let region = Region::new(reference_name, query_start..=query_end);
        let chunks = match &index {
            CoordinateIndex::Tabix(index) => index.query(0, region.interval()).unwrap(),
            CoordinateIndex::Csi(index) => index.query(0, region.interval()).unwrap(),
        };
        assert!(!chunks.is_empty());
        let duplicated_chunks = [chunks[0], chunks[0]];
        let mut stats = QueryStats::default();
        let requested_spans = spans.into_iter().collect::<HashSet<_>>();
        let matcher = CompiledMatch::new("dedup", MatchMode::Exact, false).unwrap();
        let context = QueryReadContext {
            requested_spans: &requested_spans,
            reference_ids: &reference_ids,
            configured_terms: &configured_attributes,
            case_sensitive: false,
            matcher: &matcher,
        };
        let records = read_query_chunks(
            &source,
            SortFormat::Gff,
            &duplicated_chunks,
            &context,
            &mut stats,
        )
        .unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(stats.unique_candidate_records, 2);
        assert_eq!(stats.matching_records, 2);
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_query_same_start_preserves_source_length_order() {
        let directory = tempdir().expect("should create temporary directory");
        let longer = "chr1\tsrc\tgene\t10\t30\t.\t+\t.\tName=same-start";
        let shorter = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=same-start";
        let (source, coordinate_index) = write_tbi_lines(
            directory.path(),
            "same-start",
            &["##gff-version 3", longer, shorter],
        );
        let destination = directory.path().join("same-start.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build same-start GAI");
        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open same-start GAI");
        let records = indexed.query_name("same-start").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].raw_line, longer);
        assert_eq!(records[1].raw_line, shorter);
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_feature_crossing_bgzf_block_boundary_is_queryable() {
        let directory = tempdir().expect("should create temporary directory");
        let source = directory.path().join("large.gff3.gz");
        let coordinate_index = directory.path().join("large.gff3.gz.tbi");
        let large_note = "x".repeat(70_000);
        let line = format!("chr1\tsrc\tgene\t100\t200\t.\t+\t.\tName=large;Note={large_note}");
        let mut writer = File::create(&source)
            .map(bgzf::io::Writer::new)
            .expect("should create large BGZF source");
        writeln!(writer, "##gff-version 3").unwrap();
        let start_position = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        let end_position = writer.virtual_position();
        assert!(
            end_position.compressed() > start_position.compressed(),
            "large feature should cross BGZF blocks"
        );
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
        indexer
            .add_record(
                "chr1",
                Position::try_from(100).unwrap(),
                Position::try_from(200).unwrap(),
                Chunk::new(start_position, end_position),
            )
            .unwrap();
        writer.finish().unwrap();
        let index = indexer.build();
        let mut index_writer = File::create(&coordinate_index)
            .map(tabix::io::Writer::new)
            .unwrap();
        index_writer.write_index(&index).unwrap();
        drop(index_writer);

        let destination = directory.path().join("large.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build large-feature GAI");
        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open large-feature GAI");
        let records = indexed.query_name("large").unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].start, 100);
        assert_eq!(records[0].end, 200);
        assert_eq!(
            records[0].attribute_values("Note").next(),
            Some(large_note.as_str())
        );
    }

    #[test]
    #[tracing::instrument(level = "trace", skip_all)]
    fn test_csi_build_query_and_index_fingerprints() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_csi_fixture(directory.path());
        let destination = directory.path().join("fixture-csi.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name", "Alias"], false).unwrap(),
        )
        .unwrap();
        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination).unwrap();
        let (records, stats) = indexed.query_name_with_stats("brca1").unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[2].reference_sequence_name, "chr2");
        assert_eq!(stats.requested_spans, 2);
        assert_eq!(stats.exact_interval_queries, 2);
        assert_eq!(stats.matching_records, 3);
        let metadata = indexed.metadata();
        assert_eq!(
            metadata.coordinate_index_fingerprint,
            fingerprint_file(&coordinate_index).unwrap()
        );
        assert!(metadata.reference_dictionary_fingerprint != [0; 32]);
    }
}
