use super::QueryReadContext;
use crate::*;

pub(super) fn read_bed_query_chunks(
    source: File,
    chunks: &[Chunk],
    context: &QueryReadContext<'_>,
    stats: &mut QueryStats,
) -> Result<Vec<GffRecord>> {
    let mut reader = bgzf::io::Reader::new(source);
    let mut positions = HashSet::new();
    let mut records = Vec::new();
    let mut line = Vec::new();
    let mut line_number = 0_usize;
    for chunk in chunks {
        reader.seek_to_virtual_position(chunk.start())?;
        loop {
            let source_position = reader.virtual_position();
            if source_position >= chunk.end() {
                break;
            }
            line.clear();
            let length = reader.read_until(b'\n', &mut line)?;
            if length == 0 {
                break;
            }
            line_number += 1;
            stats.bytes_read = stats
                .bytes_read
                .checked_add(length as u64)
                .ok_or(Error::InvalidCoordinate)?;
            let raw = line.strip_suffix(b"\n").unwrap_or(&line);
            let raw = raw.strip_suffix(b"\r").unwrap_or(raw);
            if raw.is_empty()
                || raw.first() == Some(&b'#')
                || raw.starts_with(b"track ")
                || raw.starts_with(b"browser ")
            {
                continue;
            }
            if !positions.insert(u64::from(source_position)) {
                continue;
            }
            stats.unique_candidate_records += 1;
            let extracted = crate::index::bed::extract_bed_record(
                raw,
                line_number,
                context.reference_ids,
                context.case_sensitive,
            )?;
            if !context.requested_spans.contains(&extracted.span) {
                continue;
            }
            let matches = extracted
                .terms
                .iter()
                .any(|value| context.matcher.matches_normalized(value));
            if !matches {
                continue;
            }
            let fields = raw.split(|byte| *byte == b'\t').collect::<Vec<_>>();
            let reference_sequence_name = String::from_utf8(fields[0].to_vec())
                .map_err(|_| Error::InvalidInput("BED record is not UTF-8".into()))?;
            let name = fields
                .get(3)
                .map(|field| String::from_utf8(field.to_vec()))
                .transpose()
                .map_err(|_| Error::InvalidInput("BED name is not UTF-8".into()))?
                .unwrap_or_default();
            let start = extracted.span.start + 1;
            let end = extracted
                .span
                .start
                .checked_add(extracted.span.length)
                .ok_or(Error::InvalidCoordinate)?;
            records.push(GffRecord {
                reference_sequence_name,
                source: ".".to_string(),
                ty: "bed".to_string(),
                start,
                end,
                score: ".".to_string(),
                strand: ".".to_string(),
                phase: ".".to_string(),
                attributes: vec![("name".to_string(), vec![name])],
                raw_line: String::from_utf8(raw.to_vec())
                    .map_err(|_| Error::InvalidInput("BED record is not UTF-8".into()))?,
            });
            stats.matching_records += 1;
        }
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use super::*;
    use crate::index::tests::write_bed_fixture;

    #[test]
    fn test_bed_name_index_builds_from_column_four() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_bed_fixture(directory.path());
        let destination = directory.path().join("fixture-bed.gai");

        // BED ignores configured GFF attributes and always indexes column 4 (`name`).
        let options = NameIndexOptions::new(["DefinitelyIgnored"], false).unwrap();
        let stats = build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build BED GAI");
        assert_eq!(stats.records_processed, 4);
        assert_eq!(stats.records_indexed, 3);
        assert_eq!(stats.distinct_terms, 2);
        assert_eq!(stats.unique_spans, 3);

        let reader = NameIndexReader::open(&destination).expect("should open BED GAI");
        assert_eq!(reader.metadata().attributes, vec!["name"]);
        assert_eq!(reader.lookup_span_ids("alpha").unwrap(), vec![0, 2]);
        assert_eq!(reader.lookup_span_ids("beta").unwrap(), vec![1]);
        assert!(
            reader
                .lookup_span_ids("DefinitelyIgnored")
                .unwrap()
                .is_empty()
        );
        assert_eq!(reader.resolve_span_id(0).unwrap().start, 0);
        assert_eq!(reader.resolve_span_id(0).unwrap().length, 10);
        assert_eq!(reader.resolve_span_id(1).unwrap().start, 20);
        assert_eq!(reader.resolve_span_id(1).unwrap().length, 5);

        let case_destination = directory.path().join("fixture-bed-case.gai");
        let case_options = NameIndexOptions::new(["Ignored"], true).unwrap();
        build_name_index(&source, &coordinate_index, &case_destination, &case_options)
            .expect("should build case-sensitive BED GAI");
        let case_reader = NameIndexReader::open(&case_destination).unwrap();
        assert_eq!(case_reader.lookup_span_ids("Alpha").unwrap(), vec![0, 2]);
        assert!(case_reader.lookup_span_ids("alpha").unwrap().is_empty());

        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open BED indexed source");
        let records = indexed.query_name("alpha").expect("should query BED name");
        assert_eq!(
            records
                .iter()
                .map(|record| record.raw_line.as_str())
                .collect::<Vec<_>>(),
            vec!["chr1\t0\t10\tAlpha", "chr2\t5\t7\tAlpha"]
        );
        assert_eq!(records[0].start, 1);
        assert_eq!(records[0].end, 10);
        assert_eq!(
            records[0].attribute_values("name").collect::<Vec<_>>(),
            vec!["Alpha"]
        );
        let prefix = indexed
            .query_name_with_mode("al", MatchMode::Prefix)
            .expect("should query BED name prefix");
        assert_eq!(prefix.len(), 2);

        let (contains, contains_stats) = indexed
            .query_name_with_mode_and_stats(" ph ", MatchMode::Contains)
            .expect("should query BED name substring");
        assert_eq!(contains.len(), 2);
        assert_eq!(contains_stats.requested_spans, 2);
        assert_eq!(contains_stats.matching_records, 2);

        let (regex, regex_stats) = indexed
            .query_name_with_mode_and_stats("^(ALPHA|BETA)$", MatchMode::Regex)
            .expect("should query BED name regex");
        assert_eq!(regex.len(), 3);
        assert_eq!(regex_stats.matching_records, 3);
        let (class_regex, _) = indexed
            .query_name_with_mode_and_stats("^ALP[A-Z]+$", MatchMode::Regex)
            .expect("should query BED regex character class");
        assert_eq!(class_regex.len(), 2);
        let (escape_regex, _) = indexed
            .query_name_with_mode_and_stats(r"\S", MatchMode::Regex)
            .expect("should preserve BED regex uppercase escape");
        assert_eq!(escape_regex.len(), 3);
        assert!(
            indexed
                .query_name_with_mode("^missing$", MatchMode::Regex)
                .unwrap()
                .is_empty()
        );
        assert!(matches!(
            indexed.query_name_with_mode("missing[", MatchMode::Regex),
            Err(Error::InvalidInput(message)) if message.contains("regex")
        ));

        let mut case_indexed = IndexedSource::open(&source, &coordinate_index, &case_destination)
            .expect("should open case-sensitive BED GAI");
        assert_eq!(
            case_indexed
                .query_name_with_mode("^Alpha$", MatchMode::Regex)
                .unwrap()
                .len(),
            2
        );
        assert!(
            case_indexed
                .query_name_with_mode("^alpha$", MatchMode::Regex)
                .unwrap()
                .is_empty()
        );
        assert!(
            case_indexed
                .query_name_with_mode("PH", MatchMode::Contains)
                .unwrap()
                .is_empty()
        );
    }
}
