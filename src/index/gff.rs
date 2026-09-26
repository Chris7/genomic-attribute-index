use super::ExtractedRecord;
use crate::*;

fn extract_gff_record(
    record: &gff::feature::RecordBuf,
    reference_ids: &HashMap<String, u32>,
    configured: &HashSet<String>,
    case_sensitive: bool,
) -> Result<ExtractedRecord> {
    let reference_sequence_name = String::from_utf8(record.reference_sequence_name().to_vec())
        .map_err(|_| Error::InvalidInput("reference sequence name is not UTF-8".into()))?;
    let reference_id = *reference_ids.get(&reference_sequence_name).ok_or_else(|| {
        Error::InvalidInput(format!(
            "GFF reference sequence {reference_sequence_name:?} is absent from coordinate index"
        ))
    })?;
    let (start, length) = gff_to_span(record.start().get() as u64, record.end().get() as u64)?;
    let mut terms = Vec::new();
    for (tag, value) in record.attributes().as_ref() {
        let tag = std::str::from_utf8(tag.as_ref())
            .map_err(|_| Error::InvalidInput("GFF attribute tag is not UTF-8".into()))?;
        if !configured.contains(tag) {
            continue;
        }
        for value in value.iter() {
            let value = String::from_utf8(<_ as AsRef<[u8]>>::as_ref(value).to_vec())
                .map_err(|_| Error::InvalidInput("GFF attribute value is not UTF-8".into()))?;
            let normalized = normalize_value(&value, case_sensitive);
            if !normalized.is_empty() {
                terms.push(normalized);
            }
        }
    }
    Ok(ExtractedRecord {
        span: SpanKey {
            reference_id,
            start,
            length,
        },
        terms,
    })
}

pub(super) fn extract_gff_line(
    raw: &[u8],
    line_number: usize,
    reference_ids: &HashMap<String, u32>,
    configured: &HashSet<String>,
    case_sensitive: bool,
) -> Result<ExtractedRecord> {
    let mut parser = gff::io::Reader::new(Cursor::new(raw));
    let mut line = gff::Line::default();
    parser.read_line(&mut line).map_err(|error| {
        Error::InvalidInput(format!("invalid GFF record at line {line_number}: {error}"))
    })?;
    let record = line
        .as_record()
        .ok_or_else(|| {
            Error::InvalidInput(format!(
                "invalid GFF record at line {line_number}: expected a feature record"
            ))
        })?
        .map_err(|error| {
            Error::InvalidInput(format!("invalid GFF record at line {line_number}: {error}"))
        })?;
    let record = gff::feature::RecordBuf::try_from_feature_record(&record).map_err(|error| {
        Error::InvalidInput(format!("invalid GFF record at line {line_number}: {error}"))
    })?;
    extract_gff_record(&record, reference_ids, configured, case_sensitive)
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use super::*;
    use crate::index::tests::{fixture_records, write_fixture};

    #[test]
    fn test_coordinate_and_normalization_boundaries() {
        assert_eq!(gff_to_span(100, 150).unwrap(), (99, 51));
        assert!(gff_to_span(0, 1).is_err());
        assert!(gff_to_span(2, 1).is_err());
        assert_eq!(normalize_value("  BRCA1\u{2003}", false), "brca1");
        assert_eq!(normalize_value("  BRCA1\u{2003}", true), "BRCA1");
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("case-sensitive.gai");
        let options = NameIndexOptions::new(["Name"], true).expect("should validate attributes");
        build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build case-sensitive GAI");
        let reader = NameIndexReader::open(destination).expect("should open GAI");
        assert!(reader.lookup_span_ids("brca1").unwrap().is_empty());
        assert_eq!(reader.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
        assert!(
            reader
                .lookup_span_ids_with_mode("br", MatchMode::Prefix)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("BR", MatchMode::Prefix)
                .unwrap(),
            vec![0, 2]
        );
        assert!(
            reader
                .lookup_span_ids_with_mode("^brca1$", MatchMode::Regex)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("^BRCA1$", MatchMode::Regex)
                .unwrap(),
            vec![0, 2]
        );
    }

    #[test]
    fn test_coordinate_index_path_and_gff_header_validation() {
        let directory = tempdir().expect("should create temporary directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let renamed_index = directory.path().join("coordinates.data");
        fs::copy(&coordinate_index, &renamed_index).expect("should copy coordinate index");
        let destination = directory.path().join("renamed.gai");
        build_name_index(
            &source,
            &renamed_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should parse a valid index with a nonstandard filename");
        let mut indexed = IndexedSource::open(&source, &renamed_index, &destination)
            .expect("should open renamed coordinate index");
        assert_eq!(indexed.query_name("brca1").unwrap().len(), 2);

        let non_gff_source = directory.path().join("bed-index.gff3.gz");
        let non_gff_index = directory.path().join("bed-index.data");
        let line = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=wrong-header";
        let mut writer = File::create(&non_gff_source)
            .map(bgzf::io::Writer::new)
            .expect("should create non-GFF source");
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::bed().build());
        let fields = line.split('\t').collect::<Vec<_>>();
        let start = Position::try_from(10).unwrap();
        let end = Position::try_from(20).unwrap();
        let start_position = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        let end_position = writer.virtual_position();
        indexer
            .add_record(
                fields[0],
                start,
                end,
                Chunk::new(start_position, end_position),
            )
            .unwrap();
        writer.finish().unwrap();
        let index = indexer.build();
        let mut index_writer = File::create(&non_gff_index)
            .map(tabix::io::Writer::new)
            .unwrap();
        index_writer.write_index(&index).unwrap();
        drop(index_writer);
        let non_gff_destination = directory.path().join("non-gff.gai");
        assert!(matches!(
            build_name_index(
                &non_gff_source,
                &non_gff_index,
                &non_gff_destination,
                &NameIndexOptions::new(["Name"], false).unwrap(),
            ),
            Err(Error::InvalidInput(message)) if message.contains("generic GFF")
        ));
    }

    #[test]
    fn test_plain_gff_uses_the_same_streaming_parser() {
        let directory = tempdir().expect("should create temp directory");
        let (_bgzf_source, coordinate_index) = write_fixture(directory.path());
        let source = directory.path().join("plain.gff3");
        let text = fixture_records().join("\n") + "\n";
        fs::write(&source, text).unwrap();
        let destination = directory.path().join("plain.gai");
        let stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name", "Alias"], false).unwrap(),
            &BuildOptions::default().with_bgzf_threads(1),
        )
        .unwrap();
        assert_eq!(stats.records_processed, 4);
        assert_eq!(stats.records_indexed, 4);
        assert_eq!(stats.distinct_terms, 4);
    }
}
