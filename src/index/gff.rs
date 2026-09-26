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
