use super::QueryReadContext;
use crate::*;

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
