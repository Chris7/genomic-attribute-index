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
