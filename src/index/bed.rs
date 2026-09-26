use super::ExtractedRecord;
use crate::*;

pub(crate) fn extract_bed_record(
    raw: &[u8],
    line_number: usize,
    reference_ids: &HashMap<String, u32>,
    case_sensitive: bool,
) -> Result<ExtractedRecord> {
    if raw.is_empty() || raw.iter().all(u8::is_ascii_whitespace) {
        return Err(Error::InvalidInput(format!(
            "invalid BED record at line {line_number}: expected at least 3 tab-separated columns"
        )));
    }

    let fields = raw.split(|byte| *byte == b'\t').collect::<Vec<_>>();
    if fields.len() < 3 {
        return Err(Error::InvalidInput(format!(
            "invalid BED record at line {line_number}: expected at least 3 tab-separated columns, found {}",
            fields.len()
        )));
    }

    let reference_sequence_name = std::str::from_utf8(fields[0]).map_err(|_| {
        Error::InvalidInput(format!(
            "invalid BED record at line {line_number}: reference sequence name is not UTF-8"
        ))
    })?;
    if reference_sequence_name.is_empty() {
        return Err(Error::InvalidInput(format!(
            "invalid BED record at line {line_number}: reference sequence name must not be empty"
        )));
    }
    let reference_id = *reference_ids.get(reference_sequence_name).ok_or_else(|| {
        Error::InvalidInput(format!(
            "BED reference sequence {reference_sequence_name:?} is absent from coordinate index"
        ))
    })?;

    let parse_coordinate = |raw: &[u8], name: &str| -> Result<u64> {
        let value = std::str::from_utf8(raw).map_err(|_| {
            Error::InvalidInput(format!(
                "invalid BED record at line {line_number}: {name} coordinate is not UTF-8"
            ))
        })?;
        value.parse::<u64>().map_err(|error| {
            Error::InvalidInput(format!(
                "invalid BED record at line {line_number}: invalid {name} coordinate {value:?}: {error}"
            ))
        })
    };

    let start = parse_coordinate(fields[1], "start")?;
    let end = parse_coordinate(fields[2], "end")?;
    let length = end.checked_sub(start).ok_or_else(|| {
        Error::InvalidInput(format!(
            "invalid BED record at line {line_number}: coordinates must satisfy start < end (got {start}..{end})"
        ))
    })?;
    if length == 0 {
        return Err(Error::InvalidInput(format!(
            "invalid BED record at line {line_number}: coordinates must satisfy start < end (got {start}..{end})"
        )));
    }
    start.checked_add(length).ok_or(Error::InvalidCoordinate)?;

    let terms = fields
        .get(3)
        .map(|name| {
            std::str::from_utf8(name)
                .map_err(|_| {
                    Error::InvalidInput(format!(
                        "invalid BED record at line {line_number}: name is not UTF-8"
                    ))
                })
                .map(|name| normalize_value(name, case_sensitive))
        })
        .transpose()?
        .filter(|name| !name.is_empty())
        .into_iter()
        .collect();

    Ok(ExtractedRecord {
        span: SpanKey {
            reference_id,
            start,
            length,
        },
        terms,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_bed_name_and_zero_based_span() {
        let reference_ids = HashMap::from([(String::from("chr1"), 7)]);
        let extracted = extract_bed_record(b"chr1\t10\t25\t Alpha ", 1, &reference_ids, false)
            .expect("should extract BED name and span");
        assert_eq!(
            extracted.span,
            SpanKey {
                reference_id: 7,
                start: 10,
                length: 15,
            }
        );
        assert_eq!(extracted.terms, vec!["alpha"]);

        let unnamed = extract_bed_record(b"chr1\t25\t30", 2, &reference_ids, false)
            .expect("should extract BED span without a name");
        assert!(unnamed.terms.is_empty());
    }
}
