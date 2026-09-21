//! Lossless coordinate sorting for GFF3 and BED text files.

use std::{
    cmp::Ordering,
    collections::{BTreeSet, HashMap, HashSet},
    fs::File,
    io::{self, BufRead, BufReader, Cursor, Write},
    path::Path,
};

use noodles::gff;

use crate::{Error, Result};

/// A supported annotation format for [`sort_file`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SortFormat {
    /// GFF or GFF3, including comment and directive lines.
    Gff,
    /// BED with at least the three required coordinate columns.
    Bed,
}

impl SortFormat {
    /// Infers a format from a path's final extension.
    pub fn from_path(path: &Path) -> Result<Self> {
        let extension = path.extension().and_then(|value| value.to_str());
        match extension {
            Some(value) if value.eq_ignore_ascii_case("gff") => Ok(Self::Gff),
            Some(value) if value.eq_ignore_ascii_case("gff3") => Ok(Self::Gff),
            Some(value) if value.eq_ignore_ascii_case("bed") => Ok(Self::Bed),
            Some(value) => Err(Error::InvalidInput(format!(
                "unsupported sort input extension .{value}; expected .gff, .gff3, or .bed"
            ))),
            None => Err(Error::InvalidInput(
                "cannot infer sort input format without a .gff, .gff3, or .bed extension".into(),
            )),
        }
    }
}

/// Sorts an annotation file inferred from its extension and writes records to
/// `writer` in lossless text order.
pub fn sort_file(input: impl AsRef<Path>, writer: impl Write) -> Result<()> {
    let input = input.as_ref();
    let format = SortFormat::from_path(input)?;
    let reader = BufReader::new(File::open(input)?);
    match format {
        SortFormat::Gff => sort_gff(reader, writer),
        SortFormat::Bed => sort_bed(reader, writer),
    }
}

/// Sorts GFF/GFF3 records by contig, start, and end position.
///
/// All lines beginning with `#` are emitted first in their original order.
/// Records sharing all three coordinate keys are topologically ordered by
/// their decoded `ID`/`Parent` attributes, with a stable source-order tie
/// break for otherwise unrelated records.
pub fn sort_gff<R: BufRead, W: Write>(mut reader: R, mut writer: W) -> Result<()> {
    let mut comments = Vec::new();
    let mut records = Vec::new();
    let mut line = Vec::new();
    let mut line_number = 0usize;

    loop {
        line.clear();
        let bytes_read = reader.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            break;
        }
        line_number += 1;
        let raw = strip_line_ending(&line);
        if raw.first() == Some(&b'#') {
            comments.push(raw.to_vec());
        } else {
            records.push(parse_gff_record(raw, line_number)?);
        }
    }

    records.sort_by(compare_gff_coordinates);
    let mut output_order = Vec::with_capacity(records.len());
    let mut offset = 0;
    while offset < records.len() {
        let end = records[offset..]
            .iter()
            .position(|record| compare_gff_coordinates(&records[offset], record) != Ordering::Equal)
            .map_or(records.len(), |relative| offset + relative);
        let order = order_gff_tie_group(&records[offset..end])?;
        for index in order {
            output_order.push(offset + index);
        }
        offset = end;
    }

    for comment in comments {
        write_line(&mut writer, &comment)?;
    }
    for index in output_order {
        write_line(&mut writer, &records[index].raw)?;
    }

    Ok(())
}

/// Sorts BED records by contig, start, and end position.
pub fn sort_bed<R: BufRead, W: Write>(mut reader: R, mut writer: W) -> Result<()> {
    let mut records = Vec::new();
    let mut line = Vec::new();
    let mut line_number = 0usize;

    loop {
        line.clear();
        let bytes_read = reader.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            break;
        }
        line_number += 1;
        let raw = strip_line_ending(&line);
        records.push(parse_bed_record(raw, line_number)?);
    }

    records.sort_by(compare_bed_records);
    for record in records {
        write_line(&mut writer, &record.raw)?;
    }
    Ok(())
}

#[derive(Debug)]
struct GffSortRecord {
    raw: Vec<u8>,
    contig: Vec<u8>,
    start: u64,
    end: u64,
    ids: Vec<Vec<u8>>,
    parents: Vec<Vec<u8>>,
    source_index: usize,
}

#[derive(Debug)]
struct BedSortRecord {
    raw: Vec<u8>,
    contig: Vec<u8>,
    start: u64,
    end: u64,
    source_index: usize,
}

fn parse_gff_record(raw: &[u8], line_number: usize) -> Result<GffSortRecord> {
    if raw.is_empty() || raw.iter().all(u8::is_ascii_whitespace) {
        return Err(invalid_line(
            "GFF",
            line_number,
            "expected a 9-column record, found a blank line",
        ));
    }

    let mut parser = gff::io::Reader::new(Cursor::new(raw));
    let mut line = gff::Line::default();
    parser.read_line(&mut line).map_err(|error| {
        invalid_line(
            "GFF",
            line_number,
            format!("could not parse record: {error}"),
        )
    })?;
    let record = line
        .as_record()
        .ok_or_else(|| invalid_line("GFF", line_number, "expected a feature record"))?
        .map_err(|error| {
            invalid_line(
                "GFF",
                line_number,
                format!("could not parse record: {error}"),
            )
        })?;

    let contig = record.reference_sequence_name().to_vec();
    if contig.is_empty() {
        return Err(invalid_line(
            "GFF",
            line_number,
            "contig column must not be empty",
        ));
    }
    let start = record
        .start()
        .map_err(|error| invalid_line("GFF", line_number, format!("invalid start: {error}")))?
        .get() as u64;
    let end = record
        .end()
        .map_err(|error| invalid_line("GFF", line_number, format!("invalid end: {error}")))?
        .get() as u64;
    if start == 0 || end == 0 || end < start {
        return Err(invalid_line(
            "GFF",
            line_number,
            format!("coordinates must satisfy 1 <= start <= end (got {start}..{end})"),
        ));
    }

    let record_buf = gff::feature::RecordBuf::try_from_feature_record(&record)
        .map_err(|error| invalid_line("GFF", line_number, format!("invalid record: {error}")))?;
    let mut ids = Vec::new();
    let mut parents = Vec::new();
    for (tag, value) in record_buf.attributes().as_ref() {
        let tag = <_ as AsRef<[u8]>>::as_ref(tag);
        if tag != b"ID" && tag != b"Parent" {
            continue;
        }
        for value in value.iter() {
            let value = <_ as AsRef<[u8]>>::as_ref(value);
            if value.is_empty() {
                continue;
            }
            if tag == b"ID" {
                ids.push(value.to_vec());
            } else {
                parents.push(value.to_vec());
            }
        }
    }

    Ok(GffSortRecord {
        raw: raw.to_vec(),
        contig,
        start,
        end,
        ids,
        parents,
        source_index: line_number,
    })
}

fn parse_bed_record(raw: &[u8], line_number: usize) -> Result<BedSortRecord> {
    if raw.is_empty() || raw.iter().all(u8::is_ascii_whitespace) {
        return Err(invalid_line(
            "BED",
            line_number,
            "expected a record with contig, start, and end columns",
        ));
    }
    let fields = raw
        .split(|byte| byte.is_ascii_whitespace())
        .filter(|field| !field.is_empty())
        .collect::<Vec<_>>();
    if fields.len() < 3 {
        return Err(invalid_line(
            "BED",
            line_number,
            format!(
                "expected at least 3 tab-separated columns, found {}",
                fields.len()
            ),
        ));
    }
    if fields[0].is_empty() {
        return Err(invalid_line(
            "BED",
            line_number,
            "contig column must not be empty",
        ));
    }
    let start = parse_coordinate(fields[1], "start", "BED", line_number)?;
    let end = parse_coordinate(fields[2], "end", "BED", line_number)?;
    if end < start {
        return Err(invalid_line(
            "BED",
            line_number,
            format!("coordinates must satisfy start <= end (got {start}..{end})"),
        ));
    }
    Ok(BedSortRecord {
        raw: raw.to_vec(),
        contig: fields[0].to_vec(),
        start,
        end,
        source_index: line_number,
    })
}

fn parse_coordinate(value: &[u8], name: &str, format: &str, line_number: usize) -> Result<u64> {
    let value = std::str::from_utf8(value).map_err(|_| {
        invalid_line(
            format,
            line_number,
            format!("{name} coordinate is not UTF-8"),
        )
    })?;
    value.parse::<u64>().map_err(|error| {
        invalid_line(
            format,
            line_number,
            format!("invalid {name} coordinate {value:?}: {error}"),
        )
    })
}

fn compare_gff_coordinates(left: &GffSortRecord, right: &GffSortRecord) -> Ordering {
    left.contig
        .cmp(&right.contig)
        .then_with(|| left.start.cmp(&right.start))
        .then_with(|| left.end.cmp(&right.end))
}

fn compare_bed_records(left: &BedSortRecord, right: &BedSortRecord) -> Ordering {
    left.contig
        .cmp(&right.contig)
        .then_with(|| left.start.cmp(&right.start))
        .then_with(|| left.end.cmp(&right.end))
        .then_with(|| left.source_index.cmp(&right.source_index))
}

fn order_gff_tie_group(records: &[GffSortRecord]) -> Result<Vec<usize>> {
    let mut id_to_records: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for (index, record) in records.iter().enumerate() {
        for id in &record.ids {
            id_to_records.entry(id).or_default().push(index);
        }
    }

    let mut children = vec![Vec::new(); records.len()];
    let mut indegree = vec![0usize; records.len()];
    let mut edges = HashSet::new();
    for (child, record) in records.iter().enumerate() {
        for parent_id in &record.parents {
            let Some(parent_records) = id_to_records.get(parent_id.as_slice()) else {
                continue;
            };
            for &parent in parent_records {
                if edges.insert((parent, child)) {
                    children[parent].push(child);
                    indegree[child] += 1;
                }
            }
        }
    }

    let mut ready = BTreeSet::new();
    for (index, (&degree, record)) in indegree.iter().zip(records).enumerate() {
        if degree == 0 {
            ready.insert((record.source_index, index));
        }
    }

    let mut order = Vec::with_capacity(records.len());
    while let Some((_, index)) = ready.pop_first() {
        order.push(index);
        for child in &children[index] {
            indegree[*child] -= 1;
            if indegree[*child] == 0 {
                ready.insert((records[*child].source_index, *child));
            }
        }
    }
    if order.len() != records.len() {
        return Err(Error::InvalidInput(format!(
            "GFF parent hierarchy contains a cycle among records on contig {:?} at {}..{}",
            String::from_utf8_lossy(&records[0].contig),
            records[0].start,
            records[0].end
        )));
    }
    Ok(order)
}

fn strip_line_ending(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

fn write_line(writer: &mut impl Write, line: &[u8]) -> io::Result<()> {
    writer.write_all(line)?;
    writer.write_all(b"\n")
}

fn invalid_line(format: &str, line_number: usize, message: impl Into<String>) -> Error {
    Error::InvalidInput(format!("{format} line {line_number}: {}", message.into()))
}
