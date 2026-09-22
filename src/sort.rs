//! Lossless coordinate sorting for GFF3 and BED text files.

use std::{
    cmp::Ordering,
    collections::{BTreeSet, HashMap, HashSet},
    fs::File,
    io::{self, BufRead, BufReader, Cursor, Write},
    path::Path,
};

use ext_sort::{ExternalSorter, ExternalSorterBuilder, LimitedBufferBuilder};
use flate2::bufread::MultiGzDecoder;
use noodles::gff;
use serde::{Deserialize, Serialize, de::DeserializeOwned};

use crate::{Error, Result};

#[cfg(not(test))]
const SORT_CHUNK_RECORDS: usize = 250_000;

#[cfg(test)]
const SORT_CHUNK_RECORDS: usize = 10;

/// A supported annotation format for [`sort_file`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SortFormat {
    /// GFF or GFF3, including comment and directive lines.
    Gff,
    /// BED with at least the three required coordinate columns.
    Bed,
}

impl SortFormat {
    /// Infers a format from a path's extension, including compressed inputs
    /// such as .gff.gz, .gff3.bgz, and .bed.bgzf.
    pub fn from_path(path: &Path) -> Result<Self> {
        let mut path = path;

        if path
            .extension()
            .and_then(|value| value.to_str())
            .is_some_and(|ext| {
                ext.eq_ignore_ascii_case("gz")
                    || ext.eq_ignore_ascii_case("bgz")
                    || ext.eq_ignore_ascii_case("bgzf")
            })
        {
            path = path.file_stem().map(Path::new).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "cannot infer sort input format from {}",
                    path.display()
                ))
            })?;
        }

        match path.extension().and_then(|value| value.to_str()) {
            Some(ext)
                if ext.eq_ignore_ascii_case("gff")
                    || ext.eq_ignore_ascii_case("gff3")
                    || ext.eq_ignore_ascii_case("gtf") =>
            {
                Ok(Self::Gff)
            }

            Some(ext) if ext.eq_ignore_ascii_case("bed") => Ok(Self::Bed),

            Some(ext) => Err(Error::InvalidInput(format!(
                "unsupported sort input extension .{ext}; expected .gff, .gff3, .gtf, or .bed \
                 (optionally followed by .gz, .bgz, or .bgzf)"
            ))),

            None => Err(Error::InvalidInput(
                "cannot infer sort input format; expected .gff, .gff3, .gtf, or .bed \
                 (optionally followed by .gz, .bgz, or .bgzf)"
                    .into(),
            )),
        }
    }
}

pub fn open_reader(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    match path.extension().and_then(|ext| ext.to_str()) {
        Some(ext)
            if ext.eq_ignore_ascii_case("gz")
                || ext.eq_ignore_ascii_case("bgz")
                || ext.eq_ignore_ascii_case("bgzf") =>
        {
            Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
        }

        _ => Ok(Box::new(reader)),
    }
}

/// Sorts an annotation file inferred from its extension and writes records to
/// `writer` in lossless text order.
pub fn sort_file(input: impl AsRef<Path>, disk_sort: bool, writer: impl Write) -> Result<()> {
    let input = input.as_ref();
    let format = SortFormat::from_path(input)?;
    let reader = open_reader(input)?;
    match format {
        SortFormat::Gff => sort_gff(reader, disk_sort, writer),
        SortFormat::Bed => sort_bed(reader, disk_sort, writer),
    }
}

/// Converts a BufRead into a stream of parsed records.
///
/// Returning Ok(None) from `parse` skips a line, which is useful for
/// comments/header lines.
fn read_records<R, T, F>(mut reader: R, mut parse: F) -> impl Iterator<Item = io::Result<T>>
where
    R: BufRead,
    F: FnMut(&[u8], usize) -> io::Result<Option<T>>,
{
    let mut line = Vec::new();
    let mut line_number = 0usize;

    std::iter::from_fn(move || {
        loop {
            line.clear();

            match reader.read_until(b'\n', &mut line) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(error) => return Some(Err(error)),
            }

            line_number += 1;
            let raw = strip_line_ending(&line);

            match parse(raw, line_number) {
                Ok(Some(record)) => return Some(Ok(record)),
                Ok(None) => continue,
                Err(error) => return Some(Err(error)),
            }
        }
    })
}

/// Sort records either in memory or externally on disk.
///
/// The input is consumed completely before this function returns, so callers
/// can safely use state borrowed by the input parser afterward.
fn sort_records<T, I>(
    records: I,
    disk_sort: bool,
    compare: fn(&T, &T) -> Ordering,
) -> Result<Box<dyn Iterator<Item = io::Result<T>>>>
where
    T: Serialize + DeserializeOwned + Send + 'static,
    I: IntoIterator<Item = io::Result<T>>,
{
    let sorter: ExternalSorter<T, io::Error, LimitedBufferBuilder> = ExternalSorterBuilder::new()
        .with_buffer(LimitedBufferBuilder::new(
            if disk_sort {
                SORT_CHUNK_RECORDS
            } else {
                usize::MAX
            },
            true,
        ))
        .build()
        .map_err(io::Error::other)?;

    let records = sorter.sort_by(records, compare).map_err(io::Error::other)?;

    Ok(Box::new(
        records.map(|record| record.map_err(io::Error::other)),
    ))
}

pub fn sort_gff<R: BufRead, W: Write>(mut reader: R, disk_sort: bool, mut writer: W) -> Result<()> {
    let mut comments = Vec::new();
    let mut fasta_line = None;
    let mut line_number = 0;

    let records = std::iter::from_fn(|| {
        if fasta_line.is_some() {
            return None;
        }

        loop {
            let mut line = Vec::new();

            match reader.read_until(b'\n', &mut line) {
                Ok(0) => return None,
                Err(err) => return Some(Err(err)),
                Ok(_) => {}
            }

            line_number += 1;

            // Work with the line minus its newline for parsing/comparison.
            let mut raw = line.as_slice();

            if let Some(stripped) = raw.strip_suffix(b"\n") {
                raw = stripped;
            }
            if let Some(stripped) = raw.strip_suffix(b"\r") {
                raw = stripped;
            }

            // This must be checked before the generic comment handling.
            if raw == b"##FASTA" {
                // Keep the original bytes, including its line ending.
                fasta_line = Some(line);
                return None;
            }

            if raw.first() == Some(&b'#') {
                comments.push(raw.to_vec());
                continue;
            }

            return Some(parse_gff_record(raw, line_number).map_err(io::Error::other));
        }
    });

    // This must completely consume `records`.
    let records = sort_records(records, disk_sort, compare_gff_coordinates)?;

    // The input iterator is finished now, so `reader` is available again.
    for comment in comments {
        write_line(&mut writer, &comment)?;
    }

    let mut group = Vec::new();

    for record in records {
        let record = record?;

        if group
            .first()
            .is_some_and(|first| compare_gff_coordinates(first, &record) != Ordering::Equal)
        {
            write_gff_group(&mut writer, &group)?;
            group.clear();
        }

        group.push(record);
    }

    if !group.is_empty() {
        write_gff_group(&mut writer, &group)?;
    }

    // Everything after ##FASTA is opaque data. Don't parse or buffer it.
    if let Some(fasta_line) = fasta_line {
        writer.write_all(&fasta_line)?;
        io::copy(&mut reader, &mut writer)?;
    }

    Ok(())
}

fn write_gff_group<W: Write>(writer: &mut W, records: &[GffSortRecord]) -> Result<()> {
    for index in order_gff_tie_group(records)? {
        write_line(writer, &records[index].raw)?;
    }

    Ok(())
}

/// Sorts BED records by contig, start, and end position.
pub fn sort_bed<R: BufRead, W: Write>(reader: R, disk_sort: bool, mut writer: W) -> Result<()> {
    let records = read_records(reader, |raw, line_number| {
        parse_bed_record(raw, line_number)
            .map(Some)
            .map_err(io::Error::other)
    });

    let records = sort_records(records, disk_sort, compare_bed_records)?;

    for record in records {
        write_line(&mut writer, &record?.raw)?;
    }

    Ok(())
}

#[derive(Debug, Serialize, Deserialize)]
struct GffSortRecord {
    raw: Vec<u8>,
    contig: String,
    start: u64,
    end: u64,
    ids: Vec<Vec<u8>>,
    parents: Vec<Vec<u8>>,
    source_index: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct BedSortRecord {
    raw: Vec<u8>,
    contig: String,
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

    let contig = record.reference_sequence_name().to_string();
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
        contig: String::from_utf8_lossy(fields[0]).to_string(),
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
    natord::compare_ignore_case(&left.contig, &right.contig)
        .then_with(|| left.start.cmp(&right.start))
        .then_with(|| left.end.cmp(&right.end))
}

fn compare_bed_records(left: &BedSortRecord, right: &BedSortRecord) -> Ordering {
    natord::compare_ignore_case(&left.contig, &right.contig)
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
            records[0].contig, records[0].start, records[0].end
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
