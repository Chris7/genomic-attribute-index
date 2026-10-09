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
use noodles::gff::{self, record::attributes::field::Value as GffAttributeValue};
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
    #[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
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
/// Opens a plain annotation file or a `.gz`, `.bgz`, or `.bgzf` file as a buffered reader.
///
/// Compressed suffixes are decoded as gzip/BGZF; all other suffixes are read as plain text.
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
pub fn open_annotation_reader(path: impl AsRef<Path>) -> io::Result<Box<dyn BufRead>> {
    let path = path.as_ref();
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
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
pub fn sort_file(input: impl AsRef<Path>, disk_sort: bool, writer: impl Write) -> Result<()> {
    let input = input.as_ref();
    let format = SortFormat::from_path(input)?;
    let reader = open_annotation_reader(input)?;
    match format {
        SortFormat::Gff => sort_gff(reader, disk_sort, writer),
        SortFormat::Bed => sort_bed(reader, disk_sort, writer),
    }
}

/// Converts a BufRead into a stream of parsed records.
///
/// Returning Ok(None) from `parse` skips a line, which is useful for
/// comments/header lines.
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
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
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
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
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
pub fn sort_gff<R: BufRead, W: Write>(mut reader: R, disk_sort: bool, mut writer: W) -> Result<()> {
    let mut comments = Vec::new();
    let mut fasta_line = None;
    let mut line_number = 0;

    let mut line = Vec::new();
    let records = std::iter::from_fn(|| {
        if fasta_line.is_some() {
            return None;
        }

        loop {
            line.clear();

            match reader.read_until(b'\n', &mut line) {
                Ok(0) => return None,
                Err(err) => return Some(Err(err)),
                Ok(_) => {}
            }

            line_number += 1;
            let raw = strip_line_ending(&line);

            // This must be checked before the generic comment handling.
            if raw == b"##FASTA" {
                // Keep the original bytes, including its line ending.
                fasta_line = Some(std::mem::take(&mut line));
                return None;
            }

            if raw.first() == Some(&b'#') {
                comments.push(raw.to_vec());
                continue;
            }

            // The sort record retains this input buffer instead of copying its bytes.
            let raw = std::mem::take(&mut line);
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
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
fn write_gff_group<W: Write>(writer: &mut W, records: &[GffSortRecord]) -> Result<()> {
    for index in order_gff_tie_group(records)? {
        write_line(writer, &records[index].raw)?;
    }

    Ok(())
}

/// Sorts BED records by contig, start, and end position.
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
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
    #[serde(with = "raw_bytes")]
    raw: Vec<u8>,
    contig: String,
    start: u64,
    end: u64,
    #[serde(with = "raw_byte_lists")]
    ids: Vec<Vec<u8>>,
    #[serde(with = "raw_byte_lists")]
    parents: Vec<Vec<u8>>,
    source_index: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct BedSortRecord {
    #[serde(with = "raw_bytes")]
    raw: Vec<u8>,
    contig: String,
    start: u64,
    end: u64,
    source_index: usize,
}

// Encode each raw row as bulk bytes instead of serializing every byte as an element in ext-sort chunks.
mod raw_bytes {
    use serde::{
        Deserializer, Serializer,
        de::{Error as DeError, Visitor},
    };

    pub(super) fn serialize<S>(raw: &[u8], serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_bytes(raw)
    }

    pub(super) fn deserialize<'de, D>(deserializer: D) -> Result<Vec<u8>, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct RawBytesVisitor;

        impl<'de> Visitor<'de> for RawBytesVisitor {
            type Value = Vec<u8>;

            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("a MessagePack binary value")
            }

            fn visit_bytes<E>(self, value: &[u8]) -> Result<Self::Value, E>
            where
                E: DeError,
            {
                Ok(value.to_vec())
            }

            fn visit_byte_buf<E>(self, value: Vec<u8>) -> Result<Self::Value, E>
            where
                E: DeError,
            {
                Ok(value)
            }
        }

        deserializer.deserialize_bytes(RawBytesVisitor)
    }
}

mod raw_byte_lists {
    use std::fmt;

    use serde::{
        Deserialize, Deserializer, Serialize, Serializer,
        de::{SeqAccess, Visitor},
        ser::SerializeSeq,
    };

    pub(super) fn serialize<S>(values: &[Vec<u8>], serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut sequence = serializer.serialize_seq(Some(values.len()))?;
        for value in values {
            sequence.serialize_element(&RawBytesRef(value))?;
        }
        sequence.end()
    }

    pub(super) fn deserialize<'de, D>(deserializer: D) -> Result<Vec<Vec<u8>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct RawByteListsVisitor;

        impl<'de> Visitor<'de> for RawByteListsVisitor {
            type Value = Vec<Vec<u8>>;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("a sequence of MessagePack binary values")
            }

            fn visit_seq<A>(self, mut sequence: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let mut values = Vec::with_capacity(sequence.size_hint().unwrap_or(0));
                while let Some(value) = sequence.next_element::<RawBytes>()? {
                    values.push(value.0);
                }
                Ok(values)
            }
        }

        deserializer.deserialize_seq(RawByteListsVisitor)
    }

    struct RawBytesRef<'a>(&'a [u8]);

    impl Serialize for RawBytesRef<'_> {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            super::raw_bytes::serialize(self.0, serializer)
        }
    }

    struct RawBytes(Vec<u8>);

    impl<'de> Deserialize<'de> for RawBytes {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>,
        {
            super::raw_bytes::deserialize(deserializer).map(Self)
        }
    }
}

#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
fn parse_gff_record(mut raw: Vec<u8>, line_number: usize) -> Result<GffSortRecord> {
    raw.truncate(strip_line_ending(&raw).len());
    if raw.is_empty() || raw.iter().all(u8::is_ascii_whitespace) {
        return Err(invalid_line(
            "GFF",
            line_number,
            "expected a 9-column record, found a blank line",
        ));
    }

    let mut parser = gff::io::Reader::new(Cursor::new(raw.as_slice()));
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

    record
        .score()
        .transpose()
        .map_err(|error| invalid_line("GFF", line_number, format!("invalid record: {error}")))?;
    record
        .strand()
        .map_err(|error| invalid_line("GFF", line_number, format!("invalid record: {error}")))?;
    record
        .phase()
        .transpose()
        .map_err(|error| invalid_line("GFF", line_number, format!("invalid record: {error}")))?;

    let mut ids = Vec::new();
    let mut parents = Vec::new();
    for result in record.attributes().iter() {
        let (tag, value) = result.map_err(|error| {
            invalid_line("GFF", line_number, format!("invalid record: {error}"))
        })?;
        let tag = <_ as AsRef<[u8]>>::as_ref(tag.as_ref());
        let is_id = tag == b"ID";
        let is_parent = tag == b"Parent";

        // RecordBuf's IndexMap keeps the last value for each percent-decoded tag.
        if is_id {
            ids.clear();
        } else if is_parent {
            parents.clear();
        }

        match value {
            GffAttributeValue::String(value) => {
                let value = <_ as AsRef<[u8]>>::as_ref(value.as_ref());
                if !value.is_empty() {
                    if is_id {
                        ids.push(value.to_vec());
                    } else if is_parent {
                        parents.push(value.to_vec());
                    }
                }
            }
            GffAttributeValue::Array(values) => {
                for value in values.iter() {
                    let value = <_ as AsRef<[u8]>>::as_ref(value.as_ref());
                    if value.is_empty() {
                        continue;
                    }
                    if is_id {
                        ids.push(value.to_vec());
                    } else if is_parent {
                        parents.push(value.to_vec());
                    }
                }
            }
        }
    }

    Ok(GffSortRecord {
        raw,
        contig,
        start,
        end,
        ids,
        parents,
        source_index: line_number,
    })
}
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
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
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
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
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
fn compare_gff_coordinates(left: &GffSortRecord, right: &GffSortRecord) -> Ordering {
    natord::compare_ignore_case(&left.contig, &right.contig)
        .then_with(|| left.start.cmp(&right.start))
        .then_with(|| left.end.cmp(&right.end))
}
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
fn compare_bed_records(left: &BedSortRecord, right: &BedSortRecord) -> Ordering {
    natord::compare_ignore_case(&left.contig, &right.contig)
        .then_with(|| left.start.cmp(&right.start))
        .then_with(|| left.end.cmp(&right.end))
        .then_with(|| left.source_index.cmp(&right.source_index))
}
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
fn order_gff_tie_group(records: &[GffSortRecord]) -> Result<Vec<usize>> {
    if records.len() == 1 {
        let record = &records[0];
        if record
            .ids
            .iter()
            .any(|id| record.parents.iter().any(|parent| parent == id))
        {
            return Err(Error::InvalidInput(format!(
                "GFF parent hierarchy contains a cycle among records on contig {:?} at {}..{}",
                record.contig, record.start, record.end
            )));
        }
        return Ok(vec![0]);
    }

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
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
fn strip_line_ending(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
fn write_line(writer: &mut impl Write, line: &[u8]) -> io::Result<()> {
    writer.write_all(line)?;
    writer.write_all(b"\n")
}
#[cfg_attr(feature = "profiling", tracing::instrument(level = "trace", skip_all))]
fn invalid_line(format: &str, line_number: usize, message: impl Into<String>) -> Error {
    Error::InvalidInput(format!("{format} line {line_number}: {}", message.into()))
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Error as IoError};

    use noodles::gff;

    use super::{parse_gff_record, sort_bed, sort_gff};

    type Relationships = (Vec<Vec<u8>>, Vec<Vec<u8>>);

    #[test]
    fn test_disk_sort_gff_spills_preserve_parent_order_and_raw_lines() {
        let mut input = String::from("##gff-version 3\n#source-order-comment\n");
        input.push_str("chr1\ts\tmRNA\t10\t20\t.\t+\t.\tID=child-escaped;Parent=pr%C3%A9%2Cfix\n");
        input.push_str("chr1\ts\tmRNA\t10\t20\t.\t+\t.\tID=child-direct;Parent=direct-λ\n");
        for index in 0..20 {
            input.push_str(&format!(
                "chr1\ts\tgene\t10\t20\t.\t+\t.\tID=independent-{index:02}\n"
            ));
        }
        input.push_str("chr1\ts\tgene\t10\t20\t.\t+\t.\tID=pr%C3%A9%2Cfix\n");
        input.push_str("chr1\ts\tgene\t10\t20\t.\t+\t.\tID=direct-λ\n");
        input.push_str("##FASTA\n>sequence\nACGT\n");

        let mut in_memory = Vec::new();
        sort_gff(Cursor::new(input.as_bytes()), false, &mut in_memory)
            .expect("in-memory GFF sort should succeed");
        let mut spilled = Vec::new();
        sort_gff(Cursor::new(input.as_bytes()), true, &mut spilled)
            .expect("multi-chunk GFF sort should succeed");
        assert_eq!(
            spilled, in_memory,
            "Spilled GFF output should match memory output"
        );

        let output = String::from_utf8(spilled).expect("GFF output should remain UTF-8");
        assert!(
            output.starts_with("##gff-version 3\n#source-order-comment\n"),
            "GFF comments should remain before sorted records"
        );
        assert!(
            output.ends_with("##FASTA\n>sequence\nACGT\n"),
            "the FASTA section should remain byte-for-byte intact"
        );

        let record_lines = output.lines().collect::<Vec<_>>();
        let line_position = |needle: &str| {
            record_lines
                .iter()
                .position(|line| line.contains(needle))
                .expect("expected feature should be present")
        };
        assert!(
            line_position("ID=pr%C3%A9%2Cfix") < line_position("ID=child-escaped"),
            "percent-decoded UTF-8 parent should precede its child across spills"
        );
        assert!(
            line_position("ID=direct-λ") < line_position("ID=child-direct"),
            "literal UTF-8 parent should precede its child across spills"
        );
        let stable_positions = (0..20)
            .map(|index| line_position(&format!("ID=independent-{index:02}")))
            .collect::<Vec<_>>();
        assert!(
            stable_positions.windows(2).all(|pair| pair[0] < pair[1]),
            "unrelated GFF ties should retain source order"
        );
    }

    #[test]
    fn test_disk_sort_gff_spills_reject_parent_cycles() {
        let self_parent = b"chr1\ts\tgene\t1\t2\t.\t+\t.\tID=self;Parent=self\n";
        let self_error = sort_gff(Cursor::new(self_parent), true, Vec::new())
            .expect_err("a self-parent record should be rejected");
        assert!(
            self_error
                .to_string()
                .contains("parent hierarchy contains a cycle"),
            "self-parent error should identify the hierarchy cycle"
        );

        let cross_parent = concat!(
            "chr1\ts\tgene\t1\t2\t.\t+\t.\tID=a;Parent=b\n",
            "chr1\ts\tgene\t1\t2\t.\t+\t.\tID=b;Parent=a\n",
        );
        let cross_error = sort_gff(Cursor::new(cross_parent.as_bytes()), true, Vec::new())
            .expect_err("a cross-parent cycle should be rejected");
        assert!(
            cross_error
                .to_string()
                .contains("parent hierarchy contains a cycle"),
            "cross-parent error should identify the hierarchy cycle"
        );
    }

    #[test]
    fn test_disk_sort_bed_spills_preserve_stable_ties_and_raw_lines() {
        let mut input = String::from("chr2\t0\t5\tother-contig\textra\n");
        for index in 0..20 {
            input.push_str(&format!("chr1\t10\t20\ttie-{index:02}\textra-{index:02}\n"));
        }
        input.push_str("chr1\t1\t30\tearly\tunaltered-extra-column\n");

        let mut in_memory = Vec::new();
        sort_bed(Cursor::new(input.as_bytes()), false, &mut in_memory)
            .expect("in-memory BED sort should succeed");
        let mut spilled = Vec::new();
        sort_bed(Cursor::new(input.as_bytes()), true, &mut spilled)
            .expect("multi-chunk BED sort should succeed");
        assert_eq!(
            spilled, in_memory,
            "Spilled BED output should match memory output"
        );

        let output = String::from_utf8(spilled).expect("BED output should remain UTF-8");
        let records = output.lines().collect::<Vec<_>>();
        assert_eq!(
            records.first(),
            Some(&"chr1\t1\t30\tearly\tunaltered-extra-column")
        );
        assert_eq!(records.last(), Some(&"chr2\t0\t5\tother-contig\textra"));
        let stable_positions = (0..20)
            .map(|index| {
                records
                    .iter()
                    .position(|line| line.contains(&format!("\ttie-{index:02}\t")))
                    .expect("expected BED record should be present")
            })
            .collect::<Vec<_>>();
        assert!(
            stable_positions.windows(2).all(|pair| pair[0] < pair[1]),
            "equal-coordinate BED records should retain source order"
        );
    }
    #[test]
    fn test_gff_sort_record_relationships_match_record_buf_duplicates() {
        let raw = b"chr1\ts\tgene\t1\t2\t.\t+\t.\tID=first;I%44=last;Parent=old;Par%65nt=new%2Ctag,second;other=a,b\n";
        let expected = record_buf_relationships(raw)
            .expect("RecordBuf should decode duplicate and array-valued attributes");
        let actual = parse_gff_record(raw.to_vec(), 1)
            .expect("the streaming parser should accept the same record");

        assert_eq!(actual.ids, expected.0);
        assert_eq!(actual.parents, expected.1);
        assert_eq!(
            expected,
            (
                vec![b"last".to_vec()],
                vec![b"new,tag".to_vec(), b"second".to_vec()],
            ),
            "last duplicate decoded tags should replace prior values"
        );
    }

    #[test]
    fn test_gff_sort_record_validation_matches_record_buf_for_unused_fields() {
        let invalid_records = [
            (
                "malformed unused attribute",
                gff_record_line(".", "+", ".", "ID=gene;unused"),
            ),
            (
                "invalid unused score",
                gff_record_line("not-a-score", "+", ".", "ID=gene"),
            ),
            (
                "invalid unused strand",
                gff_record_line(".", "x", ".", "ID=gene"),
            ),
            (
                "invalid unused phase",
                gff_record_line(".", "+", "3", "ID=gene"),
            ),
        ];

        for (description, raw) in invalid_records {
            assert!(
                record_buf_relationships(&raw).is_err(),
                "RecordBuf should reject {description}"
            );
            assert!(
                parse_gff_record(raw, 1).is_err(),
                "the streaming parser should reject {description}"
            );
        }
    }

    fn gff_record_line(score: &str, strand: &str, phase: &str, attributes: &str) -> Vec<u8> {
        format!("chr1\ts\tgene\t1\t2\t{score}\t{strand}\t{phase}\t{attributes}").into_bytes()
    }

    fn record_buf_relationships(raw: &[u8]) -> std::io::Result<Relationships> {
        let mut reader = gff::io::Reader::new(Cursor::new(raw));
        let mut line = gff::Line::default();
        reader.read_line(&mut line)?;
        let record = line
            .as_record()
            .ok_or_else(|| IoError::other("expected a feature record"))??;
        let record_buf = gff::feature::RecordBuf::try_from_feature_record(&record)?;
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

        Ok((ids, parents))
    }
}
