use std::{
    collections::HashMap,
    fs,
    io::{self, BufRead, Cursor, Write},
    path::{Path, PathBuf},
};

use noodles::{
    bed, bgzf,
    core::Position,
    csi::{
        self,
        binning_index::{
            Indexer,
            index::{
                header::{Builder as HeaderBuilder, ReferenceSequenceNames},
                reference_sequence::{bin::Chunk, index::BinnedIndex},
            },
        },
    },
    gff,
};
use tempfile::NamedTempFile;

use crate::{Error, Result, SortFormat};

const CSI_MIN_SHIFT: u8 = 14;
const CSI_DEPTH: u8 = 7;
const CSI_MAX_POSITION: u64 = (1u64 << (CSI_MIN_SHIFT + 3 * CSI_DEPTH)) - 1;

/// Compresses an already sorted GFF, GTF, or BED file to BGZF and writes a CSI beside it.
///
/// The output format is inferred from the input extension, including `.gz`, `.bgz`, and `.bgzf`
/// inputs. GTF is parsed with the GFF record model. If `coordinate_index_path` is `None`, the CSI
/// path is `<output>.csi`. Both outputs are staged as sibling temporary files and published only
/// after the source and index are complete. Existing destinations are left untouched.
///
/// # Arguments
///
/// * `input_path` - Path to the already sorted GFF, GTF, or BED source.
/// * `output_path` - Destination for the BGZF-compressed source.
/// * `coordinate_index_path` - Optional CSI destination; defaults to `<output>.csi`.
///
/// # Errors
///
/// Returns an error if the input is malformed or unsorted, a coordinate exceeds the configured
/// CSI capacity, or any input/output operation fails.
pub fn compress_file(
    input_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    coordinate_index_path: Option<&Path>,
) -> Result<()> {
    let input_path = input_path.as_ref();
    let output_path = output_path.as_ref();
    let coordinate_index_path = coordinate_index_path_for(output_path, coordinate_index_path);
    validate_output_paths(input_path, output_path, &coordinate_index_path)?;

    let format = SortFormat::from_path(input_path)?;
    let (mut output_temp, index_temp) = temporary_outputs(output_path, &coordinate_index_path)?;

    let mut bgzf_writer = bgzf::io::Writer::new(output_temp.as_file_mut());
    let index = write_bgzf_with_csi(
        crate::sort::open_reader(input_path)?,
        &mut bgzf_writer,
        format,
    )?;
    bgzf_writer.finish()?.flush()?;
    publish_outputs(
        index,
        output_temp,
        index_temp,
        output_path,
        &coordinate_index_path,
    )
}

/// Sorts an annotation file into a BGZF writer while building its CSI index from the same output
/// stream.
///
/// The format is inferred from the input extension. Sorted annotation lines are indexed as they
/// are written, so the source is read and sorted once and no sorted-text intermediate is created.
/// GFF comments and directives are skipped for indexing, and the sequence payload after
/// `##FASTA` is streamed without line buffering. The caller must finish `writer` before using or
/// serializing the returned index.
pub fn sort_bgzf_with_csi<W: Write>(
    input_path: impl AsRef<Path>,
    disk_sort: bool,
    writer: &mut bgzf::io::Writer<W>,
) -> Result<csi::Index> {
    let input_path = input_path.as_ref();
    let format = SortFormat::from_path(input_path)?;
    let mut output = SortIndexingWriter::new(writer, format);
    let sort_result = crate::sort::sort_file(input_path, disk_sort, &mut output);

    if let Some(error) = output.index_error.take() {
        return Err(error);
    }
    sort_result?;
    output.finish()
}

/// Sorts an annotation file directly to BGZF and publishes the compressed source and CSI.
///
/// Both outputs are staged as sibling temporary files and published only after sorting,
/// compression, and index serialization complete. Existing destinations are left untouched.
/// When `coordinate_index_path` is `None`, the index path is `<output>.csi`.
pub fn sort_and_compress_file(
    input_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    coordinate_index_path: Option<&Path>,
    disk_sort: bool,
) -> Result<()> {
    let input_path = input_path.as_ref();
    let output_path = output_path.as_ref();
    let coordinate_index_path = coordinate_index_path_for(output_path, coordinate_index_path);
    validate_output_paths(input_path, output_path, &coordinate_index_path)?;

    let (mut output_temp, index_temp) = temporary_outputs(output_path, &coordinate_index_path)?;
    let mut bgzf_writer = bgzf::io::Writer::new(output_temp.as_file_mut());
    let index = sort_bgzf_with_csi(input_path, disk_sort, &mut bgzf_writer)?;
    bgzf_writer.finish()?.flush()?;

    publish_outputs(
        index,
        output_temp,
        index_temp,
        output_path,
        &coordinate_index_path,
    )
}

fn coordinate_index_path_for(output_path: &Path, coordinate_index_path: Option<&Path>) -> PathBuf {
    coordinate_index_path
        .map(Path::to_path_buf)
        .unwrap_or_else(|| {
            let mut path = output_path.as_os_str().to_os_string();
            path.push(".csi");
            PathBuf::from(path)
        })
}

fn validate_output_paths(input_path: &Path, output_path: &Path, index_path: &Path) -> Result<()> {
    let input_identity = normalized_path(input_path)?;
    let output_identity = normalized_path(output_path)?;
    let index_identity = normalized_path(index_path)?;
    if input_identity == output_identity
        || input_identity == index_identity
        || output_identity == index_identity
    {
        return Err(Error::InvalidInput(
            "input, BGZF output, and CSI paths must be distinct".into(),
        ));
    }
    if output_path.exists() || index_path.exists() {
        return Err(Error::InvalidInput(
            "BGZF output and CSI destinations must not already exist".into(),
        ));
    }

    Ok(())
}

fn temporary_outputs(
    output_path: &Path,
    index_path: &Path,
) -> Result<(NamedTempFile, NamedTempFile)> {
    let output_temp = NamedTempFile::new_in(parent_directory(output_path))?;
    let index_temp = NamedTempFile::new_in(parent_directory(index_path))?;
    Ok((output_temp, index_temp))
}

fn publish_outputs(
    index: csi::Index,
    mut output_temp: NamedTempFile,
    mut index_temp: NamedTempFile,
    output_path: &Path,
    index_path: &Path,
) -> Result<()> {
    output_temp.flush()?;
    let mut csi_writer = csi::io::Writer::new(index_temp.as_file_mut());
    csi_writer.write_index(&index)?;
    csi_writer.into_inner().finish()?.flush()?;
    index_temp.flush()?;

    index_temp
        .persist_noclobber(index_path)
        .map_err(|error| Error::Io(error.error))?;
    if let Err(error) = output_temp.persist_noclobber(output_path) {
        let _ = fs::remove_file(index_path);
        return Err(Error::Io(error.error));
    }

    Ok(())
}

fn normalized_path(path: &Path) -> Result<PathBuf> {
    let absolute_path = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };

    if absolute_path.exists() {
        return fs::canonicalize(absolute_path).map_err(Error::Io);
    }

    let parent = absolute_path
        .parent()
        .ok_or_else(|| Error::InvalidInput("output path has no parent directory".into()))?;
    let file_name = absolute_path
        .file_name()
        .ok_or_else(|| Error::InvalidInput("output path has no file name".into()))?;
    Ok(fs::canonicalize(parent)?.join(file_name))
}

fn parent_directory(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

/// Writes an already sorted GFF or BED stream to BGZF while building its CSI index.
///
/// Each original line is copied unchanged, and its BGZF virtual offsets are recorded as it is
/// written. The input must be coordinate sorted. GFF directives, comments, BED headers, and any
/// sequence payload after `##FASTA` are retained; only feature records are indexed. The returned
/// CSI uses a 14-bit minimum shift and depth 7, supporting positions through 34,359,738,367. The
/// caller must finish `writer` before serializing or using the returned index so all compressed
/// blocks and the BGZF terminator have been written.
///
/// # Errors
///
/// Returns an error if the input is malformed or unsorted, a coordinate exceeds CSI capacity, or
/// reading or writing fails.
pub fn write_bgzf_with_csi<R, W>(
    mut reader: R,
    writer: &mut bgzf::io::Writer<W>,
    format: SortFormat,
) -> Result<csi::Index>
where
    R: BufRead,
    W: Write,
{
    let mut indexer = CsiIndexBuilder::new(format);
    let mut line = Vec::new();

    loop {
        let start_virtual_position = writer.virtual_position();
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }

        writer.write_all(&line)?;
        let end_virtual_position = writer.virtual_position();
        indexer.add_line(&line, start_virtual_position, end_virtual_position)?;
    }

    Ok(indexer.finish())
}

struct CsiIndexBuilder {
    format: SortFormat,
    indexer: Indexer<BinnedIndex>,
    reference_sequence_ids: HashMap<Vec<u8>, usize>,
    reference_sequence_names: Vec<Vec<u8>>,
    previous_record: Option<(usize, Position)>,
    in_fasta_payload: bool,
}

impl CsiIndexBuilder {
    fn new(format: SortFormat) -> Self {
        Self {
            format,
            indexer: Indexer::<BinnedIndex>::new(CSI_MIN_SHIFT, CSI_DEPTH),
            reference_sequence_ids: HashMap::new(),
            reference_sequence_names: Vec::new(),
            previous_record: None,
            in_fasta_payload: false,
        }
    }

    fn add_line(
        &mut self,
        line: &[u8],
        start_virtual_position: bgzf::VirtualPosition,
        end_virtual_position: bgzf::VirtualPosition,
    ) -> Result<()> {
        if self.in_fasta_payload {
            return Ok(());
        }

        let record_line = strip_line_ending(line);
        if self.format == SortFormat::Gff && record_line == b"##FASTA" {
            self.in_fasta_payload = true;
            return Ok(());
        }

        if record_line.is_empty()
            || record_line[0] == b'#'
            || record_line.starts_with(b"track ")
            || record_line.starts_with(b"browser ")
        {
            return Ok(());
        }

        let (reference_sequence_name, start, end) = match self.format {
            SortFormat::Gff => parse_gff_record(line)?,
            SortFormat::Bed => parse_bed_record(line)?,
        };

        let start_coordinate = u64::try_from(start.get()).map_err(|_| Error::InvalidCoordinate)?;
        let end_coordinate = u64::try_from(end.get()).map_err(|_| Error::InvalidCoordinate)?;
        if reference_sequence_name.is_empty() {
            return Err(Error::InvalidInput(
                "annotation reference sequence name must not be empty".into(),
            ));
        }

        if start_coordinate > CSI_MAX_POSITION || end_coordinate > CSI_MAX_POSITION {
            return Err(Error::InvalidInput(format!(
                "annotation coordinate exceeds the CSI limit of {CSI_MAX_POSITION}: {}-{}",
                start_coordinate, end_coordinate
            )));
        }

        if start > end {
            return Err(Error::InvalidInput(format!(
                "annotation start exceeds its end: {}-{}",
                start.get(),
                end.get()
            )));
        }

        let reference_sequence_id = match self.reference_sequence_ids.get(&reference_sequence_name)
        {
            Some(reference_sequence_id) => *reference_sequence_id,
            None => {
                let reference_sequence_id = self.reference_sequence_names.len();
                self.reference_sequence_ids
                    .insert(reference_sequence_name.clone(), reference_sequence_id);
                self.reference_sequence_names.push(reference_sequence_name);
                reference_sequence_id
            }
        };

        if let Some((previous_reference_sequence_id, previous_start)) = self.previous_record
            && (reference_sequence_id < previous_reference_sequence_id
                || (reference_sequence_id == previous_reference_sequence_id
                    && start < previous_start))
        {
            return Err(Error::InvalidInput(
                "annotation records must be sorted by reference sequence and start coordinate"
                    .into(),
            ));
        }

        self.indexer.add_record(
            Some((reference_sequence_id, start, end, true)),
            Chunk::new(start_virtual_position, end_virtual_position),
        )?;
        self.previous_record = Some((reference_sequence_id, start));
        Ok(())
    }

    fn finish(self) -> csi::Index {
        let header_builder = match self.format {
            SortFormat::Gff => HeaderBuilder::gff(),
            SortFormat::Bed => HeaderBuilder::bed(),
        };
        let mut names = ReferenceSequenceNames::new();
        for reference_sequence_name in &self.reference_sequence_names {
            names.insert(reference_sequence_name.as_slice().into());
        }

        let header = header_builder.set_reference_sequence_names(names).build();
        self.indexer
            .set_header(header)
            .build(self.reference_sequence_names.len())
    }
}

struct SortIndexingWriter<'a, W: Write> {
    writer: &'a mut bgzf::io::Writer<W>,
    indexer: CsiIndexBuilder,
    line: Vec<u8>,
    line_start: Option<bgzf::VirtualPosition>,
    index_error: Option<Error>,
}

impl<'a, W: Write> SortIndexingWriter<'a, W> {
    fn new(writer: &'a mut bgzf::io::Writer<W>, format: SortFormat) -> Self {
        Self {
            writer,
            indexer: CsiIndexBuilder::new(format),
            line: Vec::new(),
            line_start: None,
            index_error: None,
        }
    }

    fn process_line(&mut self, end_virtual_position: bgzf::VirtualPosition) -> io::Result<()> {
        let start_virtual_position = self
            .line_start
            .take()
            .expect("a buffered line has a virtual start position");
        if let Err(error) =
            self.indexer
                .add_line(&self.line, start_virtual_position, end_virtual_position)
        {
            self.index_error = Some(error);
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "could not index sorted annotation line",
            ));
        }
        self.line.clear();
        Ok(())
    }

    fn finish(mut self) -> Result<csi::Index> {
        if let Some(error) = self.index_error.take() {
            return Err(error);
        }
        if !self.line.is_empty()
            && let Err(error) = self.process_line(self.writer.virtual_position())
        {
            return Err(self.index_error.take().unwrap_or(Error::Io(error)));
        }
        Ok(self.indexer.finish())
    }
}

impl<W: Write> Write for SortIndexingWriter<'_, W> {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        if buffer.is_empty() {
            return Ok(0);
        }
        if self.index_error.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "could not index sorted annotation line",
            ));
        }

        let mut offset = 0;
        while offset < buffer.len() {
            if self.indexer.in_fasta_payload {
                self.writer.write_all(&buffer[offset..])?;
                return Ok(buffer.len());
            }

            let line_end = buffer[offset..]
                .iter()
                .position(|&byte| byte == b'\n')
                .map_or(buffer.len(), |index| offset + index + 1);
            let fragment = &buffer[offset..line_end];
            if self.line.is_empty() {
                self.line_start = Some(self.writer.virtual_position());
            }
            self.writer.write_all(fragment)?;
            self.line.extend_from_slice(fragment);
            offset = line_end;

            if fragment.ends_with(b"\n") {
                self.process_line(self.writer.virtual_position())?;
            }
        }

        Ok(buffer.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

fn parse_gff_record(line: &[u8]) -> Result<(Vec<u8>, Position, Position)> {
    let mut reader = gff::io::Reader::new(Cursor::new(line));
    let mut parsed_line = gff::Line::default();
    reader.read_line(&mut parsed_line)?;
    let record = parsed_line
        .as_record()
        .ok_or_else(|| Error::InvalidInput("expected a GFF feature record".into()))??;

    Ok((
        record.reference_sequence_name().to_vec(),
        record.start()?,
        record.end()?,
    ))
}

fn parse_bed_record(line: &[u8]) -> Result<(Vec<u8>, Position, Position)> {
    let mut reader = bed::io::Reader::<3, _>::new(Cursor::new(line));
    let mut record = bed::Record::<3>::default();
    reader.read_record(&mut record)?;

    let start = record.feature_start()?;
    let end = match record.feature_end().transpose()? {
        Some(end) => end,
        None if start == Position::MIN => start,
        None => {
            return Err(Error::InvalidInput(
                "BED end coordinate cannot be zero when start is nonzero".into(),
            ));
        }
    };
    let end = if end < start && start.get() - end.get() == 1 {
        start
    } else {
        end
    };

    Ok((record.reference_sequence_name().to_vec(), start, end))
}

fn strip_line_ending(mut line: &[u8]) -> &[u8] {
    if let Some(stripped) = line.strip_suffix(b"\n") {
        line = stripped;
    }
    if let Some(stripped) = line.strip_suffix(b"\r") {
        line = stripped;
    }
    line
}

#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor, Read};

    use noodles::{bgzf, core::Region, csi, csi::BinningIndex};

    use super::{SortIndexingWriter, write_bgzf_with_csi};
    use crate::SortFormat;

    fn compress_and_query(input: &[u8], format: SortFormat, region: &str) -> (Vec<u8>, Vec<u8>) {
        let mut output = Vec::new();
        let mut writer = bgzf::io::Writer::new(&mut output);
        let index = write_bgzf_with_csi(BufReader::new(input), &mut writer, format)
            .expect("should build CSI while writing BGZF");
        writer.finish().expect("should finish BGZF output");

        let mut reader = bgzf::io::Reader::new(output.as_slice());
        let mut decoded = Vec::new();
        reader
            .read_to_end(&mut decoded)
            .expect("should decompress output");

        let mut index_bytes = Vec::new();
        let mut index_writer = csi::io::Writer::new(&mut index_bytes);
        index_writer
            .write_index(&index)
            .expect("should serialize CSI");
        index_writer
            .into_inner()
            .finish()
            .expect("should finish CSI output");
        let index = csi::io::Reader::new(index_bytes.as_slice())
            .read_index()
            .expect("should read serialized CSI");

        let region = region.parse::<Region>().expect("should parse region");
        let mut indexed_reader = csi::io::IndexedReader::new(Cursor::new(output.as_slice()), index);
        let records = indexed_reader
            .query(&region)
            .expect("should query CSI")
            .map(|result| {
                result
                    .expect("should read indexed record")
                    .as_ref()
                    .as_bytes()
                    .to_vec()
            })
            .collect::<Vec<_>>()
            .concat();

        (decoded, records)
    }

    #[test]
    fn test_write_bgzf_with_csi_preserves_gff_lines_and_skips_fasta_payload() {
        let input = concat!(
            "##gff-version 3\r\n",
            "# comment\r\n",
            "chr1\tsource\tgene\t600000000\t600000010\t.\t+\t.\tName=large\r\n",
            "##FASTA\r\n",
            ">chr1\r\n",
            "ACGT\r\n",
        )
        .as_bytes();

        let (decoded, records) =
            compress_and_query(input, SortFormat::Gff, "chr1:600000000-600000010");

        assert_eq!(decoded, input);
        assert!(String::from_utf8_lossy(&records).contains("Name=large"));
        assert!(!String::from_utf8_lossy(&records).contains("ACGT"));
    }

    #[test]
    fn test_write_bgzf_with_csi_preserves_bed_headers_and_crlf() {
        let input = concat!(
            "track name=genes\r\n",
            "browser position chr1:1-10\r\n",
            "# comment\r\n",
            "chr1\t0\t10\tAlpha\r\n",
        )
        .as_bytes();

        let (decoded, records) = compress_and_query(input, SortFormat::Bed, "chr1:1-10");

        assert_eq!(decoded, input);
        assert_eq!(records, b"chr1\t0\t10\tAlpha");
    }

    #[test]
    fn test_write_bgzf_with_csi_tracks_virtual_offsets_across_bgzf_blocks() {
        let long_name = "x".repeat(70_000);
        let input = format!(
            "chr1\tsource\tgene\t1\t2\t.\t+\t.\tName={long_name}\n\
             chr1\tsource\tgene\t600000000\t600000010\t.\t+\t.\tName=after-block\n"
        );

        let (decoded, first_record) =
            compress_and_query(input.as_bytes(), SortFormat::Gff, "chr1:1-2");
        let (_, following_record) = compress_and_query(
            input.as_bytes(),
            SortFormat::Gff,
            "chr1:600000000-600000010",
        );

        assert_eq!(decoded, input.as_bytes());
        assert!(String::from_utf8_lossy(&first_record).contains(&long_name));
        assert!(!String::from_utf8_lossy(&first_record).contains("after-block"));
        assert!(String::from_utf8_lossy(&following_record).contains("Name=after-block"));
        assert!(!String::from_utf8_lossy(&following_record).contains(&long_name));
    }

    #[test]
    fn test_write_bgzf_with_csi_rejects_coordinates_beyond_index_capacity() {
        let input = b"chr1\tsource\tgene\t34359738368\t34359738368\t.\t+\t.\tName=too-large\n";
        let mut writer = bgzf::io::Writer::new(Vec::new());

        let error = write_bgzf_with_csi(
            BufReader::new(input.as_slice()),
            &mut writer,
            SortFormat::Gff,
        )
        .expect_err("should reject a position beyond the configured CSI depth");

        assert!(error.to_string().contains("CSI limit"));
    }

    #[test]
    fn test_write_bgzf_with_csi_indexes_zero_width_bed_as_a_point_candidate() {
        let input = b"chr1\t2\t2\tpoint\n";
        let mut output = Vec::new();
        let mut writer = bgzf::io::Writer::new(&mut output);
        let index = write_bgzf_with_csi(
            BufReader::new(input.as_slice()),
            &mut writer,
            SortFormat::Bed,
        )
        .expect("should build a CSI entry for the zero-width BED feature");
        writer.finish().expect("should finish BGZF output");

        let region = "chr1:3-3"
            .parse::<Region>()
            .expect("should parse point region");
        let chunks = index
            .query(0, region.interval())
            .expect("should query point candidate bins");
        assert_eq!(chunks.len(), 1);
    }

    #[test]
    fn test_compress_file_writes_default_csi_and_supports_gai_queries() {
        use std::fs;

        use crate::{BuildOptions, NameIndexOptions, build_name_index_with_options, query_index};

        let temp_dir = tempfile::tempdir().expect("should create temporary directory");
        let input_path = temp_dir.path().join("annotations.gff");
        let output_path = temp_dir.path().join("annotations.gff.bgz");
        let coordinate_index_path = temp_dir.path().join("annotations.gff.bgz.csi");
        let attribute_index_path = temp_dir.path().join("annotations.gai");
        fs::write(
            &input_path,
            b"##gff-version 3\nchr22\tsource\tgene\t600000000\t600000010\t.\t+\t.\tName=BRCA1\n",
        )
        .expect("should write input annotations");

        super::compress_file(&input_path, &output_path, None)
            .expect("should compress and index annotations");
        assert!(coordinate_index_path.exists());

        let name_options = NameIndexOptions::new(["Name"], false)
            .expect("should configure the GFF name attribute");
        build_name_index_with_options(
            &output_path,
            &coordinate_index_path,
            &attribute_index_path,
            &name_options,
            &BuildOptions::default(),
        )
        .expect("should build the attribute index from the generated CSI");
        let records = query_index(
            &output_path,
            &coordinate_index_path,
            &attribute_index_path,
            "BRCA1",
        )
        .expect("should query the generated GFF indexes");

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].reference_sequence_name, "chr22");
        assert_eq!(records[0].start, 600_000_000);
    }

    #[test]
    fn test_compress_file_writes_explicit_csi_for_bed_and_rejects_parse_failure_without_outputs() {
        use std::fs;

        use crate::{BuildOptions, NameIndexOptions, build_name_index_with_options, query_index};

        let temp_dir = tempfile::tempdir().expect("should create temporary directory");
        let input_path = temp_dir.path().join("annotations.bed");
        let output_path = temp_dir.path().join("annotations.bed.bgzf");
        let coordinate_index_path = temp_dir.path().join("indexes/custom.index");
        let attribute_index_path = temp_dir.path().join("annotations.bed.gai");
        fs::create_dir_all(coordinate_index_path.parent().expect("should have parent"))
            .expect("should create index directory");
        fs::write(&input_path, b"chr2\t600000000\t600000010\tMYO2B\n")
            .expect("should write BED annotations");

        super::compress_file(&input_path, &output_path, Some(&coordinate_index_path))
            .expect("should compress BED annotations with an explicit index path");
        assert!(output_path.exists());
        assert!(coordinate_index_path.exists());

        build_name_index_with_options(
            &output_path,
            &coordinate_index_path,
            &attribute_index_path,
            &NameIndexOptions::bed(false),
            &BuildOptions::default(),
        )
        .expect("should build BED attribute index from generated CSI");
        let matching_records = query_index(
            &output_path,
            &coordinate_index_path,
            &attribute_index_path,
            "MYO2B",
        )
        .expect("should query BED attribute and coordinate indexes");
        assert_eq!(matching_records.len(), 1);
        assert_eq!(matching_records[0].start, 600_000_001);

        let index = csi::fs::read(&coordinate_index_path).expect("should read generated CSI");
        let region = "chr2:600000001-600000010"
            .parse()
            .expect("should parse query region");
        let records = csi::io::IndexedReader::new(
            fs::File::open(&output_path).expect("should open BGZF source"),
            index,
        )
        .query(&region)
        .expect("should query BED CSI")
        .map(|record| {
            record
                .expect("should read indexed BED record")
                .as_ref()
                .as_bytes()
                .to_vec()
        })
        .collect::<Vec<_>>();
        assert_eq!(records.concat(), b"chr2\t600000000\t600000010\tMYO2B");

        let invalid_input_path = temp_dir.path().join("invalid.gff");
        let invalid_output_path = temp_dir.path().join("invalid.gff.bgz");
        fs::write(&invalid_input_path, b"not a GFF record\n").expect("should write malformed GFF");

        assert!(super::compress_file(&invalid_input_path, &invalid_output_path, None).is_err());
        assert!(!invalid_output_path.exists());
        assert!(!temp_dir.path().join("invalid.gff.bgz.csi").exists());
    }

    #[test]
    fn test_sort_indexing_writer_tracks_fragmented_lines_and_final_unterminated_record() {
        use std::io::Write;

        let long_value = "x".repeat(70_000);
        let first = format!("chr1\tsource\tgene\t1\t2\t.\t+\t.\tName={long_value}\n");
        let final_record = b"chr1\tsource\tgene\t600000000\t600000010\t.\t+\t.\tName=final";
        let mut output = Vec::new();
        let mut writer = bgzf::io::Writer::new(&mut output);
        let mut indexing_writer = SortIndexingWriter::new(&mut writer, SortFormat::Gff);

        indexing_writer
            .write_all(&first.as_bytes()[..17])
            .expect("should accept first fragmented record prefix");
        indexing_writer
            .write_all(&first.as_bytes()[17..65_537])
            .expect("should accept first fragmented record middle");
        indexing_writer
            .write_all(&first.as_bytes()[65_537..])
            .expect("should accept first fragmented record suffix");
        indexing_writer
            .write_all(&final_record[..31])
            .expect("should accept final fragmented record prefix");
        indexing_writer
            .write_all(&final_record[31..])
            .expect("should accept final fragmented record suffix");

        let index = indexing_writer
            .finish()
            .expect("should index final record without a newline");
        writer.finish().expect("should finish BGZF output");

        let region = "chr1:600000000-600000010"
            .parse::<Region>()
            .expect("should parse final record region");
        let records = csi::io::IndexedReader::new(Cursor::new(output), index)
            .query(&region)
            .expect("should query final record")
            .map(|record| {
                record
                    .expect("should read final record")
                    .as_ref()
                    .as_bytes()
                    .to_vec()
            })
            .collect::<Vec<_>>();
        assert_eq!(records.concat(), final_record);
    }
}
