use std::{
    fs,
    io::{Cursor, Read},
    path::{Path, PathBuf},
    process::Command,
};

use flate2::read::MultiGzDecoder;
use gai::{
    BuildOptions, NameIndexOptions, build_name_index_with_options, query_index, sort_bed,
    sort_file, sort_gff,
};
use noodles::{
    bgzf,
    core::Region,
    csi::{self, BinningIndex},
};
use tempfile::tempdir;

fn lines(output: Vec<u8>) -> Vec<String> {
    String::from_utf8(output)
        .expect("sort output should be UTF-8")
        .lines()
        .map(str::to_owned)
        .collect()
}

fn decode_bgzf(input: &[u8]) -> Vec<u8> {
    let mut reader = bgzf::io::Reader::new(Cursor::new(input));
    let mut output = Vec::new();
    reader
        .read_to_end(&mut output)
        .expect("BGZF output should decode");
    output
}

fn query_csi(source: &[u8], index_bytes: &[u8], region: &str) -> Vec<u8> {
    let index = csi::io::Reader::new(Cursor::new(index_bytes))
        .read_index()
        .expect("should read CSI index");
    let region = region.parse::<Region>().expect("should parse query region");
    let mut reader = csi::io::IndexedReader::new(Cursor::new(source), index);
    reader
        .query(&region)
        .expect("should query CSI")
        .map(|result| {
            result
                .expect("should read indexed annotation")
                .as_ref()
                .as_bytes()
                .to_vec()
        })
        .collect::<Vec<_>>()
        .concat()
}

#[test]
fn gff_preserves_headers_and_orders_coordinates_and_hierarchy() {
    let input = concat!(
        "##gff-version 3\n",
        "chr2\ts\tgene\t1\t2\t.\t+\t.\tID=chr2\n",
        "chr1\ts\texon\t10\t20\t.\t+\t.\tID=child;Parent=transcript%201,unrelated\n",
        "#source-order-comment\n",
        "chr1\ts\tgene\t10\t20\t.\t+\t.\tID=gene%201\n",
        "GLP1\ts\ttranscript\t10\t20\t.\t+\t.\tID=capped\n",
        "chr1\ts\tgene\t10\t19\t.\t+\t.\tID=shorter\n",
        "chr1\ts\tgene\t10\t20\t.\t+\t.\tID=unrelated\n",
        "chr1\ts\ttranscript\t10\t20\t.\t+\t.\tID=transcript%201;Parent=gene%201\n",
    );
    let mut output = Vec::new();
    sort_gff(Cursor::new(input.as_bytes()), false, &mut output).expect("GFF should sort");
    let output = lines(output);

    assert_eq!(
        output,
        [
            "##gff-version 3",
            "#source-order-comment",
            "chr1\ts\tgene\t10\t19\t.\t+\t.\tID=shorter",
            "chr1\ts\tgene\t10\t20\t.\t+\t.\tID=gene%201",
            "chr1\ts\tgene\t10\t20\t.\t+\t.\tID=unrelated",
            "chr1\ts\ttranscript\t10\t20\t.\t+\t.\tID=transcript%201;Parent=gene%201",
            "chr1\ts\texon\t10\t20\t.\t+\t.\tID=child;Parent=transcript%201,unrelated",
            "chr2\ts\tgene\t1\t2\t.\t+\t.\tID=chr2",
            "GLP1\ts\ttranscript\t10\t20\t.\t+\t.\tID=capped",
        ]
    );
}

#[test]
fn gff_rejects_parent_cycles() {
    let input = concat!(
        "chr1\ts\tgene\t1\t2\t.\t+\t.\tID=a;Parent=b\n",
        "chr1\ts\tgene\t1\t2\t.\t+\t.\tID=b;Parent=a\n",
    );
    let error =
        sort_gff(Cursor::new(input.as_bytes()), false, Vec::new()).expect_err("cycle should fail");
    assert!(
        error
            .to_string()
            .contains("parent hierarchy contains a cycle")
    );
}

#[test]
fn gff_parent_hierarchy_does_not_override_coordinate_order() {
    let input = concat!(
        "chr1\ts\texon\t10\t20\t.\t+\t.\tID=child;Parent=gene\n",
        "chr1\ts\tgene\t30\t40\t.\t+\t.\tID=gene\n",
    );
    let mut output = Vec::new();
    sort_gff(Cursor::new(input.as_bytes()), false, &mut output).expect("GFF should sort");
    assert_eq!(
        lines(output),
        [
            "chr1\ts\texon\t10\t20\t.\t+\t.\tID=child;Parent=gene",
            "chr1\ts\tgene\t30\t40\t.\t+\t.\tID=gene",
        ]
    );
}

#[test]
fn bed_sorts_by_contig_start_and_end_with_stable_ties() {
    let input = concat!(
        "chr2\t0\t5\ttwo\n",
        "chr1\t10\t30\tlong\n",
        "chr1\t2\t20\twide\n",
        "chr1\t2\t10\tfirst\n",
        "chr1\t2\t10\tsecond\n",
    );
    let mut output = Vec::new();
    sort_bed(Cursor::new(input.as_bytes()), false, &mut output).expect("BED should sort");
    assert_eq!(
        lines(output),
        [
            "chr1\t2\t10\tfirst",
            "chr1\t2\t10\tsecond",
            "chr1\t2\t20\twide",
            "chr1\t10\t30\tlong",
            "chr2\t0\t5\ttwo",
        ]
    );
}

#[test]
fn malformed_rows_are_actionable() {
    let gff_error = sort_gff(
        Cursor::new(b"chr1\ts\tgene\tbad\t2\t.\t+\t.\t."),
        false,
        Vec::new(),
    )
    .expect_err("invalid GFF coordinates should fail");
    assert!(gff_error.to_string().contains("GFF line 1"));
    assert!(gff_error.to_string().contains("invalid start"));

    let bed_error = sort_bed(Cursor::new(b"chr1\tnope\t3\n"), false, Vec::new())
        .expect_err("invalid BED coordinates should fail");
    assert!(bed_error.to_string().contains("BED line 1"));
    assert!(bed_error.to_string().contains("invalid start"));
}

#[test]
fn real_ecoli_bed_fixture_can_be_reordered_and_sorted() {
    let fixture =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/sorting/ecoli_k12_mg1655.bed.gz");
    let mut decoder =
        MultiGzDecoder::new(fs::File::open(&fixture).expect("real BED fixture exists"));
    let mut source = String::new();
    decoder
        .read_to_string(&mut source)
        .expect("real BED fixture should decompress");
    let mut fixture_lines = source.lines().collect::<Vec<_>>();
    assert!(fixture_lines.len() > 100);
    let last = fixture_lines.len() - 1;
    fixture_lines.swap(0, last);

    let directory = tempdir().expect("should create temporary directory");
    let input = directory.path().join("ecoli-reordered.bed");
    fs::write(&input, fixture_lines.join("\n") + "\n").expect("should write reordered BED");
    let mut output = Vec::new();
    sort_file(&input, false, &mut output).expect("fixture should sort");
    let output = lines(output);
    assert_eq!(output.len(), fixture_lines.len());
    assert_eq!(
        output
            .first()
            .unwrap()
            .split('\t')
            .take(3)
            .collect::<Vec<_>>(),
        ["NC_000913", "189", "255"]
    );
    assert!(output.windows(2).all(|pair| {
        let left = pair[0].split('\t').collect::<Vec<_>>();
        let right = pair[1].split('\t').collect::<Vec<_>>();
        (
            left[0],
            left[1].parse::<u64>().unwrap(),
            left[2].parse::<u64>().unwrap(),
        ) <= (
            right[0],
            right[1].parse::<u64>().unwrap(),
            right[2].parse::<u64>().unwrap(),
        )
    }));
}

#[test]
fn real_ecoli_gff_fixture_sorts_without_changing_record_count() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/sorting/ecoli_unsorted.gff");
    let source = fs::read_to_string(&fixture).expect("supplied GFF fixture should exist");
    let expected_records = source.lines().filter(|line| !line.starts_with('#')).count();
    let expected_comments = source.lines().filter(|line| line.starts_with('#')).count();

    let mut output = Vec::new();
    sort_file(&fixture, false, &mut output).expect("supplied GFF fixture should sort");
    let output = lines(output);
    assert_eq!(
        output.iter().filter(|line| !line.starts_with('#')).count(),
        expected_records
    );
    assert_eq!(
        output.iter().filter(|line| line.starts_with('#')).count(),
        expected_comments
    );

    let mut records = Vec::new();
    let mut saw_record = false;
    for line in &output {
        if line.starts_with('#') {
            assert!(!saw_record, "GFF comments must remain before records");
        } else {
            saw_record = true;
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(fields.len(), 9);
            records.push((
                fields[0],
                fields[3].parse::<u64>().expect("sorted GFF start"),
                fields[4].parse::<u64>().expect("sorted GFF end"),
            ));
        }
    }
    assert!(records.windows(2).all(|pair| pair[0] <= pair[1]));
}

#[test]
fn disk_sort_gff() {
    for (fixture_path, expected_path) in ["fixtures/random.gff.gz", "fixtures/random_fasta.gff.gz"]
        .iter()
        .zip([
            "fixtures/random_sorted.gff.gz",
            "fixtures/random_fasta_sorted.gff.gz",
        ])
    {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join(fixture_path);
        let mut output = Vec::new();
        sort_file(&fixture, true, &mut output).expect("supplied GFF fixture should sort");
        let output = lines(output);
        let mut reader = MultiGzDecoder::new(
            fs::File::open(Path::new(env!("CARGO_MANIFEST_DIR")).join(expected_path))
                .expect("should open file"),
        );
        let mut expected = Vec::new();
        reader.read_to_end(&mut expected).expect("should read");
        let expected = lines(expected);
        assert_eq!(
            expected
                .iter()
                .filter(|line| !line.starts_with('#'))
                .collect::<Vec<&String>>(),
            output
                .iter()
                .filter(|line| !line.starts_with('#'))
                .collect::<Vec<&String>>(),
            "Sorted using disk does not match"
        );
    }
}

#[test]
fn cli_sort_help_stdout_and_extension_errors() {
    let directory = tempdir().expect("should create temporary directory");
    let input = directory.path().join("input.GFF3");
    fs::write(
        &input,
        "#header\nchr2\ts\tgene\t1\t2\t.\t+\t.\t.\nchr1\ts\tgene\t1\t2\t.\t+\t.\t.\n",
    )
    .expect("should write input");
    let binary = env!("CARGO_BIN_EXE_gai");
    let sorted = Command::new(binary)
        .args(["sort", input.to_str().unwrap()])
        .output()
        .expect("should run gai sort");
    assert!(sorted.status.success(), "stderr: {:?}", sorted.stderr);
    assert_eq!(
        String::from_utf8(sorted.stdout).unwrap(),
        "#header\nchr1\ts\tgene\t1\t2\t.\t+\t.\t.\nchr2\ts\tgene\t1\t2\t.\t+\t.\t.\n"
    );

    let top_help = Command::new(binary)
        .arg("--help")
        .output()
        .expect("should show top-level help");
    assert!(top_help.status.success());
    assert!(String::from_utf8_lossy(&top_help.stdout).contains("sort"));

    let help = Command::new(binary)
        .args(["sort", "--help"])
        .output()
        .expect("should show sort help");
    assert!(help.status.success());
    let help = String::from_utf8_lossy(&help.stdout);
    assert!(help.contains("Input .gff, .gff3, .gtf, or .bed"));
    assert!(help.contains("--compress"));
    assert!(help.contains("--output"));
    assert!(help.contains("--coordinate-index"));
    assert!(help.contains("Write sorted output as BGZF"));

    let unsupported = directory.path().join("input.txt");
    fs::write(&unsupported, b"not an annotation\n").expect("should write unsupported input");
    let rejected = Command::new(binary)
        .args(["sort", unsupported.to_str().unwrap()])
        .output()
        .expect("should reject unsupported extension");
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("unsupported sort input extension"));
}

#[test]
fn cli_sort_compresses_gff_gtf_and_bed_output_as_bgzf() {
    let directory = tempdir().expect("should create temporary directory");
    let fasta_sequence = "A".repeat(70 * 1024);
    let long_annotation = "x".repeat(70 * 1024);
    let mut gff = String::from(
        "##gff-version 3\n\
#source-order-comment\n\
chr1\ts\texon\t10\t20\t.\t+\t.\tID=exon1;Parent=transcript1\n\
chr1\ts\tgene\t10\t20\t.\t+\t.\tID=gene1;Name=gene-one\n\
chr1\ts\ttranscript\t10\t20\t.\t+\t.\tID=transcript1;Parent=gene1\n\
",
    );
    gff.push_str(&format!(
        "chr2\ts\tgene\t300000000\t300000010\t.\t+\t.\tID=long;Name=longmatch;Note={long_annotation}\n"
    ));
    gff.push_str(
        "chr2\ts\tgene\t600000000\t600000010\t.\t+\t.\tID=after;Name=afterblock\n##FASTA\n>chr1\n",
    );
    gff.push_str(&fasta_sequence);
    gff.push('\n');

    // The .gtf path uses the sorter's existing GFF attribute parser.
    let gtf = concat!(
        "#gtf-version 2.2\n",
        "chr2\ts\texon\t1\t5\t.\t+\t.\tgene_id=g2;transcript_id=t2\n",
        "chr1\ts\texon\t2\t8\t.\t+\t.\tgene_id=g1;transcript_id=t1\n",
    );
    let bed = concat!(
        "chr2\t0\t5\ttwo\n",
        "chr1\t10\t30\tlong\n",
        "chr1\t2\t20\twide\n",
        "chr1\t2\t10\tfirst\n",
        "chr1\t2\t10\tsecond\n",
        "chr3\t2\t2\tpoint\n",
    );
    let cases = [
        ("stream.gff3", gff.into_bytes()),
        ("stream.gtf", gtf.as_bytes().to_vec()),
        ("stream.bed", bed.as_bytes().to_vec()),
    ];
    let binary = env!("CARGO_BIN_EXE_gai");

    for (file_name, input_bytes) in cases {
        let input = directory.path().join(file_name);
        fs::write(&input, input_bytes).expect("should write sort input");
        let run_sort = || {
            Command::new(binary)
                .args(["sort"])
                .arg(&input)
                .output()
                .expect("should run gai sort")
        };

        let plain = run_sort();
        assert!(plain.status.success(), "stderr: {:?}", plain.stderr);
        assert!(!plain.stdout.starts_with(&[0x1f, 0x8b]));

        let output_path = directory.path().join(format!("{file_name}.bgzf"));
        let explicit_index = file_name == "stream.bed";
        let index_path = if explicit_index {
            let index_directory = directory.path().join("indexes");
            fs::create_dir_all(&index_directory).expect("should create custom index directory");
            index_directory.join("stream.bed.csi")
        } else {
            PathBuf::from(format!("{}.csi", output_path.display()))
        };
        let mut compressed_command = Command::new(binary);
        compressed_command
            .args(["sort"])
            .arg(&input)
            .args(["--compress", "--output"])
            .arg(&output_path);
        if explicit_index {
            compressed_command
                .arg("--coordinate-index")
                .arg(&index_path);
        }
        let compressed = compressed_command
            .output()
            .expect("should run gai sort --compress");
        assert!(
            compressed.status.success(),
            "stderr: {:?}",
            compressed.stderr
        );
        assert!(compressed.stdout.is_empty());
        let compressed_bytes = fs::read(&output_path).expect("should read BGZF output");
        let index_bytes = fs::read(&index_path).expect("should read generated CSI");
        assert_eq!(decode_bgzf(&compressed_bytes), plain.stdout);

        let (region, expected_record) = match file_name {
            "stream.gff3" => ("chr1:10-20", "ID=gene1"),
            "stream.gtf" => ("chr1:2-8", "gene_id=g1"),
            "stream.bed" => ("chr1:11-30", "long"),
            _ => unreachable!(),
        };
        let queried_records = query_csi(&compressed_bytes, &index_bytes, region);
        assert!(
            String::from_utf8_lossy(&queried_records).contains(expected_record),
            "CSI query {region} returned {:?}",
            String::from_utf8_lossy(&queried_records)
        );

        if file_name == "stream.gff3" {
            assert!(plain.stdout.len() > 64 * 1024);
            let expected_prefix = concat!(
                "##gff-version 3\n",
                "#source-order-comment\n",
                "chr1\ts\tgene\t10\t20\t.\t+\t.\tID=gene1;Name=gene-one\n",
                "chr1\ts\ttranscript\t10\t20\t.\t+\t.\tID=transcript1;Parent=gene1\n",
                "chr1\ts\texon\t10\t20\t.\t+\t.\tID=exon1;Parent=transcript1\n",
                "chr2\ts\tgene\t300000000\t300000010\t.\t+\t.\tID=long;Name=longmatch;Note=",
            );
            assert!(plain.stdout.starts_with(expected_prefix.as_bytes()));
            let expected_fasta_tail = format!("##FASTA\n>chr1\n{fasta_sequence}\n");
            assert!(plain.stdout.ends_with(expected_fasta_tail.as_bytes()));

            let bgzf_eof = bgzf::io::Writer::new(Vec::new())
                .finish()
                .expect("should create a BGZF EOF marker");
            assert!(compressed_bytes.ends_with(&bgzf_eof));
            assert!(
                String::from_utf8_lossy(&query_csi(
                    &compressed_bytes,
                    &index_bytes,
                    "chr2:300000000-300000010"
                ))
                .contains("Name=longmatch")
            );
            assert!(
                String::from_utf8_lossy(&query_csi(
                    &compressed_bytes,
                    &index_bytes,
                    "chr2:600000000-600000010"
                ))
                .contains("Name=afterblock")
            );

            let disk_compressed_index = directory.path().join("disk-sorted.csi");
            let disk_compressed = Command::new(binary)
                .args(["sort"])
                .arg(&input)
                .args(["--disk-sort", "--compress", "--coordinate-index"])
                .arg(&disk_compressed_index)
                .output()
                .expect("should sort compressed stdout and CSI with --disk-sort");
            assert!(
                disk_compressed.status.success(),
                "stderr: {:?}",
                disk_compressed.stderr
            );
            assert_eq!(decode_bgzf(&disk_compressed.stdout), plain.stdout);
            let stdout_index = fs::read(&disk_compressed_index).expect("should read stdout CSI");
            assert!(
                String::from_utf8_lossy(&query_csi(
                    &disk_compressed.stdout,
                    &stdout_index,
                    "chr2:600000000-600000010"
                ))
                .contains("Name=afterblock")
            );

            let attribute_index_path = directory.path().join("stream.gff3.gai");
            let name_options =
                NameIndexOptions::new(["Name"], false).expect("should configure GFF Name indexing");
            build_name_index_with_options(
                &output_path,
                &index_path,
                &attribute_index_path,
                &name_options,
                &BuildOptions::default(),
            )
            .expect("should build GAI using sort-generated CSI");
            let matching_records = query_index(
                &output_path,
                &index_path,
                &attribute_index_path,
                "longmatch",
            )
            .expect("should query the generated name and coordinate indexes");
            assert_eq!(matching_records.len(), 1);
            assert_eq!(matching_records[0].start, 300_000_000);

            if cfg!(feature = "profiling") {
                let profile_index = directory.path().join("profile-sorted.csi");
                let profiled = Command::new(binary)
                    .args(["profile", "sort"])
                    .arg(&input)
                    .args(["--disk-sort", "--compress", "--coordinate-index"])
                    .arg(&profile_index)
                    .output()
                    .expect("should run profiled sort with CSI output");
                assert!(profiled.status.success(), "stderr: {:?}", profiled.stderr);
                assert_eq!(decode_bgzf(&profiled.stdout), plain.stdout);
            }
        } else if file_name == "stream.bed" {
            let index = csi::io::Reader::new(Cursor::new(&index_bytes))
                .read_index()
                .expect("should read BED CSI");
            let point_region = "chr3:3-3"
                .parse::<Region>()
                .expect("should parse BED point region");
            assert!(
                !index
                    .query(2, point_region.interval())
                    .expect("should query BED point candidate bins")
                    .is_empty(),
                "zero-width BED intervals should be indexed as point candidates"
            );

            let stdout_index_path = directory.path().join("stdout-sorted.bed.csi");
            let stdout_compressed = Command::new(binary)
                .args(["sort"])
                .arg(&input)
                .args(["--disk-sort", "--compress", "--coordinate-index"])
                .arg(&stdout_index_path)
                .output()
                .expect("should sort BED to compressed stdout and CSI");
            assert!(
                stdout_compressed.status.success(),
                "stderr: {:?}",
                stdout_compressed.stderr
            );
            assert_eq!(decode_bgzf(&stdout_compressed.stdout), plain.stdout);
            let stdout_index_bytes =
                fs::read(&stdout_index_path).expect("should read stdout BED CSI");
            let stdout_records =
                query_csi(&stdout_compressed.stdout, &stdout_index_bytes, "chr1:11-30");
            assert!(String::from_utf8_lossy(&stdout_records).contains("long"));
        }
    }
}

#[test]
fn cli_sort_compress_requires_destination_and_preserves_final_paths_on_failure() {
    let directory = tempdir().expect("should create temporary directory");
    let binary = env!("CARGO_BIN_EXE_gai");
    let valid_input = directory.path().join("valid.gff3");
    fs::write(&valid_input, b"chr1\ts\tgene\t1\t2\t.\t+\t.\tName=ok\n")
        .expect("should write valid GFF input");

    let missing_destination = Command::new(binary)
        .args(["sort"])
        .arg(&valid_input)
        .arg("--compress")
        .output()
        .expect("should reject compressed stdout without a CSI destination");
    assert!(!missing_destination.status.success());
    assert!(missing_destination.stdout.is_empty());
    assert!(
        String::from_utf8_lossy(&missing_destination.stderr)
            .contains("--compress requires --output or --coordinate-index")
    );

    let invalid_input = directory.path().join("invalid.gff3");
    fs::write(&invalid_input, b"not a GFF record\n").expect("should write invalid input");
    let failed_output = directory.path().join("failed.gff3.bgzf");
    let failed = Command::new(binary)
        .args(["sort"])
        .arg(&invalid_input)
        .args(["--compress", "--output"])
        .arg(&failed_output)
        .output()
        .expect("should reject malformed GFF");
    assert!(!failed.status.success());
    assert!(!failed_output.exists());
    assert!(!PathBuf::from(format!("{}.csi", failed_output.display())).exists());

    let existing_output = directory.path().join("existing.gff3.bgzf");
    let existing_index = PathBuf::from(format!("{}.csi", existing_output.display()));
    fs::write(&existing_output, b"preserve compressed output")
        .expect("should write existing output sentinel");
    let rejected_output = Command::new(binary)
        .args(["sort"])
        .arg(&valid_input)
        .args(["--compress", "--output"])
        .arg(&existing_output)
        .output()
        .expect("should reject existing output destination");
    assert!(!rejected_output.status.success());
    assert_eq!(
        fs::read(&existing_output).unwrap(),
        b"preserve compressed output"
    );
    assert!(!existing_index.exists());

    let output_for_existing_index = directory.path().join("index-exists.gff3.bgzf");
    let existing_index_path = directory.path().join("custom.csi");
    fs::write(&existing_index_path, b"preserve CSI").expect("should write CSI sentinel");
    let rejected_index = Command::new(binary)
        .args(["sort"])
        .arg(&valid_input)
        .args(["--compress", "--output"])
        .arg(&output_for_existing_index)
        .arg("--coordinate-index")
        .arg(&existing_index_path)
        .output()
        .expect("should reject existing CSI destination");
    assert!(!rejected_index.status.success());
    assert!(!output_for_existing_index.exists());
    assert_eq!(fs::read(&existing_index_path).unwrap(), b"preserve CSI");

    let rejected_stdout_index = Command::new(binary)
        .args(["sort"])
        .arg(&valid_input)
        .args(["--compress", "--coordinate-index"])
        .arg(&existing_index_path)
        .output()
        .expect("should reject existing stdout CSI destination");
    assert!(!rejected_stdout_index.status.success());
    assert!(rejected_stdout_index.stdout.is_empty());
    assert_eq!(fs::read(&existing_index_path).unwrap(), b"preserve CSI");

    let empty_input = directory.path().join("empty.gff3");
    let empty_output = directory.path().join("empty.gff3.bgzf");
    fs::write(&empty_input, b"").expect("should write empty input");
    let empty = Command::new(binary)
        .args(["sort"])
        .arg(&empty_input)
        .args(["--compress", "--output"])
        .arg(&empty_output)
        .output()
        .expect("should compress empty source");
    assert!(empty.status.success(), "stderr: {:?}", empty.stderr);
    assert!(decode_bgzf(&fs::read(&empty_output).unwrap()).is_empty());
    assert!(PathBuf::from(format!("{}.csi", empty_output.display())).exists());
}
