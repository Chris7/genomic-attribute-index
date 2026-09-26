use std::{
    fs,
    io::{Cursor, Read},
    path::Path,
    process::Command,
};

use flate2::read::MultiGzDecoder;
use gai::{sort_bed, sort_file, sort_gff};
use tempfile::tempdir;

fn lines(output: Vec<u8>) -> Vec<String> {
    String::from_utf8(output)
        .expect("sort output should be UTF-8")
        .lines()
        .map(str::to_owned)
        .collect()
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
    assert!(String::from_utf8_lossy(&help.stdout).contains("Input .gff, .gff3, or .bed"));

    let unsupported = directory.path().join("input.txt");
    fs::write(&unsupported, b"not an annotation\n").expect("should write unsupported input");
    let rejected = Command::new(binary)
        .args(["sort", unsupported.to_str().unwrap()])
        .output()
        .expect("should reject unsupported extension");
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("unsupported sort input extension"));
}
