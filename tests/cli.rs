use std::{fs::File, io::Write, path::Path, process::Command};

use noodles::{
    bgzf, core::Position, csi::binning_index::index::reference_sequence::bin::Chunk, tabix,
};
use tempfile::tempdir;

fn write_fixture(directory: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let source_path = directory.join("cli.gff3.gz");
    let index_path = directory.join("cli.gff3.gz.tbi");
    let lines = [
        "##gff-version 3",
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1;Alias=BRCC1",
        "chr1\tsrc\tgene\t100\t110\t.\t-\t.\tName=Other%20Gene",
    ];
    let mut writer = File::create(&source_path)
        .map(bgzf::io::Writer::new)
        .expect("should create BGZF source");
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(noodles::csi::binning_index::index::header::Builder::gff().build());
    for line in lines {
        if line.starts_with('#') {
            writeln!(writer, "{line}").expect("should write directive");
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let start = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
        let end = Position::try_from(fields[4].parse::<usize>().unwrap()).unwrap();
        let start_position = writer.virtual_position();
        writeln!(writer, "{line}").expect("should write record");
        let end_position = writer.virtual_position();
        indexer
            .add_record(
                fields[0],
                start,
                end,
                Chunk::new(start_position, end_position),
            )
            .expect("should index record");
    }
    writer.finish().expect("should finish BGZF");
    let index = indexer.build();
    let mut index_writer = File::create(&index_path)
        .map(tabix::io::Writer::new)
        .expect("should create TBI");
    index_writer.write_index(&index).expect("should write TBI");
    drop(index_writer);
    (source_path, index_path)
}

#[test]
fn cli_index_query_and_inspect() {
    let directory = tempdir().expect("should create temp directory");
    let (source, coordinate_index) = write_fixture(directory.path());
    let destination = directory.path().join("cli.gai");
    let binary = env!("CARGO_BIN_EXE_gai");
    let source = source.to_string_lossy().into_owned();
    let coordinate_index = coordinate_index.to_string_lossy().into_owned();
    let destination_string = destination.to_string_lossy().into_owned();

    let indexed = Command::new(binary)
        .args([
            "build-index",
            &source,
            "--attribute",
            "Name",
            "--coordinate-index",
            &coordinate_index,
            "--output",
            &destination_string,
        ])
        .output()
        .expect("should run gai build-index");
    assert!(indexed.status.success(), "stderr: {:?}", indexed.stderr);
    assert!(destination.exists());
    assert!(String::from_utf8_lossy(&indexed.stdout).contains("wrote"));
    let indexing_stderr = String::from_utf8_lossy(&indexed.stderr);
    assert!(indexing_stderr.contains("Scan"));
    assert!(indexing_stderr.contains("Complete"));

    let queried = Command::new(binary)
        .args([
            "query-index",
            &source,
            "BRCA1",
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &destination_string,
        ])
        .output()
        .expect("should run gai query-index");
    assert!(queried.status.success(), "stderr: {:?}", queried.stderr);
    assert_eq!(
        String::from_utf8_lossy(&queried.stdout).trim(),
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1;Alias=BRCC1"
    );

    let prefix = Command::new(binary)
        .args([
            "query-index",
            &source,
            "BRCA",
            "--match",
            "prefix",
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &destination_string,
        ])
        .output()
        .expect("should run prefix query");
    assert!(prefix.status.success(), "stderr: {:?}", prefix.stderr);
    assert_eq!(
        String::from_utf8_lossy(&prefix.stdout).trim(),
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1;Alias=BRCC1"
    );

    let invalid_match = Command::new(binary)
        .args([
            "query-index",
            &source,
            "BRCA",
            "--match",
            "substring",
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &destination_string,
        ])
        .output()
        .expect("should reject invalid match mode");
    assert!(!invalid_match.status.success());
    assert!(String::from_utf8_lossy(&invalid_match.stderr).contains("exact"));

    let inspected = Command::new(binary)
        .args(["inspect-index", &destination_string])
        .output()
        .expect("should run gai inspect-index");
    assert!(inspected.status.success(), "stderr: {:?}", inspected.stderr);
    let inspection = String::from_utf8_lossy(&inspected.stdout);
    assert!(inspection.contains("zero-based half-open"));
    assert!(inspection.contains("FST bytes:"));
    assert!(inspection.contains("postings compression"));
    assert!(inspection.contains("starts compression"));
    assert!(inspection.contains("lengths compression"));

    let missing_attribute = Command::new(binary)
        .args([
            "build-index",
            &source,
            "--coordinate-index",
            &coordinate_index,
        ])
        .output()
        .expect("should run invalid build command");
    assert!(!missing_attribute.status.success());
    assert!(String::from_utf8_lossy(&missing_attribute.stderr).contains("attribute"));

    let help = Command::new(binary)
        .arg("--help")
        .output()
        .expect("should run gai help");
    assert!(help.status.success());
    let help = String::from_utf8_lossy(&help.stdout);
    assert!(help.contains("Usage: gai"));
    assert!(help.contains("build-index"));
    assert!(help.contains("query-index"));
    assert!(help.contains("inspect-index"));
    assert!(!help.contains("gff"));

    let old_surface = Command::new(binary)
        .args(["gff", "query-name", &source, "BRCA1"])
        .output()
        .expect("should run removed command check");
    assert!(!old_surface.status.success());
}
