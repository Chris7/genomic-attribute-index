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
    let destination = directory.path().join("cli.gni");
    let binary = env!("CARGO_BIN_EXE_gen");
    let source = source.to_string_lossy().into_owned();
    let coordinate_index = coordinate_index.to_string_lossy().into_owned();
    let destination_string = destination.to_string_lossy().into_owned();

    let indexed = Command::new(binary)
        .args([
            "gff",
            "index-names",
            &source,
            "--attribute",
            "Name",
            "--coordinate-index",
            &coordinate_index,
            "--output",
            &destination_string,
        ])
        .output()
        .expect("should run gen index-names");
    assert!(indexed.status.success(), "stderr: {:?}", indexed.stderr);
    assert!(destination.exists());
    assert!(String::from_utf8_lossy(&indexed.stdout).contains("wrote"));
    let indexing_stderr = String::from_utf8_lossy(&indexed.stderr);
    assert!(indexing_stderr.contains("Scan"));
    assert!(indexing_stderr.contains("Complete"));

    let queried = Command::new(binary)
        .args([
            "gff",
            "query-name",
            &source,
            "BRCA1",
            "--coordinate-index",
            &coordinate_index,
            "--gni",
            &destination_string,
        ])
        .output()
        .expect("should run gen query-name");
    assert!(queried.status.success(), "stderr: {:?}", queried.stderr);
    assert_eq!(
        String::from_utf8_lossy(&queried.stdout).trim(),
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1;Alias=BRCC1"
    );

    let inspected = Command::new(binary)
        .args(["gff", "inspect-name-index", &destination_string])
        .output()
        .expect("should run gen inspect-name-index");
    assert!(inspected.status.success(), "stderr: {:?}", inspected.stderr);
    let inspection = String::from_utf8_lossy(&inspected.stdout);
    assert!(inspection.contains("zero-based half-open"));
    assert!(inspection.contains("FST bytes:"));
    assert!(inspection.contains("postings compression"));
    assert!(inspection.contains("span compression"));

    let missing_attribute = Command::new(binary)
        .args([
            "gff",
            "index-names",
            &source,
            "--coordinate-index",
            &coordinate_index,
        ])
        .output()
        .expect("should run invalid gen command");
    assert!(!missing_attribute.status.success());
    assert!(String::from_utf8_lossy(&missing_attribute.stderr).contains("attribute"));
}
