use std::{fs::File, io::Write, path::Path, process::Command};

use gai::{
    MatchMode, NameIndexOptions, build_name_index, inspect_index, open_index, query_index,
    query_index_with_mode,
};
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

#[derive(Clone, Copy)]
enum FixtureFormat {
    Gff,
    Bed,
}

fn write_multi_contig_fixture(
    directory: &Path,
    source_name: &str,
    format: FixtureFormat,
    lines: &[&str],
) -> (std::path::PathBuf, std::path::PathBuf) {
    let source_path = directory.join(source_name);
    let index_path = std::path::PathBuf::from(format!("{}.tbi", source_path.display()));
    let mut writer = File::create(&source_path)
        .map(bgzf::io::Writer::new)
        .expect("should create BGZF source");
    let mut indexer = tabix::index::Indexer::default();
    let is_gff = matches!(format, FixtureFormat::Gff);
    let header = if is_gff {
        noodles::csi::binning_index::index::header::Builder::gff().build()
    } else {
        noodles::csi::binning_index::index::header::Builder::bed().build()
    };
    indexer.set_header(header);

    for line in lines {
        if line.starts_with('#') {
            writeln!(writer, "{line}").expect("should write directive");
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let (reference, start, end) = if is_gff {
            (
                fields[0],
                fields[3].parse::<usize>().unwrap(),
                fields[4].parse::<usize>().unwrap(),
            )
        } else {
            let start = fields[1].parse::<usize>().unwrap();
            (fields[0], start + 1, fields[2].parse::<usize>().unwrap())
        };
        let start = Position::try_from(start).unwrap();
        let end = Position::try_from(end).unwrap();
        let start_position = writer.virtual_position();
        writeln!(writer, "{line}").expect("should write record");
        let end_position = writer.virtual_position();
        indexer
            .add_record(
                reference,
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
    let alias_destination = directory.path().join("cli-alias.gai");
    let alias_destination_string = alias_destination.to_string_lossy().into_owned();

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

    let alias_built = Command::new(binary)
        .args([
            "build",
            &source,
            "--attribute",
            "Name",
            "--coordinate-index",
            &coordinate_index,
            "--output",
            &alias_destination_string,
        ])
        .output()
        .expect("should run gai build alias");
    assert!(
        alias_built.status.success(),
        "stderr: {:?}",
        alias_built.stderr
    );
    assert!(alias_destination.exists());

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

    let alias_queried = Command::new(binary)
        .args([
            "query",
            &source,
            "BRCA1",
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &alias_destination_string,
        ])
        .output()
        .expect("should run gai query alias");
    assert!(
        alias_queried.status.success(),
        "stderr: {:?}",
        alias_queried.stderr
    );
    assert_eq!(alias_queried.stdout, queried.stdout);

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

    let contains = Command::new(binary)
        .args([
            "query-index",
            &source,
            "RCA",
            "--match",
            "contains",
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &destination_string,
        ])
        .output()
        .expect("should run contains query");
    assert!(contains.status.success(), "stderr: {:?}", contains.stderr);
    assert_eq!(
        String::from_utf8_lossy(&contains.stdout).trim(),
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1;Alias=BRCC1"
    );

    let regex = Command::new(binary)
        .args([
            "query-index",
            &source,
            "^(BRCA[0-9]|Other Gene)$",
            "--match",
            "regex",
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &destination_string,
        ])
        .output()
        .expect("should run regex query");
    assert!(regex.status.success(), "stderr: {:?}", regex.stderr);
    assert_eq!(
        String::from_utf8_lossy(&regex.stdout).trim(),
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1;Alias=BRCC1\nchr1\tsrc\tgene\t100\t110\t.\t-\t.\tName=Other%20Gene"
    );

    let uppercase_escape = Command::new(binary)
        .args([
            "query-index",
            &source,
            r"\S",
            "--match",
            "regex",
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &destination_string,
        ])
        .output()
        .expect("should preserve uppercase regex escapes");
    assert!(
        uppercase_escape.status.success(),
        "stderr: {:?}",
        uppercase_escape.stderr
    );
    assert_eq!(
        String::from_utf8_lossy(&uppercase_escape.stdout)
            .lines()
            .count(),
        2
    );

    let invalid_regex = Command::new(binary)
        .args([
            "query-index",
            &source,
            "no-match[",
            "--match",
            "regex",
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &destination_string,
        ])
        .output()
        .expect("should reject invalid regex even with no matching term");
    assert!(!invalid_regex.status.success());
    assert!(String::from_utf8_lossy(&invalid_regex.stderr).contains("invalid regex"));

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
    assert!(!String::from_utf8_lossy(&inspected.stderr).contains("time.busy"));

    let alias_inspected = Command::new(binary)
        .args(["inspect", &alias_destination_string])
        .output()
        .expect("should run gai inspect alias");
    assert!(
        alias_inspected.status.success(),
        "stderr: {:?}",
        alias_inspected.stderr
    );
    assert_eq!(alias_inspected.stdout, inspected.stdout);

    if cfg!(feature = "profiling") {
        let profiled = Command::new(binary)
            .args(["profile", "inspect-index", &destination_string])
            .output()
            .expect("should run profiled inspect-index");
        assert!(profiled.status.success(), "stderr: {:?}", profiled.stderr);
        assert_eq!(profiled.stdout, inspected.stdout);
        let profiled_alias = Command::new(binary)
            .args([
                "profile",
                "query",
                &source,
                "BRCA1",
                "--coordinate-index",
                &coordinate_index,
                "--gai",
                &alias_destination_string,
            ])
            .output()
            .expect("should run profiled query alias");
        assert!(
            profiled_alias.status.success(),
            "stderr: {:?}",
            profiled_alias.stderr
        );
        assert_eq!(profiled_alias.stdout, queried.stdout);
        let profiling_stderr = String::from_utf8_lossy(&profiled.stderr);
        assert!(profiling_stderr.contains("Profile results"));
        assert!(profiling_stderr.contains("Total (ms)"));
        assert!(profiling_stderr.contains("Calls"));
        assert!(!profiling_stderr.contains("time.busy"));
        assert!(profiling_stderr.lines().any(|line| {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            fields.first() == Some(&"gai::compression_ratio") && fields.get(1) == Some(&"3")
        }));
    }

    if cfg!(all(unix, feature = "profiling")) {
        let sampled = Command::new(binary)
            .args(["profile", "--sample", "inspect-index", &destination_string])
            .output()
            .expect("should run sample-profiled inspect-index");
        assert!(sampled.status.success(), "stderr: {:?}", sampled.stderr);
        assert_eq!(sampled.stdout, inspected.stdout);
        let sample_stderr = String::from_utf8_lossy(&sampled.stderr);
        assert!(sample_stderr.contains("Sampling profile results"));
        assert!(sample_stderr.contains("Samples"));
        assert!(sample_stderr.contains("Time (ms)"));
        assert!(sample_stderr.contains("Pct"));
        assert!(!sample_stderr.contains("Report {"));
        assert!(!sample_stderr.contains("time.busy"));
    }

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
    assert!(help.contains("[alias: build]"));
    assert!(help.contains("[alias: query]"));
    assert!(help.contains("[alias: inspect]"));
    assert_eq!(
        help.contains("profile"),
        cfg!(feature = "profiling"),
        "profile command should only be included with the profiling feature"
    );
    assert!(!help.contains("gff"));

    let old_surface = Command::new(binary)
        .args(["gff", "query-name", &source, "BRCA1"])
        .output()
        .expect("should run removed command check");
    assert!(!old_surface.status.success());
}

#[test]
fn cli_query_filters_exact_repeatable_contig_names_for_gff_and_bed() {
    let directory = tempdir().expect("should create temp directory");
    let binary = env!("CARGO_BIN_EXE_gai");
    let gff_records = [
        "##gff-version 3",
        "chr2\tsrc\tgene\t10\t20\t.\t+\t.\tName=Shared",
        "chr1\tsrc\tgene\t30\t40\t.\t+\t.\tName=Shared",
        "chr1\tsrc\tgene\t70\t80\t.\t+\t.\tName=Shared",
        "chr1\tsrc\tgene\t90\t100\t.\t+\t.\tName=SharedLong",
        "chr3\tsrc\tgene\t110\t120\t.\t+\t.\tName=Shared",
    ];
    let bed_records = [
        "chr2\t0\t10\tShared",
        "chr1\t20\t30\tShared",
        "chr1\t60\t70\tShared",
        "chr1\t80\t90\tSharedLong",
        "chr3\t100\t110\tShared",
    ];
    let cases = [
        (FixtureFormat::Gff, "multi-contig.gff3.gz", &gff_records[..]),
        (FixtureFormat::Bed, "multi-contig.bed.gz", &bed_records[..]),
    ];

    let query_help = Command::new(binary)
        .args(["query-index", "--help"])
        .output()
        .expect("should show query-index help");
    assert!(query_help.status.success());
    let query_help = String::from_utf8_lossy(&query_help.stdout);
    assert!(query_help.contains("--contig <CONTIG>"));
    assert!(query_help.contains("case-insensitive contig name"));
    assert!(query_help.contains("may be repeated"));

    for (format, source_name, records) in cases {
        let (source, coordinate_index) =
            write_multi_contig_fixture(directory.path(), source_name, format, records);
        let destination = directory.path().join(format!("{source_name}.gai"));
        let source = source.to_string_lossy().into_owned();
        let coordinate_index = coordinate_index.to_string_lossy().into_owned();
        let destination = destination.to_string_lossy().into_owned();
        let attribute = if matches!(format, FixtureFormat::Gff) {
            "Name"
        } else {
            "ignored"
        };
        let indexed = Command::new(binary)
            .args([
                "build-index",
                &source,
                "--attribute",
                attribute,
                "--coordinate-index",
                &coordinate_index,
                "--output",
                &destination,
            ])
            .output()
            .expect("should build multi-contig GAI");
        assert!(indexed.status.success(), "stderr: {:?}", indexed.stderr);

        let query =
            |prefix: &[&str], command: &str, term: &str, mode: Option<&str>, contigs: &[&str]| {
                let mut arguments = prefix
                    .iter()
                    .map(|argument| (*argument).to_owned())
                    .collect::<Vec<_>>();
                arguments.extend([
                    command.to_owned(),
                    source.clone(),
                    term.to_owned(),
                    "--coordinate-index".to_owned(),
                    coordinate_index.clone(),
                    "--gai".to_owned(),
                    destination.clone(),
                ]);
                if let Some(mode) = mode {
                    arguments.push("--match".to_owned());
                    arguments.push(mode.to_owned());
                }
                for contig in contigs {
                    arguments.push("--contig".to_owned());
                    arguments.push((*contig).to_owned());
                }
                Command::new(binary)
                    .args(arguments)
                    .output()
                    .expect("should query multi-contig index")
            };
        let assert_records = |output: &std::process::Output, expected: &[&str]| {
            assert!(output.status.success(), "stderr: {:?}", output.stderr);
            let actual = String::from_utf8_lossy(&output.stdout)
                .lines()
                .map(str::to_owned)
                .collect::<Vec<_>>();
            let expected = expected
                .iter()
                .map(|record| (*record).to_owned())
                .collect::<Vec<_>>();
            assert_eq!(actual, expected);
        };
        let expected_all = records
            .iter()
            .filter(|record| record.ends_with("=Shared") || record.ends_with("\tShared"))
            .copied()
            .collect::<Vec<_>>();
        let expected_chr1 = expected_all
            .iter()
            .filter(|record| record.starts_with("chr1\t"))
            .copied()
            .collect::<Vec<_>>();
        let expected_chr1_prefix = records
            .iter()
            .filter(|record| record.starts_with("chr1\t"))
            .copied()
            .collect::<Vec<_>>();
        let expected_multiple_contigs = expected_all
            .iter()
            .filter(|record| record.starts_with("chr1\t") || record.starts_with("chr2\t"))
            .copied()
            .collect::<Vec<_>>();

        let mut indexed_api =
            open_index(&source, &coordinate_index, &destination).expect("should open GAI API");
        let exact_api = indexed_api
            .query_on_contigs("SHARED", &["CHR1", "cHr2"])
            .expect("should run exact default contig query");
        assert_eq!(
            exact_api
                .iter()
                .map(|record| record.raw_line.as_str())
                .collect::<Vec<_>>(),
            expected_multiple_contigs
        );

        let prefix_api = indexed_api
            .query_on_contigs_with_mode("SHA", &["cHr1"], MatchMode::Prefix)
            .expect("should run prefix contig query");
        assert_eq!(
            prefix_api
                .iter()
                .map(|record| record.raw_line.as_str())
                .collect::<Vec<_>>(),
            expected_chr1_prefix
        );

        let (exact_with_stats, exact_stats) = indexed_api
            .query_on_contigs_with_stats("SHARED", &["CHR1", "cHr2"])
            .expect("should run exact default contig query with stats");
        assert_eq!(
            exact_with_stats
                .iter()
                .map(|record| record.raw_line.as_str())
                .collect::<Vec<_>>(),
            expected_multiple_contigs
        );
        assert_eq!(
            exact_stats.requested_spans,
            expected_multiple_contigs.len() as u64
        );
        assert_eq!(
            exact_stats.matching_records,
            expected_multiple_contigs.len() as u64
        );

        let (prefix_with_stats, prefix_stats) = indexed_api
            .query_on_contigs_with_mode_and_stats("SHA", &["cHr1"], MatchMode::Prefix)
            .expect("should run prefix contig query with stats");
        assert_eq!(
            prefix_with_stats
                .iter()
                .map(|record| record.raw_line.as_str())
                .collect::<Vec<_>>(),
            expected_chr1_prefix
        );
        assert_eq!(
            prefix_stats.requested_spans,
            expected_chr1_prefix.len() as u64
        );
        assert_eq!(
            prefix_stats.matching_records,
            expected_chr1_prefix.len() as u64
        );

        let all = query(&[], "query-index", "SHARED", None, &[]);
        assert_records(&all, &expected_all);

        let chr1 = query(&[], "query", "SHARED", None, &["cHr1"]);
        assert_records(&chr1, &expected_chr1);

        let repeated = query(
            &[],
            "query-index",
            "SHARED",
            None,
            &["CHR1", "cHr2", "cHr1"],
        );
        assert_records(&repeated, &expected_multiple_contigs);

        let unknown = query(&[], "query-index", "SHARED", None, &["missing"]);
        assert_records(&unknown, &[]);

        let mixed_case = query(&[], "query-index", "SHARED", None, &["cHR1"]);
        assert_records(&mixed_case, &expected_chr1);

        let prefix = query(&[], "query-index", "SHA", Some("prefix"), &["cHr1"]);
        assert_records(&prefix, &expected_chr1_prefix);

        if cfg!(feature = "profiling") {
            let profiled = query(&["profile"], "query", "SHARED", None, &["chr2"]);
            let expected_chr2 = expected_all
                .iter()
                .filter(|record| record.starts_with("chr2\t"))
                .copied()
                .collect::<Vec<_>>();
            assert_records(&profiled, &expected_chr2);
        }
    }
}

#[test]
fn rust_querying_api_gff_helpers_and_reusable_source() {
    let directory = tempdir().expect("should create temp directory");
    let (source, coordinate_index) = write_fixture(directory.path());
    let destination = directory.path().join("api.gai");
    let options = NameIndexOptions::new(["Name"], false).expect("should configure Name");
    build_name_index(&source, &coordinate_index, &destination, &options)
        .expect("should build GAI fixture");

    let metadata = inspect_index(&destination).expect("should inspect GAI");
    assert_eq!(metadata.attributes, vec!["Name".to_string()]);
    assert_eq!(metadata.term_count, 2);

    let exact = query_index(&source, &coordinate_index, &destination, "BRCA1")
        .expect("one-shot exact query should work");
    assert_eq!(exact.len(), 1);
    assert_eq!(
        exact[0].raw_line,
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1;Alias=BRCC1"
    );

    let prefix = query_index_with_mode(
        &source,
        &coordinate_index,
        &destination,
        "BRCA",
        MatchMode::Prefix,
    )
    .expect("one-shot prefix query should work");
    assert_eq!(prefix, exact);

    let mut indexed =
        open_index(&source, &coordinate_index, &destination).expect("should open indexed source");
    assert_eq!(indexed.metadata().term_count, metadata.term_count);
    assert_eq!(
        indexed.query("BRCA1").expect("exact query should work"),
        exact
    );
    assert_eq!(
        indexed
            .query_with_mode("RCA", MatchMode::Contains)
            .expect("contains query should work"),
        exact
    );

    let (exact_with_stats, exact_stats) = indexed
        .query_with_stats("BRCA1")
        .expect("exact query stats should work");
    assert_eq!(exact_with_stats, exact);
    assert_eq!(exact_stats.matching_records, 1);
    let (prefix_with_stats, prefix_stats) = indexed
        .query_with_mode_and_stats("BRCA", MatchMode::Prefix)
        .expect("prefix query stats should work");
    assert_eq!(prefix_with_stats, exact);
    assert_eq!(prefix_stats.matching_records, 1);

    assert!(
        indexed
            .query("missing")
            .expect("unknown terms are empty")
            .is_empty()
    );
    assert!(matches!(
        indexed.query_with_mode("missing[", MatchMode::Regex),
        Err(gai::Error::InvalidInput(_))
    ));

    let stale_source = directory.path().join("stale.gff3.gz");
    std::fs::copy(&source, &stale_source).expect("should copy source fixture");
    std::fs::write(&stale_source, b"changed source").expect("should stale source fixture");
    assert!(matches!(
        open_index(&stale_source, &coordinate_index, &destination),
        Err(gai::Error::Stale(message)) if message.contains("source fingerprint")
    ));
    assert!(matches!(
        inspect_index(directory.path().join("missing.gai")),
        Err(gai::Error::Io(_))
    ));
}

#[test]
fn cli_bed_index_query_supports_all_match_modes() {
    let directory = tempdir().expect("should create temp directory");
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/sorting");
    let source = fixture_dir.join("ecoli_k12_mg1655.bed.gz");
    let coordinate_index = fixture_dir.join("ecoli_k12_mg1655.bed.gz.tbi");
    let destination = directory.path().join("fixture-bed.gai");
    let binary = env!("CARGO_BIN_EXE_gai");
    let source = source.to_string_lossy().into_owned();
    let coordinate_index = coordinate_index.to_string_lossy().into_owned();
    let destination = destination.to_string_lossy().into_owned();

    let indexed = Command::new(binary)
        .args([
            "build-index",
            &source,
            "--attribute",
            "ignored",
            "--coordinate-index",
            &coordinate_index,
            "--output",
            &destination,
        ])
        .output()
        .expect("should build BED name index");
    assert!(indexed.status.success(), "stderr: {:?}", indexed.stderr);

    let exact_api = gai::query_index(&source, &coordinate_index, &destination, "thrL")
        .expect("one-shot BED exact query should work");
    assert_eq!(exact_api.len(), 1);
    let mut indexed_api = gai::open_index(&source, &coordinate_index, &destination)
        .expect("should open BED indexed source");
    let (prefix_api, prefix_stats) = indexed_api
        .query_with_mode_and_stats("thr", MatchMode::Prefix)
        .expect("reusable BED prefix query should work");
    assert!(prefix_api.len() >= 4);
    assert_eq!(prefix_stats.matching_records as usize, prefix_api.len());

    let query = |term: &str, mode: Option<&str>| {
        let mut command = Command::new(binary);
        command.args([
            "query-index",
            &source,
            term,
            "--coordinate-index",
            &coordinate_index,
            "--gai",
            &destination,
        ]);
        if let Some(mode) = mode {
            command.args(["--match", mode]);
        }
        command.output().expect("should query BED names")
    };

    let exact = query("thrL", None);
    assert!(exact.status.success(), "stderr: {:?}", exact.stderr);
    assert_eq!(String::from_utf8_lossy(&exact.stdout).lines().count(), 1);
    assert!(String::from_utf8_lossy(&exact.stdout).contains("\tthrL\t"));

    let prefix = query("thr", Some("prefix"));
    assert!(prefix.status.success(), "stderr: {:?}", prefix.stderr);
    assert!(String::from_utf8_lossy(&prefix.stdout).lines().count() >= 4);

    let contains = query("hr", Some("contains"));
    assert!(contains.status.success(), "stderr: {:?}", contains.stderr);
    assert!(String::from_utf8_lossy(&contains.stdout).contains("\tthrL\t"));

    let regex = query("^thr[ABCL]$", Some("regex"));
    assert!(regex.status.success(), "stderr: {:?}", regex.stderr);
    assert_eq!(String::from_utf8_lossy(&regex.stdout).lines().count(), 4);
}
