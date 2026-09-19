//! Generates the small deterministic BGZF/TBI/CSI fixture used by the Python tests.

use std::{env, fs::File, io::Write, path::PathBuf};

use noodles::{
    bgzf,
    core::Position,
    csi::{
        self,
        binning_index::{
            Indexer,
            index::reference_sequence::{bin::Chunk, index::BinnedIndex},
        },
    },
    tabix,
};

const RECORDS: [&str; 7] = [
    "##gff-version 3",
    "##sequence-region chr1 1 1000",
    "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=Alpha;Alias=Beta,Gamma",
    "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=Alpha;ID=identical",
    "chr1\tsrc\tgene\t30\t40\t.\t-\t.\tName=Other%20Gene",
    "chr2\tsrc\tgene\t5\t15\t.\t+\t.\tAlias=Alpha",
    "chr2\tsrc\tgene\t5\t25\t.\t+\t.\tName=Alpha",
];

fn main() {
    let directory = env::args()
        .nth(1)
        .map(PathBuf::from)
        .expect("usage: python_fixture OUTPUT_DIRECTORY");
    std::fs::create_dir_all(&directory).expect("should create fixture directory");

    let source_path = directory.join("fixture.gff3.gz");
    let tbi_path = directory.join("fixture.gff3.gz.tbi");
    let csi_path = directory.join("fixture.gff3.gz.csi");
    let mut writer = File::create(&source_path)
        .map(bgzf::io::Writer::new)
        .expect("should create BGZF source");
    let mut tbi_indexer = tabix::index::Indexer::default();
    tbi_indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
    let mut csi_names = csi::binning_index::index::header::ReferenceSequenceNames::new();
    csi_names.insert("chr1".into());
    csi_names.insert("chr2".into());
    let csi_header = csi::binning_index::index::Header::builder()
        .set_format(csi::binning_index::index::header::Format::Generic(
            csi::binning_index::index::header::format::CoordinateSystem::Gff,
        ))
        .set_reference_sequence_names(csi_names)
        .build();
    let mut csi_indexer = Indexer::<BinnedIndex>::default().set_header(csi_header);

    for line in RECORDS {
        if line.starts_with('#') {
            writeln!(writer, "{line}").expect("should write directive");
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let reference_id = if fields[0] == "chr1" { 0 } else { 1 };
        let start = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
        let end = Position::try_from(fields[4].parse::<usize>().unwrap()).unwrap();
        let start_position = writer.virtual_position();
        writeln!(writer, "{line}").expect("should write record");
        let end_position = writer.virtual_position();
        let chunk = Chunk::new(start_position, end_position);
        tbi_indexer
            .add_record(fields[0], start, end, chunk)
            .expect("should index TBI record");
        csi_indexer
            .add_record(Some((reference_id, start, end, true)), chunk)
            .expect("should index CSI record");
    }
    writer.finish().expect("should finish BGZF source");

    let tbi = tbi_indexer.build();
    let mut tbi_writer =
        tabix::io::Writer::new(File::create(tbi_path).expect("should create TBI output"));
    tbi_writer.write_index(&tbi).expect("should write TBI");
    tbi_writer.try_finish().expect("should finish TBI");

    let csi = csi_indexer.build(2);
    let mut csi_writer =
        csi::io::Writer::new(File::create(csi_path).expect("should create CSI output"));
    csi_writer.write_index(&csi).expect("should write CSI");
}
