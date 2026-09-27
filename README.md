# GAI: Genomic Attribute Index

GAI is an index for fast lookup of GFF and bed attributes. It is designed to complement
a tabix index to quickly extract matching records.

```text
annotations.gff3.gz
annotations.gff3.gz.tbi       # or .csi
annotations.gff3.gz.gai       # Attribute index
```

GAI can be installed from pypi via `pip install genomic-attribute-index` or for Rust users 
via `cargo install genomic-attribute-index`.

## Design Choices

Compliance of GFF files is notoriously bad. Thus, we have very little validation beyond:
 * The GFF must be indexable by tabix
 * The GFF record must be parseable by noodles (our GFF parser of choice here)

GFF files need to be sorted to be indexed by tabix. However, sorting a GFF is a pain. Thus,
we have a sort function to sort entries correctly. Even better, if a GFF file does make use of
Parent/Child tags, the parent entries will be sorted above the child entries to make the lives of
downstream renderers and processors easier. The sort function will similarly sort bed files.

Although multiple attributes can be indexed, the values will be collapsed. For example, 
`name=BRCA,alias=BRCA` will be stored as a single `BRCA` match. The logic is to extract the relevant
portions of the GFF/bed as quickly as possible, with any downstream filtering the responsibility
of the application.

We have designed for speed and compressibility of the attribute index. Because a GFF can be so feature
rich, it makes little sense to have an index whose size is similar to a compressed GFF. This influences our
decisions around things such as attributes can be looked up exactly, by prefix, by literal
substring, or with a regular expression.

## Example commands

```console
$ gai build-index annotations.gff3.gz \
    --attribute Name --attribute Alias --attribute gene_name
$ gai query-index annotations.gff3.gz BRCA1
$ gai query-index annotations.gff3.gz BRCA --match prefix
$ gai query-index annotations.gff3.gz RCA --match contains
$ gai query-index annotations.gff3.gz '^BRCA[0-9]+$' --match regex
$ gai inspect-index annotations.gff3.gz.gai
$ cargo run --features profiling -- profile query-index annotations.gff3.gz BRCA1
$ cargo run --features profiling -- profile --sample query-index annotations.gff3.gz BRCA1
$ gai sort annotations.gff3 > annotations.sorted.gff3
$ gai sort annotations.bed > annotations.sorted.bed
```

## General arguments

A coordinate index is discovered from an unambiguous sibling `.tbi` or `.csi`, but
can be explicitly referenced by the `--coordinate-index` flag. Query output is lossless source
record text on stdout, while build progress and phase timings go to stderr.
Profiling support is opt-in and is not included in default builds. Build with
the Cargo `profiling` feature and use `gai profile <command> ...` to print an
aggregated tracing summary to stderr, grouped by call stack and sorted by total
time. Each row reports call count, inclusive total duration, and average duration.
Command results remain on stdout. Add `--sample` after `profile` to use a 100 Hz
CPU sampler instead; it groups sampled application call stacks and reports sample
counts, estimated time, and percentage to stderr without tracing.
Sampling is currently supported on Unix platforms; function tracing is available
on all supported platforms.

## Building an index

An index is built via the `build-index` command. `--attribute` identifies which attribute values
to extract and index. It can be repeated to extract multiple attributes By default, the index is created
as a .gai file with the same prefix as input. `--output` can be used to save to an explicit path.

For GFF files, values are parsed, percent-decoded, and split into
valid array values. Normalization trims surrounding Unicode whitespace and,
by default, lowercases ASCII letters. `--case-sensitive` disables only the
lowercasing step. Punctuation, identifier versions, and other biological
normalizations are preserved.

For advanced control `--memory-budget`, `--compression-threads`, and
`--bgzf-threads` can be used to constrain the parallelization and memory utility of index building.

## Querying

Queries use exact normalized value matching by default. `--match prefix` matches values beginning
with the query, and `--match contains` matches a literal substring anywhere in a value. Regex mode
uses an unanchored Unicode-aware regular expression search; `^` and `$` can anchor a full-value
match. Regex syntax is preserved, so escapes such as `\S` and `\D` retain their meaning. In a
case-insensitive index, regex matching uses the regex engine's Unicode case folding, while indexed
values continue to use GAI's ASCII-only lowercasing normalization. Query boundary whitespace is
trimmed, and an empty query returns no results. Contains and regex queries scan the distinct indexed
terms incrementally; exact and prefix queries retain their direct FST lookup paths. Existing indexes
work with the new modes without rebuilding.

## Sorting

`gai sort` infers GFF/GFF3 or BED from the input extension and writes sorted
records to stdout. GFF comment and directive lines remain first in source order;
feature records sort by contig, start, and end, with `ID`/`Parent` hierarchy
putting parents before children when all three coordinates tie. BED records use
the same contig/start/end ordering.

For vey large files, `--disk-sort` can be passed to spill to disk for sorting where
an in-memory sort is not possible.

For GFF files containing sequences, sorting and parsing stops when the `##FASTA` tag
is encountered. The records are then sorted and the sorted GFF is emited including the
trailing fasta contents.

 and reports phase progress to stderr while keeping indexed
records on stdout. Library callers can use `BuildOptions` and its optional
progress callback without any library-level stderr output.

## Benchmarks

Run [`scripts/benchmark_index.sh`](scripts/benchmark_index.sh) with
`./scripts/benchmark_index.sh > benchmark-results.md` to reproduce these tables. The input files
are compressed BGZF files; reported sizes use decimal MB. Indexes use the default case-insensitive
matching. Timings are single-run wall-clock measurements. Index build timings exclude compilation;
query timings include the complete CLI invocation, including index opening, source fingerprint
validation, and captured output.

### Index build results

| Index / file | Attributes indexed | Raw compressed file size (MB) | Index size (MB) | Build time |
| :--- | :--- | ---: | ---: | ---: |
| gencode_sorted.gff.gz | gene_name | 83.53 | 4.17 | 8.622s |
| gencode_sorted.gff.gz | gene_name, transcript_name | 83.53 | 7.28 | 10.896s |
| gencode_sorted.gff.gz | ID, gene_name, transcript_name | 83.53 | 22.01 | 15.899s |
| gencode_v46.bed.gz | name | 10.34 | 3.27 | 0.437s |

### Query results

| Index queried | Attributes indexed | Query type | Query | Records matched | Query time |
| :--- | :--- | :--- | :--- | ---: | ---: |
| gencode_sorted.gff.gz | gene_name | exact | `brca1` | 1436 | 0.088s |
| gencode_sorted.gff.gz | gene_name | contains | `orf` | 16797 | 1.205s |
| gencode_sorted.gff.gz | gene_name | prefix | `brca` | 2312 | 0.120s |
| gencode_sorted.gff.gz | gene_name | regex | `c\d+orf` | 16033 | 1.153s |
| gencode_sorted.gff.gz | gene_name, transcript_name | exact | `brca1` | 1436 | 0.110s |
| gencode_sorted.gff.gz | gene_name, transcript_name | contains | `orf` | 16797 | 1.223s |
| gencode_sorted.gff.gz | gene_name, transcript_name | prefix | `brca` | 2312 | 0.107s |
| gencode_sorted.gff.gz | gene_name, transcript_name | regex | `c\d+orf` | 16033 | 1.192s |
| gencode_sorted.gff.gz | ID, gene_name, transcript_name | exact | `brca1` | 1436 | 0.099s |
| gencode_sorted.gff.gz | ID, gene_name, transcript_name | contains | `orf` | 16797 | 1.415s |
| gencode_sorted.gff.gz | ID, gene_name, transcript_name | prefix | `brca` | 2312 | 0.108s |
| gencode_sorted.gff.gz | ID, gene_name, transcript_name | regex | `c\d+orf` | 16033 | 1.413s |
| gencode_v46.bed.gz | name | exact | `ENST00000607096.1` | 1 | 0.039s |
| gencode_v46.bed.gz | name | contains | `607096` | 1 | 0.077s |
| gencode_v46.bed.gz | name | prefix | `ENST000006070` | 50 | 0.105s |
| gencode_v46.bed.gz | name | regex | `^ENST00000607096\.1$` | 1 | 0.067s |


## Python API

```python
from pathlib import Path
import gai

source = Path("annotations.gff3.gz")
tbi = Path("annotations.gff3.gz.tbi")
stats = gai.build_index(source, tbi, Path("annotations.gff3.gz.gai"), ["Name", "Alias"])
indexed = gai.open_index(source, tbi, Path("annotations.gff3.gz.gai"))
for record in indexed.query("BRCA1"):  # exact is the default
    print(record.reference_sequence_name, record.start, record.attributes)
for record in indexed.query("BRCA", match="prefix"):
    print(record.reference_sequence_name, record.start, record.attributes)
for record in indexed.query("RCA", match="contains"):
    print(record.reference_sequence_name, record.start, record.attributes)
for record in indexed.query(r"^BRCA[0-9]+$", match="regex"):
    print(record.reference_sequence_name, record.start, record.attributes)
print(indexed.metadata().term_count, stats.records_processed)
```
