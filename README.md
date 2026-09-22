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
decisions around things such as attributes can be looked up by prefix or exact matches only.

## Example commands

```console
$ gai build-index annotations.gff3.gz \
    --attribute Name --attribute Alias --attribute gene_name
$ gai query-index annotations.gff3.gz BRCA1
$ gai query-index annotations.gff3.gz BRCA --match prefix
$ gai inspect-index annotations.gff3.gz.gai
$ gai sort annotations.gff3 > annotations.sorted.gff3
$ gai sort annotations.bed > annotations.sorted.bed
```

## General arguments

A coordinate index is discovered from an unambiguous sibling `.tbi` or `.csi`, but
can be explicitly referenced by the `--coordinate-index` flag. Query output is lossless GFF3
record text on stdout, while build progress and phase timings go to stderr.

## Building an index

An index is built via the `build-index` command. `--attribute` identifies which attribute values
to extract and index. It can be repeated to extract multiple attributes By default, the index is streamed
to stdout and `--output` can be used to save to an explicit path.

For GFF files, values are parsed, percent-decoded, and split into
valid array values. Normalization trims surrounding Unicode whitespace and,
by default, lowercases ASCII letters. `--case-sensitive` disables only the
lowercasing step. Punctuation, identifier versions, and other biological
normalizations are preserved.

For advanced control `--memory-budget`, `--compression-threads`, and
`--bgzf-threads` can be used to constrain the parallelization and memory utility of index building.

## Querying

Queries use exact normalized value matching by default. Passing `--match prefix` can be used to match based
on prefix. No substring or fuzzy matching is currently available.

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
print(indexed.metadata().term_count, stats.records_processed)
```
