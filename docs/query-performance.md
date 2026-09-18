# GNI indexed-query performance

`benches/query-performance.rs` measures a real indexed lookup against the
checked-in `fixtures/gencode_sorted.gff.gz`, its `.tbi`, and its `.gni`. It
does not modify the fixture. The benchmark is a custom ignored target, so it
is not run by ordinary tests:

```console
cargo bench --bench query-performance --all-features
```

The measurement below was captured on x86-64 Linux with Rust 1.98,
noodles-bgzf 0.43.0, and the release profile in this repository. The fixture
contains 1,436 matching `gene_name=brca1` records spread across 160 exact
spans (75 reference/start groups). The open phase includes the GNI mmap and
the source, coordinate-index, and reference-dictionary fingerprint checks.

| phase | elapsed | result |
| --- | ---: | ---: |
| open and fingerprint | 76.108 ms | reader opened successfully |
| indexed query | 21.857 ms | 1,436 records, 65,701 records/s |

The query instrumentation for that run was:

| metric | value |
| --- | ---: |
| requested spans | 160 |
| span blocks decoded | 1 |
| exact TBI interval queries | 160 |
| raw BGZF chunks returned | 3,360 |
| merged BGZF chunks read | 21 |
| unique candidate records parsed | 6,984 |
| matching records | 1,436 |
| uncompressed candidate bytes read | 3,318,100 |

Each span query is exact; no disjoint span is widened into a bounding
interval. Chunks are sorted and unioned by virtual offset, then read through
one BGZF reader with one reusable noodles GFF line/record parser. Candidates
reached through overlapping chunks are deduplicated by their source virtual
position, not by text, so identical source lines at distinct positions remain
distinct results. Span blocks are decoded once per query even when many
postings refer to rows in the same block.

For context, a prior release binary using the old per-span/full-source-order
query path measured about 4.79 s for the same CLI lookup on this machine. That
number is a machine-local before/after observation rather than a stable
benchmark baseline; the reproducible target above reports the current path's
phase and amplification data. The current indexed query is also substantially
below the approximately 2.03 s local `gunzip | grep` comparison previously
used for this fixture.
