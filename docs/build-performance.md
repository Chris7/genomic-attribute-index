# GNI build-performance benchmark

`benches/build-performance.rs` creates a deterministic 100,000-record,
four-reference BGZF GFF3/TBI fixture. It compares the same GNI 1.0 builder in
single-pass, no-spill mode (1 GiB scan budget) and bounded spill/merge mode
(1 MiB scan budget). This is not a historical nested-tree baseline. Both runs
use four BGZF and four block-compression workers. The
benchmark is ignored by ordinary tests and runs with:

Spill compaction is multi-pass with a fixed 64-run fan-in, so the merge keeps
file-descriptor use bounded even when a small memory budget produces many runs.

```console
cargo bench --bench build-performance --all-features
```

Reference run: Rust 1.98, x86-64 Linux, zstd 1.5.7, Rayon 1.12.0. The fixture
contained 6,973,032 uncompressed GFF bytes and 1,141,726 compressed BGZF
bytes; deterministic fixture generation took 460.925 ms.

| configuration | elapsed | records/sec | peak working-set proxy | spill | scan | merge | postings | spans | serialize |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| single-pass no-spill baseline (1 GiB) | 194.134 ms | 515,108 | 13,800,000 B | 47.443 ms | 93.595 ms | 4.256 ms | 41.470 ms | 5.212 ms | 0.149 ms |
| bounded spill (1 MiB) | 178.057 ms | 561,618 | 1,048,662 B | 28.982 ms | 87.276 ms | 5.721 ms | 50.532 ms | 4.669 ms | 0.085 ms |

The bounded run was 1.09x the baseline throughput and used 0.08x the measured
working-set proxy. The two output GNI files were byte-identical. The working
set is the collector's bounded term/span observation estimate; it excludes the
final serialized output and OS page cache, so it is a reproducible bound proxy,
not a process RSS measurement.

## Real GENCODE fixture

An additional release build used the repository's
`fixtures/gencode_sorted.gff.gz` (3,766,032 records, 83,531,736 compressed
source bytes) with `gene_name`. The builder completed its measured library
work in 22.571 s (166,850 records/s), with phase timings of scan 20.565 s,
spill 1.504 s, merge 0.070 s, postings 0.348 s, spans 0.070 s, and serialize
0.005 s. Its bounded working-set proxy was 67,108,927 bytes. The resulting
4,200,460-byte GNI was byte-identical to the checked-in reference GNI; the
existing fixture was not overwritten.
