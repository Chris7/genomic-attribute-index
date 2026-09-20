# GAI build-performance benchmark

`benches/build-performance.rs` creates a deterministic 100,000-record,
four-reference BGZF GFF3/TBI fixture. It compares the same default GAI 1.0 builder in
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
bytes; deterministic fixture generation took 457.456 ms.

| configuration | elapsed | records/sec | peak working-set proxy | spill | scan | merge | postings | spans | serialize |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| single-pass no-spill baseline (1 GiB) | 192.551 ms | 519,344 | 13,800,000 B | 50.202 ms | 93.659 ms | 4.224 ms | 41.345 ms | 1.332 ms | 0.129 ms |
| bounded spill (1 MiB) | 174.406 ms | 573,375 | 1,048,662 B | 29.105 ms | 86.930 ms | 5.505 ms | 50.623 ms | 1.324 ms | 0.102 ms |

The bounded run was 1.09x the baseline throughput and used 0.08x the measured
working-set proxy. The two output GAI files were byte-identical. The working
set is the collector's bounded term/span observation estimate; it excludes the
final serialized output and OS page cache, so it is a reproducible bound proxy,
not a process RSS measurement.

## Real GENCODE fixture

An additional release build used the repository's
`fixtures/gencode_sorted.gff.gz` (3,766,032 records, 83,531,736 compressed
source bytes) with `gene_name`. The split delta-start/separate-length
writer completed its measured library work in 21.960 s (171,494 records/s),
with phase timings of scan 19.886 s, spill 1.480 s, merge 0.083 s, postings
0.386 s, spans 0.112 s, and serialize 0.004 s. Its bounded working-set proxy
was 67,108,927 bytes. The resulting 4,170,123-byte GAI was byte-identical to a
second release build (SHA-256 `9fef0a0313f89a8a72a8396dbe7aab9d5d9b7a89e940b775aac142b8e451184c`).
For size context, the pre-split checked-in predecessor was 4,200,460 bytes and
the earlier split-start experiment was 4,806,892 bytes; those layouts are
intentionally not accepted by the current reader.
