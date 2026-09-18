# GNI compression benchmark

This report is produced from the deterministic synthetic workload in
`benches/compression.rs`. It contains 7,200 coordinate-sorted records on three
references, 7,172 globally unique spans (including intentionally identical
spans shared by multiple records), structured and pseudo-random terms,
repeated aliases, clustered and widely separated starts, and short and long
spans.
The benchmark is ignored by default because it is intended for comparative
measurements rather than a test gate.

Reference run: Rust 1.98, x86-64 Linux, zstd 1.5.7, compression level 3.
The numbers below are the captured output of the same algorithms with the
workload seed fixed to `0x47534e49`.

| candidate | bytes | ratio vs fixed triples |
| --- | ---: | ---: |
| interleaved fixed-width triples (unique spans) | 143,440 | 1.00x |
| interleaved unsigned varints (all records) | 43,613 | 3.29x |
| row-oriented unsigned varints (unique spans) | 43,442 | 3.30x |
| columnar start delta-varints (reference-reset) | 23,040 | 6.23x |
| columnar frame-of-reference (reference-local) | 35,860 | 4.00x |
| columnar FOR + zstd blocks | 30,640 | 4.68x |

The selected default is 4,096 spans per reference-specific block. The same
benchmark also builds a BGZF GFF3/TBI/GNI fixture and reports these measured
values (the GFF byte count is uncompressed):

| block size | structural bytes | zstd bytes | span blocks | total GNI bytes | avg lookup bytes | avg lookup latency |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,024 | 23,192 | 20,781 | 8 | 133,199 | 34,142 | 236.975332 ms |
| 4,096 | 23,118 | 20,672 | 4 | 132,882 | 40,086 | 234.376279 ms |
| 16,384 | 23,100 | 20,844 | 3 | 133,002 | 46,672 | 231.810042 ms |

| end-to-end metric | measured value |
| --- | ---: |
| sampled query terms (all known) | 58 |
| uncompressed GFF bytes | 585,385 |
| BGZF GFF bytes | 129,342 |
| TBI bytes | 15,669 |
| selected GNI bytes | 132,882 |
| FST bytes | 90,713 |
| postings data bytes | 20,740 |
| span-table data bytes | 20,672 |
| fixed block-directory bytes | 240 |
| bytes per term | 18.20 |
| bytes per posting | 9.25 |
| bytes per unique span | 18.53 |
| average known-term lookup bytes decompressed (4,096 blocks) | 40,086 |
| average known-term lookup latency (4,096 blocks) | 234.376279 ms |

The builder chooses frame-of-reference streams only when their exact
structural stream is smaller, then applies zstd only when the compressed block
is smaller than its structural payload. Delta starts reset at each reference,
and all sampled query terms are selected from the generated workload (there
are no unknown-term zero-I/O samples). The 1,024-row target reduces lookup
amplification by decoding smaller blocks at a modest directory and total-size
cost; 16,384 rows slightly reduces directory overhead but increases decoded
bytes. 4,096 rows is the selected balance in this run. These measurements do
not impose an arbitrary GNI-to-GFF ratio; real datasets should be benchmarked
with their own term and coordinate distributions. Lookup latency includes the
required TBI query and exact-span candidate parsing.

The report includes the format candidates requested by the design. The
benchmark constructs its own deterministic BGZF GFF3 and TBI fixture, so it
does not require samtools or any external indexing executable.
