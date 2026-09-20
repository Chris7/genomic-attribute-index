# GAI compression benchmark

This report is produced from the deterministic synthetic workload in
`benches/compression.rs`. It contains 7,200 coordinate-sorted records on three
references, 7,172 globally unique spans (including intentionally identical
spans shared by multiple records), structured and pseudo-random terms,
repeated aliases, clustered and widely separated starts, and short and long
spans.
The benchmark is ignored by default because it is intended for comparative
measurements rather than a test gate.

Reference run: Rust 1.98, x86-64 Linux, zstd 1.5.7, compression level 3.
The current output uses GAI 1.0 split span blocks: canonical delta-varint starts
and independently encoded varint/FOR lengths. The three block-size runs below
use the real writer/reader path; the candidate table is an apples-to-apples
structural comparison only.
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
benchmark also builds a BGZF GFF3/TBI/GAI fixture and reports these measured
values (the GFF byte count is uncompressed):

| block size | starts raw | starts zstd | lengths raw | lengths zstd | span blocks | total GAI bytes | avg lookup bytes | avg lookup latency |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,024 | 11,012 | 9,350 | 12,180 | 10,371 | 8 | 132,339 | 34,142 | 16.783544 ms |
| 4,096 | 11,018 | 9,530 | 12,100 | 10,214 | 4 | 132,074 | 40,086 | 16.868991 ms |
| 16,384 | 11,020 | 9,681 | 12,080 | 10,316 | 3 | 132,255 | 46,672 | 17.213468 ms |

| end-to-end metric | measured value |
| --- | ---: |
| sampled query terms (all known) | 58 |
| uncompressed GFF bytes | 585,385 |
| BGZF GFF bytes | 129,342 |
| TBI bytes | 15,669 |
| selected GAI bytes | 132,074 |
| FST bytes | 90,713 |
| postings data bytes | 20,740 |
| span component data bytes (starts/lengths) | 9,530 / 10,214 |
| span component raw bytes (starts/lengths) | 11,018 / 12,100 |
| fixed block-directory bytes | 320 |
| bytes per term | 18.09 |
| bytes per posting | 9.19 |
| bytes per unique span | 18.42 |
| average known-term lookup bytes decompressed (4,096 blocks) | 40,086 |
| average known-term lookup latency (4,096 blocks) | 16.868991 ms |

The builder chooses FOR lengths only when their exact structural stream is
smaller than varints, then applies zstd independently to starts and lengths
only when each compressed component is smaller than its raw payload. The
1,024-row target reduces lookup amplification at a modest directory and total
size cost; 16,384 rows reduces directory overhead but increases decoded bytes.
4,096 rows remains the selected balance for this workload. Lookup latency
includes the required TBI query and exact-span candidate parsing.

### Split versus combined span payloads

The selected 4,096-row run was also compared with a measurement-only combined
layout. The candidate concatenates each block's delta-start stream and its
adaptive length payload, compresses that complete payload once, and uses a
52-byte combined block-directory entry plus one fewer 40-byte section-directory
entry. This includes the existing 20-byte length header and all projected
directory overhead; it does not change the shipped split layout. On this
workload, split components were 23,118 raw bytes and 19,744 compressed bytes;
the combined candidate was 23,118 raw and 20,632 compressed bytes (4 zstd
blocks), projecting to 132,842 total GAI bytes versus 132,074 split bytes
(+768 bytes, +0.58%). Both layouts decode the same 40,086 average bytes for a
known-term lookup because starts and lengths are both required; only the split
reader latency (16.868991 ms average) was measured. The result supports
retaining independent compression/checksums: combined compression did not
recover enough redundancy to offset its block-directory and section overhead.

## Real GENCODE comparison

On the checked-in `fixtures/gencode_sorted.gff.gz` with `gene_name`, the
pre-split predecessor GAI was 4,200,460 bytes and the earlier split-start
experiment was 4,806,892 bytes. The current split delta-varint
layout is 4,170,123 bytes (SHA-256
`9fef0a0313f89a8a72a8396dbe7aab9d5d9b7a89e940b775aac142b8e451184c`),
0.72% smaller than the pre-split predecessor and 13.25% smaller than the
earlier split-start experiment. It is an intentional incompatible replacement, not a dual-format
reader. Component measurements are:

| component | raw bytes | compressed bytes |
| --- | ---: | ---: |
| starts (791 delta-varint blocks) | 1,805,223 | 1,376,301 |
| lengths (791 varint/FOR blocks) | 2,248,450 | 2,011,278 |
| span directory (72 bytes/block) | 56,952 | 56,952 |

The release build processed 3,766,032 records in 21.960 s with one BGZF and
one compression worker (scan 19.886 s, spill 1.480 s, merge 0.083 s,
postings 0.386 s, spans 0.112 s, serialize 0.004 s). A cold-ish exact
`gene_name=BRCA1` query returned 1,436 records in 18.587 ms; prefix `BRCA`
returned 2,312 records. The exact query
requested 160 spans, decoded one span block, issued 160 exact TBI intervals,
unioned 3,360 raw chunks into 21 reads, parsed 6,984 unique candidates, and
read 3,318,100 uncompressed candidate bytes.

The same measurement-only combined-span projection used for the synthetic
workload produced 4,053,673 raw and 3,432,876 compressed span payload bytes
(544 zstd blocks), a 41,132-byte combined block directory, and a projected
4,199,560-byte GAI. That is 29,437 bytes (+0.71%) larger than the current
split file. The BRCA1 lookup decoded 77,680 GAI postings/span bytes in either
projection; combined query latency was not measured because no combined reader
is shipped. The real-fixture result reinforces the synthetic recommendation to
retain independently compressed and checksummed starts and lengths.

The report includes the format candidates requested by the design. The
benchmark constructs its own deterministic BGZF GFF3 and TBI fixture, so it
does not require samtools or any external indexing executable.
