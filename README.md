# GAI: Genomic Attribute Index

GAI means **Genomic Attribute Index**. It is a versioned custom binary format for
indexed lookup of explicitly configured GFF3 attribute values. It complements a
BGZF GFF3 file and its TBI or CSI interval index:

```text
annotations.gff3.gz
annotations.gff3.gz.tbi       # or .csi
annotations.gff3.gz.gai
```

`ID` is not special. It is optional, is not required to build an index, and is
searchable only when passed explicitly with `--attribute ID`. GAI does not
group records into semantic features, infer relationships, create hashes, or
store BGZF virtual offsets.

## Building and querying with Rust

```console
$ gai build-index annotations.gff3.gz \
    --attribute Name --attribute Alias --attribute gene_name
$ gai query-index annotations.gff3.gz BRCA1
$ gai query-index annotations.gff3.gz BRCA --match prefix
$ gai inspect-index annotations.gff3.gz.gai
$ gai sort annotations.gff3 > annotations.sorted.gff3
$ gai sort annotations.bed > annotations.sorted.bed
```

The `gai` binary has four top-level subcommands: `build-index`, `query-index`,
`inspect-index`, and `sort`. A coordinate index is discovered from an
unambiguous sibling `.tbi` or `.csi`; pass `--coordinate-index` when both are
present or when the index has a nonstandard name. Query output is lossless GFF3
record text on stdout, while build progress and phase timings go to stderr.

`--attribute` is repeatable and required. The builder deduplicates repeated
names while preserving their first deterministic occurrence. `--coordinate-index`
and `--output` select explicit paths; otherwise an unambiguous `.tbi` or `.csi`
and `<input>.gai` are discovered.

Queries use exact normalized value matching by default. Pass the typed
`--match prefix` option to stream every configured attribute value beginning
with the normalized query from the FST; `--match exact` is equivalent to the
default. No substring or fuzzy matching is provided.

`gai sort` infers GFF/GFF3 or BED from the input extension and writes sorted
records to stdout. GFF comment and directive lines remain first in source order;
feature records sort by contig, start, and end, with `ID`/`Parent` hierarchy
putting parents before children when all three coordinates tie. BED records use
the same contig/start/end ordering.

The build scan uses one reusable noodles GFF3 parser over either BGZF or plain
input, stops feature extraction at `##FASTA`, and hashes the exact complete
source bytes (including the FASTA tail). Term/span observations are sorted into
same-directory temporary runs once the configurable working-set budget is
reached, then merged deterministically; temporary runs are removed on success
or failure. Multi-pass compaction caps each merge at 64 open runs, keeping file
descriptor use bounded even for very large inputs. The CLI accepts `--memory-budget`, `--compression-threads`, and
`--bgzf-threads`, and reports phase progress to stderr while keeping indexed
records on stdout. Library callers can use `BuildOptions` and its optional
progress callback without any library-level stderr output.

Values are parsed by the noodles GFF3 parser, percent-decoded, and split into
valid array values. Normalization trims surrounding Unicode whitespace and,
by default, lowercases ASCII letters. `--case-sensitive` disables only the
lowercasing step. Punctuation, identifier versions, and other biological
normalizations are preserved.

## Coordinates and retrieval

GFF3 coordinates are one-based and inclusive. At the parser boundary GAI
performs:

```rust
let start = gff_start.checked_sub(1)?;
let length = gff_end.checked_sub(start)?;
```

Thus `100..150` becomes `[99, 150)` with length `51`. GAI stores only
`reference_id`, `start`, and `length`; `end` is reconstructed with checked
addition. A query resolves all posted spans in one batch, issues one exact
TBI/CSI interval query per span, unions overlapping BGZF chunks, and reads the
merged chunks once. It then retains only records whose decoded
`(reference_id, start, length)` exactly equals one requested span and whose
configured attribute value matches. Postings are coordinate ordered, while
TBI/CSI preserves source order for records sharing one start (including
different lengths), so overlapping spans cannot duplicate a record and
byte-identical records are still returned independently. Disjoint spans are
never replaced with bounding intervals.

## Binary layout (GAI 1.0)

All integers are little-endian. Fixed fields use exactly their stated `u8`,
`u16`, `u32`, or `u64` width; variable integers are unsigned LEB128. Checksums
are CRC-32. GAI 1.0 uses the split start/length layout described below (there
is no legacy combined-span layout). Readers require major `1`, minor `0`, and
reject newer versions, unknown flags, nonzero reserved bytes, invalid counts,
or unsafe ranges. The magic is `GAI\x01`.

The fixed 256-byte header is:

| byte range | field |
| ---: | --- |
| 0..4 | magic `GAI\x01` |
| 4..6 | major version `u16` (`1`) |
| 6..8 | minor version `u16` (`0`) |
| 8..12 | flags `u32` (bit 0: independent zstd blocks) |
| 12 | byte order (`1` = little-endian) |
| 13 | coordinates (`1` = zero-based half-open `start + length`) |
| 14 | normalization (`0` = ASCII lowercase, `1` = case-sensitive) |
| 15 | reserved `u8`, zero |
| 16..20 | header size `u32` (`256`) |
| 20..24 | section-directory entry size `u32` (`40`) |
| 24..28 | section count `u32` (`7`) |
| 28..32 | reserved `u32`, zero |
| 32..40 | term count `u64` |
| 40..48 | globally unique span count `u64` |
| 48..56 | term-to-span posting count `u64` |
| 56..64 | postings block count `u64` |
| 64..72 | span block count `u64` |
| 72..104 | source GFF SHA-256 |
| 104..136 | TBI/CSI SHA-256 |
| 136..168 | reference-dictionary SHA-256 |
| 168..172 | configured attribute count `u32` |
| 172..176 | target rows per reference-local span block `u32` |
| 176..184 | absolute section-directory offset `u64` (`256`) |
| 184..192 | section-directory byte length `u64` |
| 192..200 | complete file size `u64` |
| 200..204 | coordinate-index reference count `u32` |
| 204..256 | reserved bytes, all zero |

The directory follows the header. Every entry is exactly 40 bytes: `kind u32`,
`flags u32` (zero), absolute file `offset u64`, section `length u64`,
`item_count u64`, section CRC-32 `u32`, and reserved `u32` (zero), in that
order. Section offsets are absolute, non-overlapping file ranges and section
checksums cover complete payloads. The seven section kinds are 1 configured
attributes, 2 immutable `fst::Map` terms, 3 postings directory, 4 postings
data, 5 span directory, 6 starts data, and 7 lengths data. Item counts are
attribute count, term count, postings-block count, postings-data byte length,
span-block count, starts-data byte length, and lengths-data byte length,
respectively; all are checked against the payload.

The configured-attribute payload is a `u32` count followed by that many
`u32` UTF-8 byte length plus bytes. Names are nonempty and deduplicated. The
term payload is a directly mappable `fst::Map`. Each FST value is a packed
locator: postings block ID in bits 63..32 and record byte offset in bits 31..0.
Each 32-byte postings-directory entry contains, in order, compressed offset
relative to postings data (`u64`), compressed length (`u32`), uncompressed
length (`u32`), CRC-32 (`u32`), compression (`u8`, 0 raw or 1 zstd), three
reserved zero bytes, and a reserved zero `u64`. A postings record is a varint
count, one absolute first span ID varint, then `count - 1` strictly positive
delta varints.

Span IDs are globally deduplicated `(reference_id, start, length)` tuples in
reference/start/length order. Each reference-local block is represented by a
72-byte span-directory entry:

| offset | width | field |
| ---: | ---: | --- |
| 0 | 8 | first global span ID `u64` |
| 8 | 4 | row count `u32` |
| 12 | 4 | reference ID `u32` |
| 16 | 8 | first absolute start `u64` |
| 24 | 8 | starts-data offset `u64` |
| 32 | 4 | starts compressed length `u32` |
| 36 | 4 | starts uncompressed length `u32` |
| 40 | 4 | starts CRC-32 `u32` |
| 44 | 8 | lengths-data offset `u64` |
| 52 | 4 | lengths compressed length `u32` |
| 56 | 4 | lengths uncompressed length `u32` |
| 60 | 4 | lengths CRC-32 `u32` |
| 64 | 1 | start encoding (`1` = canonical delta-varint) |
| 65 | 1 | length encoding (`0` = varint, `2` = FOR) |
| 66 | 1 | starts compression (`0` raw, `1` zstd) |
| 67 | 1 | lengths compression (`0` raw, `1` zstd) |
| 68 | 4 | reserved `u32`, zero |

Starts and lengths are physically separate and independently compressed and
checksummed. A starts payload has no header: the directory's `span_count` and
`first_start` are its metadata. It contains exactly `span_count - 1` canonical
unsigned LEB128 deltas after the implicit first row, so equal starts encode as
zero and a single-row payload is empty. Deltas are checked-added to
`first_start`; nonminimal varints, truncation, trailing bytes, and overflow are
rejected. The decoder bounds the row count and payload before allocating and
requires exact end-of-stream consumption.

Each lengths payload begins with a 20-byte header: row count `u32`, encoding
`u8`, bit width `u8`, reserved `u16` (zero), FOR base `u64`, and stream length
`u32`; the stream follows immediately with no padding. Varint lengths are one
positive absolute value per row. FOR lengths are packed low-bit-first values
with a base and width (at most 63 bits). Length row counts, exact stream size,
padding, positive values, and checked `start + length` are validated.

All block offsets are relative to their own data section, never BGZF virtual
offsets. Component CRCs cover uncompressed payloads and section CRCs cover
complete sections. Zstandard is selected independently for each component only
when smaller. Opening the index validates all fixed directories and counts;
lookup decompresses only referenced posting/start/length blocks. `gai
inspect-index` reports the three independent data sizes, uncompressed/compressed
ratios, and selected encodings. `NameIndexReader::open_mmap` keeps the FST and
fixed directories backed by an operating-system read-only mapping.

Opening a GAI validates the magic, version, all section ranges, fixed-width
directories, section checksums, and cross-section counts. Block checksums,
compression sizes, varints, bit widths, and checked coordinate arithmetic are
validated as each block is accessed, so a lookup never decompresses unrelated
blocks. `IndexedGff::open` additionally checks all three fingerprints and
rejects stale source/index pairs. Writes use a same-directory temporary file,
`fsync`, and atomic rename.

`gai inspect-index` reports selected span encodings, the block target, section
sizes, uncompressed/compressed block totals and their ratios, and independent
zstd-block counts. `NameIndexReader::open_mmap` is available when
callers want the fixed directories and FST backed directly by an operating-
system read-only mapping.

## Compression measurements

The ignored benchmark in `benches/compression.rs` generates deterministic
structured and random terms, shared spans, clustered and dispersed starts,
and short and long intervals. It compares fixed-width triples, interleaved
varints, row-oriented and columnar streams, delta-varints, frame-of-reference
packing, and optional zstd. Run it with:

```console
cargo bench --bench compression
```

The current reference run (Rust 1.98, x86-64 Linux, zstd level 3) selected
4,096-span blocks, canonical delta-varint starts, adaptive varint/FOR lengths,
and optional zstd-per-component. The report compares total index size and measured lookup
amplification/latency for 1,024, 4,096, and 16,384-row alternatives in
[`docs/compression-benchmark.md`](docs/compression-benchmark.md); it is a
measurement guide rather than a promised GAI/GFF size ratio.

The build-performance benchmark uses a deterministic 100,000-record fixture
and compares the current single-pass, no-spill pipeline with a 1 MiB
spill/merge budget. It reports phase timings, records per second, and the
bounded working-set proxy:

```console
cargo bench --bench build-performance
```

The captured reference output is in
[`docs/build-performance.md`](docs/build-performance.md).

The indexed query benchmark uses the checked-in GENCODE fixture and reports
the exact interval/chunk amplification and candidate counts exposed by
`query_name_with_stats` (and its mode-aware companion
`query_name_with_mode_and_stats`):

```console
cargo bench --bench query-performance --all-features
```

See [`docs/query-performance.md`](docs/query-performance.md) for a captured
release run and the comparison with the earlier repeated-scan implementation.

## Build performance

The builder parallelizes independent posting and reference-local span block
compression with an indexed Rayon pool. Collection order and all final bytes
remain deterministic across worker counts. BGZF decompression can likewise use
the noodles multithreaded reader; plain input and query retrieval stay
single-stream. Exact source and coordinate-index fingerprints are retained.

## Current limitations

The public reader currently uses a small per-call decode path rather than an
optional LRU cache. Query retrieval expects BGZF source data when using a
TBI/CSI index.

## Python package

The optional `genomic-attribute-index` Python package (imported as `gai`) is a PyO3/ABI3 extension built from
`python/`. It returns structured `GffRecord`, `IndexMetadata`, `BuildStats`,
and `QueryStats` objects, and maps stale/corrupt/input/I/O failures to typed
exceptions. Long Rust build, open, and query operations release the GIL.

Install a published wheel or build an editable checkout:

```console
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install genomic-attribute-index  # published wheel

# or, from a checkout:
python -m pip install maturin pytest
maturin develop --manifest-path python/Cargo.toml --features abi3,extension-module
pytest -q python/tests
```

Example API usage:

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

The package's deterministic fixture helper creates both TBI and CSI test
indexes. To run the same installed-extension suite locally, use
`maturin develop` followed by `pytest -q python/tests`; the tests do not mock
the extension.

## Docker development

The repository pins Rust nightly-2026-06-26 in `rust-toolchain.toml`; the root
Dockerfile uses a reproducible Rust 1.98 base and installs that same pinned
toolchain, then runs the locked all-target test suite and starts a non-root
shell:

```console
docker build --build-arg RUST_VERSION=1.98.0 \
  --build-arg RUST_TOOLCHAIN=nightly-2026-06-26 -t gai:rust-dev .
docker run --rm -it gai:rust-dev
```

The Python development image builds the ABI3 wheel and runs pytest by default.
Its build context is the repository root so the Python crate's path dependency
on the Rust library resolves correctly:

```console
docker build -f python/Dockerfile --build-arg RUST_VERSION=1.98.0 \
  --build-arg RUST_TOOLCHAIN=nightly-2026-06-26 -t gai:python-dev .
docker run --rm gai:python-dev
```

## CI and releases

`.github/workflows/ci.yml` runs Rust formatting, clippy, all-target tests,
documentation, release builds, `cargo package`, and a Python 3.10--3.14
wheel/test matrix. `docker.yml` validates both Dockerfiles on pushes and pull
requests. Cargo and Python dependency caches are keyed by the checked-in lock
files, and build/package jobs use locked resolution.

`release.yml` is intentionally gated: a published GitHub release runs the
verification and artifact jobs, while `workflow_dispatch` defaults to
build-only and requires the `publish` boolean to be enabled for publication.
The `vX.Y.Z` tag must match the root Cargo version, the Python crate version,
and `python/pyproject.toml`. The crates.io job uses the `crates-io-auth-action`
OIDC exchange and the `crates-io` environment; the PyPI job uses trusted
publishing with the `pypi` environment and `id-token: write`. Configure those
two repository environments with the corresponding crates.io and PyPI trusted
publisher policies before enabling publication. No package is published from
pull requests.
