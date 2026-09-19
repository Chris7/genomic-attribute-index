# GNI: GFF Name Index

GNI means **GFF Name Index**. It is a versioned custom binary format for
exact lookup of explicitly configured GFF3 attribute values. It complements a
BGZF GFF3 file and its TBI or CSI interval index:

```text
annotations.gff3.gz
annotations.gff3.gz.tbi       # or .csi
annotations.gff3.gz.gni
```

`ID` is not special. It is optional, is not required to build an index, and is
searchable only when passed explicitly with `--attribute ID`. GNI does not
group records into semantic features, infer relationships, create hashes, or
store BGZF virtual offsets.

## Building and querying with Rust

```console
$ gni build-index annotations.gff3.gz \
    --attribute Name --attribute Alias --attribute gene_name
$ gni query-index annotations.gff3.gz BRCA1
$ gni inspect-index annotations.gff3.gz.gni
```

The `gni` binary has exactly three top-level subcommands: `build-index`,
`query-index`, and `inspect-index`. A coordinate index is discovered from an
unambiguous sibling `.tbi` or `.csi`; pass `--coordinate-index` when both are
present or when the index has a nonstandard name. Query output is lossless GFF3
record text on stdout, while build progress and phase timings go to stderr.

`--attribute` is repeatable and required. The builder deduplicates repeated
names while preserving their first deterministic occurrence. `--coordinate-index`
and `--output` select explicit paths; otherwise an unambiguous `.tbi` or `.csi`
and `<input>.gni` are discovered.

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

GFF3 coordinates are one-based and inclusive. At the parser boundary GNI
performs:

```rust
let start = gff_start.checked_sub(1)?;
let length = gff_end.checked_sub(start)?;
```

Thus `100..150` becomes `[99, 150)` with length `51`. GNI stores only
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

## Binary layout (GNI 1.0)

All integers are little-endian. Fields identified as `u8`, `u16`, `u32`, or
`u64` have exactly that fixed width; variable integer fields are unsigned
LEB128 varints. Checksums are CRC-32. The first 256 bytes are the fixed GNI
1.0 header:

| byte range | field |
| ---: | --- |
| 0..4 | magic `GNI\x01` |
| 4..6 | major version `u16` |
| 6..8 | minor version `u16` |
| 8..12 | flags `u32` (only bit 0, independent zstd blocks, is defined) |
| 12 | byte order (`1` = little-endian) |
| 13 | coordinate convention (`1` = zero-based half-open `start + length`) |
| 14 | normalization (`0` = ASCII lowercase, `1` = case-sensitive) |
| 15 | reserved `u8`, zero |
| 16..20 | header size `u32` (256) |
| 20..24 | section-directory entry size `u32` (40) |
| 24..28 | section count `u32` (currently 6) |
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
| 172..176 | target spans per span block `u32` |
| 176..184 | absolute section-directory offset `u64` (256) |
| 184..192 | section-directory byte length `u64` |
| 192..200 | complete file size `u64` |
| 200..204 | coordinate-index reference count `u32` |
| 204..256 | reserved bytes, all zero |

The reader accepts the exact major version and any minor version less than or
equal to the implementation's supported minor version. A newer major or minor
version, an unknown flag, a nonzero reserved field, or nonzero reserved header
byte is rejected. The three fingerprints bind a GNI to the exact source,
coordinate index, and reference dictionary.

Immediately after the header is the section directory. Each entry is exactly
40 bytes, in little-endian field order:

| offset within entry | width | field |
| ---: | ---: | --- |
| 0 | 4 | section kind `u32` |
| 4 | 4 | section flags `u32` (zero) |
| 8 | 8 | absolute file offset `u64` |
| 16 | 8 | section length `u64` |
| 24 | 8 | section item count `u64` |
| 32 | 4 | section CRC-32 `u32` |
| 36 | 4 | reserved `u32`, zero |

Section offsets are absolute GNI file offsets and sections are non-overlapping.
The section CRC covers the complete section payload. The six current sections
are, in order: configured attributes, immutable `fst::Map` term dictionary,
postings directory, postings data, span directory, and span data. The `item
count` is the number of attributes, terms, postings blocks, or span blocks for
the first, second, third, and fifth sections; for the fourth and sixth
(postings-data and span-data) sections it is the complete payload byte count.

The configured-attribute section starts with a `u32` count, followed by that
many entries of `u32` UTF-8 byte length and exactly that many UTF-8 bytes. Names
are nonempty and deduplicated. The term section is a directly mappable
`fst::Map`; each FST value is a packed `u64` locator with the postings block
ID in bits 63..32 and the byte offset of the posting record in bits 31..0.

Each postings-directory entry is exactly 32 bytes:

| offset within entry | width | field |
| ---: | ---: | --- |
| 0 | 8 | compressed-block offset relative to postings-data section `u64` |
| 8 | 4 | compressed length `u32` |
| 12 | 4 | uncompressed length `u32` |
| 16 | 4 | CRC-32 of the uncompressed block `u32` |
| 20 | 1 | compression (`0` = raw, `1` = zstd) |
| 21 | 3 | reserved bytes, zero |
| 24 | 8 | reserved tail `u64`, zero |

Posting data blocks are independently bounded and checksummed; a block may
contain multiple records but a record never crosses a block boundary. A
posting record is an unsigned-LEB128 count, one unsigned-LEB128 absolute first
span ID, and then `count - 1` unsigned-LEB128 strictly positive deltas. Span
IDs are reconstructed cumulatively and must remain below the global span count.

Span IDs are globally deduplicated `(reference_id, start, length)` tuples,
assigned in reference/start/length order. Span blocks contain one reference
only and are capped by the header's target row count. The GNI 1.0 encoder and
decoder use a deliberate 52-byte span-directory entry; this adjusted width
includes the encoding, compression, and reserved metadata:

| offset within entry | width | field |
| ---: | ---: | --- |
| 0 | 8 | first global span ID `u64` |
| 8 | 4 | row count `u32` |
| 12 | 4 | reference ID `u32` |
| 16 | 8 | first start `u64` |
| 24 | 8 | compressed-block offset relative to span-data section `u64` |
| 32 | 4 | compressed length `u32` |
| 36 | 4 | uncompressed length `u32` |
| 40 | 4 | CRC-32 of the uncompressed block `u32` |
| 44 | 1 | start encoding (`0`, absolute varint; `1`, delta varint; `2`, FOR) |
| 45 | 1 | length encoding (`0`, varint; `1`, delta varint; `2`, FOR) |
| 46 | 1 | compression (`0` = raw, `1` = zstd) |
| 47 | 1 | reserved `u8`, zero |
| 48 | 4 | reserved tail `u32`, zero |

The span-block payload begins with a 20-byte header: `u32` start-stream byte
length, `u32` length-stream byte length, `u64` length FOR base, `u8` start
encoding, `u8` length encoding, `u8` start FOR bit width, and `u8` length FOR
bit width. The start stream immediately follows that header, then the length
stream; no padding is inserted. For the encoder's delta-varint start encoding
(`1`), the first start is omitted because it is `first_start` in the fixed
directory entry, and the stream contains `row_count - 1` deltas. For start FOR
encoding (`2`), the stream likewise contains only `row_count - 1` values,
packed least-significant-bit first, as offsets from `first_start` with an
implicit base of zero. The decoder also understands absolute start varints
(`0`) for forward/corruption compatibility, but the current encoder does not
emit them. Length varints (`0`) contain one absolute value per row; length
delta varints (`1`) contain one first value followed by deltas; length FOR
(`2`) contains one packed value per row using the `length_base`. Non-FOR
streams have zero bit width, and `length_base` is zero unless length FOR is
selected. FOR widths are bounded to 63 bits.

All block offsets are relative to their own GNI data section, never BGZF virtual
offsets into the source file. Block CRCs cover uncompressed payloads, while
section CRCs cover complete sections. Zstandard is selected independently per
block only when it is smaller than the raw payload.

Opening a GNI validates the magic, version, all section ranges, fixed-width
directories, section checksums, and cross-section counts. Block checksums,
compression sizes, varints, bit widths, and checked coordinate arithmetic are
validated as each block is accessed, so a lookup never decompresses unrelated
blocks. `IndexedGff::open` additionally checks all three fingerprints and
rejects stale source/index pairs. Writes use a same-directory temporary file,
`fsync`, and atomic rename.

`gni inspect-index` reports selected span encodings, the block target,
section sizes, uncompressed/compressed block totals and their ratios, and
independent zstd-block counts. `NameIndexReader::open_mmap` is available when
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
4,096-span blocks, delta-varint starts, ordinary-varint lengths, and optional
zstd-per-block. The report compares total index size and measured lookup
amplification/latency for 1,024, 4,096, and 16,384-row alternatives in
[`docs/compression-benchmark.md`](docs/compression-benchmark.md); it is a
measurement guide rather than a promised GNI/GFF size ratio.

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
`query_name_with_stats`:

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

The optional `genomic-attribute-index` Python package (imported as `gni`) is a PyO3/ABI3 extension built from
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
import gni

source = Path("annotations.gff3.gz")
tbi = Path("annotations.gff3.gz.tbi")
stats = gni.build_index(source, tbi, Path("annotations.gff3.gz.gni"), ["Name", "Alias"])
indexed = gni.open_index(source, tbi, Path("annotations.gff3.gz.gni"))
for record in indexed.query("BRCA1"):
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
  --build-arg RUST_TOOLCHAIN=nightly-2026-06-26 -t gni:rust-dev .
docker run --rm -it gni:rust-dev
```

The Python development image builds the ABI3 wheel and runs pytest by default.
Its build context is the repository root so the Python crate's path dependency
on the Rust library resolves correctly:

```console
docker build -f python/Dockerfile --build-arg RUST_VERSION=1.98.0 \
  --build-arg RUST_TOOLCHAIN=nightly-2026-06-26 -t gni:python-dev .
docker run --rm gni:python-dev
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
