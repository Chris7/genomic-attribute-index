# GAI index benchmark script

`python/benchmarks/benchmark_index.py` measures GAI index builds and indexed
queries for sorted GFF/GFF3 or BED files with an existing TBI or CSI coordinate
index. It reports the input file and type, decoded input size, indexed
attributes, resulting `.gai` size, build time, query mode and text, result
count, and query time. BED always indexes its fourth-column name as `name`;
GFF/GFF3 attribute sets can be varied to compare index size and query time.

For representative performance measurements, build the release Python
extension and run the script from the repository root:

```console
maturin develop --release --manifest-path python/Cargo.toml
python3 python/benchmarks/benchmark_index.py fixtures/gencode_sorted.gff.gz \
  --attribute-set gene_name \
  --attribute-set gene_name,hgnc_id \
  --query-mode exact=brca1 \
  --query-mode prefix=brca \
  --query-mode contains=rca \
  --query-mode 'regex=^brca[0-9]+$' \
  --output benchmark-results.md
```

The script discovers an unambiguous sibling `.tbi` or `.csi` file. Pass
`--coordinate-index path` when the coordinate index is elsewhere or both
formats exist. Each attribute set builds a temporary GAI, so the report
compares each configuration without reusing an index built for a different
set. The input and coordinate index are read only; temporary attribute
indexes are removed when the run finishes.

For BED input, omit `--attribute-set`; GAI indexes its fourth-column name
field. For example:

```console
python3 python/benchmarks/benchmark_index.py fixtures/sorting/gencode_v46.bed.gz \
  --query-mode exact=ENST00000607096.1 \
  --query-mode prefix=ENST000006070 \
  --query-mode contains=607096 \
  --query-mode 'regex=^ENST00000607096\.1$'
```

Use `--query-mode MODE=TEXT` to provide representative input for each mode.
The supported modes are `exact`, `prefix`, `contains` (literal substring),
and `regex` (unanchored Unicode-aware search). Alternatively, repeat `--query
TEXT` to run each query text with all four modes; this is useful when measuring
the modes against the same input text. The report includes the query text so
results remain interpretable.

By default, the script performs one build per attribute set, one query warmup,
and five timed query calls. Build times are medians across `--build-repeats`
builds. Query times are medians across `--repeats` calls after `--warmups`
calls. Build timing covers the Python `build_index` call. Query timing covers
`IndexedSource.query` and excludes the one-time index open and fingerprint
validation. Use `--format csv` for a single machine-readable CSV stream; the
CSV includes the repeat and warmup counts.
