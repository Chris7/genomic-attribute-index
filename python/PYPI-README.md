# `genomic-attribute-index` Python bindings

The `gai` package wraps the Rust Genomic Attribute Index compressor, builder, and indexed reader.
It provides BGZF compression with a CSI built in the same pass, deterministic GAI construction,
TBI/CSI-backed exact, prefix, contains, or regex attribute queries, structured records and metadata,
and typed exceptions for stale or corrupt inputs.

```python
from pathlib import Path
import gai

sorted_source = Path("annotations.sorted.gff3")
bgzf_source = Path("annotations.sorted.gff3.gz")
gai.compress(sorted_source, bgzf_source)  # writes annotations.sorted.gff3.gz.csi
csi = Path("annotations.sorted.gff3.gz.csi")
gai.build_index(bgzf_source, csi, Path("annotations.sorted.gff3.gz.gai"), ["Name"])
records = gai.query_index(
    bgzf_source, csi, Path("annotations.sorted.gff3.gz.gai"), "BRCA1"
)
prefix_records = gai.query_index(
    bgzf_source, csi, Path("annotations.sorted.gff3.gz.gai"), "BRCA", match="prefix"
)
contains_records = gai.query_index(
    bgzf_source, csi, Path("annotations.sorted.gff3.gz.gai"), "RCA", match="contains"
)
regex_records = gai.query_index(
    bgzf_source, csi, Path("annotations.sorted.gff3.gz.gai"), r"^BRCA[0-9]+$", match="regex"
)
```

See the repository README for the complete API and development instructions.
