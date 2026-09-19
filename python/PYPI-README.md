# `genomic-attribute-index` Python bindings

The `gni` package wraps the Rust GFF Name Index builder and indexed reader.
It provides deterministic GNI construction, TBI/CSI-backed exact attribute
queries, structured records and metadata, and typed exceptions for stale or
corrupt inputs.

```python
from pathlib import Path
import gni

source = Path("annotations.gff3.gz")
tbi = Path("annotations.gff3.gz.tbi")
gni.build_index(source, tbi, Path("annotations.gff3.gz.gni"), ["Name"])
records = gni.query_index(source, tbi, Path("annotations.gff3.gz.gni"), "BRCA1")
```

See the repository README for the complete API and development instructions.
