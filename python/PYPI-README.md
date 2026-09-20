# `genomic-attribute-index` Python bindings

The `gai` package wraps the Rust Genomic Attribute Index builder and indexed reader.
It provides deterministic GAI construction, TBI/CSI-backed exact or prefix
attribute queries, structured records and metadata, and typed exceptions for
stale or corrupt inputs.

```python
from pathlib import Path
import gai

source = Path("annotations.gff3.gz")
tbi = Path("annotations.gff3.gz.tbi")
gai.build_index(source, tbi, Path("annotations.gff3.gz.gai"), ["Name"])
records = gai.query_index(source, tbi, Path("annotations.gff3.gz.gai"), "BRCA1")
prefix_records = gai.query_index(
    source, tbi, Path("annotations.gff3.gz.gai"), "BRCA", match="prefix"
)
```

See the repository README for the complete API and development instructions.
