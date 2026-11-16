"""FlowCellar: quality control for Oxford Nanopore sequencing runs."""
from flowcellar._version import __version__

#: Version of the bundle schema this package writes (src/flowcellar/model/schemas/).
SCHEMA_VERSION = "0.1.1"
#: Layout version of a library folder (the LIBRARY marker file).
LAYOUT_VERSION = 1

__all__ = ["__version__", "SCHEMA_VERSION", "LAYOUT_VERSION"]
