"""ont-guide-aligner: count CRISPR sgRNAs in Oxford Nanopore reads, one guide per read.

Use it from the command line (`ont-guide-aligner --help`). The Python modules are internal and
may change between releases.
"""

from importlib.metadata import PackageNotFoundError, version as _version

try:
    __version__ = _version("ont-guide-aligner")  # set once, in Cargo.toml
except PackageNotFoundError:  # running from an unbuilt source checkout
    __version__ = "unknown"

__all__ = ["__version__"]
