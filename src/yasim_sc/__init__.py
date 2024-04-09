"""
yasim_sc -- Yet Another Single Cell RNA-Seq SIMulator

This is a simulator designed to simulate scRNA-Seq data.
However, what it can only do now is to generate barcodes,
and make LLRGs work in parallel.
It cannot simulate 10x Genomics data.

It is recommended for users to provide their own count matrix.

.. versionadded:: 0.1.0
"""

description = __doc__.splitlines()[1]
try:
    import labw_utils

    _labw_utils_version = labw_utils.__version__

except ImportError as e:
    raise e

from labw_utils import UnmetDependenciesError

try:
    import yasim
except ImportError as e:
    raise UnmetDependenciesError("yasim") from e
