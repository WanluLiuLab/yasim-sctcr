"""
yasim_sctcr -- Yet Another Single-Cell TCR-Seq Simulator

.. versionadded:: 0.1.0
"""

__version__ = "0.1.1"

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

# Pandas and Numpy are required by yasim.
