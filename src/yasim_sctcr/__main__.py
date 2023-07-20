"""
TODO: docs

.. versionadded:: 0.1.0
"""

from labw_utils.commonutils.libfrontend import setup_frontend

from yasim_sctcr import __version__
from yasim_sctcr import description

if __name__ == '__main__':
    setup_frontend(
        f"{__package__}._main",
        description,
        __version__
    )
