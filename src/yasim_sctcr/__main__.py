from labw_utils.commonutils.libfrontend import setup_frontend

from yasim import __version__
from yasim_sctcr import description

if __name__ == '__main__':
    setup_frontend(
        "yasim_sctcr._main",
        description,
        __version__
    )
