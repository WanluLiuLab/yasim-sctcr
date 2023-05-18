from labw_utils.commonutils.libfrontend import setup_frontend

from yasim import \
    __version__ as yasim_ver
from yasim_sc import description

if __name__ == '__main__':
    setup_frontend(
        "yasim_sc._main",
        description,
        yasim_ver,
    )
