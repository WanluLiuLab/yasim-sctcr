__all__ = (
    "get_sample_data_path",
)

import os

_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
_SAMPLE_DATA_PATH = os.path.join(_DIR_PATH, "sample_data")


def get_sample_data_path(name: str) -> str:
    return os.path.join(_SAMPLE_DATA_PATH, name)
