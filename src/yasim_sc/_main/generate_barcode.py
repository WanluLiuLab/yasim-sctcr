"""
generate_barcode.py -- Generate barcodes with give length and number.

.. versionadded:: 0.1.0
"""

__all__ = ("create_parser", "main")

import argparse
import random
from labw_utils.typing_importer import List

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

logger = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_sc generate_barcode", description=__doc__.splitlines()[1]
    )
    parser.add_argument("-n", "--num_cells", required=True, help="Number of Cells", nargs="?", type=int, action="store")
    parser.add_argument(
        "-l", "--length", required=False, default=17, help="Length of barcode", nargs="?", type=int, action="store"
    )
    parser.add_argument(
        "-o", "--out", required=True, help="Output path to barcode TXT", nargs="?", type=str, action="store"
    )
    return parser


def generate(barcode_length: int) -> str:
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    _generated_barcodes = set()
    while True:
        new_barcode = "".join(map(lambda _: random.choice("ATCG"), range(barcode_length)))
        if new_barcode not in _generated_barcodes:
            _generated_barcodes.add(new_barcode)
            return new_barcode


def main(args: List[str]):
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    args = create_parser().parse_args(args)
    with get_writer(args.out) as writer:
        for _ in tqdm(range(args.num_cells), desc="Generating"):
            writer.write(generate(args.length) + "\n")
