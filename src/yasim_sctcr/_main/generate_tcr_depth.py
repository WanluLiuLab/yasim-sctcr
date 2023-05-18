"""
generate_tcr_depth.py -- Generate scTCR depth TSV without clonal expansion
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse
from labw_utils.typing_importer import List

from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from yasim.helper.depth_io import write_depth
from yasim.helper.frontend import patch_frontend_argument_parser


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_sctcr generate_tcr_depth",
        description=__doc__.splitlines()[1]
    )
    parser = patch_frontend_argument_parser(parser, "-b")
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Path to output depth TSV",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-d',
        '--depth',
        required=True,
        help="Simulated depth",
        nargs='?',
        type=int,
        action='store'
    )
    return parser


def main(args: List[str]) -> int:
    args = create_parser().parse_args(args)
    depth_db = {}
    for barcode in get_tqdm_line_reader(args.barcodes):
        depth_db[f"{barcode}:A"] = args.depth
        depth_db[f"{barcode}:B"] = args.depth

    write_depth(depth_db, args.out, "TRANSCRIPT_ID")
    return 0
