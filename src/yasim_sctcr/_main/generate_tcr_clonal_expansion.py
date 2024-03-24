"""
generate_tcr_clonal_expansion.py -- Generate ground-truth TCR Contigs

.. versionadded:: 0.1.1
"""

__all__ = ("main", "create_parser")

import argparse

import numpy as np
import numpy.typing as npt
import pandas as pd

from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.appender import load_table_appender_class, TableAppenderConfig
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.mlutils.ndarray_helper import describe
from labw_utils.typing_importer import List
from yasim.helper.frontend import patch_frontend_argument_parser

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    TODO: docs

    .. versionadded:: 0.1.1
    """
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_sctcr generate_tcr_clonal_expansion",
        description=__doc__.splitlines()[1],
    )
    parser = patch_frontend_argument_parser(parser, "-b")
    parser.add_argument(
        "-i",
        "--src_tcr_stats_tsv",
        help="Path to simulated TCR status TSV.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--dst_nt_fasta",
        help="Output nucleotide FASTA for each barcode",
        required=True,
    )
    parser.add_argument(
        "--alpha",
        required=False,
        help="Zipf's Coefficient, larger for larger differences",
        nargs="?",
        type=int,
        action="store",
        default=1,
    )
    return parser


def zipf_yuzj_y4p(nspecies: int, ntotal: int, alpha: float) -> npt.NDArray:
    """
    From the following R code:

    zipf <- function(nspecies, ntotal, alpha = 1) {
        final_num <- seq(1, nspecies, 1)^(-alpha)
        final_num <- final_num - min(final_num)
        final_num <- as.integer(
            final_num * (ntotal - nspecies) / sum(final_num)
        ) + 1
        return(final_num)
    }
    """
    final_num = np.arange(1, nspecies + 1, dtype=float) ** (-alpha)
    final_num -= np.min(final_num)
    final_num = np.array(final_num * (ntotal - nspecies) / np.sum(final_num) + 1, dtype=int)
    return final_num.tolist()


def main(args: List[str]):
    """
    TODO: docs

    .. versionadded:: 0.1.1
    """
    argv = create_parser().parse_args(args)
    df = pd.read_csv(argv.src_tcr_stats_tsv, sep="\t", quotechar="'")[["UUID", "ALPHA_NT", "BETA_NT"]]
    n_tcr = len(df)
    alpha = argv.alpha
    barcodes = list(get_tqdm_line_reader(argv.barcodes))
    n_t_cell = len(barcodes)
    _lh.info("Generating %d T-Cells from %s TCRs.", n_t_cell, n_tcr)

    clon_size = np.sort(zipf_yuzj_y4p(n_tcr, n_t_cell, alpha))[::-1]
    while clon_size.sum() < n_t_cell:
        clon_size[0] += 1
    while clon_size.sum() > n_t_cell:
        clon_size[0] -= 1
    _lh.info("Generated: " + describe(clon_size))

    barcode_uuid_map = []
    for ii, i in enumerate(clon_size):
        barcode_uuid_map.extend(([ii] * i.item()))

    with FastaWriter(argv.dst_nt_fasta) as faw, load_table_appender_class("TSVTableAppender")(
        filename=argv.dst_nt_fasta + ".stats",
        header=["UUID", "barcode"],
        tac=TableAppenderConfig(buffer_size=1024),
    ) as appender:
        for barcode_id, tcr_id in enumerate(barcode_uuid_map):
            tcr_uuid = df.iloc[tcr_id, 0]
            barcode = barcodes[barcode_id]
            appender.append((tcr_uuid, barcode))
            faw.write(FastaRecord(f"{barcode}_A", df.iloc[tcr_id, 1]))
            faw.write(FastaRecord(f"{barcode}_B", df.iloc[tcr_id, 2]))
