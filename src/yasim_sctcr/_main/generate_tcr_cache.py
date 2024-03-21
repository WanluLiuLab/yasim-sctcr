"""
generate_tcr_cache.py -- Generation of TCR Cache.

.. versionadded:: 0.1.0
"""

__all__ = ("main", "create_parser")

import argparse
import itertools
import json
import sys

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List
from yasim_sctcr.helper.tcr import align

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_sctcr generate_tcr_cache", description=__doc__.splitlines()[1]
    )
    parser.add_argument(
        "--tcr_cdna_fa_path",
        required=True,
        help="TCR cDNA FASTA path. The sequence name should conform TR[AB][VDJC].+ regular expression.",
        type=str,
    )
    parser.add_argument(
        "--tcr_pep_fa_path",
        required=True,
        help="TCR Peptide FASTA path. The sequence name should be same as `tcr_cdna_fa_path`.",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output TCR Cache",
        nargs="?",
        type=str,
        action="store",
    )
    return parser


def create_tcr_cache(
    tcr_cdna_fa_path: str,
    tcr_pep_fa_path: str,
    tcr_cache_path: str,
):
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    tcr_cdna_fav = FastaViewFactory(tcr_cdna_fa_path, show_tqdm=True, read_into_memory=True)
    tcr_pep_fav = FastaViewFactory(tcr_pep_fa_path, show_tqdm=True, read_into_memory=True)
    tcr_genelist = {
        "trav_names": [],
        "traj_names": [],
        "trbv_names": [],
        "trbj_names": [],
        "trbc_names": [],
        "trac_names": [],
    }

    for chr_name in tcr_cdna_fav.chr_names:
        chr_name_head = chr_name[0:4].lower()
        try:
            tcr_genelist[f"{chr_name_head}_names"].append(chr_name)
        except KeyError:
            _lh.error(f"TCR Gene Name {chr_name} invalid!")
            continue

    tcrs = {}
    for gene_name in tqdm(list(itertools.chain(*tcr_genelist.values())), desc="Generating TCR Cache"):
        real_nt_seq = tcr_cdna_fav.sequence(gene_name)
        if real_nt_seq == "":
            _lh.warning("Gene %s have no NT sequence!", gene_name)
            continue
        try:
            real_aa_seq = tcr_pep_fav.sequence(gene_name)
        except KeyError:
            _lh.warning("Gene %s have no AA sequence!", gene_name)
            continue
        tcrs[gene_name] = align(ref_nt_seq=real_nt_seq, imgt_aa_seq=real_aa_seq)
    with get_writer(tcr_cache_path) as writer:
        json.dump(tcrs, writer)

    with get_writer(tcr_cache_path + ".nt.fa") as writer:
        for tcr_name, tcr_tt in tcrs.items():
            writer.write(f">{tcr_name}\n{''.join(list(zip(*tcr_tt))[0])}\n")

    with get_writer(tcr_cache_path + ".aa.fa") as writer:
        for tcr_name, tcr_tt in tcrs.items():
            writer.write(f">{tcr_name}\n{''.join(list(zip(*tcr_tt))[1])}\n")


def main(args: List[str]) -> int:
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    args = create_parser().parse_args(args)
    create_tcr_cache(
        tcr_cdna_fa_path=args.tcr_cdna_fa_path,
        tcr_pep_fa_path=args.tcr_pep_fa_path,
        tcr_cache_path=args.out,
    )
    return 0
