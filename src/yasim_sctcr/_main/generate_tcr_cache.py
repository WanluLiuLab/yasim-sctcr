"""
generate_tcr_cache.py -- Generation of TCR Cache.

.. versionadded:: 0.1.0
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse
import itertools
import json
from labw_utils.typing_importer import List

from labw_utils.bioutils.algorithm.sequence import is_valid_chrname
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_reader, get_writer
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.frontend import patch_frontend_argument_parser
from yasim_sctcr.helper.tcr import align
from yasim_sctcr._main import get_sample_data_path

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_sctcr generate_tcr_cache",
        description=__doc__.splitlines()[1]
    )
    parser.add_argument(
        '--tcr_genelist_path',
        required=True,
        help=f"TCR Gene List JSON. The IMGT version for human is {get_sample_data_path('tcr_gene_list.min.json')}",
        nargs='?',
        type=str,
        action='store'
    )
    parser = patch_frontend_argument_parser(parser, "-f")
    parser.add_argument(
        '-g',
        '--gtf',
        required=True,
        help="Path to reference genome GTF",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--tcr_aa_table_path',
        required=True,
        help=f"TCR AA Table Path. The IMGT version for human is {get_sample_data_path('tcr_aa_table.min.json')}",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Output TCR Cache",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def create_tcr_cache(
        ref_fa_path: str,
        ref_gtf_path: str,
        tcr_genelist_path: str,
        tcr_aa_table_path: str,
        tcr_cache_path: str
):
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    with get_reader(tcr_genelist_path) as reader:
        tcr_genelist = json.load(reader)
    with get_reader(tcr_aa_table_path) as reader:
        tcr_aa_table = json.load(reader)
    ref_fasta_view = FastaViewFactory(
        ref_fa_path,
        show_tqdm=True,
        read_into_memory=False
    )
    ref_gene_view = GeneViewFactory.from_file(ref_gtf_path)
    tcrs = {}
    for gene_name in tqdm(
            list(itertools.chain(*tcr_genelist.values())),
            desc="Generating TCR Cache"
    ):
        try:
            gene = ref_gene_view.get_gene(gene_name)
        except KeyError:
            _lh.warning("Gene %s not found!", gene_name)
            continue
        valid_transcript = []
        for transcript in gene.iter_transcripts():
            if not is_valid_chrname(transcript.seqname):
                continue
            valid_transcript.append(transcript)
        if len(valid_transcript) != 1:
            _lh.warning(
                "Gene %s have %d != 1 valid transcript (%s)!",
                gene_name,
                len(valid_transcript),
                str(list(map(lambda x: x.transcript_id, valid_transcript)))
            )
            continue
        real_nt_seq = valid_transcript[0].cdna_sequence(ref_fasta_view.sequence)
        if real_nt_seq == "":
            _lh.warning(
                "Gene %s have no NT sequence!",
                gene_name
            )
            continue
        try:
            real_aa_seq = tcr_aa_table[gene_name]
        except KeyError:
            _lh.warning(
                "Gene %s have no AA sequence!",
                gene_name
            )
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
        ref_fa_path=args.fasta,
        ref_gtf_path=args.gtf,
        tcr_genelist_path=args.tcr_genelist_path,
        tcr_aa_table_path=args.tcr_aa_table_path,
        tcr_cache_path=args.out
    )
    return 0
