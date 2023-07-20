"""
sample_pcg.py -- Sample Protein-Coding genes for scRNA-Seq Simulation

.. versionadded:: 0.1.0
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse
import random

from labw_utils.bioutils.parser.gtf import GtfIterator, GtfIteratorWriter
from labw_utils.bioutils.record.feature import Feature
from labw_utils.typing_importer import List, Iterable

import pandas as pd

from labw_utils.bioutils.algorithm.sequence import is_valid_chrname
from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.datastructure.gv.gene import DumbGene
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    parser = ArgumentParserWithEnhancedFormatHelp(prog="python -m yasim_scripts sample_pcg",
                                                  description=__doc__.splitlines()[1])
    parser.add_argument(
        '-i',
        '--ncbi_dataset',
        required=True,
        help="The `ncbi_dataset.tsv` retrived from <https://www.ncbi.nlm.nih.gov/labs/data-hub/gene/table/taxon/9606/>",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-g',
        '--gtf',
        required=True,
        help="Path to genome annotation in GTF format to be sampled from",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Path to output GTF",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--num_genes_to_sample',
        required=False,
        help="Number of genes to be sampled",
        nargs='?',
        type=int,
        default=2000,
        action='store'
    )
    parser.add_argument(
        '--gene_name_attribute',
        required=False,
        help="Attribute name that shows Entrez symbol. Should be `gene_id` or `gene_name`.",
        nargs='?',
        type=str,
        default="gene_id",
        action='store'
    )
    return parser


def gtfi(gtf_path: str, gene_name_attribute: str) -> Iterable[Feature]:
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    nw = 0
    nr = 0
    with GtfIterator(gtf_path) as _gtfi:
        for gtf_record in _gtfi:
            nr += 1
            gene_id = gtf_record.attribute_get(gene_name_attribute)
            if (
                    gene_id is not None and
                    not gene_id.startswith("TR")
                    and is_valid_chrname(gtf_record.seqname)
            ):
                yield gtf_record
                nw += 1
    _lh.info(
        "GTF: Loaded %d gtf records out of %d (%.2f%%)",
        nw, nr, nw / nr * 100
    )


def main(args: List[str]):
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    args = create_parser().parse_args(args)
    df = pd.read_table(args.ncbi_dataset)
    pcg_df = (
        df.
        query("`Gene Type` == 'PROTEIN_CODING'").
        query("`Gene Group Method` == 'NCBI Ortholog'")
    )
    ncbi_pcg_gene_ids = set(pcg_df["Symbol"])
    _lh.info(
        "NCBI: Filtered protein-coding %d genes out of %d (%.2f%%)",
        len(ncbi_pcg_gene_ids), len(df), len(pcg_df) / len(df) * 100
    )

    nw = 0
    input_gtf_recoeds = list(gtfi(args.gtf, args.gene_name_attribute))
    gv = DiploidGeneTree.from_feature_iterator(
        input_gtf_recoeds,
        gene_implementation=DumbGene
    )
    all_existing_gene_ids = set(gv.gene_ids)
    _lh.info("GTF: Load %d genes", len(all_existing_gene_ids))

    both_existing_gene_ids = all_existing_gene_ids.intersection(ncbi_pcg_gene_ids)
    if both_existing_gene_ids:
        _lh.info("INTERSECT: Load %d genes", len(both_existing_gene_ids))
    else:
        _lh.error("INTERSECT: Intersection between NCBI and GTF becomes empty!")
        return 1
    rdg = random.SystemRandom()
    sampled_gene_ids = rdg.choices(list(both_existing_gene_ids), k=args.num_genes_to_sample)

    with GtfIteratorWriter(args.out) as gtfw:
        for gene in gv.gene_values:
            if gene.gene_id in sampled_gene_ids:
                gtfw.write(gene)
                nw += 1
                for transcript in gene.transcript_values:
                    gtfw.write(transcript)
                    nw += 1
                    for exon in transcript.exons:
                        gtfw.write(exon)
                        nw += 1

    _lh.info(
        "Filtered %d gtf records out of %d (%.2f%%)",
        nw, len(input_gtf_recoeds), nw / len(input_gtf_recoeds) * 100
    )
