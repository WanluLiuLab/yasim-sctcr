"""
scaffold.py -- Generate paired scRNA-Seq and scTCR-Seq data from configuration files or options

.. versionadded:: 0.1.1
"""

import argparse
import os
import random
import re
from collections import defaultdict

import pandas as pd

from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
from labw_utils.mlutils.ndarray_helper import describe
from labw_utils.typing_importer import List
from yasim_sc.helper.rna_seq import generate_barcodes

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    TODO: docs

    .. versionadded:: 0.1.1
    """
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_sctcr scaffold",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "--transcript_gene_mapping",
        help="A headless two-column TSV whose first column are transcript IDs and second column are gene IDs.",
        required=True,
    )
    parser.add_argument(
        "--src_sc_data",
        help="A TSV or Apache Parquet file, one cell per column, with gene names on 'FEATURE' column. "
        "The column names should contain cell type information "
        "that allows distinguish between T-cells and other cells using regular expression. "
        "Notice that if you want to read scRNA-Seq files in Apache Parquet format, "
        "you need to install Apache Arrow or FastParquet."
        "If you're using anndata, please convert it using `convert_anndata`.",
        required=True,
    )
    parser.add_argument(
        "--n_genes",
        help="Number of genes to sample.",
        default=500,
        type=int,
        required=False,
    )
    parser.add_argument(
        "--mean_depth",
        help="Mean sequencing depth.",
        default=1.0,
        type=float,
        required=False,
    )
    parser.add_argument(
        "--out",
        help="Path to output directory.",
        required=True,
    )
    parser.add_argument(
        "--t_cell_regex",
        help="Regular expression for T-cells name",
        required=True,
    )
    return parser


def main(args: List[str]) -> None:
    argv = create_parser().parse_args(args)
    out_dir = argv.out
    rm_rf(out_dir)
    os.makedirs(out_dir)
    g_t_map = defaultdict(lambda: [])
    t_g_map = {}

    for tup in pd.read_csv(argv.transcript_gene_mapping, sep="\t", names=["transcript_id", "gene_id"]).itertuples(
        index=False
    ):
        g_t_map[tup.gene_id].append(tup.transcript_id)
        t_g_map[tup.transcript_id] = tup.gene_id
    _lh.info(
        "Got %d genes with %d transcripts",
        len(g_t_map),
        len(t_g_map),
    )
    _lh.info("Generating gene expression")
    if argv.src_sc_data.endswith(".parquet"):
        df = pd.read_parquet(argv.src_sc_data)
    elif argv.src_sc_data.endswith(".tsv"):
        df = pd.read_csv(argv.src_sc_data, sep="\t")
    else:
        _lh.error("Unknown scRNA-Seq format.")
        exit(1)
    df = df.set_index("FEATURE")
    _lh.info("Read sample with %d cells and %d genes", len(df.columns), len(df))
    df = df.loc[list(map(lambda x: x in g_t_map.keys(), df.index)), :]
    _lh.info("%d genes left after filtering unknown genes", len(df))
    df = df.sample(n=argv.n_genes)
    _lh.info("Sampled to 500 genes. Current data distribution: %s", describe(df.to_numpy()))
    df /= df.to_numpy().mean() * argv.mean_depth
    _lh.info("Re-scale to %.2f. Current data distribution: %s", argv.mean_depth, describe(df.to_numpy()))
    odf = df.reset_index()
    try:
        odf.to_parquet(os.path.join(out_dir, "sim_dw_sampled.parquet"), index=False)
    except (ImportError, AttributeError):
        odf.to_csv(os.path.join(out_dir, "sim_dw_sampled.parquet"), index=False, sep="\t")
    gene_ids_selected = df.index
    # In this step, we will select 1 isoform per gene.
    rng = random.SystemRandom()
    g_t_map_sel = {
        gene_id: rng.choice(transcript_ids)
        for gene_id, transcript_ids in g_t_map.items()
        if gene_id in gene_ids_selected
    }
    df.index = [g_t_map_sel[gene_id] for gene_id in df.index]

    col_names_t = []
    col_names_n = []
    t_cell_regex = re.compile(argv.t_cell_regex)
    for col_name in df.columns:
        (col_names_t if t_cell_regex.match(col_name) != None else col_names_n).append(col_name)

    num_t_cells = len(col_names_t)
    num_n_cells = len(col_names_n)
    _lh.info("Got %d T-cells vs. %d normal cells", num_t_cells, num_n_cells)

    barcodes = list(generate_barcodes(15, num_t_cells + num_n_cells))
    barcodes_t = barcodes[:num_t_cells]
    barcodes_n = barcodes[num_t_cells:]

    with get_writer(os.path.join(out_dir, "t_cell_bc.txt")) as w:
        for barcode in barcodes[:num_t_cells]:
            w.write(barcode)
            w.write("\n")
    with get_writer(os.path.join(out_dir, "n_cell_bc.txt")) as w:
        for barcode in barcodes[num_t_cells:]:
            w.write(barcode)
            w.write("\n")
    depth_d = os.path.join(out_dir, "scGEX.depth.d")
    os.makedirs(depth_d)
    for barcode, col_name in zip(barcodes_t, col_names_t):
        df[[col_name]].rename(columns={col_name: "DEPTH"}).to_csv(
            os.path.join(depth_d, f"{barcode}.tsv"),
            index_label="TRANSCRIPT_ID",
            sep="\t",
        )
    for barcode, col_name in zip(barcodes_n, col_names_n):
        df[[col_name]].rename(columns={col_name: "DEPTH"}).to_csv(
            os.path.join(depth_d, f"{barcode}.tsv"),
            index_label="TRANSCRIPT_ID",
            sep="\t",
        )
