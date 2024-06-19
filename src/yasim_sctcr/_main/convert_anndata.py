"""
convert_anndata -- Convert Scanpy AnnData object to TSV/Apache Parquet for ``scaffold``

.. versionadded :: 1.0.0

"""

__all__ = ("main", "create_parser")

import argparse
import os
import sys

from labw_utils import UnmetDependenciesError
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

if os.environ.get("LABW_UTILS_UNDER_PYTEST", None) is not None:
    import pytest

    try:
        fp = pytest.importorskip("fastparquet")
    except AttributeError:
        # Error in Numba under PyPy 3.7
        pytest.skip(allow_module_level=True)
else:
    pytest = None
    try:
        import anndata as ad
    except ImportError as e:
        raise UnmetDependenciesError("anndata") from e

import pandas as pd

from labw_utils.typing_importer import List

_lh = get_logger(__name__)

SUPPORTING_EXTENSIONS = (
    "h5ad",
    "csv",
    "txt",
    "tab",
    "data",
    "loom",
    "mtx",
    "gz",
    "npz",
    "zarr",
)


def read_anndata(filename: str) -> ad.AnnData:
    """
    Read various AnnData supporting formats.

    .. author:: WU JX
    """
    file_format = filename.split(".")[-1].lower()

    if file_format == "h5ad":
        adata = ad.read_h5ad(filename)
    elif file_format == "csv":
        adata = ad.read_csv(filename)
    elif file_format in ["txt", "tab", "data"]:
        adata = ad.read_text(filename)
    elif file_format == "loom":
        adata = ad.read_loom(filename)
    elif file_format == "mtx":
        adata = ad.read_mtx(filename)
    elif file_format in ["gz", "npz"]:
        adata = ad.read_umi_tools(filename)
    elif file_format == "zarr":
        adata = ad.read_zarr(filename)
    else:
        _lh.error(
            "Error: Unsupported file format '%s'. Should be in one of '%s'",
            file_format,
            " ".join(SUPPORTING_EXTENSIONS),
        )
        sys.exit(1)
    return adata


def create_dataframe_from_anndata(adata: ad.AnnData, column_name: str) -> pd.DataFrame:
    """
    Create a dataframe from an AnnData object.

    .. author:: WU JX

    :param adata: AnnData object
    :param column_name: Column name for cell type annotation defined in ``adata.obs``.
    """
    if not hasattr(adata, "obs") or not hasattr(adata, "var"):
        raise ValueError("adata must be an anndata object")

    adata.obs["annotation"] = adata.obs[column_name].astype(str)

    gene_ids = adata.var_names

    adata.obs["annotation_count"] = adata.obs.groupby("annotation").cumcount() + 1
    annotations = adata.obs["annotation"] + adata.obs["annotation_count"].astype(str)

    df = pd.DataFrame(adata.X.T, index=gene_ids, columns=annotations)
    return df


def create_parser() -> argparse.ArgumentParser:
    """
    TODO: docs

    .. versionadded:: 0.1.0
    """
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_sctcr convert_anndata.",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "--adata",
        required=True,
        help="scRNA-Seq file readable by anndata",
        type=str,
    )
    parser.add_argument(
        "--annotation_column_name",
        required=True,
        help="Column name for cell type annotation defined in ``adata.obs``.",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output file, should be in Apache Parquet or TSV format. "
        "Notice that if you want to read scRNA-Seq files in Apache Parquet format, "
        "you need to install Apache Arrow or FastParquet.",
        type=str,
    )
    return parser


def main(args: List[str]) -> int:
    argv = create_parser().parse_args(args)
    _lh.info("Reading anndata...")
    adata = ad.read(argv.adata)
    _lh.info("Read %d cells for %d genes", *adata.shape)
    df = create_dataframe_from_anndata(adata, argv.annotation_column_name)
    if argv.out.endswith(".parquet"):
        df.to_parquet(argv.src_sc_data, index=False)
    elif argv.out.endswith(".tsv"):
        df.to_csv(argv.out, sep="\t", index=False)
    else:
        _lh.error("Unknown output scRNA-Seq format.")
        return 1
    return 0
