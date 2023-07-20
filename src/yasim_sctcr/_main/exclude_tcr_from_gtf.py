"""
exclude_tcr_from_gtf.py -- Filter out gene names started with TR.

For GENCODE Reference Only.

.. versionadded:: 0.1.0
"""

from labw_utils.typing_importer import List

from labw_utils.bioutils.parser.gtf import GtfIterator, GtfIteratorWriter


def main(args: List[str]):
    for arg in args:
        with GtfIteratorWriter(arg + ".tcr_filtered.gtf") as writer:
            for record in GtfIterator(arg, show_tqdm=True):
                if not record.attribute_get("gene_name", "TRBV").startswith("TR"):
                    writer.write_feature(record)
