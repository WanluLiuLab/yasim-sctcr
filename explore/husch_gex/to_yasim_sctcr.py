import os

import pandas as pd

from labw_utils.commonutils.lwio.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
from yasim_sc.helper.rna_seq import generate_barcodes

sample_names = [
    "HU_0196_Kidney_GSE109564",
    "HU_0043_Blood_10x",
    "HU_0148_Decidua_EBI",
    "HU_0223_Muscle_GSE134355",
    "HU_0125_Cerebrospinal-Fluid_GSE134577",
]

MANE_SUMMARY_PAH = "MANE.GRCh38.v1.3.summary.txt"

if __name__ == "__main__":
    mane_pd = pd.read_csv(MANE_SUMMARY_PAH, sep="\t")[["symbol", "RefSeq_nuc"]].set_index("symbol")
    for sample_name in sample_names:
        print(sample_name)
        df = (
            mane_pd.join(pd.read_parquet(f"{sample_name}_sim.parquet").set_index("FEATURE"), how="inner")
            .sample(n=500)
            .set_index("RefSeq_nuc")
        )
        dfm = df.to_numpy()
        df /= dfm.mean() * 0.2
        df.reset_index().rename(columns={"RefSeq_nuc": "FEATURE"}).to_parquet(f"{sample_name}_sim_dw_sampled.parquet")
        rm_rf(f"{sample_name}.sim.d")
        os.makedirs(f"{sample_name}.sim.d")

        col_names_t = []
        col_names_n = []
        for col_name in df.columns:
            (col_names_t if col_name.split("_")[0] in {"CD4T", "CD8T"} else col_names_n).append(col_name)

        num_t_cells = len(col_names_t)
        num_n_cells = len(col_names_n)

        barcodes = list(generate_barcodes(15, num_t_cells + num_n_cells))
        barcodes_t = barcodes[:num_t_cells]
        barcodes_n = barcodes[num_t_cells:]

        with get_writer(os.path.join(f"{sample_name}.sim.d", "t_cell_bc.txt")) as w:
            for barcode in barcodes[:num_t_cells]:
                w.write(barcode)
                w.write("\n")
        with get_writer(os.path.join(f"{sample_name}.sim.d", "n_cell_bc.txt")) as w:
            for barcode in barcodes[num_t_cells:]:
                w.write(barcode)
                w.write("\n")
        os.makedirs(os.path.join(f"{sample_name}.sim.d", "depth.d"))
        for barcode, col_name in zip(barcodes_t, col_names_t):
            df[[col_name]].rename(columns={col_name: "DEPTH"}).to_csv(
                os.path.join(f"{sample_name}.sim.d", "depth.d", f"{barcode}.tsv"),
                index_label="TRANSCRIPT_ID",
                sep="\t",
            )
        for barcode, col_name in zip(barcodes_n, col_names_n):
            df[[col_name]].rename(columns={col_name: "DEPTH"}).to_csv(
                os.path.join(f"{sample_name}.sim.d", "depth.d", f"{barcode}.tsv"),
                index_label="TRANSCRIPT_ID",
                sep="\t",
            )
