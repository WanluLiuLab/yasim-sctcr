import concurrent.futures
import glob
import os
import shutil
import subprocess

import pandas as pd

from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf

ART_PATH = "/home/yuzj/Documents/pbsim3_modern/build/art"
SALMON_PATH = "/opt/yuzj/conda/envs/yasim-salmon/bin/salmon"
ART_PROFILE_PATH = "/home/yuzj/Documents/pbsim3_modern/art/Illumina_profiles"
MANE_FA_PATH = "/home/yuzj/Documents/pbsim3_modern/raw_data/MANE.GRCh38.v1.3.refseq_rna.fna"
MANE_SALMON_IDX_PATH = "/home/yuzj/Documents/pbsim3_modern/raw_data/MANE.GRCh38.v1.3.refseq_rna.salmon_idx"
MANE_MAPPING_PATH = "/home/yuzj/Documents/pbsim3_modern/raw_data/MANE_salmon_genemap.tsv"
MANE_SUMMARY_PAH = "/home/yuzj/Documents/pbsim3_modern/raw_data/MANE.GRCh38.v1.3.summary.txt"


def art_salmon(colname: str, sample_name: str) -> None:
    depth_path = os.path.join(f"{sample_name}.sim.d", "depth.d", f"{colname}.tsv")
    out_prefix = os.path.join(f"{sample_name}.sim.d", "gex.d", f"{colname}")
    art_out_prefix = out_prefix + "_art"
    salmon_out_dir = out_prefix + "_salmon.d"
    art_out_dir = art_out_prefix + ".d"
    subprocess.Popen(
        [
            ART_PATH,
            "--qual_file_1",
            os.path.join(ART_PROFILE_PATH, "HiSeq2500L125R1.txt"),
            "--qual_file_2",
            os.path.join(ART_PROFILE_PATH, "HiSeq2500L125R2.txt"),
            "--seq_file",
            MANE_FA_PATH,
            "--out_file_prefix",
            art_out_dir,
            "--read_len",
            "125",
            "--fcov",
            depth_path,
            "--is_pe",
            "--parallel",
            "-1",
            "--pe_frag_dist_std_dev",
            "50",
            "--pe_frag_dist_mean",
            "500",
            "--no_sam",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).wait()
    fq1_path = f"{art_out_prefix}_1.fq"
    fq2_path = f"{art_out_prefix}_2.fq"
    with open(fq1_path, "w") as faw1, open(fq2_path, "w") as faw2:
        for fn in glob.glob(f"{art_out_dir}/*_1.fq"):
            shutil.copyfileobj(open(fn), faw1)
        for fn in glob.glob(f"{art_out_dir}/*_2.fq"):
            shutil.copyfileobj(open(fn), faw2)
    shutil.rmtree(art_out_dir)
    subprocess.Popen(
        [
            SALMON_PATH,
            "--no-version-check",
            "quant",
            "-i",
            MANE_SALMON_IDX_PATH,
            "-l",
            "IU",
            "-1",
            fq1_path,
            "-2",
            fq2_path,
            "-o",
            salmon_out_dir,
            "-p",
            "1",
            "--geneMap",
            MANE_MAPPING_PATH,
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).wait()


def run_sample(_sample_name: str):
    barcodes = [
        *get_tqdm_line_reader(os.path.join(f"{_sample_name}.sim.d", "n_cell_bc.txt")),
        *get_tqdm_line_reader(os.path.join(f"{_sample_name}.sim.d", "t_cell_bc.txt")),
    ]
    rm_rf(os.path.join(f"{_sample_name}.sim.d", "gex.d"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=40) as pool:
        for barcode in barcodes:
            pool.submit(art_salmon, barcode, _sample_name)
    odf = pd.read_csv(MANE_SUMMARY_PAH, sep="\t")[["symbol"]].set_index("symbol")

    for colname in barcodes:
        salmon_csv_path = os.path.join(
            f"{_sample_name}.sim.d",
            "gex.d",
            f"{colname}_salmon.d",
            "quant.genes.sf",
        )
        if not os.path.exists(salmon_csv_path):
            continue
        odf = odf.join(
            pd.read_csv(
                salmon_csv_path,
                sep="\t",
            )[
                ["Name", "NumReads"]
            ].set_index("Name"),
            how="left",
        ).rename(columns={"NumReads": colname})
    odf.reset_index().rename(columns={"index": "FEATURE"}).to_parquet(f"{_sample_name}.art_salmon.parquet")


if __name__ == "__main__":
    for _sample_name in [
        "HU_0196_Kidney_GSE109564",
        "HU_0043_Blood_10x",
        "HU_0148_Decidua_EBI",
        "HU_0223_Muscle_GSE134355",
        "HU_0125_Cerebrospinal-Fluid_GSE134577",
    ]:
        run_sample(_sample_name)
