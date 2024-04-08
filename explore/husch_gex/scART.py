import concurrent.futures
import glob
import os
import shutil

import pandas as pd

from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.parallel_helper import easyexec
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf

ART_PATH = "/slurm/home/yrd/liulab/yuzhejian/art"
SALMON_PATH = "/slurm/home/yrd/liulab/yuzhejian/conda/envs/yasim-salmon/bin/salmon"
ART_PROFILE_PATH = "/slurm/home/yrd/liulab/yuzhejian/Illumina_profiles"
MANE_FA_PATH = "ref/MANE.GRCh38.v1.3.ensembl_rna.fna"
MANE_SALMON_IDX_PATH = "ref/MANE.GRCh38.v1.3.refseq_rna.salmon_idx"
MANE_MAPPING_PATH = "ref/hgnc_salmon_genemap.tsv"
NJOBS = 10


def art_salmon(colname: str, sample_name: str) -> None:
    depth_path = os.path.join(f"{sample_name}.sim.d", "depth.d", f"{colname}.tsv")
    out_prefix = os.path.join(f"{sample_name}.sim.d", "gex.d", f"{colname}")
    art_out_prefix = out_prefix + "_art"
    salmon_out_dir = out_prefix + "_salmon.d"
    art_out_dir = art_out_prefix + ".d"
    print(f"ART {colname}")
    easyexec(
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
    )
    fq1_path = f"{art_out_prefix}_1.fq"
    fq2_path = f"{art_out_prefix}_2.fq"
    with open(fq1_path, "w") as faw1, open(fq2_path, "w") as faw2:
        for fn in glob.glob(f"{art_out_dir}/*_1.fq"):
            shutil.copyfileobj(open(fn), faw1)
        for fn in glob.glob(f"{art_out_dir}/*_2.fq"):
            shutil.copyfileobj(open(fn), faw2)
    shutil.rmtree(art_out_dir)
    print(f"SALMON {colname}")
    easyexec(
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
    )
    print(f"FIN {colname}")


def run_sample(_sample_name: str):
    barcodes = [
        *get_tqdm_line_reader(os.path.join(f"{_sample_name}.sim.d", "n_cell_bc.txt")),
        *get_tqdm_line_reader(os.path.join(f"{_sample_name}.sim.d", "t_cell_bc.txt")),
    ]
    rm_rf(os.path.join(f"{_sample_name}.sim.d", "gex.d"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=NJOBS) as pool:
        for barcode in barcodes:
            pool.submit(art_salmon, barcode, _sample_name)
    odf = pd.read_csv(MANE_MAPPING_PATH, sep="\t", names=["ensembl", "symbol"])[["symbol"]].set_index("symbol")

    for colname in barcodes:
        salmon_dir_path = os.path.join(f"{_sample_name}.sim.d", "gex.d", f"{colname}_salmon.d")
        salmon_csv_path = os.path.join(salmon_dir_path, "quant.genes.sf")
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
        rm_rf(salmon_dir_path)
    odf.reset_index().rename(columns={"symbol": "FEATURE"}).to_parquet(f"parquets/{_sample_name}.art_salmon.parquet")


if __name__ == "__main__":
    if os.path.exists(MANE_SALMON_IDX_PATH):
        easyexec([SALMON_PATH, "index", "-t", MANE_FA_PATH, "-i", MANE_SALMON_IDX_PATH, "-p", str(NJOBS)])
    for _sample_name in [
        # "HU_0196_Kidney_GSE109564",
        # "HU_0043_Blood_10x",
        # "HU_0148_Decidua_EBI",
        # "HU_0223_Muscle_GSE134355",
        # "HU_0125_Cerebrospinal-Fluid_GSE134577",
    ]:
        rm_rf(os.path.join(f"{_sample_name}.sim.d", "gex.d"))
        run_sample(_sample_name)
