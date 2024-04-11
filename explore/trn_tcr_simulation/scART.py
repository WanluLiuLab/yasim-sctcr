import concurrent.futures
import glob
import os
import shutil
import sys

import pandas as pd

from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.parallel_helper import easyexec
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf

ART_PATH = "/slurm/home/yrd/liulab/yuzhejian/art"
ART_PROFILE_PATH = "/slurm/home/yrd/liulab/yuzhejian/Illumina_profiles"
MANE_FA_PATH = "ref/ens.sel_genes.cdna.fa"
MANE_MAPPING_PATH = "ref/ens.trans_gene_map.tsv"
NJOBS = 40


RLEN_ERRORPROFILE_STRATEGY1 = [
    (50, "EmpR50R1.txt"),
    (100, "HiSeq2kL100R1.txt"),
    (150, "HiSeq2500L150R1.txt"),
    (250, "MiSeqv3L250R1.txt"),
]

def art_salmon(
        colname: str, 
        depth_path: str,
        read_len: int,
        err_profile: str,
        output_dir_name: str,
    ) -> None:
    art_out_prefix = output_dir_name + "_art"
    art_out_dir = art_out_prefix + ".d"
    print(f"ART {colname}", file=sys.stderr)
    easyexec(
        [
            ART_PATH,
            "--qual_file_1",
            os.path.join(ART_PROFILE_PATH, err_profile),
            "--seq_file",
            MANE_FA_PATH,
            "--out_file_prefix",
            art_out_dir,
            "--read_len",
            str(read_len),
            "--fcov",
            depth_path,
            "--is_amplicon",
            "--parallel",
            "-1",
            "--no_sam",
        ],
        log_path=art_out_prefix+".log"
    )
    fq1_path = f"{art_out_prefix}.fq"
    with open(fq1_path, "w") as faw1:
        for fn in glob.glob(f"{art_out_dir}/*.fq"):
            shutil.copyfileobj(open(fn), faw1)

    shutil.rmtree(art_out_dir)


def run_strategy_1():
    _sample_name = "HU_0043_Blood_10x"
    barcodes = [
        *get_tqdm_line_reader(os.path.join("sim", f"{_sample_name}.sim.d", "n_cell_bc.txt")),
        *get_tqdm_line_reader(os.path.join("sim",f"{_sample_name}.sim.d", "t_cell_bc.txt")),
    ]
    os.makedirs(os.path.join("sim", f"{_sample_name}.sim.d", "art_gex_diff_rlen"), exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=NJOBS) as pool:
        for read_len, error_profile in RLEN_ERRORPROFILE_STRATEGY1:
            output_dir_name = os.path.join("sim", f"{_sample_name}.sim.d", "art_gex_diff_rlen", f"sim_gex_rlen{read_len}")
            for barcode in barcodes:
                pool.submit(
                    art_salmon,
                    barcode,
                    os.path.join("sim", f"{_sample_name}.sim.d", "scGEX.depth.d", f"{barcode}.tsv"),
                    read_len,
                    error_profile,
                    os.path.join(output_dir_name, f"{barcode}.d"),
                )

def run_strategy_2():
    for t_cell_num in (100, 500, 1000):
        for replicate_num in range(1, 11):
            with concurrent.futures.ThreadPoolExecutor(max_workers=NJOBS) as pool:
                _sample_name = f"HU_0043_Blood_10x_sim_tcell_only_ncells{t_cell_num}_rep{replicate_num}.d"
                barcodes = [
                    *get_tqdm_line_reader(os.path.join("sim", _sample_name, "n_cell_bc.txt")),
                    *get_tqdm_line_reader(os.path.join("sim", _sample_name, "t_cell_bc.txt")),
                ]
                output_dir_name = os.path.join("sim", _sample_name, "art_sim_gex_rlen150.d")
                for barcode in barcodes:
                    pool.submit(
                        art_salmon,
                        barcode,
                        os.path.join("sim", _sample_name, "scGEX.depth.d", f"{barcode}.tsv"),
                        150,
                        "HiSeq2500L150R1.txt",
                        os.path.join(output_dir_name, f"{barcode}.d"),
                    )


if __name__ == "__main__":
    # run_strategy_1()
    run_strategy_2()

