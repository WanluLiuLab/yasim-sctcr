import glob
import os
import re

import pandas as pd

from labw_utils.commonutils.stdlib_helper.parallel_helper import parallel_map
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf

SAMPLE_NAME_REGEX = re.compile("_Volumes_rsch_HuarcBackup_(.+)_all_contig_annotations.json.[0-9]+.parquet")


def run(sample: str):
    flist = glob.glob(f"./pq/_Volumes_rsch_HuarcBackup_{sample}_all_contig_annotations.json.*.parquet")
    merged_lists = []
    for tup in pd.read_parquet(flist).itertuples(index=False):
        if tup.stop_codon_pos is None:
            continue
        retd = {}
        for gene_name in ("v", "d", "j", "c"):
            if tup.__getattribute__(f"{gene_name}") is not None:
                retd[gene_name] = (
                    tup.__getattribute__(f"{gene_name}_start")
                    <= tup.stop_codon_pos
                    <= tup.__getattribute__(f"{gene_name}_end")
                )
            else:
                retd[gene_name] = None
        merged_lists.append(retd)
    pd.DataFrame(merged_lists).to_parquet(os.path.join("pq_sc", f"{sample}.parquet"))


if __name__ == "__main__":
    rm_rf("pq_sc")
    os.makedirs("pq_sc", exist_ok=True)
    flist = glob.glob("./pq/*.parquet")
    samples = set(SAMPLE_NAME_REGEX.search(fname).group(1) for fname in flist)
    parallel_map(
        run,
        samples,
        n_jobs=30,
        backend="loky",
    )
    pd.read_parquet(glob.glob("./pq_sc/*.parquet")).to_parquet("merged_pq_sc.parquet")
