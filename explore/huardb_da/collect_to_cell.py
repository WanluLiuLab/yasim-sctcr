import glob
import os
import re
import hashlib
import pandas as pd

from labw_utils.commonutils.stdlib_helper.parallel_helper import parallel_map
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf

SAMPLE_NAME_REGEX = re.compile("_Volumes_rsch_HuarcBackup_(.+)_all_contig_annotations.json.[0-9]+.parquet")


def run(sample:str):
    flist = glob.glob(f"./pq/_Volumes_rsch_HuarcBackup_{sample}_all_contig_annotations.json.*.parquet")
    merged_lists = []
    for group_key, grouped_df in pd.read_parquet(flist, columns=["chain", "v", "j", "sample_level_barcode", "nt"]).groupby(
            "sample_level_barcode"):
        grouped_df = grouped_df.sort_values(by=["chain"])
        if list(grouped_df["chain"]) != ["TRA", "TRB"]:
            continue
        else:
            tra_chain, trb_chain = grouped_df.itertuples(index=False)
        merged_lists.append({
            "sample_level_barcode": group_key,
            "traj": tra_chain.j,
            "trav": tra_chain.v,
            "trbj": trb_chain.j,
            "trbv": trb_chain.v,
            "sample": sample,
            "nt_blake2b": hashlib.blake2b(bytes(tra_chain.nt+"_"+trb_chain.nt, encoding="UTF-8")).hexdigest()
        })
    pd.DataFrame(merged_lists).to_parquet(os.path.join("pq_cell", f"{sample}.parquet"))


if __name__ == "__main__":
    rm_rf("pq_cell")
    os.makedirs("pq_cell", exist_ok=True)
    flist = glob.glob("./pq/*.parquet")
    samples = set(
        SAMPLE_NAME_REGEX.search(fname).group(1)
        for fname in flist
    )
    parallel_map(
        run,
        samples,
        n_jobs=30,
        backend="loky",
    )
    pd.read_parquet(glob.glob("./pq_cell/*.parquet")).to_parquet("merged_pq_cell.parquet")
