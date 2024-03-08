import glob
import json
import os

import pandas as pd

from labw_utils.commonutils.stdlib_helper.parallel_helper import parallel_map
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf


def run(fn: str):
    part_id = 0
    sample_name = (
        fn.replace("raw_data/_Volumes_rsch_HuarcBackup_", "")
        .replace("_all_contig_annotations.json", "")
        .replace("_all_contig_annotations.json", "")
    )
    print(f"SAMPLE: {sample_name} START")
    with open(fn) as reader:
        try:
            fn_ser = json.load(reader)
        except json.decoder.JSONDecodeError:
            print(f"SAMPLE: {sample_name} READ ERR")
            return
    retl = []
    for fn_ser_item in fn_ser:
        if not fn_ser_item["full_length"]:
            continue
        retd = {
            "high_confidence": fn_ser_item["high_confidence"],
            "productive": fn_ser_item["productive"],
            "filtered": fn_ser_item["filtered"],
            "umi_count": fn_ser_item["umi_count"],
            "nt": fn_ser_item["sequence"],
            "quals": fn_ser_item["quals"],
            "chain": None,
            "v": None,
            "d": None,
            "j": None,
            "c": None,
            "v_start": None,
            "v_end": None,
            "d_start": None,
            "d_end": None,
            "j_start": None,
            "j_end": None,
            "c_start": None,
            "c_end": None,
            "cdr3_aa": fn_ser_item["cdr3"],
            "cdr3_nt": fn_ser_item["cdr3_seq"],
            "cdr3_start": fn_ser_item["cdr3_start"],
            "cdr3_stop": fn_ser_item["cdr3_stop"],
            "sample": sample_name
        }
        for annotation in fn_ser_item["annotations"]:
            retd["chain"] = annotation["feature"]["chain"]
            segment_name = annotation["feature"]["gene_name"][3].lower()
            if segment_name not in {"v", "d", "j", "c"}:
                continue
            retd[segment_name] = annotation["feature"]["gene_name"]
            retd[segment_name + "_start"] = annotation["contig_match_start"]
            retd[segment_name + "_end"] = annotation["contig_match_end"]
        retl.append(retd)
        if len(retl) == (1 << 14):
            pd.DataFrame(retl).to_parquet(os.path.join("pq", os.path.basename(fn) + f".{part_id}.parquet"))
            print(f"SAVED PARQUET: {sample_name}:{part_id}")
            retl.clear()
            part_id += 1
    if retl:
        pd.DataFrame(retl).to_parquet(os.path.join("pq", os.path.basename(fn) + f".{part_id}.parquet"))
        retl.clear()
        print(f"SAVED PARQUET: {sample_name}:{part_id}")

    print(f"SAMPLE: {sample_name} FIN")


if __name__ == "__main__":
    rm_rf("pq")
    os.makedirs("pq", exist_ok=True)
    parallel_map(run, glob.glob(
        os.path.join("raw_data", "*_all_contig_annotations.json")), n_jobs=20, backend="multiprocessing")
