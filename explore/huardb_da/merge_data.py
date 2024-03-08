import glob
from collections import defaultdict

import pandas as pd

from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_appender
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf

if __name__ == "__main__":
    for filename in tqdm(glob.glob("./pq/*.parquet")):
        sample_name_id = filename.replace("./pq/_Volumes_rsch_HuarcBackup_", "").replace("_all_contig_annotations.json", "").replace(".parquet", "")
        gene_seqs = defaultdict(lambda: [])
        intermediate_seqs = defaultdict(lambda: [])
        for item in pd.read_parquet(filename).itertuples(index=False):
            # for gene_name in "vdjc":
            #     if item.__getattribute__(gene_name) is not None:
            #         gene_symbol = item.__getattribute__(gene_name).replace("/", "")
            #         if gene_symbol.startswith("IG"):
            #             continue
            #         seq = item.nt[int(item.__getattribute__(gene_name+"_start")):int(item.__getattribute__(gene_name+"_end"))]
            #         gene_seqs[item.__getattribute__(gene_name).replace("/", "")].append(seq)
            if item.v_end is not None and item.j_start is not None and item.v_end < item.j_start and item.v.startswith("TR"):
                intermediate_seqs[
                    item.v.replace("/", "") + "_" + item.j.replace("/", "")
                ].append(item.nt[int(item.v_end):int(item.j_start)])
        # for gene_name, seqs in gene_seqs.items():
        #     with get_appender(f"out4msa_samples.nt.fa.d/{gene_name}.fa") as faw:
        #         for seq_id, seq in enumerate(seqs):
        #             faw.write(str(FastaRecord(f"{sample_name_id}.{seq_id}", seq)))
        #             faw.write("\n")
        for gene_name, seqs in intermediate_seqs.items():
            with get_appender(f"out4msa_samples_d.nt.fa.d/{gene_name}.fa") as faw:
                for seq_id, seq in enumerate(seqs):
                    faw.write(str(FastaRecord(f"{sample_name_id}.{seq_id}", seq)))
                    faw.write("\n")
    for filename in tqdm(glob.glob("./out4msa_samples_d.nt.fa.d/*.fa")):
        if len(FastaViewFactory(filename, read_into_memory=True, show_tqdm=False).chr_names) <= 100:
            rm_rf(filename)
    # for filename in tqdm(glob.glob("./out4msa_samples.nt.fa.d/*.fa")):
    #     if len(FastaViewFactory(filename, read_into_memory=True, show_tqdm=False).chr_names) <= 100:
    #         rm_rf(filename)
