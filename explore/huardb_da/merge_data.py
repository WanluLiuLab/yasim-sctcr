import glob
from collections import defaultdict

import pandas as pd

from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.parser.fastq import FastqWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.bioutils.record.fastq import FastqRecord
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio.safe_io import get_appender, get_writer

if __name__ == "__main__":
    pd.read_parquet(glob.glob("./pq/*.parquet")).to_parquet("merged.parquet")
    with FastaWriter("all_contigs.fa") as contigs_faw, get_writer(
        "all_contigs.cellranger.bed", is_binary=False
    ) as contig_annotw, FastqWriter("all_contigs.fq") as contigs_fqw:
        for filename in tqdm(glob.glob("./pq/*.parquet")):
            sample_name_id = (
                filename.replace("./pq/_Volumes_rsch_HuarcBackup_", "")
                .replace("_all_contig_annotations.json", "")
                .replace(".parquet", "")
            )
            gene_seqs = defaultdict(lambda: [])
            gene_quals = defaultdict(lambda: [])
            for item_id, item in enumerate(pd.read_parquet(filename).itertuples(index=False)):
                fasta_header = f"{sample_name_id}.{item_id}"
                contigs_faw.write(FastaRecord(fasta_header, item.nt))
                contigs_fqw.write(FastqRecord(fasta_header, item.nt, item.quals))
                for segment_name in "vdjc":
                    gene_name = item.__getattribute__(segment_name)
                    if item.__getattribute__(segment_name) is not None:
                        gene_start = int(item.__getattribute__(segment_name + "_start"))
                        gene_end = int(item.__getattribute__(segment_name + "_end"))
                        contig_annotw.write(
                            "\t".join(
                                [
                                    fasta_header,
                                    str(gene_start),
                                    str(gene_end),
                                    gene_name,
                                    ".",
                                    "+",
                                ]
                            )
                        )
                        contig_annotw.write("\n")
                        seq = item.nt[gene_start:gene_end]
                        qual = item.quals[gene_start:gene_end]
                        gene_seqs[gene_name].append(seq)
                        gene_quals[gene_name].append(qual)

            for gene_name, seqs in gene_seqs.items():
                with get_appender(f"out4msa_samples.nt.fa.d/{gene_name}.fa") as faw, get_appender(
                    f"out4msa_samples.nt.fq.d/{gene_name}.fq"
                ) as fqw:
                    for seq_id, seq in enumerate(seqs):
                        faw.write(str(FastaRecord(f"{sample_name_id}.{seq_id}", seq)))
                        faw.write("\n")
                        fqw.write(str(FastqRecord(f"{sample_name_id}.{seq_id}", seq, gene_quals[gene_name][seq_id])))
                        fqw.write("\n")
