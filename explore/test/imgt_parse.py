import glob
import re
from collections import defaultdict

from labw_utils.bioutils.parser.fasta import FastaIterator, FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.importer.tqdm_importer import tqdm

if __name__ == "__main__":
    gene_seqs_aa = defaultdict(lambda: defaultdict(lambda: ""))
    gene_seqs_nt = defaultdict(lambda: defaultdict(lambda: ""))
    with FastaIterator("IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP", full_header=True) as fai:
        for record in fai:
            split_header = record.seq_id.split("|")
            if not split_header[2].startswith("Homo sapiens"):
                continue
            gene_name = split_header[1].split("*")[0]
            gene_name = gene_name.replace("/", "")
            if not gene_name.startswith("TR"):
                continue
            seqname = f"IMGT-{split_header[1]}-{split_header[0]}"
            gene_seqs_aa[gene_name][seqname] += record.sequence
    with FastaIterator("IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP", full_header=True) as fai:
        for record in fai:
            split_header = record.seq_id.split("|")
            if not split_header[2].startswith("Homo sapiens"):
                continue
            gene_name = split_header[1].split("*")[0]
            gene_name = gene_name.replace("/", "")
            if not gene_name.startswith("TR"):
                continue
            seqname = f"IMGT-{split_header[1]}-{split_header[0]}"
            gene_seqs_nt[gene_name][seqname] += record.sequence
    with FastaIterator("Homo_sapiens.GRCh38.cdna.all.fa", full_header=True) as fai:
        for record in fai:
            split_header = record.seq_id.split(" ")
            accn = split_header[0]
            gene_symbol = None
            gene_biotype = None
            for k_v in split_header[1:]:
                kvs = k_v.split(":")
                if kvs[0] == "gene_biotype":
                    gene_biotype = kvs[1]
                if kvs[0] == "gene_symbol":
                    gene_symbol = kvs[1]
            if gene_biotype is None or not gene_biotype.startswith("TR"):
                continue
            if gene_symbol is None:
                continue
            seqname = f"ENSG-{accn}"
            gene_seqs_nt[gene_symbol][seqname] += record.sequence
    with FastaIterator("Homo_sapiens.GRCh38.pep.all.fa", full_header=True) as fai:
        for record in fai:
            split_header = record.seq_id.split(" ")
            accn = split_header[0]
            gene_symbol = None
            gene_biotype = None
            for k_v in split_header[1:]:
                kvs = k_v.split(":")
                if kvs[0] == "gene_biotype":
                    gene_biotype = kvs[1]
                if kvs[0] == "gene_symbol":
                    gene_symbol = kvs[1]
            if gene_biotype is None or not gene_biotype.startswith("TR"):
                continue
            if gene_symbol is None:
                continue
            seqname = f"ENSG-{accn}"
            gene_seqs_aa[gene_symbol][seqname] += record.sequence

    UNIPROT_GENE_NAME_REGEX = re.compile(r" GN=(\S+) ")
    with FastaIterator("uniprot_sprot.fasta", full_header=True) as fai:
        for record in fai:
            if " OX=9606 " not in record.seq_id:
                continue
            match_result = UNIPROT_GENE_NAME_REGEX.search(record.seq_id)
            if match_result is not None:
                gene_name = match_result.group(1)
                if not gene_name.startswith("TR"):
                    continue
            else:
                continue
            accn = "-".join(record.seq_id.split(" ")[0].split("|")[1:2])
            seqname = f"SP-{accn}"
            gene_seqs_aa[gene_name][seqname] += record.sequence

    with FastaIterator("igblast.fa", full_header=False) as fai:
        for record in fai:
            split_header = record.seq_id.split("*")
            gene_name = split_header[0].replace("/", "")
            if not gene_name.startswith("TR"):
                continue
            seqname = f"IGBLAST-{record.seq_id}"
            gene_seqs_aa[gene_name][seqname] += record.sequence

    for gene_name, seqs in gene_seqs_nt.items():
        if len(set(map(lambda x: x.split("-")[0], seqs.keys()))) > 1:
            with FastaWriter(f"out4msa.nt.fa.d/{gene_name}.fa") as faw:
                for seq_id, seq in seqs.items():
                    faw.write(FastaRecord(seq_id, seq))
    for gene_name, seqs in gene_seqs_aa.items():
        if len(set(map(lambda x: x.split("-")[0], seqs.keys()))) > 1:
            with FastaWriter(f"out4msa.aa.fa.d/{gene_name}.fa") as faw:
                for seq_id, seq in seqs.items():
                    faw.write(FastaRecord(seq_id, seq))

