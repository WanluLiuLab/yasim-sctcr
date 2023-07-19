---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Tutorial on Simulation of scTCR-Seq Data

This is a tutorial on how to generate scTCR-Seq data from reference genome.

```{warning}
This version of scTCR-Seq simulator **does not** support clonal expansion or generation of raw FASTQ reads from manufacturers like 10xGenomics, etc.
```

Download of real statistical data:

- [cdr3_deletion_table.min.json.xz](../data/cdr3_deletion_table.min.json.xz): Human CDR3 deletion table.
- [cdr3_insertion_table.min.json.xz](../data/cdr3_insertion_table.min.json.xz): Human CDR3 insertion table.
- [tcr_cache.min.json.xz](../data/tcr_cache.min.json.xz): TCR Cache built against [UCSC hg38](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz) reference genome sequence and [UCSC NCBIRefSeq](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz) reference genome annotation.
- [tcr_genelist.min.json.xz](../data/tcr_genelist.min.json.xz): Human TCR gene list curated from IMGT.
- [IMGT_Protein_Display.min.json.xz](../data/IMGT_Protein_Display.min.json.xz): Human TCR reference Amino-Acid (AA) Sequence curated from IMGT.

+++

## Preparation

Following code retrieves reference genome sequence and annotations used in this example. Since TCRs are only located at chromosome 14 and 7, only those 2 references are used.

```{code-cell}
:tags: [skip-execution]

%%bash
axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
gzip -dcf hg38.ncbiRefSeq.gtf.gz \
    | grep -e '^chr7\s' -e '^chr14\s' \
    > hg38.ncbiRefSeq_chr7_14.gtf

axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

Generate barcodes for 10 cells.

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim_sc generate_barcode -n 10 -o barcode.txt
```

Generate TCR depth. Following code generates a depth file with 20 depth for each barcode.

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim_sctcr generate_tcr_depth \
    -b barcode.txt \
    -o tcr_depth.tsv \
    -d 20
```

Generates TCR cache. The TCR cache is a JSON containing TCR Nucleotide (NT) to AA alignment information.

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim_sctcr generate_tcr_cache \
    --tcr_genelist_path ../data/tcr_genelist.min.json.xz \
    -o tcr_cache.json \
    --tcr_aa_table_path ../data/IMGT_Protein_Display.min.json.xz \
    -f hg38.fa \
    -g hg38.ncbiRefSeq_chr7_14.gtf
```

This generates following files:

- `tcr_cache.json`, the generated TCR cache. User may not use this file.
- `tcr_cache.json.aa.fa` and `tcr_cache.json.nt.fa`, NT and AA sequence of generated TCR cache. User may not use this file.

+++

Simulate ground-truth TCR contigs. This step may report `finish with N failures` -- don't worry! None of your cells would lost.

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim_sctcr rearrange_tcr \
    --tcr_genelist_path ../data/tcr_genelist.min.json.xz \
    --tcr_cache_path tcr_cache.json \
    --cdr3_deletion_table_path ../data/cdr3_deletion_table.min.json.xz \
    --cdr3_insertion_table_path ../data/cdr3_insertion_table.min.json.xz \
    -b barcode.txt \
    -o sim_tcr
```

This generates following files:

- `sim_tcr.stats.tsv`, the ground-truth TCR contig statistics, with following columns:
  - `UUID`, the barcode, which is previously in UUID format.
  - `TRAV` \& `TRAJ` \& `TRBV` \& `TRBJ`, gene name of corresponding TCR segments.
  - `ACDR3_AA` \& `BCDR3_AA`, AA sequence of CDR3 on corresponding TCR segments.
  - `ACDR3_NT` \& `BCDR3_NT`, NT sequence of CDR3 on corresponding TCR segments.
  - `ALPHA_AA` \& `BETA_AA`, AA sequence of corresponding TCR chain.
  - `ALPHA_NT` \& `BETA_NT`, NT sequence of corresponding TCR chain.
- `sim_tcr.aa.fa` and `sim_tcr.nt.fa`, ground-truth contig in AA and NT with seqname `{barcode}:A` for TCR alpha chain and `{barcode}:B` for TCR beta chain.

+++

Strip the FASTA before performing calling YASIM RNA-Seq interface, and perform single-end bulk RNA-Seq using ART.

```{code-cell}
:tags: [skip-execution]

%%bash
python -m labw_utils.bioutils split_fasta sim_tcr.nt.fa
python -m yasim art \
    -F sim_tcr.nt.fa.d \
    -o sim_tcr_50 \
    --sequencer_name GA2 \
    --read_length 50 \
    -d tcr_depth.tsv \
    -e art_illumina \
    -j 20
```

Generates:

- `sim_tcr_50.d`: The target directory. See FASTQ files inside.
- `sim_tcr_50.fq`: Useless merged file. DO NOT USE.
- `sim_tcr_50.fq.stats`: Real simulation statistics.

```{code-cell}
!ls -lFh
```
