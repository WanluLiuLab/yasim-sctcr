# `yasim_sctcr` -- Yet Another SIMulator for Single-Cell T-Cell Receptor Sequencing (scTCR-Seq)

**Markdown compatibility guide** This file is written in [Myst-flavored Markdown](https://myst-parser.readthedocs.io/), and may show errors on the default landing page of PYPI or Git Hostings. You can correctly preview it on generated Sphinx documentation or [Visual Studio Code](https://code.visualstudio.com) with [ExecutableBookProject.myst-highlight](https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight) plugin.

---

## Introduction

Single-Cel T-Cell Receptor (TCR) Sequencing (scTCR-Seq) is an important method in studying the diversity and dynamics of T-cell populations in organisms. This software provides an easy way to simulated Next-Generation Sequencing (NGS)-based scTCR-Seq using Illumina sequencer simulator.

**Significance** This software can mimic realistic TCR recombination events. Currently it uses human TCR V/J frequency data and CDR3 statistics from [hUARdb](https://huarc.net).

**Limitations** This software **CANNOT** simulate 10X genomics reads, realistic TCR expression differences and T-Cell clonal expansion status. The TCR it generated are not biologically tested.

## Installation

Since this software had not been released onto PYPI, you can only install it by building from source yourself. Retrive its source code using:

% Source code publication TODO

You need Python interpreter (CPython implementation) >= 3.7, latest PYPA [`build`](https://pypa-build.readthedocs.io), and [`setuptools`](https://setuptools.pypa.io/) to build this software. You are recommended to build the software in a virtual environment provided by [`virtualenv`](https://virtualenv.pypa.io), etc.

Build and install the simulator using:

```shell
cd yasim-sctcr
python3 -m build
pip install dist/yasim-sctcr-XXX-py3-none-any.whl
```

Apart from above instructions, you should also install [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) which is a general-purpose NGS DNA-Seq simulator and is available from [Conda](https://anaconda.org/bioconda/art) and [APT](https://packages.debian.org/stable/art-nextgen-simulation-tools). Tested versions are `2.5.8 (June 6, 2016)`.

## Minima Working Example (MWE)

This MWE retrives human data and test `yasim_sctcr`. We assume you're inside the root directory of this repository (i.e., directory of this Readme).

```shell
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
gzip -dcf hg38.ncbiRefSeq.gtf.gz \
    | grep -e '^chr7\s' -e '^chr14\s' \
    > hg38.ncbiRefSeq_chr7_14.gtf

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

python -m yasim_sc generate_barcode -n 10 -o barcode.txt

python -m yasim_sctcr generate_tcr_depth \
    -b barcode.txt \
    -o tcr_depth.tsv \
    -d 20
python -m yasim_sctcr generate_tcr_cache \
    --tcr_genelist_path data/tcr_genelist.min.json.xz \
    -o tcr_cache.json \
    --tcr_aa_table_path data/IMGT_Protein_Display.min.json.xz \
    -f hg38.fa \
    -g hg38.ncbiRefSeq_chr7_14.gtf
python -m yasim_sctcr rearrange_tcr \
    --tcr_genelist_path data/tcr_genelist.min.json.xz \
    --tcr_cache_path tcr_cache.json \
    --cdr3_deletion_table_path data/cdr3_deletion_table.min.json.xz \
    --cdr3_insertion_table_path data/cdr3_insertion_table.min.json.xz \
    -b barcode.txt \
    -o sim_tcr
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
