# Quickstart

**Simulate scRNA-Seq and scTCR-Seq data using YASIM-scTCR**

[Single Cell Immune Profiling](https://www.10xgenomics.com/support/single-cell-immune-profiling) by [10X Genomics](https://www.10xgenomics.com/) allows sense of both TCR/BCR sequence and gene expression data using 5' sequencing on a single-cell basis. This tutorial tells you how to simulate such data using YASIM-scTCR.

Make sure you had correctly set up the environment as-is specified in {doc}`Readme </_root/Readme>`.

```{note}
Modern GNU/Linux is assumed for this tutorial. If you're working under Microsoft Windows, you may use [WSL](https://learn.microsoft.com/en-us/windows/wsl/).

Basic knowledge on Shell programming and Python are assumed for this tutorial.

You should also know the basis of [regular expressions](https://en.wikipedia.org/wiki/Regular_expression) and its [Python implementation](https://docs.python.org/3/library/re.html).
```

Before simulation, you may download the reference data and simulated scRNA-Seq data from [Zenodo](TODO).

## Step 1. Generating YASIM-sc Scaffold

```{note}
**YASIM scRNA-Seq Format**

The scRNA-Seq data should be provided in a tabular format with one cell per column and one gene per row. The first column should be named `FEATURES` and contains gene name, with other column names started with the cell type information that allows separation of T cells and other (normal) cells using [Python regular expression](https://docs.python.org/3/library/re.html).

We would recommend using [Apache Parquet](https://parquet.apache.org) format for faster I/O. However, it will need [arrow](https://cran.r-project.org/web/packages/arrow/index.html) R package and [pyarrow](https://pypi.org/project/pyarrow) Python package. If it is not possible, you may also write the data in TSV format.
```


YASIM single cell scaffold allows clear representation of YASIM-generated files for paired scRNA-Seq and scTCR-Seq data from cell type annotated scRNA-Seq data. It would down-sample genes and normalize expression data. The scaffold for the sample selected above can be generated using the following command:

```shell
python -m yasim_sctcr scaffold \
    --transcript_gene_mapping ref/ens.trans_gene_map.tsv \
    --src_sc_data HU_0043_Blood_10x_sim.parquet \
    --out HU_0043_Blood_10x_sim.sim.d \
    --t_cell_regex 'CD[48]T'
```

**Generates:**

- `HU_0043_Blood_10x_sim.sim.d`, the generated directory for simulation.
  - `scGEX.depth.d`, a directory of depth files for each cell.
  - `n_cell_bc.txt`, barcodes for normal cells, one barcode per line..
  - `t_cell_bc.txt`, barcodes for T cells, one barcode per line.
  - `sim_dw_sampled.parquet` or `sim_dw_sampled.tsv`, down sampled normalized gene counts for each cell.

## Step 2. Simulate scTCR-Seq Data

The simulation of scTCR-Seq data mainly includes genesis of TCR contig sequences, clonal expansion and invocation of LLRGs.

### Simulate Ground-Truth TCR Contigs

Download the following real statistical data:

- [cdr3_deletion_table.min.json.xz](../../../data/cdr3_deletion_table.min.json.xz): Human CDR3 deletion table.
- [cdr3_insertion_table.min.json.xz](../../../data/cdr3_insertion_table.min.json.xz): Human CDR3 insertion table.

This version of YASIM-scTCR separates generation of TCR clonotypes and TCR sequences for distinct T-cells, allowing the simulation of clonal expansion. To do this, the number of TCR contigs should be less than 1/2 T-cells. The following example generates 100 TCRs:

```shell
python -m yasim_sctcr rearrange_tcr \
    --tcr_cache_path tcr_cache.json \
    --cdr3_deletion_table_path cdr3_deletion_table.json \
    --cdr3_insertion_table_path cdr3_insertion_table.json \
    --usage_bias_json usage_bias.json \
    -n 100 \
    -o HU_0043_Blood_10x.sim.d/sim_tcr
```

You may also generate unexpanded T-cells (i.e., number of T-cells equals to number of TCRs):

```shell
python -m yasim_sctcr rearrange_tcr \
    --tcr_cache_path tcr_cache.json \
    --cdr3_deletion_table_path cdr3_deletion_table.json \
    --cdr3_insertion_table_path cdr3_insertion_table.json \
    --usage_bias_json usage_bias.json \
    -n "$(cat HU_0043_Blood_10x.sim.d/t_cell_bc.txt | wc -l)" \
    -o HU_0043_Blood_10x.sim.d/sim_tcr
```

This step may report `finish with N failures` -- don't worry! None of your cells would lose.

**Generates:**

- `sim_tcr.stats.tsv`, the ground-truth TCR contig statistics, with following columns:
  - `UUID`, the unique identifier of each TCR.
  - `TRAV` \& `TRAJ` \& `TRAC` \& `TRBV` \& `TRBJ` \& `TRBC`, gene name of corresponding TCR segments.
  - `ACDR3_AA` \& `BCDR3_AA`, AA sequence of CDR3 on corresponding TCR segments.
  - `ACDR3_NT` \& `BCDR3_NT`, NT sequence of CDR3 on corresponding TCR segments.
  - `ALPHA_AA` \& `BETA_AA`, AA sequence of the corresponding TCR chain.
  - `ALPHA_NT` \& `BETA_NT`, NT sequence of the corresponding TCR chain.
- `sim_tcr.nt.fa`, ground-truth contig in nucleotide with seqname `{UUID}_A` for TCR alpha chain and `{UUID}_B` for TCR beta chain.

### Simulation of T-Cell Clonal Expansion

TCR clonal expansion was found to follow a Zipf's distribution. Here we would generate it. A clonal expansion distribution is generated, and clonal-expanded TCR contigs will be allocated to each T-cell (represented by distinct barcode).

```shell
python -m yasim_sctcr generate_tcr_clonal_expansion \
    -b HU_0043_Blood_10x.sim.d/t_cell_bc.txt \
    --src_tcr_stats_tsv HU_0043_Blood_10x.sim.d/sim_tcr.stats.tsv \
    --dst_nt_fasta HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa \
    --alpha 1
```

**Generates**:

- `HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa`: Nucleotide sequences for TCR of each T-cell with seqname `{BARCODE}_A` for TCR alpha chain and `{BARCODE}_B` for TCR beta chain.

### Generation of TCR Depths

Before introducing machine errors, we will add depth information to each T-cell. The following code generates a depth file with 20 (uniform distribution) as sequencing depth for each T-cell.

```shell
python -m yasim_sctcr generate_tcr_depth \
    -b HU_0043_Blood_10x.sim.d/t_cell_bc.txt \
    -o HU_0043_Blood_10x.sim.d/scTCR.depth.tsv \
    -d 20
```

**Generates:**

- `HU_0043_Blood_10x.sim.d/scTCR.depth.tsv`: A YASIM Isoform-level-depth compatible depth file. Here, TRA and TRB contig of a cell is regarded as one isoform.

### Invocation of LLRGs

The following example would perform single-end scTCR-Seq using ART. To simulate full-length scRNA-Seq (i.e., SMART-Seq, SMART-Seq2, etc.), do the fiollowing:

```shell
# Split the FASTA before performing calling YASIM RNA-Seq interface
python -m labw_utils.bioutils split_fasta HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa
python -m yasim art \
    -F HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa.d \
    -o HU_0043_Blood_10x/sim_tcr_50 \
    --sequencer_name GA2 \
    --read_length 50 \
    -d HU_0043_Blood_10x.sim.d/scTCR.depth.tsv \
    -e art_illumina \
    --preserve_intermediate_files \
    -j 20
```

You may simulate 3' amplified scTCR-Seq with following step:

```shell
python -m labw_utils.bioutils split_fasta HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa
python -m yasim art \
    -F HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa.d \
    -o HU_0043_Blood_10x/sim_tcr_50 \
    --sequencer_name GA2 \
    --read_length 50 \
    -d HU_0043_Blood_10x.sim.d/scTCR.depth.tsv \
    -e art_illumina \
    --preserve_intermediate_files \
    -j 20 \
    --amplicon
```

You may simulate 5' amplified scTCR-Seq with following step:

```shell
seqkit seq --reverse --complement \
    <HU_0043_Blood_10x.sim.d/sim_t_cell.nt.fa \
    >HU_0043_Blood_10x.sim.d/sim_t_cell.rc.nt.fa
python -m labw_utils.bioutils split_fasta HU_0043_Blood_10x.sim.d/sim_t_cell.rc.nt.fa
python -m yasim art \
    -F HU_0043_Blood_10x.sim.d/sim_t_cell.rc.nt.fa.d \
    -o HU_0043_Blood_10x/sim_tcr_50 \
    --sequencer_name GA2 \
    --read_length 50 \
    -d HU_0043_Blood_10x.sim.d/scTCR.depth.tsv \
    -e art_illumina \
    --preserve_intermediate_files \
    -j 20 \
    --amplicon
```

**Generates:**

- `HU_0043_Blood_10x.sim.d/sim_tcr_50.d`: The target directory where reads generated from TRA and TRB sequences of all cells are separated into single file.
- `HU_0043_Blood_10x.sim.d/sim_tcr_50.fq`: A merged file which contains all simulated reads.
- `HU_0043_Blood_10x.sim.d/sim_tcr_50.fq.stats`: Real simulation statistics.

## Step 3. Simulate scRNA-Seq Data

Since scRNA-Seq depth information had already been contained in the scaffold, we only need to invoke LLRG. The following example would perform scRNA-Seq using ART. The following example simulated full-length scRNA-Seq, with 5' or 3' amplified scRNA-Seq available using the method described above.

```shell
python -m labw_utils.bioutils split_fasta ens.sel_genes.cdna.fa
python -m yasim_sc art \
    -F ens.sel_genes.cdna.fa.d \
    -o HU_0043_Blood_10x/sim_notcr_50.d \
    -d HU_0043_Blood_10x/scGEX.depth.d \
    -e art_illumina \
    -j 20
```

Generates:

- `HU_0043_Blood_10x/sim_notcr_50.d`: The target directory where reads generated from all cells is separated into single files.

