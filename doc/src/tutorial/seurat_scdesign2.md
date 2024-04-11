# Simulate Paired 5' scRNA-Seq and scTCR-Seq data using Seurat, scDesign2, and YASIM-scTCR

[Single Cell Immune Profiling](https://www.10xgenomics.com/support/single-cell-immune-profiling) by [10X Genomics](https://www.10xgenomics.com/) allows sense of both TCR/BCR sequence and gene expression data using 5' sequencing on a single-cell basis. This tutorial tells you how to simulate such data using YASIM-scTCR.

Make sure you had correctly set up the environment as-is specified in {doc}`Readme </_root/Readme>`.

```{note}
Basic knowledge on Shell programming, R (especially [Tidyverse](https://www.tidyverse.org/) series and [Seurat](https://satijalab.org/seurat/)) and Python are assumed for this tutorial.

You should also know the basis of [regular expressions](https://en.wikipedia.org/wiki/Regular_expression) and its [Python implementation](https://docs.python.org/3/library/re.html).
```

## Step 0. Preparation

### Making the Reference Transcriptome

We already know that YASIM LLRGs works on transcript-level, while scRNA-Seq data are usually provided on a gene-level basis. So, it is important to convert gene IDs to its **representative** transcripts. Here we will use protein-coding genes (PCGs) in Ensembl release 97 that have HGNC symbol and use MANE-selected transcripts as representative transcripts. If a gene has multiple MANE-selected transcripts, we will randomly decide one.

Firstly, we need a mapping between HGNC gene symbols and transcript IDs. The following R script could perform this:

```r
library("tidyverse")
library("biomaRt")

# You may use other Ensembl release version if you like.
ensembl97 <- useMart(
    host = 'https://jul2019.archive.ensembl.org',
    biomart = 'ENSEMBL_MART_ENSEMBL',
    dataset = 'hsapiens_gene_ensembl'
)
# Get a list of usable attributes using following command:
# attrs <- listAttributes(ensembl97) %>% as_tibble()

# Gene ID-Transcript ID mapping.
gene_trans <- getBM(
    attributes = c("ensembl_gene_id_version", "ensembl_transcript_id_version"),
    mart = ensembl97
)
# Ensembl transcript ID-MANE selected RefSeq ID mapping.
mane_trans <- getBM(
    attributes = c("transcript_mane_select", "ensembl_transcript_id_version"),
    mart = ensembl97
) %>%
    dplyr::filter(transcript_mane_select != "")
# Ensembl gene ID-HGNC symbol mapping.
hgnc_gene <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id_version", "gene_biotype"),
    mart = ensembl97
) %>%
    dplyr::filter(hgnc_symbol != "", gene_biotype == "protein_coding")

selected_gene_trans <- gene_trans %>%
    dplyr::inner_join(mane_trans, by="ensembl_transcript_id_version") %>%
    dplyr::select(!transcript_mane_select) %>%
    dplyr::inner_join(hgnc_gene, by="ensembl_gene_id_version") %>%
    dplyr::select(!c(ensembl_gene_id_version, gene_biotype))
readr::write_tsv(selected_gene_trans, "ens.trans_gene_map.tsv", col_names=FALSE)
```

**Generates:**

- `ens.trans_gene_map.tsv`, which is a Tab-Separated Value (TSV) whose first column is Ensembl transcript ID and second column is HGNC symbol of corresponding gene. [Salmon](https://salmon.readthedocs.io/en/latest/index.html) could use this file to perform gene-level quantification.

Now we would download Ensembl transcript cDNAs and subset them to exclude unselected transcripts.

```shell
wget https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
# Use sektk to select selected transcripts.
seqtk subseq \
    Homo_sapiens.GRCh38.cdna.all.fa \
    <(cut -f 1 ens.trans_gene_map.tsv) \
    >ens.sel_genes.cdna.fa
```

**Generates:**

- `ens.sel_genes.cdna.fa`, a FASTA of all selected transcripts that is considerably smaller than the original one.

### Constructing the TCR Cache

The TCR cache is a **reusable** JSON containing TCR nucleotide to amino acid alignment information, which is required for simulation of TCR rearrangement. The generated file can be reused for all studies using the same reference TCR sequence.

Before generating the TCR cache, we would download TCR sequences in cDNA and amino acid from Ensembl. The following code generates FASTA for peptides and cDNA sequences for TRAV/TRBV genes used in this simulator from Ensembl. The simulator uses data built from huARdb v1 which uses Ensembl 97 for reference, so we use this version in our tutorial.

```shell
# The cDNA reference had already been downloaded.
wget https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
for name in cdna pep; do
    seqkit grep \
        --by-name \
        --use-regexp \
        -p 'chromosome:GRCh38:(7|14):' \
        Homo_sapiens.GRCh38."${name}".all.fa |
        seqkit grep \
            --by-name \
            --use-regexp \
            -p 'gene_biotype:TR_[VDJC]_gene' \
            /dev/stdin |
        seqkit grep \
            --by-name \
            --use-regexp \
            -p 'gene_symbol:TR[AB]' \
            /dev/stdin | 
        sed -E 's;^>.+ gene_symbol:(\S+) .+$;>\1;' \
        >ens."${name}".fa
    samtools faidx ens."${name}".fa
done
```

The above command will select genes that:

- Are on chromosome 7 or 14, which excludes pesudogene `TRBV20OR9-2` (see it on [ENSG00000205274](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000205274;r=9:33617845-33618508;t=ENST00000379435;redirect=no) or [HGNC:12197](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:12197)).
- Have `TR_V_gene`/`TR_D_gene`/`TR_J_gene`/`TR_C_gene` as their [biotype](https://www.ensembl.org/info/genome/genebuild/biotypes.html).
- Gene symbol matches `TR[AB]` regular expression.

Since there will be only 1 transcript for all TCR gene segments in Ensembl release 97, we replaced FASTA header for each record with the corresponding gene symbol. **This may not be true in other reference genomes or further releases of Esnembl.**

It will finally generate `ens.cdna.fa` for cDNA sequences and `ens.pep.fa` for peptide sequences and their indices in FAI format if no error occurs. We may generate a TCR cache using:

```shell
python -m yasim_sctcr generate_tcr_cache \
    --tcr_cdna_fa_path ens.cdna.fa \
    --tcr_pep_fa_path ens.pep.fa \
    -o tcr_cache.json
```

**Generating:**

- `tcr_cache.json`, the generated TCR cache. User may not use this file.

## Step 1. Learn from scRNA-Seq Data

Although we can directly use raw scRNA-Seq count matrix, count-level scRNA-Seq simulators like [scDesign2](https://github.com/JSB-UCLA/scDesign2) could provide more flexible options in selecting genes, cell types, and generate any number of cells. Here we will use the [`HU_0043_Blood_10x`](http://husch.comp-genomics.org/#/detail/HU_0043_Blood_10x) data from [HUSCH](http://husch.comp-genomics.org/) database as example.

### Downloading Raw Data

We will firstly download its expression data and cell type annotations. Download the data using:

```shell
wget https://biostorage.s3.ap-northeast-2.amazonaws.com/HUSCH/HUSCH_data/HU_0043_Blood_10x/HU_0043_Blood_10x_gene_count.h5
wget https://biostorage.s3.ap-northeast-2.amazonaws.com/HUSCH/HUSCH_data/HU_0043_Blood_10x/HU_0043_Blood_10x_meta.txt
```

Then read it into Seurat:

```r
library("Seurat")
library("tidyverse")

# Decalre sample name for convenience
sample_name <- "HU_0043_Blood_10x"

# Create seurat object
so <- CreateSeuratObject(Read10X_h5(sprintf("%s_gene_count.h5", sample_name)))
# Read cell type annotation
sannot <- readr::read_tsv(
    sprintf("%s_meta.txt", sample_name),
    show_col_types = FALSE
)

# Add annotation to metadata
so <- AddMetaData(
    object = so,
    metadata = sannot$Celltype,
    col.name = "Celltype"
)
```

### Quality Control

We assume that the gene symbols used by HUSCH are HGNC-compatible. For reduction of computation power and consistency across reference genomes, we will remove unknown or unselected genes:

```r
# Read HGNC genes retrived in previous steps
hgnc_genes <- unique(
    (
        readr::read_tsv(
            "ens.trans_gene_map.tsv", 
            col_names = c("TRANSCRIPT_ID", "GENE_ID")
        )
    )$GENE_ID
)

# Subset seurat object to include selected genes only.
so <- subset(so, features = names(which(so[["RNA"]]@features %in% hgnc_genes)))
```

### Fitting and Simulation using scDesign2

Firstly, we will allow scDesign2 to fit the scRNA-Seq data. We use `"poisson"` for `marginal` parameter instead of defaults since default parameter generates several errors.

```r
library("scDesign2")

# Extract count matrux
scm <- as.matrix(so[["RNA"]]$counts)
# Rename the column names to cell types as-is required by scDesign2
colnames(scm) <- so$Celltype
# Fit the data
copula_result <- fit_model_scDesign2(scm, unique(so$Celltype), marginal = "poisson")
```

Now we simulate 500 new cells using the same proportion of each cell type:

```r
n_cells <- 500
sim_count_copula <- simulate_count_scDesign2(
    copula_result,
    n_cells,
    cell_type_prop = table(so$Celltype) / length(so$Celltype)
)
# scDesign2 may forget gene names. Use this to add them back.
row.names(sim_count_copula) <- row.names(scm)
```

It is also possible to perform a T-cell-only simulation:

```r
# Construction of "cell_type_prop" that sets probability of cells except CD4T/CD8T to zero.
all_cell_types <- unique(so$Celltype)
t_cells_names <- grep("CD[48]T", so$Celltype, value = TRUE)
other_cell_names <- so$Celltype[ ! so$Celltype %in% unique(t_cells_names)]
cell_type_prop <- c(table(t_cells_names) / length(t_cells_names), table(other_cell_names) * 0)
cell_type_prop <- cell_type_prop[sort(names(cell_type_prop), index.return=TRUE)$ix]

sim_count_copula <- simulate_count_scDesign2(
    copula_result,
    n_cells,
    cell_type_prop = cell_type_prop
)
row.names(sim_count_copula) <- row.names(scm)
```

The `sim_count_copula` object now contains the simulated data.

### Writing the Simulated Data to Disk

The on-disk simulated data format can be used by YASIM-scTCR `scaffold` command introduced below.

```r
library("arrow")

# Convert to dataframe. Notice the column names now are cell types and are duplicated.
sim_count_copula_df <- as.data.frame(sim_count_copula)

# Normalize column names. The column names is now cell type and n-th cell separated by `-`
colnames(sim_count_copula_df) <- paste(
    colnames(sim_count_copula_df),
    as.character(seq_len(ncol(sim_count_copula_df))),
    sep="_"
)

# Write to disk using Apache Parquet
arrow::write_parquet(
    tibble::as_tibble(sim_count_copula_df, rownames = "FEATURE"),
    sprintf("%s_sim.parquet", sample_name)
)
```

```{note}
**YASIM scRNA-Seq Format**

The scRNA-Seq data should be provided in a tabular format with one cell per column and one gene per row. The first column should be named `FEATURES` and contains gene name, with other column names started with the cell type information that allows separation of T cells and other (normal) cells using [Python regular expression](https://docs.python.org/3/library/re.html).

We would recommend using [Apache Parquet](https://parquet.apache.org) format for faster I/O. However, it will need [arrow](https://cran.r-project.org/web/packages/arrow/index.html) R package and [pyarrow](https://pypi.org/project/pyarrow) Python package. If it is not possible, you may also write the data in TSV format.
```

## Step 2. Generating YASIM-sc Scaffold

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

## Step 3. Simulate scTCR-Seq Data

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

The following example would perform single-end scTCR-Seq using ART.

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
    -j 20 \
    --amplicon
```

Please note that the `--amplicon` parameter (which is directly passed to ART) is crucial; without which the generated data will not be 5' amplified.

**Generates:**

- `HU_0043_Blood_10x.sim.d/sim_tcr_50.d`: The target directory where reads generated from TRA and TRB sequences of all cells are separated into single file.
- `HU_0043_Blood_10x.sim.d/sim_tcr_50.fq`: A merged file which contains all simulated reads.
- `HU_0043_Blood_10x.sim.d/sim_tcr_50.fq.stats`: Real simulation statistics.

## Step 4. Simulate scRNA-Seq Data

Since scRNA-Seq depth information had already been contained in the scaffold, we only need to invoke LLRG. The following example would perform scRNA-Seq using ART.

```shell
python -m labw_utils.bioutils split_fasta ens.sel_genes.cdna.fa
python -m yasim_sc art \
    -F ens.sel_genes.cdna.fa.d \
    -o HU_0043_Blood_10x/sim_notcr_50.d \
    -d HU_0043_Blood_10x/scGEX.depth.d \
    -e art_illumina \
    -j 20 \
    --amplicon
```

Generates:

- `HU_0043_Blood_10x/sim_notcr_50.d`: The target directory where reads generated from all cells is separated into single files.

