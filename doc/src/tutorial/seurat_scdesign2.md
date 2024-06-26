# Simulating scRNA-Seq Data using Seurat and scDesign2

Although we can directly use raw scRNA-Seq count matrix, count-level scRNA-Seq simulators like [scDesign2](https://github.com/JSB-UCLA/scDesign2) could provide more flexible options in selecting genes, cell types, and generate any number of cells. Here we will use the [`HU_0043_Blood_10x`](http://husch.comp-genomics.org/#/detail/HU_0043_Blood_10x) data from [HUSCH](http://husch.comp-genomics.org/) database as example.

```{note}
This tutorial does not require GNU/Linux. You may execute it in anywhere with [R](https://www.r-project.org/) and reqired packages.

Basic knowledge on R (especially [Tidyverse](https://www.tidyverse.org/) series and [Seurat](https://satijalab.org/seurat/)) are assumed for this tutorial.
```

## Downloading Raw Data

We will firstly download its expression data and cell type annotations. For example, following code downloads the data using GNU WGet:

```shell
wget https://biostorage.s3.ap-northeast-2.amazonaws.com/HUSCH/HUSCH_data/HU_0043_Blood_10x/HU_0043_Blood_10x_gene_count.h5
wget https://biostorage.s3.ap-northeast-2.amazonaws.com/HUSCH/HUSCH_data/HU_0043_Blood_10x/HU_0043_Blood_10x_meta.txt
```

````{tip}
Also available through AWS CLI.

```shell
aws s3 cp --no-sign-request s3://biostorage/HUSCH/HUSCH_data/HU_0043_Blood_10x/HU_0043_Blood_10x_gene_count.h5 sample_gene_count.h5
aws s3 cp --no-sign-request s3://biostorage/HUSCH/HUSCH_data/HU_0043_Blood_10x/HU_0043_Blood_10x_meta.txt sample_meta.txt
```

````

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

## Quality Control

We assume that the gene symbols used by HUSCH are HGNC-compatible, and assume you've got `ens.trans_gene_map.tsv` from [Zenodo](https://doi.org/10.5281/zenodo.12513698). For reduction of computation power and consistency across reference genomes, we will remove unknown or unselected genes:

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

## Fitting and Simulation using scDesign2

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

## Writing the Simulated Data to Disk

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
