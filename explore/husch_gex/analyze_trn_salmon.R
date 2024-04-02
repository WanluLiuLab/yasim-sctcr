library(Seurat)
library(tidyverse)

mut_table <- function(x) {
    x <- as.data.frame(x)
    row.names(x) <- x$FEATURE
    x$FEATURE <- NULL
    colnames(x) <- as.character(seq_len(ncol(x)))
    return(x)
}

sample_name <- "HU_0196_Kidney_GSE109564"


sos <- list()
sos[[sample_name]] <- list(
    real = Seurat::CreateSeuratObject(
        arrow::read_parquet(sprintf("%s_real.parquet", sample_name), show_col_types = FALSE) %>% mut_table(),
        project = "real",
        min.cells = 3,
        min.features = 200
    ),
    sim = Seurat::CreateSeuratObject(
        arrow::read_parquet(sprintf("%s_sim.parquet", sample_name), show_col_types = FALSE) %>% mut_table(),
        project = "sim",
        min.cells = 3,
        min.features = 200
    ),
    art = Seurat::CreateSeuratObject(
        arrow::read_parquet("/home/yuzj/Documents/pbsim3_modern/HU_0196_Kidney_GSE109564.art_salmon.parquet", show_col_types = FALSE) %>% mut_table(),
        project = "art",
        min.cells = 3,
        min.features = 200
    )
)
sos[[sample_name]][["merged"]] <- merge(
    merge(
        sos[[sample_name]][["real"]],
        sos[[sample_name]][["sim"]]
    ),
    sos[[sample_name]][["art"]]
)

sos[[sample_name]][["merged_processed"]] <- sos[[sample_name]][["merged"]] %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(reduction = "pca", resolution = 1)

pdf(sprintf("%s_data_withart.pdf", sample_name), width = 6, height = 5)
print(DimPlot(sos[[sample_name]][["merged_processed"]],
              reduction = "umap",
              group.by = "orig.ident",
              pt.size = 0.05
))
dev.off()

