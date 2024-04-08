library(Seurat)
library(tidyverse)

mut_table <- function(x) {
    x <- as.data.frame(x)
    if (length(x$symbol) != 0) {
        x <- x[!duplicated(x$symbol), ]
        row.names(x) <- x$symbol
        x$symbol <- NULL
    }
    else {
        row.names(x) <- x$FEATURE
        x$FEATURE <- NULL
    }
    colnames(x) <- as.character(seq_len(ncol(x)))
    print(x[1:10, 1:10])
    return(x)
}

sample_names <- c(
    "HU_0043_Blood_10x",
    "HU_0196_Kidney_GSE109564",
    "HU_0148_Decidua_EBI",
    "HU_0223_Muscle_GSE134355",
    "HU_0125_Cerebrospinal-Fluid_GSE134577"
)
sos <- list()
for (sample_name in sample_names) {
    print(sample_name)
    sos[[sample_name]] <- list(
        real = Seurat::CreateSeuratObject(
            arrow::read_parquet(sprintf("%s_real.parquet", sample_name)) %>% mut_table(),
            project = "real"
        ),
        sim = Seurat::CreateSeuratObject(
            arrow::read_parquet(sprintf("%s_sim.parquet", sample_name)) %>% mut_table(),
            project = "sim"
        ),
        art = Seurat::CreateSeuratObject(
            arrow::read_parquet(sprintf("%s.art_salmon.parquet", sample_name)) %>% mut_table(),
            project = "art"
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
        FindVariableFeatures(selection.method = "mean.var.plot", nfeatures = 2000) %>%
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
}