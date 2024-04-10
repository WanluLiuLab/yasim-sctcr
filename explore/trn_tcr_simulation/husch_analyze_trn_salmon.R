library(Seurat)
library(tidyverse)

mut_table <- function(x, suffix) {
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
    colnames(x) <- paste(colnames(x), suffix, sep="-")
    return(x)
}

sample_names <- c(
    "HU_0043_Blood_10x"
    # "HU_0196_Kidney_GSE109564",
    # "HU_0148_Decidua_EBI",
    # "HU_0223_Muscle_GSE134355",
    # "HU_0125_Cerebrospinal-Fluid_GSE134577"
)

for (sample_name in sample_names) {
    print(sample_name)
    so <- list(
        real = Seurat::CreateSeuratObject(
            arrow::read_parquet(sprintf("parquets/%s_real.parquet", sample_name)) %>% mut_table("real"),
            project = "real"
        ) %>% AddMetaData(metadata = "real", col.name = "dtype"),
        sim = Seurat::CreateSeuratObject(
            arrow::read_parquet(sprintf("parquets/%s_sim.parquet", sample_name)) %>% mut_table("sim"),
            project = "sim"
        ) %>% AddMetaData(metadata = "sim", col.name = "dtype"),
        art = Seurat::CreateSeuratObject(
            arrow::read_parquet(sprintf("%s.sim.d/sim_dw_sampled.parquet", sample_name)) %>% mut_table("art"),
            project = "art"
        )%>% AddMetaData(metadata = "art", col.name = "dtype")
    )
    sel_features <- rownames(so[["art"]][["RNA"]]@features)

    so[["merged"]] <- merge(
        merge(so[["real"]],  so[["sim"]]),
        so[["art"]]
    ) %>%
        subset(features=sel_features)

    so[["merged_processed"]] <- so[["merged"]] %>%
        NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 500) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(npcs = 30, verbose = FALSE) %>%
        RunUMAP(reduction = "pca", dims = 1:30) %>%
        FindNeighbors(reduction = "pca", dims = 1:30) %>%
        FindClusters(reduction = "pca", resolution = 1)

    pdf(sprintf("figs/%s_data_withart.pdf", sample_name), width = 16, height = 5)

    print(DimPlot(so[["merged_processed"]],
                  reduction = "umap",
                  group.by = "orig.ident",
                  split.by = "dtype",
                  pt.size = 0.05
    ))
    dev.off()
}