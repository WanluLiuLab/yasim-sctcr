library(Seurat)
library(tidyverse)

mut_table <- function(x) {
    x <- as.data.frame(x)
    row.names(x) <- x$FEATURE
    x$FEATURE <- NULL
    return(x)
}

sos <- list()
sample_names <- c(
    "HU_0043_Blood_10x"
    # "HU_0196_Kidney_GSE109564",
    # "HU_0148_Decidua_EBI",
    # "HU_0223_Muscle_GSE134355",
    # "HU_0125_Cerebrospinal-Fluid_GSE134577"
)
for (sample_name in sample_names) {
    sos[[sample_name]] <- list(
        real = Seurat::CreateSeuratObject(
            arrow::read_parquet(sprintf("parquets/%s_real.parquet", sample_name), show_col_types = FALSE) %>%
                mut_table(),
            project = "real",
            min.cells = 3,
            min.features = 200
        ) %>%
            AddMetaData("real", col.name = "dtype"),
        sim = Seurat::CreateSeuratObject(
            arrow::read_parquet(sprintf("parquets/%s_sim.parquet", sample_name), show_col_types = FALSE) %>%
                mut_table(),
            project = "sim",
            min.cells = 3,
            min.features = 200
        ) %>%
            AddMetaData("sim", col.name = "dtype")
    )
    sos[[sample_name]][["merged"]] <- merge(
        sos[[sample_name]][["sim"]],
        sos[[sample_name]][["real"]]
    )
}

for (sample_name in sample_names) {
    so <- sos[[sample_name]][["merged"]] %>%
        NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(npcs = 30, verbose = FALSE) %>%
        RunUMAP(reduction = "pca", dims = 1:30)
    umap_data <- so@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        dplyr::mutate(
            cell_name = so@meta.data$dtype,
            cell_type = so@meta.data$orig.ident %in% c("CD4T", "CD8T")
        )
    p <- ggplot(umap_data) +
        geom_point(aes(x = umap_1, y = umap_2, color = cell_name), size = 0.5, alpha = 0.1) +
        facet_grid(. ~ cell_type) +
        scale_color_manual(values = c("#1d3557", "#ff8fab")) +
        theme_bw()
    ggsave(sprintf("figs/%s_data.pdf", sample_name), p, width = 11, height = 5)
}

