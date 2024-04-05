library("Seurat")
library("tidyverse")
library("UpSetR")

hgnc_genes <- (readr::read_tsv("raw_data/hgnc_complete_set.txt", show_col_types = FALSE) %>%
                   dplyr::filter(locus_group == "protein-coding gene") %>%
                   dplyr::filter(!is.na(mane_select)) %>%
                   dplyr::select(symbol))$symbol
sos <- list()

sample_names <- c(
    "HU_0043_Blood_10x", 
    "HU_0196_Kidney_GSE109564", 
    "HU_0148_Decidua_EBI",
    "HU_0223_Muscle_GSE134355",
    "HU_0125_Cerebrospinal-Fluid_GSE134577"
)
# Synthesize shared gene set
orig_genes_table <- data.frame(
    sample_name="HGNC",
    gene_name=hgnc_genes
)

for (sample_name in sample_names) {
    print(sample_name)
    so <- sprintf("raw_data/%s_gene_count.h5", sample_name) %>%
        Seurat::Read10X_h5() %>%
        Seurat::CreateSeuratObject()
    sannot <- readr::read_tsv(
        sprintf("raw_data/%s_meta.txt", sample_name),
        show_col_types = FALSE
    )
    so <- AddMetaData(
        object = so,
        metadata = sannot$Celltype,
        col.name = 'Celltype'
    )
    rm(sannot)
    sos[[sample_name]] <- so
    orig_genes_table <- rbind(
        orig_genes_table,
        data.frame(
            sample_name=sample_name,
            gene_name=rownames(so[["RNA"]]@features)
        )
    )
}
orig_genes_table_wide <- orig_genes_table %>%
    dplyr::mutate(value=1) %>%
    tidyr::pivot_wider(
        id_cols=gene_name,
        names_from = sample_name, 
        values_fill = 0) %>%
    as.data.frame()
row.names(orig_genes_table_wide) <- orig_genes_table_wide$gene_name
orig_genes_table_wide$gene_name <- NULL
pdf("genes_share.pdf", width=10, height=10)
UpSetR::upset(orig_genes_table_wide, order.by = "freq")
dev.off()

sel_genes <- names(which(apply(orig_genes_table_wide, 1, prod) == 1))

cell_level_qc_df <- data.frame()

for (sample_name in sample_names) {
    so <- sos[[sample_name]]
    so_sel <- subset(so, features = names(which(so[["RNA"]]@features %in% hgnc_genes)))
    sos[[sample_name]] <- so_sel
    scm_c <- as.matrix(so_sel[["RNA"]]$counts)
    cell_level_qc_df <- rbind(
        cell_level_qc_df,
        data.frame(
            sample_name=sample_name,
            nFeature=apply(scm_c, 2, function(x){sum(x>0)})
        )
    )
}
print("Plotting QC...")
cell_level_qc_df_long <- cell_level_qc_df %>%
    tidyr::pivot_longer(c("nFeature"))
p <- ggplot(cell_level_qc_df_long) +
    geom_violin(aes(y=sample_name, x=value)) +
    facet_wrap(name~., scales="free") +
    theme_bw()
ggsave("ncount_nfeature.pdf",p, width=10, height=5)


for (sample_name in sample_names) {
    sos[[sample_name]] <- sos[[sample_name]] %>%
        NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(npcs = 30, verbose = FALSE) %>%
        RunUMAP(reduction = "pca", dims = 1:30) %>%
        FindNeighbors(reduction = "pca", dims = 1:30) %>%
        FindClusters(reduction = "pca", resolution = 1)
    umap_data <- sos[[sample_name]]@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        dplyr::mutate(is_t_cell=sos[[sample_name]]$Celltype %in% c("CD4T", "CD8T"))
    p <- ggplot(umap_data) +
        geom_point(aes(x=umap_1, y=umap_2, color=is_t_cell), size=0.5, alpha=0.1) +
        theme_bw()
    ggsave(sprintf("%s_ctype.pdf", sample_name), p, width = 6,height = 5)
}

saveRDS(object = sos, file = "sos_post_qc.rds")
