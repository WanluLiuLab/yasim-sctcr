options(error = function() traceback(20))
library("scDesign2")
library("Seurat")
library("tidyverse")

hgnc_genes <- (readr::read_tsv("raw_data/hgnc_complete_set.txt", show_col_types = FALSE) %>%
    dplyr::filter(locus_group == "protein-coding gene") %>%
    dplyr::filter(!is.na(mane_select)) %>%
    dplyr::select(symbol))$symbol

# Synthesize shared gene set
for (sample_name in c("sample", "HU_0133_Cerebral-Cortex_GSE134355", "HU_0196_Kidney_GSE109564", "HU_0281_Blood_GSE138867")) {
    so <- Seurat::CreateSeuratObject(Seurat::Read10X_h5(sprintf("raw_data/%s_gene_count.h5", sample_name)))
    hgnc_genes <- intersect(hgnc_genes, rownames(so[["RNA"]]@features))
}
message(sprintf("Getting %d genes at last", length(hgnc_genes)))
gc()

w_sample <- function(x, this_sample_name, data_type){
    write.csv(x, sprintf("%s_%s.csv", this_sample_name, data_type))

    x <- x %>% as.data.frame()
    colnames(x) <- as.character(seq_len(ncol(x)))
    x %>%
        tibble::as_tibble(rownames = "FEATURE") %>%
        arrow::write_parquet(sprintf("%s_%s.parquet", this_sample_name, data_type))
}


for (sample_name in c("sample", "HU_0133_Cerebral-Cortex_GSE134355", "HU_0196_Kidney_GSE109564", "HU_0281_Blood_GSE138867")) {
    print(sample_name)
    so <- Seurat::CreateSeuratObject(Seurat::Read10X_h5(sprintf("raw_data/%s_gene_count.h5", sample_name)))
    sannot <- readr::read_tsv(sprintf("raw_data/%s_meta.txt", sample_name), show_col_types = FALSE)
    scm <- as.matrix(so[["RNA"]]$counts)
    colnames(scm) <- sannot$Celltype
    scm <- scm[rownames(scm) %in% hgnc_genes, ]
    w_sample(scm, sample_name, "real")

    copula_result <- fit_model_scDesign2(scm, unique(sannot$Celltype), marginal = "poisson")
    sim_count_copula <- simulate_count_scDesign2(
        copula_result,
        2000,
        sim_method = 'copula',
        cell_type_prop = table(sannot$Celltype) / length(sannot$Celltype)
    )
    row.names(sim_count_copula) <- row.names(scm)
    w_sample(sim_count_copula, sample_name, "sim")
    gc()
}
