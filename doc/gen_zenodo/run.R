library("tidyverse")
library("arrow")
library("biomaRt")
library("Seurat")
library("scDesign2")

ensembl97 <- useMart(
    host = 'https://jul2019.archive.ensembl.org',
    biomart = 'ENSEMBL_MART_ENSEMBL',
    dataset = 'hsapiens_gene_ensembl'
)

gene_trans <- getBM(
    attributes = c("ensembl_gene_id_version", "ensembl_transcript_id_version"),
    mart = ensembl97
)
mane_trans <- getBM(
    attributes = c("transcript_mane_select", "ensembl_transcript_id_version"),
    mart = ensembl97
) %>%
    dplyr::filter(transcript_mane_select != "")
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
readr::write_tsv(selected_gene_trans, "ref/ens.trans_gene_map.tsv",col_names=FALSE)

sample_name <- "HU_0043_Blood_10x"
sannot <- readr::read_tsv(
    sprintf("raw_data/%s_meta.txt", sample_name),
    show_col_types = FALSE
)
so <- sprintf("raw_data/%s_gene_count.h5", sample_name) %>%
    Read10X_h5() %>%
    CreateSeuratObject() %>%
    AddMetaData(metadata = sannot$Celltype, col.name = "Celltype") %>%
    subset(features = names(which(so[["RNA"]]@features %in% selected_gene_trans$hgnc_symbol)))

scm <- as.matrix(so[["RNA"]]$counts)
colnames(scm) <- so$Celltype
copula_result <- fit_model_scDesign2(scm, unique(so$Celltype), marginal = "poisson")
n_cells <- 500
sim_count_copula <- simulate_count_scDesign2(
    copula_result,
    n_cells,
    cell_type_prop = table(so$Celltype) / length(so$Celltype)
)
row.names(sim_count_copula) <- row.names(scm)

sim_count_copula_df <- as.data.frame(sim_count_copula)
colnames(sim_count_copula_df) <- paste(
    colnames(sim_count_copula_df),
    as.character(seq_len(ncol(sim_count_copula_df))),
    sep="_"
)
arrow::write_parquet(
    tibble::as_tibble(sim_count_copula_df, rownames = "FEATURE"),
    sprintf("ref/%s_sim.parquet", sample_name)
)
