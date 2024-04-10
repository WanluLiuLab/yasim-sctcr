library("scDesign2")
library("Seurat")
library("tidyverse")

sample_names <- c(
    "HU_0043_Blood_10x"
    # "HU_0196_Kidney_GSE109564",
    # "HU_0148_Decidua_EBI",
    # "HU_0223_Muscle_GSE134355",
    # "HU_0125_Cerebrospinal-Fluid_GSE134577"
)

w_sample <- function(x, this_sample_name, data_type){
    x <- x %>% as.data.frame()
    colnames(x) <- paste(colnames(x), as.character(seq_len(ncol(x))), sep="_")
    x %>%
        tibble::as_tibble(rownames = "FEATURE") %>%
        arrow::write_parquet(sprintf("parquets/%s_%s.parquet", this_sample_name, data_type))
}
sos_post_qc <- readRDS("./sos_post_qc.rds")

for (sample_name in sample_names) {
    message(sample_name)
    so <- sos_post_qc[[sample_name]]
    scm <- as.matrix(so[["RNA"]]$counts)
    colnames(scm) <- so$Celltype
    # w_sample(scm, sample_name, "real")

    n_cells <- 500
    all_cell_types <- unique(so$Celltype)
    t_cells_names <- grep("CD[48]T", so$Celltype, value = TRUE)
    other_cell_names <- so$Celltype[ ! so$Celltype %in% unique(t_cells_names)]
    cell_type_prop <- c(table(t_cells_names) / length(t_cells_names), table(other_cell_names) * 0)
    cell_type_prop <- cell_type_prop[sort(names(cell_type_prop), index.return=TRUE)$ix]

    message("Learning...")
    copula_result <- fit_model_scDesign2(scm, unique(so$Celltype), marginal = "poisson")
    # message(sprintf("Simulating mixed sample for %d cells...", n_cells))
    # sim_count_copula <- simulate_count_scDesign2(
    #     copula_result,
    #     n_cells,
    #     cell_type_prop = table(so$Celltype) / length(so$Celltype)
    # )
    # row.names(sim_count_copula) <- row.names(scm)
    # w_sample(sim_count_copula, sample_name, "sim")
    for (n_t_cells in c(100, 500, 1000)){
        for (replicate_num in 1:10){
            message(sprintf("Simulating T-Cell only sample for %d cells rep %d...", n_t_cells, replicate_num))

            sim_count_copula <- simulate_count_scDesign2(
                copula_result,
                n_t_cells,
                cell_type_prop = cell_type_prop
            )
            row.names(sim_count_copula) <- row.names(scm)
            w_sample(sim_count_copula, sample_name, sprintf("sim_tcell_only_ncells%d_rep%d", n_t_cells, replicate_num))
        }
    }
    gc()
}
