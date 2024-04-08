library("scDesign2")
library("Seurat")
library("tidyverse")

sample_names <- c(
    "HU_0196_Kidney_GSE109564",
    "HU_0043_Blood_10x",  
    "HU_0148_Decidua_EBI",
    "HU_0223_Muscle_GSE134355",
    "HU_0125_Cerebrospinal-Fluid_GSE134577"
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
    print(sample_name)
    so <- sos_post_qc[[sample_name]]
    scm <- as.matrix(so[["RNA"]]$counts)
    colnames(scm) <- so$Celltype
    w_sample(scm, sample_name, "real")

    copula_result <- fit_model_scDesign2(scm, unique(so$Celltype), marginal = "poisson")
    sim_count_copula <- simulate_count_scDesign2(
        copula_result,
        500,
        cell_type_prop = table(so$Celltype) / length(so$Celltype)
    )
    row.names(sim_count_copula) <- row.names(scm)
    w_sample(sim_count_copula, sample_name, "sim")
    gc()
}
