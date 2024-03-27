library("SPARSim")
library("Seurat")


for (data_name in c("C12_PBMC", "CPInc_NC1", "nc_CT1", "U4_PBMC")) {
    df <- CreateSeuratObject(
        Seurat::Read10X(sprintf("%s/filtered_feature_bc_matrix", data_name))
    )
    df <- subset(df, nFeature_RNA > 200 & nCount_RNA > 500)
    cm <- as.matrix(df[["RNA"]]$counts)
    ncm <- SPARSim::scran_normalization(cm)
    params <- SPARSim::SPARSim_estimate_parameter_from_data(
        cm,
        ncm,
        list(a = colnames(cm))
    )
    names(params[[1]]$lib_size) <- colnames(cm)
    sim_result <- SPARSim_simulation(params)
    write.csv(cm, sprintf("%s.cnt_mtx.csv", data_name))
    write.csv(sim_result$count_matrix, sprintf("%s.sim_mtx.csv", data_name))
}
