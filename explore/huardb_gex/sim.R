library("SPARSim")
library("Seurat")

df <-CreateSeuratObject(
    Seurat::Read10X("C12_PBMC/filtered_feature_bc_matrix")
)
df <- subset(df, nFeature_RNA > 200 & nCount_RNA > 500)

cm  <- as.matrix(df[["RNA"]]$counts)
ncm <- SPARSim::scran_normalization(cm)

params <- SPARSim::SPARSim_estimate_parameter_from_data(
    cm,
    ncm,
    list(a=colnames(cm))
)
names(params[[1]]$lib_size) <- colnames(cm)

sim_result <- SPARSim_simulation(params)

T_10X_like_dataset <- SPARSim_simulation(T_10X_param_preset)
