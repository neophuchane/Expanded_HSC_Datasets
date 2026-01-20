library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(dplyr)

# =============================================================================
# LOAD FROM RAW FILES - NO OLD R OBJECTS
# =============================================================================

# STEP 1: Load AGM data from CSV files
message("Loading AGM CSV files...")

agm_files <- list.files(pattern = "GSM.*\\.csv\\.gz", full.names = TRUE)
message(paste("Found", length(agm_files), "AGM files"))

agm_objects <- lapply(agm_files, function(file) {
  message(paste("Reading:", basename(file)))
  data <- read.csv(file, row.names = 1)
  data_matrix <- as.matrix(data)
  data_sparse <- as(data_matrix, "sparseMatrix")
  
  # Extract sample name from filename
  sample_name <- gsub("GSM[0-9]+_Aorta-", "", basename(file))
  sample_name <- gsub(" \\(1\\)\\.csv\\.gz", "", sample_name)
  sample_name <- paste0("agm-", sample_name)
  
  obj <- CreateSeuratObject(counts = data_sparse, project = sample_name)
  obj$orig.ident <- sample_name
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  
  message(paste("Loaded:", sample_name))
  return(obj)
})

# STEP 2: Load Zheng .h5 files
message("\nLoading Zheng .h5 files...")
h5_files <- list.files(path = "Zheng", pattern = "\\.h5$", full.names = TRUE)
message(paste("Found", length(h5_files), ".h5 files"))

zheng_objects <- lapply(h5_files, function(file) {
  message(paste("Processing:", basename(file)))
  data <- Read10X_h5(file)
  obj <- CreateSeuratObject(counts = data, project = basename(file))
  obj$orig.ident <- basename(file)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  return(obj)
})

# STEP 3: Load Somarin MTX files
message("\nLoading Somarin fetal liver datasets...")

fl1_data <- Read10X(data.dir = "Somarin/FL1_hpc")
if(is.list(fl1_data)) fl1_data <- fl1_data$`Gene Expression`
fl1_obj <- CreateSeuratObject(counts = fl1_data, project = "FL1_CS22")
fl1_obj$orig.ident <- "FL1_hpc_CS22"
fl1_obj <- NormalizeData(fl1_obj)
fl1_obj <- FindVariableFeatures(fl1_obj)
fl1_obj <- ScaleData(fl1_obj)

fl2_data <- Read10X(data.dir = "Somarin/FL2_hpc")
if(is.list(fl2_data)) fl2_data <- fl2_data$`Gene Expression`
fl2_obj <- CreateSeuratObject(counts = fl2_data, project = "FL2_CS22")
fl2_obj$orig.ident <- "FL2_hpc_CS22"
fl2_obj <- NormalizeData(fl2_obj)
fl2_obj <- FindVariableFeatures(fl2_obj)
fl2_obj <- ScaleData(fl2_obj)

fl3_data <- Read10X(data.dir = "Somarin/FL_CS16_and_W9")
if(is.list(fl3_data)) fl3_data <- fl3_data$`Gene Expression`
fl3_obj <- CreateSeuratObject(counts = fl3_data, project = "FL_CS16_W9")
fl3_obj$orig.ident <- "FL_CS16_W9"
fl3_obj <- NormalizeData(fl3_obj)
fl3_obj <- FindVariableFeatures(fl3_obj)
fl3_obj <- ScaleData(fl3_obj)

# STEP 4: Merge everything
message("\nMerging all datasets...")
all_objects <- c(agm_objects, zheng_objects, list(fl1_obj, fl2_obj, fl3_obj))

combined <- merge(x = all_objects[[1]],
                  y = all_objects[-1],
                  add.cell.ids = paste0("dataset", 1:length(all_objects)))

message("Re-normalizing combined data...")
combined <- NormalizeData(combined)

# STEP 5: Define nascent HSC genes
Nascent_HSC <- c("RUNX1","MLLT3","HOXA9","MECOM","HLF","SPINK2","MYB","GFI1",
                 "STAT5A","ZBTB16","HOPX","GATA2","GBP4","ITGA2B","KCNK17",
                 "SVOPL","C2orf88","SELP","CD82","ITGA4","GP9","TMEM163",
                 "RAB27B","SMIM24","GMPR","PDLIM1","ALDH1A1","NRGN","CCDC173",
                 "CXCL3","CYTL1","PRSS57","ANGPT1","CD34","PECAM1","CDH5",
                 "ECSCR","CALCRL","PROCR","ESAM","TIE1","EMCN")

nascent_genes_available <- Nascent_HSC[Nascent_HSC %in% rownames(combined)]
message(paste("Using", length(nascent_genes_available), "nascent HSC genes"))

# STEP 6: Define groups
agm_samples <- grep("agm-", unique(combined$orig.ident), value = TRUE)
somarin_samples <- c("FL1_hpc_CS22", "FL2_hpc_CS22", "FL_CS16_W9")
bm_samples <- grep("PCW.*BM.*\\.h5", unique(combined$orig.ident), value = TRUE)
spleen_samples <- grep("hSP.*\\.h5", unique(combined$orig.ident), value = TRUE)

# STEP 7: Process each stage - WITH DUPLICATE REMOVAL
process_stage <- function(seurat_obj, stage_name, resolution = 0.2) {
  message(paste("\n========== Processing:", stage_name, "=========="))
  
  # CRITICAL: re-normalize on the subset only
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  
  seurat_obj <- RunPCA(seurat_obj, npcs = 15, verbose = FALSE)
  
  pdf("elbow_plots.pdf", width = 10, height = 6)
  
  # duplicate removal (optional, but OK)
  pca_embeddings <- Embeddings(seurat_obj, "pca")[,1:15]
  duplicates <- duplicated(pca_embeddings)
  if (any(duplicates)) {
    message(paste("Removing", sum(duplicates), "duplicate cells"))
    seurat_obj <- seurat_obj[, !duplicates]
  }
  
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:15)
  seurat_obj$developmental_stage <- stage_name
  
  message(paste("Final cells:", ncol(seurat_obj)))
  return(seurat_obj)
}

message("\n=== COMPLETE ===")
library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(dplyr)

# =============================================================================
# LOAD FROM RAW FILES - NO OLD R OBJECTS
# =============================================================================

# STEP 1: Load AGM data from CSV files
message("Loading AGM CSV files...")

agm_files <- list.files(pattern = "GSM.*\\.csv\\.gz", full.names = TRUE)
message(paste("Found", length(agm_files), "AGM files"))

agm_objects <- lapply(agm_files, function(file) {
  message(paste("Reading:", basename(file)))
  data <- read.csv(file, row.names = 1)
  data_matrix <- as.matrix(data)
  data_sparse <- as(data_matrix, "sparseMatrix")
  
  # Extract sample name from filename
  sample_name <- gsub("GSM[0-9]+_Aorta-", "", basename(file))
  sample_name <- gsub(" \\(1\\)\\.csv\\.gz", "", sample_name)
  sample_name <- paste0("agm-", sample_name)
  
  obj <- CreateSeuratObject(counts = data_sparse, project = sample_name)
  obj$orig.ident <- sample_name
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  
  message(paste("Loaded:", sample_name))
  return(obj)
})

# STEP 2: Load Zheng .h5 files
message("\nLoading Zheng .h5 files...")
h5_files <- list.files(path = "Zheng", pattern = "\\.h5$", full.names = TRUE)
message(paste("Found", length(h5_files), ".h5 files"))

zheng_objects <- lapply(h5_files, function(file) {
  message(paste("Processing:", basename(file)))
  data <- Read10X_h5(file)
  obj <- CreateSeuratObject(counts = data, project = basename(file))
  obj$orig.ident <- basename(file)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  return(obj)
})

# STEP 3: Load Somarin MTX files
message("\nLoading Somarin fetal liver datasets...")

fl1_data <- Read10X(data.dir = "Somarin/FL1_hpc")
if(is.list(fl1_data)) fl1_data <- fl1_data$`Gene Expression`
fl1_obj <- CreateSeuratObject(counts = fl1_data, project = "FL1_CS22")
fl1_obj$orig.ident <- "FL1_hpc_CS22"
fl1_obj <- NormalizeData(fl1_obj)
fl1_obj <- FindVariableFeatures(fl1_obj)
fl1_obj <- ScaleData(fl1_obj)

fl2_data <- Read10X(data.dir = "Somarin/FL2_hpc")
if(is.list(fl2_data)) fl2_data <- fl2_data$`Gene Expression`
fl2_obj <- CreateSeuratObject(counts = fl2_data, project = "FL2_CS22")
fl2_obj$orig.ident <- "FL2_hpc_CS22"
fl2_obj <- NormalizeData(fl2_obj)
fl2_obj <- FindVariableFeatures(fl2_obj)
fl2_obj <- ScaleData(fl2_obj)

fl3_data <- Read10X(data.dir = "Somarin/FL_CS16_and_W9")
if(is.list(fl3_data)) fl3_data <- fl3_data$`Gene Expression`
fl3_obj <- CreateSeuratObject(counts = fl3_data, project = "FL_CS16_W9")
fl3_obj$orig.ident <- "FL_CS16_W9"
fl3_obj <- NormalizeData(fl3_obj)
fl3_obj <- FindVariableFeatures(fl3_obj)
fl3_obj <- ScaleData(fl3_obj)

# STEP 4: Merge everything
message("\nMerging all datasets...")
all_objects <- c(agm_objects, zheng_objects, list(fl1_obj, fl2_obj, fl3_obj))

combined <- merge(x = all_objects[[1]],
                  y = all_objects[-1],
                  add.cell.ids = paste0("dataset", 1:length(all_objects)))

message("Re-normalizing combined data...")
combined <- NormalizeData(combined)

# STEP 5: Define nascent HSC genes
Nascent_HSC <- c("RUNX1","MLLT3","HOXA9","MECOM","HLF","SPINK2","MYB","GFI1",
                 "STAT5A","ZBTB16","HOPX","GATA2","GBP4","ITGA2B","KCNK17",
                 "SVOPL","C2orf88","SELP","CD82","ITGA4","GP9","TMEM163",
                 "RAB27B","SMIM24","GMPR","PDLIM1","ALDH1A1","NRGN","CCDC173",
                 "CXCL3","CYTL1","PRSS57","ANGPT1","CD34","PECAM1","CDH5",
                 "ECSCR","CALCRL","PROCR","ESAM","TIE1","EMCN")

nascent_genes_available <- Nascent_HSC[Nascent_HSC %in% rownames(combined)]
message(paste("Using", length(nascent_genes_available), "nascent HSC genes"))

# STEP 6: Define groups
agm_samples <- grep("agm-", unique(combined$orig.ident), value = TRUE)
somarin_samples <- c("FL1_hpc_CS22", "FL2_hpc_CS22", "FL_CS16_W9")
bm_samples <- grep("PCW.*BM.*\\.h5", unique(combined$orig.ident), value = TRUE)
spleen_samples <- grep("hSP.*\\.h5", unique(combined$orig.ident), value = TRUE)

# STEP 7: Process each stage - WITH DUPLICATE REMOVAL
process_stage <- function(seurat_obj, stage_name, resolution = 0.2) {
  message(paste("\n========== Processing:", stage_name, "=========="))
  
  # CRITICAL: re-normalize on the subset only
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  
  seurat_obj <- RunPCA(seurat_obj, npcs = 15, verbose = FALSE)
  
  pdf(paste0(stage_name, "_elbow_plot.pdf"), width = 10, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 15))
  dev.off()
  
  # duplicate removal (optional, but OK)
  pca_embeddings <- Embeddings(seurat_obj, "pca")[,1:15]
  duplicates <- duplicated(pca_embeddings)
  if (any(duplicates)) {
    message(paste("Removing", sum(duplicates), "duplicate cells"))
    seurat_obj <- seurat_obj[, !duplicates]
  }
  
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:15)
  seurat_obj$developmental_stage <- stage_name
  
  message(paste("Final cells:", ncol(seurat_obj)))
  return(seurat_obj)
}

message("\n=== STAGE 1: AGM ===")
agm <- subset(combined, subset = orig.ident %in% agm_samples)
agm <- process_stage(agm, "AGM_Week4-5", resolution = 0.17)

message("\n=== STAGE 2: SOMARIN ===")
somarin <- subset(combined, subset = orig.ident %in% somarin_samples)
somarin <- process_stage(somarin, "FetalLiver_Week8-9", resolution = 0.2)

message("\n=== STAGE 3: BONE MARROW ===")
bone_marrow <- subset(combined, subset = orig.ident %in% bm_samples)
bone_marrow <- process_stage(bone_marrow, "BoneMarrow_Week10-14", resolution = 0.25)

message("\n=== STAGE 4: SPLEEN ===")
spleen <- subset(combined, subset = orig.ident %in% spleen_samples)
spleen <- process_stage(spleen, "Spleen_Week12-14", resolution = 0.25)

# STEP 8: Generate scorecards
message("\n=== GENERATING SCORECARDS ===")

pdf("Nascent_HSC_three_datasets.pdf", width = 16, height = 12)

p1 <- DimPlot(agm, reduction = "tsne", label = TRUE, pt.size = 1) +
  ggtitle("AGM (Week 4-5) - NASCENT HSCs")
print(p1)

p1_score <- DotPlot(agm, features = nascent_genes_available,
                    cols = c("grey90", "red3"),
                    group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("AGM - Nascent HSC Signature (EXPECT HIGH)")
print(p1_score)

p2 <- DimPlot(somarin, reduction = "tsne", label = TRUE, pt.size = 1) +
  ggtitle("Fetal Liver (Week 8-9)")
print(p2)

p2_score <- DotPlot(somarin, features = nascent_genes_available,
                    cols = c("grey90", "red3"),
                    group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Fetal Liver - Nascent HSC Signature")
print(p2_score)

p3 <- DimPlot(bone_marrow, reduction = "tsne", label = TRUE, pt.size = 0.8) +
  ggtitle("Bone Marrow (Week 10-14)")
print(p3)

p3_score <- DotPlot(bone_marrow, features = nascent_genes_available,
                    cols = c("grey90", "red3"),
                    group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Bone Marrow - Nascent HSC Signature (EXPECT LOW)")
print(p3_score)

p4 <- DimPlot(spleen, reduction = "tsne", label = TRUE, pt.size = 0.8) +
  ggtitle("Spleen (Week 12-14)")
print(p4)

p4_score <- DotPlot(spleen, features = nascent_genes_available,
                    cols = c("grey90", "red3"),
                    group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Spleen - Nascent HSC Signature (EXPECT LOW)")
print(p4_score)

dev.off()

message("\n=== Finding markers ===")
agm_markers <- FindAllMarkers(agm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(agm_markers, "AGM_markers.csv", row.names = FALSE)

somarin_markers <- FindAllMarkers(somarin, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(somarin_markers, "Somarin_markers.csv", row.names = FALSE)

bm_markers <- FindAllMarkers(bone_marrow, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(bm_markers, "BoneMarrow_markers.csv", row.names = FALSE)

spleen_markers <- FindAllMarkers(spleen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(spleen_markers, "Spleen_markers.csv", row.names = FALSE)

saveRDS(agm, "agm_processed.rds")
saveRDS(somarin, "somarin_processed.rds")
saveRDS(bone_marrow, "bone_marrow_processed.rds")
saveRDS(spleen, "spleen_processed.rds")

message("\n=== COMPLETE ===")
library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(dplyr)

# =============================================================================
# LOAD FROM RAW FILES - NO OLD R OBJECTS
# =============================================================================

# STEP 1: Load AGM data from CSV files
message("Loading AGM CSV files...")

agm_files <- list.files(pattern = "GSM.*\\.csv\\.gz", full.names = TRUE)
message(paste("Found", length(agm_files), "AGM files"))

agm_objects <- lapply(agm_files, function(file) {
  message(paste("Reading:", basename(file)))
  data <- read.csv(file, row.names = 1)
  data_matrix <- as.matrix(data)
  data_sparse <- as(data_matrix, "sparseMatrix")
  
  # Extract sample name from filename
  sample_name <- gsub("GSM[0-9]+_Aorta-", "", basename(file))
  sample_name <- gsub(" \\(1\\)\\.csv\\.gz", "", sample_name)
  sample_name <- paste0("agm-", sample_name)
  
  obj <- CreateSeuratObject(counts = data_sparse, project = sample_name)
  obj$orig.ident <- sample_name
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  
  message(paste("Loaded:", sample_name))
  return(obj)
})

# STEP 2: Load Zheng .h5 files
message("\nLoading Zheng .h5 files...")
h5_files <- list.files(path = "Zheng", pattern = "\\.h5$", full.names = TRUE)
message(paste("Found", length(h5_files), ".h5 files"))

zheng_objects <- lapply(h5_files, function(file) {
  message(paste("Processing:", basename(file)))
  data <- Read10X_h5(file)
  obj <- CreateSeuratObject(counts = data, project = basename(file))
  obj$orig.ident <- basename(file)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  return(obj)
})

# STEP 3: Load Somarin MTX files
message("\nLoading Somarin fetal liver datasets...")

fl1_data <- Read10X(data.dir = "Somarin/FL1_hpc")
if(is.list(fl1_data)) fl1_data <- fl1_data$`Gene Expression`
fl1_obj <- CreateSeuratObject(counts = fl1_data, project = "FL1_CS22")
fl1_obj$orig.ident <- "FL1_hpc_CS22"
fl1_obj <- NormalizeData(fl1_obj)
fl1_obj <- FindVariableFeatures(fl1_obj)
fl1_obj <- ScaleData(fl1_obj)

fl2_data <- Read10X(data.dir = "Somarin/FL2_hpc")
if(is.list(fl2_data)) fl2_data <- fl2_data$`Gene Expression`
fl2_obj <- CreateSeuratObject(counts = fl2_data, project = "FL2_CS22")
fl2_obj$orig.ident <- "FL2_hpc_CS22"
fl2_obj <- NormalizeData(fl2_obj)
fl2_obj <- FindVariableFeatures(fl2_obj)
fl2_obj <- ScaleData(fl2_obj)

fl3_data <- Read10X(data.dir = "Somarin/FL_CS16_and_W9")
if(is.list(fl3_data)) fl3_data <- fl3_data$`Gene Expression`
fl3_obj <- CreateSeuratObject(counts = fl3_data, project = "FL_CS16_W9")
fl3_obj$orig.ident <- "FL_CS16_W9"
fl3_obj <- NormalizeData(fl3_obj)
fl3_obj <- FindVariableFeatures(fl3_obj)
fl3_obj <- ScaleData(fl3_obj)

# STEP 4: Merge everything
message("\nMerging all datasets...")
all_objects <- c(agm_objects, zheng_objects, list(fl1_obj, fl2_obj, fl3_obj))

combined <- merge(x = all_objects[[1]],
                  y = all_objects[-1],
                  add.cell.ids = paste0("dataset", 1:length(all_objects)))

message("Re-normalizing combined data...")
combined <- NormalizeData(combined)

# STEP 6: Define groups
agm_samples <- grep("agm-", unique(combined$orig.ident), value = TRUE)
somarin_samples <- c("FL1_hpc_CS22", "FL2_hpc_CS22", "FL_CS16_W9")
bm_samples <- grep("PCW.*BM.*\\.h5", unique(combined$orig.ident), value = TRUE)
spleen_samples <- grep("hSP.*\\.h5", unique(combined$orig.ident), value = TRUE)

# STEP 7: Process each stage - ONLY ELBOW PLOTS
process_stage_elbow <- function(seurat_obj, stage_name) {
  message(paste("\n========== Processing:", stage_name, "=========="))

  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 15, verbose = FALSE)

  pdf(paste0(stage_name, "_elbow_plot.pdf"), width = 10, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 15))
  dev.off()

  message(paste("Elbow plot saved:", paste0(stage_name, "_elbow_plot.pdf")))
}

message("\n=== STAGE 1: AGM ===")
agm <- subset(combined, subset = orig.ident %in% agm_samples)
process_stage_elbow(agm, "AGM_Week4-5")

message("\n=== STAGE 2: SOMARIN ===")
somarin <- subset(combined, subset = orig.ident %in% somarin_samples)
process_stage_elbow(somarin, "FetalLiver_Week8-9")

message("\n=== STAGE 3: BONE MARROW ===")
bone_marrow <- subset(combined, subset = orig.ident %in% bm_samples)
process_stage_elbow(bone_marrow, "BoneMarrow_Week10-14")

message("\n=== STAGE 4: SPLEEN ===")
spleen <- subset(combined, subset = orig.ident %in% spleen_samples)
process_stage_elbow(spleen, "Spleen_Week12-14")

message("\n=== COMPLETE - Elbow plots generated ===")
library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(dplyr)

# =============================================================================
# LOAD FROM RAW FILES - NO OLD R OBJECTS
# =============================================================================

# STEP 1: Load AGM data from CSV files
message("Loading AGM CSV files...")

agm_files <- list.files(pattern = "GSM.*\\.csv\\.gz", full.names = TRUE)
message(paste("Found", length(agm_files), "AGM files"))

agm_objects <- lapply(agm_files, function(file) {
  message(paste("Reading:", basename(file)))
  data <- read.csv(file, row.names = 1)
  data_matrix <- as.matrix(data)
  data_sparse <- as(data_matrix, "sparseMatrix")
  
  # Extract sample name from filename
  sample_name <- gsub("GSM[0-9]+_Aorta-", "", basename(file))
  sample_name <- gsub(" \\(1\\)\\.csv\\.gz", "", sample_name)
  sample_name <- paste0("agm-", sample_name)
  
  obj <- CreateSeuratObject(counts = data_sparse, project = sample_name)
  obj$orig.ident <- sample_name
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  
  message(paste("Loaded:", sample_name))
  return(obj)
})

# STEP 2: Load Zheng .h5 files
message("\nLoading Zheng .h5 files...")
h5_files <- list.files(path = "Zheng", pattern = "\\.h5$", full.names = TRUE)
message(paste("Found", length(h5_files), ".h5 files"))

zheng_objects <- lapply(h5_files, function(file) {
  message(paste("Processing:", basename(file)))
  data <- Read10X_h5(file)
  obj <- CreateSeuratObject(counts = data, project = basename(file))
  obj$orig.ident <- basename(file)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  return(obj)
})

# STEP 3: Load Somarin MTX files
message("\nLoading Somarin fetal liver datasets...")

fl1_data <- Read10X(data.dir = "Somarin/FL1_hpc")
if(is.list(fl1_data)) fl1_data <- fl1_data$`Gene Expression`
fl1_obj <- CreateSeuratObject(counts = fl1_data, project = "FL1_CS22")
fl1_obj$orig.ident <- "FL1_hpc_CS22"
fl1_obj <- NormalizeData(fl1_obj)
fl1_obj <- FindVariableFeatures(fl1_obj)
fl1_obj <- ScaleData(fl1_obj)

fl2_data <- Read10X(data.dir = "Somarin/FL2_hpc")
if(is.list(fl2_data)) fl2_data <- fl2_data$`Gene Expression`
fl2_obj <- CreateSeuratObject(counts = fl2_data, project = "FL2_CS22")
fl2_obj$orig.ident <- "FL2_hpc_CS22"
fl2_obj <- NormalizeData(fl2_obj)
fl2_obj <- FindVariableFeatures(fl2_obj)
fl2_obj <- ScaleData(fl2_obj)

fl3_data <- Read10X(data.dir = "Somarin/FL_CS16_and_W9")
if(is.list(fl3_data)) fl3_data <- fl3_data$`Gene Expression`
fl3_obj <- CreateSeuratObject(counts = fl3_data, project = "FL_CS16_W9")
fl3_obj$orig.ident <- "FL_CS16_W9"
fl3_obj <- NormalizeData(fl3_obj)
fl3_obj <- FindVariableFeatures(fl3_obj)
fl3_obj <- ScaleData(fl3_obj)

# STEP 4: Merge everything
message("\nMerging all datasets...")
all_objects <- c(agm_objects, zheng_objects, list(fl1_obj, fl2_obj, fl3_obj))

combined <- merge(x = all_objects[[1]],
                  y = all_objects[-1],
                  add.cell.ids = paste0("dataset", 1:length(all_objects)))

message("Re-normalizing combined data...")
combined <- NormalizeData(combined)

# STEP 6: Define groups
agm_samples <- grep("agm-", unique(combined$orig.ident), value = TRUE)
somarin_samples <- c("FL1_hpc_CS22", "FL2_hpc_CS22", "FL_CS16_W9")
bm_samples <- grep("PCW.*BM.*\\.h5", unique(combined$orig.ident), value = TRUE)
spleen_samples <- grep("hSP.*\\.h5", unique(combined$orig.ident), value = TRUE)

# STEP 7: Process each stage - ONLY ELBOW PLOTS
process_stage_elbow <- function(seurat_obj, stage_name) {
  message(paste("\n========== Processing:", stage_name, "=========="))

  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 15, verbose = FALSE)

  pdf(paste0(stage_name, "_elbow_plot.pdf"), width = 10, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 15))
  dev.off()

  message(paste("Elbow plot saved:", paste0(stage_name, "_elbow_plot.pdf")))
}

message("\n=== STAGE 1: AGM ===")
agm <- subset(combined, subset = orig.ident %in% agm_samples)
process_stage_elbow(agm, "AGM_Week4-5")

message("\n=== STAGE 2: SOMARIN ===")
somarin <- subset(combined, subset = orig.ident %in% somarin_samples)
process_stage_elbow(somarin, "FetalLiver_Week8-9")

message("\n=== STAGE 3: BONE MARROW ===")
bone_marrow <- subset(combined, subset = orig.ident %in% bm_samples)
process_stage_elbow(bone_marrow, "BoneMarrow_Week10-14")

message("\n=== STAGE 4: SPLEEN ===")
spleen <- subset(combined, subset = orig.ident %in% spleen_samples)
process_stage_elbow(spleen, "Spleen_Week12-14")

message("\n=== COMPLETE - Elbow plots generated ===")
library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(dplyr)

# =============================================================================
# LOAD FROM RAW FILES - NO OLD R OBJECTS
# =============================================================================

# STEP 1: Load AGM data from CSV files
message("Loading AGM CSV files...")

agm_files <- list.files(pattern = "GSM.*\\.csv\\.gz", full.names = TRUE)
message(paste("Found", length(agm_files), "AGM files"))

agm_objects <- lapply(agm_files, function(file) {
  message(paste("Reading:", basename(file)))
  data <- read.csv(file, row.names = 1)
  data_matrix <- as.matrix(data)
  data_sparse <- as(data_matrix, "sparseMatrix")
  
  # Extract sample name from filename
  sample_name <- gsub("GSM[0-9]+_Aorta-", "", basename(file))
  sample_name <- gsub(" \\(1\\)\\.csv\\.gz", "", sample_name)
  sample_name <- paste0("agm-", sample_name)
  
  obj <- CreateSeuratObject(counts = data_sparse, project = sample_name)
  obj$orig.ident <- sample_name
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  
  message(paste("Loaded:", sample_name))
  return(obj)
})

# STEP 2: Load Zheng .h5 files
message("\nLoading Zheng .h5 files...")
h5_files <- list.files(path = "Zheng", pattern = "\\.h5$", full.names = TRUE)
message(paste("Found", length(h5_files), ".h5 files"))

zheng_objects <- lapply(h5_files, function(file) {
  message(paste("Processing:", basename(file)))
  data <- Read10X_h5(file)
  obj <- CreateSeuratObject(counts = data, project = basename(file))
  obj$orig.ident <- basename(file)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  return(obj)
})

# STEP 3: Load Somarin MTX files
message("\nLoading Somarin fetal liver datasets...")

fl1_data <- Read10X(data.dir = "Somarin/FL1_hpc")
if(is.list(fl1_data)) fl1_data <- fl1_data$`Gene Expression`
fl1_obj <- CreateSeuratObject(counts = fl1_data, project = "FL1_CS22")
fl1_obj$orig.ident <- "FL1_hpc_CS22"
fl1_obj <- NormalizeData(fl1_obj)
fl1_obj <- FindVariableFeatures(fl1_obj)
fl1_obj <- ScaleData(fl1_obj)

fl2_data <- Read10X(data.dir = "Somarin/FL2_hpc")
if(is.list(fl2_data)) fl2_data <- fl2_data$`Gene Expression`
fl2_obj <- CreateSeuratObject(counts = fl2_data, project = "FL2_CS22")
fl2_obj$orig.ident <- "FL2_hpc_CS22"
fl2_obj <- NormalizeData(fl2_obj)
fl2_obj <- FindVariableFeatures(fl2_obj)
fl2_obj <- ScaleData(fl2_obj)

fl3_data <- Read10X(data.dir = "Somarin/FL_CS16_and_W9")
if(is.list(fl3_data)) fl3_data <- fl3_data$`Gene Expression`
fl3_obj <- CreateSeuratObject(counts = fl3_data, project = "FL_CS16_W9")
fl3_obj$orig.ident <- "FL_CS16_W9"
fl3_obj <- NormalizeData(fl3_obj)
fl3_obj <- FindVariableFeatures(fl3_obj)
fl3_obj <- ScaleData(fl3_obj)

# STEP 4: Merge everything
message("\nMerging all datasets...")
all_objects <- c(agm_objects, zheng_objects, list(fl1_obj, fl2_obj, fl3_obj))

combined <- merge(x = all_objects[[1]],
                  y = all_objects[-1],
                  add.cell.ids = paste0("dataset", 1:length(all_objects)))

message("Re-normalizing combined data...")
combined <- NormalizeData(combined)

# STEP 6: Define groups
agm_samples <- grep("agm-", unique(combined$orig.ident), value = TRUE)
somarin_samples <- c("FL1_hpc_CS22", "FL2_hpc_CS22", "FL_CS16_W9")
bm_samples <- grep("PCW.*BM.*\\.h5", unique(combined$orig.ident), value = TRUE)
spleen_samples <- grep("hSP.*\\.h5", unique(combined$orig.ident), value = TRUE)

# STEP 7: Process each stage - ONLY ELBOW PLOTS
process_stage_elbow <- function(seurat_obj, stage_name) {
  message(paste("\n========== Processing:", stage_name, "=========="))

  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 15, verbose = FALSE)

  pdf(paste0(stage_name, "_elbow_plot.pdf"), width = 10, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 15))
  dev.off()

  message(paste("Elbow plot saved:", paste0(stage_name, "_elbow_plot.pdf")))
}

message("\n=== STAGE 1: AGM ===")
agm <- subset(combined, subset = orig.ident %in% agm_samples)
process_stage_elbow(agm, "AGM_Week4-5")

message("\n=== STAGE 2: SOMARIN ===")
somarin <- subset(combined, subset = orig.ident %in% somarin_samples)
process_stage_elbow(somarin, "FetalLiver_Week8-9")

message("\n=== STAGE 3: BONE MARROW ===")
bone_marrow <- subset(combined, subset = orig.ident %in% bm_samples)
process_stage_elbow(bone_marrow, "BoneMarrow_Week10-14")

message("\n=== STAGE 4: SPLEEN ===")
spleen <- subset(combined, subset = orig.ident %in% spleen_samples)
process_stage_elbow(spleen, "Spleen_Week12-14")

message("\n=== COMPLETE - Elbow plots generated ===")
library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(dplyr)

# =============================================================================
# LOAD FROM RAW FILES - NO OLD R OBJECTS
# =============================================================================

# STEP 1: Load AGM data from CSV files
message("Loading AGM CSV files...")

agm_files <- list.files(pattern = "GSM.*\\.csv\\.gz", full.names = TRUE)
message(paste("Found", length(agm_files), "AGM files"))

agm_objects <- lapply(agm_files, function(file) {
  message(paste("Reading:", basename(file)))
  data <- read.csv(file, row.names = 1)
  data_matrix <- as.matrix(data)
  data_sparse <- as(data_matrix, "sparseMatrix")
  
  # Extract sample name from filename
  sample_name <- gsub("GSM[0-9]+_Aorta-", "", basename(file))
  sample_name <- gsub(" \\(1\\)\\.csv\\.gz", "", sample_name)
  sample_name <- paste0("agm-", sample_name)
  
  obj <- CreateSeuratObject(counts = data_sparse, project = sample_name)
  obj$orig.ident <- sample_name
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  
  message(paste("Loaded:", sample_name))
  return(obj)
})

# STEP 2: Load Zheng .h5 files
message("\nLoading Zheng .h5 files...")
h5_files <- list.files(path = "Zheng", pattern = "\\.h5$", full.names = TRUE)
message(paste("Found", length(h5_files), ".h5 files"))

zheng_objects <- lapply(h5_files, function(file) {
  message(paste("Processing:", basename(file)))
  data <- Read10X_h5(file)
  obj <- CreateSeuratObject(counts = data, project = basename(file))
  obj$orig.ident <- basename(file)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  return(obj)
})

# STEP 3: Load Somarin MTX files
message("\nLoading Somarin fetal liver datasets...")

fl1_data <- Read10X(data.dir = "Somarin/FL1_hpc")
if(is.list(fl1_data)) fl1_data <- fl1_data$`Gene Expression`
fl1_obj <- CreateSeuratObject(counts = fl1_data, project = "FL1_CS22")
fl1_obj$orig.ident <- "FL1_hpc_CS22"
fl1_obj <- NormalizeData(fl1_obj)
fl1_obj <- FindVariableFeatures(fl1_obj)
fl1_obj <- ScaleData(fl1_obj)

fl2_data <- Read10X(data.dir = "Somarin/FL2_hpc")
if(is.list(fl2_data)) fl2_data <- fl2_data$`Gene Expression`
fl2_obj <- CreateSeuratObject(counts = fl2_data, project = "FL2_CS22")
fl2_obj$orig.ident <- "FL2_hpc_CS22"
fl2_obj <- NormalizeData(fl2_obj)
fl2_obj <- FindVariableFeatures(fl2_obj)
fl2_obj <- ScaleData(fl2_obj)

fl3_data <- Read10X(data.dir = "Somarin/FL_CS16_and_W9")
if(is.list(fl3_data)) fl3_data <- fl3_data$`Gene Expression`
fl3_obj <- CreateSeuratObject(counts = fl3_data, project = "FL_CS16_W9")
fl3_obj$orig.ident <- "FL_CS16_W9"
fl3_obj <- NormalizeData(fl3_obj)
fl3_obj <- FindVariableFeatures(fl3_obj)
fl3_obj <- ScaleData(fl3_obj)

# STEP 4: Merge everything
message("\nMerging all datasets...")
all_objects <- c(agm_objects, zheng_objects, list(fl1_obj, fl2_obj, fl3_obj))

combined <- merge(x = all_objects[[1]],
                  y = all_objects[-1],
                  add.cell.ids = paste0("dataset", 1:length(all_objects)))

message("Re-normalizing combined data...")
combined <- NormalizeData(combined)

# STEP 6: Define groups
agm_samples <- grep("agm-", unique(combined$orig.ident), value = TRUE)
somarin_samples <- c("FL1_hpc_CS22", "FL2_hpc_CS22", "FL_CS16_W9")
bm_samples <- grep("PCW.*BM.*\\.h5", unique(combined$orig.ident), value = TRUE)
spleen_samples <- grep("hSP.*\\.h5", unique(combined$orig.ident), value = TRUE)

# STEP 7: Process each stage - ONLY ELBOW PLOTS
process_stage_elbow <- function(seurat_obj, stage_name) {
  message(paste("\n========== Processing:", stage_name, "=========="))
  
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 15, verbose = FALSE)
  
  pdf(paste0(stage_name, "_elbow_plot.pdf"), width = 10, height = 6)
  print(ElbowPlot(seurat_obj, ndims = 15))
  dev.off()
  
  message(paste("Elbow plot saved:", paste0(stage_name, "_elbow_plot.pdf")))
}

message("\n=== STAGE 1: AGM ===")
agm <- subset(combined, subset = orig.ident %in% agm_samples)
process_stage_elbow(agm, "AGM_Week4-5")

message("\n=== STAGE 2: SOMARIN ===")
somarin <- subset(combined, subset = orig.ident %in% somarin_samples)
process_stage_elbow(somarin, "FetalLiver_Week8-9")

message("\n=== STAGE 3: BONE MARROW ===")
bone_marrow <- subset(combined, subset = orig.ident %in% bm_samples)
process_stage_elbow(bone_marrow, "BoneMarrow_Week10-14")

message("\n=== STAGE 4: SPLEEN ===")
spleen <- subset(combined, subset = orig.ident %in% spleen_samples)
process_stage_elbow(spleen, "Spleen_Week12-14")

message("\n=== COMPLETE - Elbow plots generated ===")
