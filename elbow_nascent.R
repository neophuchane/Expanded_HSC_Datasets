  library(Seurat)
  library(Matrix)
  library(hdf5r)
  library(ggplot2)
  library(dplyr)

  # STEP 1: Load AGM data from CSV files
  message("Loading AGM CSV files...")
  agm_files <- list.files(pattern = "GSM.*\\.csv\\.gz", full.names = TRUE)
  message(paste("Found", length(agm_files), "AGM files"))

  agm_objects <- lapply(agm_files, function(file) {
    message(paste("Reading:", basename(file)))
    data <- read.csv(file, row.names = 1)
    data_matrix <- as.matrix(data)
    data_sparse <- as(data_matrix, "sparseMatrix")
    sample_name <- gsub("GSM[0-9]+_Aorta-", "", basename(file))
    sample_name <- gsub(" \\(1\\)\\.csv\\.gz", "", sample_name)
    sample_name <- paste0("agm-", sample_name)
    obj <- CreateSeuratObject(counts = data_sparse, project = sample_name)
    obj$orig.ident <- sample_name
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    return(obj)
  })

  # STEP 2: Load Zheng .h5 files
  message("\nLoading Zheng .h5 files...")
  h5_files <- list.files(path = "Zheng", pattern = "\\.h5$", full.names = TRUE)
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
  combined <- merge(x = all_objects[[1]], y = all_objects[-1],
                    add.cell.ids = paste0("dataset", 1:length(all_objects)))
  combined <- NormalizeData(combined)

  # STEP 5: Define groups
  agm_samples <- grep("agm-", unique(combined$orig.ident), value = TRUE)
  somarin_samples <- c("FL1_hpc_CS22", "FL2_hpc_CS22", "FL_CS16_W9")
  bm_samples <- grep("PCW.*BM.*\\.h5", unique(combined$orig.ident), value = TRUE)
  spleen_samples <- grep("hSP.*\\.h5", unique(combined$orig.ident), value = TRUE)

  # STEP 6: Generate elbow plots only
  process_stage_elbow <- function(seurat_obj, stage_name) {
    message(paste("\n=== Processing:", stage_name, "==="))
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, npcs = 15, verbose = FALSE)

    pdf(paste0(stage_name, "_elbow_plot.pdf"), width = 10, height = 6)
    print(ElbowPlot(seurat_obj, ndims = 15))
    dev.off()
    message(paste("Saved:", paste0(stage_name, "_elbow_plot.pdf")))
  }

  agm <- subset(combined, subset = orig.ident %in% agm_samples)
  process_stage_elbow(agm, "AGM_Week4-5")

  somarin <- subset(combined, subset = orig.ident %in% somarin_samples)
  process_stage_elbow(somarin, "FetalLiver_Week8-9")

  bone_marrow <- subset(combined, subset = orig.ident %in% bm_samples)
  process_stage_elbow(bone_marrow, "BoneMarrow_Week10-14")

  spleen <- subset(combined, subset = orig.ident %in% spleen_samples)
  process_stage_elbow(spleen, "Spleen_Week12-14")

  message("\n=== COMPLETE ===")
