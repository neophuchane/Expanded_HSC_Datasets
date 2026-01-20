#Mature Scorecard


library(Matrix)
library(Seurat)
library(SeuratObject)
library(hdf5r)
library(ggplot2)

# Checkpoint 1: Load existing scorecard
message("Loading existing Seurat object...")
load("~/hsc_project/seurat_object_23.Rdata")

# FIX: Remove Graph/neighbor slots BEFORE updating to prevent errors
message("Removing problematic Graph/neighbor slots...")
sample@graphs <- list()
sample@neighbors <- list()
sample@commands <- list()
message("Object loaded")

# Checkpoint 2: Define HSC gene list
message("Defining HSC maturation gene list...")
HSC_maturation <- c("CDH5","MEIS2","RUNX1T1","ESAM","PLVAP","SELP","ITGA2B",
                    "GP9","RAB27B","GMPR","GFI1","GBP4","MECOM","IGFBP2","HMGA2",
                    "LIN28B","IL3RA","IL6R","IL11RA","IFNAR1","IFNAR2","IFNGR1",
                    "CSF1R","CSF2RA","MKI67","TOP2A","AURKB","HLA-A","HLA-B",
                    "HLA-C","HLA-E","B2M","HLA-DMA","HLA-DPB1","HLA-DRA",
                    "HLA-DQA1","HLA-DPA1","HLA-DQB1","MLLT3","HLF","MALAT1",
                    "MSI2","EVI2B","SOCS2","HEMGN","HOPX","SPINK2","CD52","SELL","PROM1")
message("Gene list defined.")

# Checkpoint 3: Load new .h5 files
message("Locating new .h5 files...")
h5_files <- list.files(path = "~/hsc_project/Zheng",
                       pattern = "\\.h5$",
                       full.names = TRUE)
message(paste("Found", length(h5_files), ".h5 files."))

# Checkpoint 4: Load and preprocess .h5 datasets
message("Loading and preprocessing .h5 datasets...")
new_objects <- lapply(h5_files, function(file) {
  message(paste("Processing:", basename(file)))
  data <- Read10X_h5(file)
  obj <- CreateSeuratObject(counts = data, project = basename(file))
  obj$orig.ident <- basename(file)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj@graphs <- list()
  obj@neighbors <- list()
  message(paste("Finished:", basename(file)))
  return(obj)
})
message("All .h5 datasets loaded and preprocessed.")

# NEW: Checkpoint 5: Load Somarin fetal liver datasets (MTX format)
message("Loading Somarin fetal liver datasets...")

# FL1
message("Loading FL1_CS22...")
fl1_data <- Read10X(data.dir = "~/hsc_project/Somarin/FL1_hpc")
if(is.list(fl1_data)) fl1_data <- fl1_data$`Gene Expression`
fl1_obj <- CreateSeuratObject(counts = fl1_data, project = "FL1_CS22_FetalLiver")
fl1_obj$orig.ident <- "FL1_hpc_CS22"
fl1_obj <- NormalizeData(fl1_obj)
fl1_obj <- FindVariableFeatures(fl1_obj)
fl1_obj <- ScaleData(fl1_obj)
fl1_obj@graphs <- list()
fl1_obj@neighbors <- list()

# FL2
message("Loading FL2_hpc...")
fl2_data <- Read10X(data.dir = "~/hsc_project/Somarin/FL2_hpc")
if(is.list(fl2_data)) fl2_data <- fl2_data$`Gene Expression`
fl2_obj <- CreateSeuratObject(counts = fl2_data, project = "FL2_hpc_FetalLiver")
fl2_obj$orig.ident <- "FL2_hpc_CS22"
fl2_obj <- NormalizeData(fl2_obj)
fl2_obj <- FindVariableFeatures(fl2_obj)
fl2_obj <- ScaleData(fl2_obj)
fl2_obj@graphs <- list()
fl2_obj@neighbors <- list()

# FL3
message("Loading FL_CS16_W9...")
fl3_data <- Read10X(data.dir = "~/hsc_project/Somarin/FL_CS16_and_W9")
if(is.list(fl3_data)) fl3_data <- fl3_data$`Gene Expression`
fl3_obj <- CreateSeuratObject(counts = fl3_data, project = "FL_CS16_W9_FetalLiver")
fl3_obj$orig.ident <- "FL_CS16_W9"
fl3_obj <- NormalizeData(fl3_obj)
fl3_obj <- FindVariableFeatures(fl3_obj)
fl3_obj <- ScaleData(fl3_obj)
fl3_obj@graphs <- list()
fl3_obj@neighbors <- list()

message("All Somarin datasets loaded and preprocessed.")

# Add FL objects to the list
new_objects <- c(new_objects, list(fl1_obj, fl2_obj, fl3_obj))

# Checkpoint 6: Merge datasets
message("Merging all datasets...")
combined <- merge(x = sample,
                  y = new_objects,
                  add.cell.ids = c("original", paste0("new", 1:length(new_objects))))
message("Merge complete.")

# Checkpoint 7: Re-normalize the combined object
message("Re-normalizing combined data...")
combined <- NormalizeData(combined)
message("Re-normalization complete.")

# Checkpoint 7.5: Check current sample names
message("===========================================")
message("Current Sample Names (orig.ident values):")
message("===========================================")
sample_names <- unique(combined$orig.ident)
for(i in 1:length(sample_names)) {
  message(paste(i, ":", sample_names[i]))
}
message("===========================================")
message("Names for Developmental Order")
message("===========================================")

# Checkpoint 7.6: Set developmental order
message("Setting developmental order for samples...")

# UPDATED to match actual sample names in data
developmental_order <- c(
  # AGM - Weeks 4-5
  "agm-4wk-658",
  "agm-5wk-555",
  "agm-5wk-575",
  
  # Fetal Liver - Week 6
  "liver-6wk-563",
  
  # Fetal Liver - Week 8 (CS22)
  "liver-8wk-553",
  "FL1_hpc_CS22",
  "FL2_hpc_CS22",
  
  # Fetal Liver - Week 6-9 (mixed, positioned here due to W9 cells)
  "FL_CS16_W9",
  
  # Fetal Liver - Week 11
  "liver-11wk-569",
  
  # Fetal Liver - Week 15 (ADDED: was missing)
  "liver-15wk-101",
  
  # Bone Marrow - Weeks 10-14
  "PCW10_BM_1_10x.h5",
  "PCW10_BM_2_10x.h5",
  "PCW11_BM_1_10x.h5",
  "PCW11_BM_2_10x.h5",
  "PCW12_BM_1_10x.h5",
  "PCW12_BM_2_10x.h5",
  "PCW13_BM_1_10x.h5",
  "PCW13_BM_2_10x.h5",
  "PCW14_BM_10x.h5",
  
  # Spleen - Weeks 12-14
  "hSP_12w_1.h5",
  "hSP_12w_2.h5",
  "hSP_13w_1.h5",
  "hSP_13w_2.h5",
  "hSP_14w.h5",
  
  # Cord Blood - Week 40
  "cb-40wk-201",
  
  # Mouse Bone Marrow - E16.5-E18.5 (for comparison)
  "mBM_E16.5_1.h5",
  "mBM_E16.5_2.h5",
  "mBM_E17.5.h5",
  "mBM_E18.5.h5"
)



# Try to set the factor order
message("Attempting to set developmental order...")
if(all(developmental_order %in% unique(combined$orig.ident))) {
  combined$orig.ident <- factor(combined$orig.ident, levels = developmental_order)
  message("✓ Developmental order set successfully!")
} else {
  message("⚠ WARNING: Some sample names don't match.")
  message("Missing from combined data:")
  missing <- setdiff(developmental_order, unique(combined$orig.ident))
  for(m in missing) message(paste("  -", m))
  message("Extra samples in combined data (not in developmental_order):")
  extra <- setdiff(unique(combined$orig.ident), developmental_order)
  for(e in extra) message(paste("  -", e))
  message("")
  message(">>> ACTION REQUIRED: Update the developmental_order list in the script!")
  message(">>> Then re-run the script.")
}

# Filter gene list to only include genes present in data
message("Checking which genes are present in combined dataset...")
available_genes <- HSC_maturation[HSC_maturation %in% rownames(combined)]
missing_genes <- HSC_maturation[!HSC_maturation %in% rownames(combined)]
if(length(missing_genes) > 0) {
  message(paste("WARNING: Missing", length(missing_genes), "genes:", paste(missing_genes, collapse=", ")))
}
message(paste("Using", length(available_genes), "out of", length(HSC_maturation), "genes"))

# Checkpoint 8: Generate expanded HSC scorecard
message("Generating expanded HSC maturation scorecard...")
pdf("HSC_maturation_scorecard_expanded.pdf", width = 20, height = 12)
p <- DotPlot(combined,
             features = available_genes,
             cols = c("grey90","red3"),
             group.by = "orig.ident") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  ggtitle("HSC Maturation Scorecard - Developmental Timeline")
print(p)
dev.off()
message("Scorecard plot saved as PDF.")

# Save the combined object
message("Saving combined object...")
saveRDS(combined, file = "combined_expanded_maturation_scorecard.rds")
message("Combined object saved.")
message("Analysis complete!")