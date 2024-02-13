#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 2022/02/22
# Description: CAR T project colors and themes
#==============================================================================#

library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(nord)
library(plyr)
library(circlize)
library(googlesheets4)

#==============================================================================
# Colors and themes
#==============================================================================

# Colors for plotting
# Define colors for each level of categorical variables

# Clusters
product_clusters <- as.factor(c(0, seq(1:13)))
product_clusters_woutsmall <- as.factor(c(0, seq(1:11)))
product_cd8_clusters <- as.factor(c(0, seq(1:7)))
leukPBMC_clusters <- as.factor(c(0, seq(1:18)))
leukPBMC_Tcell_clusters <- as.factor(c(0, seq(1:7)))
tumor_clusters <- as.factor(c(0, seq(1:20)))
tumor_15_clusters <- as.factor(c(0, seq(1:14)))
leukPBMC_myeloid_clusters <- as.factor(c(0, seq(1:15)))
matched_csf_clusters <- as.factor(c(0, seq(1:12)))
matched_pbmc_clusters <- as.factor(c(0, seq(1:14)))
csf_pbmc_clusters <- as.factor(c(0, seq(1:10)))
csf_clusters <- as.factor(c(0, seq(1:12)))
pbmc_clusters <- as.factor(c(0, seq(1:15)))
upn109_clusters <- as.factor(c(0, seq(1:15)))

upn109_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(upn109_clusters))
names(upn109_cluster_col) <- levels(upn109_clusters)

csf_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(csf_clusters))
names(csf_cluster_col) <- levels(csf_clusters)

pbmc_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(pbmc_clusters))
names(pbmc_cluster_col) <- levels(pbmc_clusters)

product_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(product_clusters))
names(product_cluster_col) <- levels(product_clusters)

product_cluster_col_woutsmall <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(product_clusters_woutsmall))
names(product_cluster_col_woutsmall) <- levels(product_clusters_woutsmall)

product_cd8_cluster_col <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols <- length(product_cd8_clusters))
names(product_cd8_cluster_col) <- levels(product_cd8_clusters)

leukPBMC_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(length(leukPBMC_clusters))
names(leukPBMC_cluster_col) <- levels(leukPBMC_clusters)

leukPBMC_Tcell_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(length(leukPBMC_Tcell_clusters))
names(leukPBMC_Tcell_cluster_col) <- levels(leukPBMC_Tcell_clusters)

tumor_cluster_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(tumor_clusters))
names(tumor_cluster_col) <- levels(tumor_clusters)

tumor_15_cluster_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(tumor_15_clusters))
names(tumor_15_cluster_col) <- levels(tumor_15_clusters)

leukPBMC_myeloid_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(length(leukPBMC_myeloid_clusters))
names(leukPBMC_myeloid_cluster_col) <- levels(leukPBMC_myeloid_clusters)

matched_csf_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(length(matched_csf_clusters))
names(matched_csf_cluster_col) <- levels(matched_csf_clusters)

matched_pbmc_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(length(matched_pbmc_clusters))
names(matched_pbmc_cluster_col) <- levels(matched_pbmc_clusters)

csf_pbmc_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(length(csf_pbmc_clusters))
names(csf_pbmc_cluster_col) <- levels(csf_pbmc_clusters)

annotation_level_1 <- c("Fibroblasts", "T cells", "Myeloid", "Non-immune", "NA")
annotation_level_1_col <- colorRampPalette(brewer.pal(10, "Paired"))(length(annotation_level_1))
names(annotation_level_1_col) <- annotation_level_1

# CD3_high_low
CD3_high_low_col <- c("High" = brewer.pal(3, "Set1")[1],
                      "Low" = brewer.pal(3, "Set1")[2],
                      "NA" = "azure3")


# CD4/CD8
cd4_cd8_col <- c("CD8" = brewer.pal(3, "Set1")[1],
                 "CD4" = brewer.pal(3, "Set1")[2])

# Response
response_col <- c("NA" = nord("silver_mine", 5)[3],
                  "Progression Disease (PD)" = nord("aurora", 5)[1],
                  "Partial Response (PR)" = nord("aurora", 5)[3],
                  "Stable Disease (SD)" = nord("aurora", 5)[5],
                  "Complete Response (CR)" = nord("aurora", 5)[4])

# Batch
batch_col <- colorRampPalette(brewer.pal(10, "Spectral"))(19)
names(batch_col) <- paste0("Batch", seq(1,19))
colScale_batch <- scale_fill_manual(name = "Batch", values = batch_col)

simple_batch_col <- colorRampPalette(brewer.pal(10, "Spectral"))(19)
names(simple_batch_col) <- seq(1,19)

# Diagnosis
diagnosis_col <- c(colorRampPalette(brewer.pal(10, "Spectral"))(10), "#666666")
diagnosis <- c("Glioblastoma, NOS", "Ependymoma, NOS", "Astrocytoma, anaplastic", 
               "Oligodendroglioma, anaplastic", "Primitive neuroectodermal tumor",                            
               "Astrocytoma, NOS", "Fibrillary astrocytoma", "Pleomorphic xanthoastrocytoma",
               "Oligodendroglioma, NOS", "Ependymoma, anaplastic", "NA")

names(diagnosis_col) <- diagnosis
colScale_diagnosis <- scale_fill_manual(name = "Diagnosis", values = diagnosis_col)


# Tumor UPN color
tumor_UPN_col <- colorRampPalette(brewer.pal(11, "Spectral"))(44)
names(tumor_UPN_col) <- sort(c("109", "243", "185", "277", "241", "315", "303", "282", "275", "289", "265", "125", "181", "234", "117", "232", "191", "129", "146", "157", "208", "145", "214", "237", "240", "218", "350", "228", "248", "266", "131", "239", "301", "215", "223", "224", "230", "260", "213", "193", "149", "141", "124", "122"))

csf_pbmc_UPN_col <- colorRampPalette(brewer.pal(12, "Paired"))(22)
names(csf_pbmc_UPN_col) <- c("109", "266", "265", "210", "254", "315", "350", "235", "239", "275", "338", "260", "282", "409", "303", "444", "362", "226", "289", "240", "221", "230")

# Cell cycle
cell_cycle_col <- c("G1" = nord("lumina", 5)[1],
                    "G2M" = nord("lumina", 5)[3],
                    "S" = nord("lumina", 5)[5])

# Manufacture type
#manufacture_col <- c("TCM" = nord("mountain_forms", 5)[2],
#                     "Tnmem" = nord("mountain_forms", 5)[3])
manufacture_col <- c("Tnmem" = "#2774B7",
                     "TCM" = "#F7A080")

# Antibody panel
abpanel_col <- c("TS-C_99328" = "bisque4",
                 "TS-C_99814" = "bisque2",
                 "TS-C_399905" = "bisque",
                 "nonCITEseq" = "snow1")

# TCR clone size
tcr_clonetype_col <- c("Single (0 < X <= 1)" = "lightyellow",
                       "Small (1 < X <= 5)" = "yellow2",
                       "Medium (5 < X <= 20)" = "orange",
                       "Large (20 < X <= 100)" = "red",
                       "Hyperexpanded (100 < X <= 500)" = "darkred")

# Canonical T cell markers
gs4_deauth()
markers  <- gs4_get("https://docs.google.com/spreadsheets/d/186PIcU3Jb0Xm3IuZgENEChluykvwsK4JnHF8FbeJASA/edit?usp=sharing")
sheet_names(markers)
canonical_markers <- read_sheet(markers, sheet = "T cells")
head(canonical_markers)

marker_col <- nord("afternoon_prarie", length(unique(canonical_markers$RNA)))
names(marker_col) <- unique(canonical_markers$RNA)

col <- list(product_cluster = product_cluster_col,
            product_cd8_cluster = product_cd8_cluster_col,
            leukPBMCs_cluster = leukPBMC_cluster_col,
            leukPBMC_Tcells_cluster = leukPBMC_Tcell_cluster_col,
            leukPBMC_myeloid_cluster = leukPBMC_myeloid_cluster_col,
            cd4_cd8_cluster = cd4_cd8_col,
            response = response_col,
            batch = simple_batch_col,
            cell_cycle_phase = cell_cycle_col,
            manufacture = manufacture_col,
            ab_panel = abpanel_col,
            marker_type = marker_col)

# Heatmap colors
col_fun = colorRamp2(c(-2, 0, 2), c("cadetblue4", "white", "coral2"))

# Theme for plotting
my_theme <- theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  #legend.position = "none",
                  plot.title = element_text(size = 12))

manuscript_theme <- theme(axis.text = element_text(size = 6),
                          axis.title = element_text(size = 7),
                          plot.title = element_text(size = 7),
                          legend.text = element_text(size = 6),
                          legend.title = element_text(size = 7),
                          strip.text.x = element_text(size = 7))

slide_theme <- theme(axis.text = element_text(size = 9),
                     axis.title = element_text(size = 10),
                     plot.title = element_text(size = 10),
                     legend.text = element_text(size = 9),
                     legend.title = element_text(size = 9))
