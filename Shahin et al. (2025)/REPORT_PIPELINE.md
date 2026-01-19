# 🧬 Systems Immunology of Malaria  
## Replication of Inter-Ethnic Immune Variation Using Single-Cell Transcriptomics

---

## 📌 Overview

This project is a **replication and systems-level interpretation** of the study:

**Shahin et al. (2025)**  
*Single-cell transcriptomics reveals inter-ethnic variation in immune response to Falciparum malaria*  
**The American Journal of Human Genetics**

Using publicly available scRNA-seq data, this analysis investigates **why Fulani children experience less severe malaria than Mossi children**, despite comparable parasite exposure.

The core hypothesis is that **baseline immune programming shaped by ethnicity**—rather than infection alone—determines disease outcomes.

---

## 🔬 Study Design

- **71,784 high-quality immune cells**
- **126 children** (Fulani and Mossi)
- Infected vs non-infected states
- Endemic malaria setting (Burkina Faso)
- Integration of **transcriptomics and genetics**
- Data source: GEO `GSE273781`, `GSE273785`

---

## 💡 Key Biological Insights

### 1️⃣ Ethnicity Dominates Immune Variation
Ethnicity drives **10–40× more gene expression changes** than malaria infection itself, indicating that immune outcomes are largely pre-configured before exposure.

---

### 2️⃣ Trained Tolerance in Fulani Monocytes
Fulani monocytes suppress inflammatory pathways (TNF-α, IL-6, NF-κB) while enhancing parasite clearance via phagocytosis (CD36).

**Outcome:** Efficient parasite removal without immunopathology.

---

### 3️⃣ Compartmentalized Immune Strategy

| Cell Type | Dominant Function |
|---------|------------------|
| Monocytes | Silent parasite clearance |
| B cells | Antibody production |
| CD4 T cells | Dampened activation |
| CD8 / NK | Preserved interferon signaling |

This challenges vaccine strategies that broadly amplify T-cell activation.

---

### 4️⃣ Population-Specific Genetic Regulation
Variants such as **CD36 (rs1049654)** and **MT2A (rs35402964)** exhibit **gene-by-ethnicity interactions**, meaning the same variant behaves differently across populations.

---

### 5️⃣ Metallothioneins as System-Wide Protectors
Pan-cellular upregulation of **MT1X** and **MT2A** supports:
- Antioxidant defense
- Zinc homeostasis
- Anti-inflammatory signaling
- Enhanced oxidative burst

---

## 🧪 Translational Implications

- Predictive biomarkers for severe malaria
- Population-aware vaccine design
- Precision nutrition (e.g. zinc modulation)
- Trained-tolerance-based immunotherapies

---

# 🧬 Analysis Pipeline (R / Seurat)
**All code below is part of a single, reproducible workflow**

---

```r
# =============================================================================
# Publication-Quality Figure Generation for scRNA-seq Analysis
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)

# Load processed Seurat object
seurat_combined <- readRDS("MoFu_unfiltered.rds")

# =============================================================================
# FIGURE 1: UMAP with Cell Type Annotations
# =============================================================================

cell_colors <- c(
  "B intermediate" = "#FF8C00",
  "B memory" = "#FF6500",
  "B naive" = "#FFD700",
  "CD14 Mono" = "#D2691E",
  "CD16 Mono" = "#CD853F",
  "CD4 Naive" = "#00BFFF",
  "CD4 TCM" = "#1E90FF",
  "CD4 TEM" = "#0066CC",
  "CD8 Naive" = "#00FF7F",
  "CD8 TCM" = "#00C957",
  "CD8 TEM" = "#008B45",
  "NK" = "#DA70D6",
  "pDC" = "#9932CC",
  "Platelet" = "#FF69B4"
)

p1 <- DimPlot(
  seurat_combined,
  reduction = "azimuth_umap",
  group.by = "predicted.celltype.l2",
  cols = cell_colors,
  pt.size = 0.1,
  label = TRUE,
  repel = TRUE
) + theme_classic()

ggsave("Fig1_UMAP_CellTypes.png", p1, width = 10, height = 8, dpi = 300)
ggsave("Fig1_UMAP_CellTypes.pdf", p1, width = 10, height = 8)

# =============================================================================
# FIGURE 2: UMAP Split by Infection Status
# =============================================================================

p2 <- DimPlot(
  seurat_combined,
  reduction = "azimuth_umap",
  split.by = "Infection",
  group.by = "Ethnicity",
  pt.size = 0.1,
  ncol = 2
) + theme_classic()

ggsave("Fig2_UMAP_Infection_Split.png", p2, width = 12, height = 6, dpi = 300)
ggsave("Fig2_UMAP_Infection_Split.pdf", p2, width = 12, height = 6)

# =============================================================================
# FIGURE 3: PCA of Monocyte Subsets
# =============================================================================

mono_subset <- subset(
  seurat_combined,
  subset = predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono")
)

mono_subset <- FindVariableFeatures(mono_subset, nfeatures = 2000)
mono_subset <- ScaleData(mono_subset)
mono_subset <- RunPCA(mono_subset, npcs = 30)

p3 <- DimPlot(
  mono_subset,
  reduction = "pca",
  group.by = "Ethnicity",
  shape.by = "predicted.celltype.l2"
) + theme_classic()

ggsave("Fig3_PCA_Monocytes.png", p3, width = 8, height = 6, dpi = 300)
ggsave("Fig3_PCA_Monocytes.pdf", p3, width = 8, height = 6)

# =============================================================================
# FIGURE 4: Differential Expression – CD14+ and CD16+ Monocytes
# =============================================================================

cd14_subset <- subset(seurat_combined, subset = predicted.celltype.l2 == "CD14 Mono")
Idents(cd14_subset) <- "Infection"
cd14_degs <- FindMarkers(cd14_subset, ident.1 = "Inf", ident.2 = "NI")

cd16_subset <- subset(seurat_combined, subset = predicted.celltype.l2 == "CD16 Mono")
Idents(cd16_subset) <- "Infection"
cd16_degs <- FindMarkers(cd16_subset, ident.1 = "Inf", ident.2 = "NI")

# =============================================================================
# FIGURE 5: Heatmap of Top DEGs (CD14+ Monocytes)
# =============================================================================

top_genes <- rownames(cd14_degs)[1:20]
cd14_subset <- ScaleData(cd14_subset, features = top_genes)

expr_mat <- as.matrix(cd14_subset[["RNA"]]@scale.data[top_genes, ])

col_anno <- HeatmapAnnotation(
  Infection = cd14_subset$Infection,
  col = list(Infection = c("NI" = "blue", "Inf" = "red"))
)

png("Fig4_Heatmap_CD14_Mono.png", width = 3000, height = 3000, res = 300)
Heatmap(
  expr_mat,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = col_anno,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_title = "CD14+ Monocytes: Infected vs Non-infected"
)
dev.off()

# =============================================================================
# FIGURE 6: Venn Diagram – Upregulated Genes
# =============================================================================

cd14_up <- rownames(subset(cd14_degs, avg_log2FC > 0 & p_val_adj < 0.05))
cd16_up <- rownames(subset(cd16_degs, avg_log2FC > 0 & p_val_adj < 0.05))

png("Fig5_Venn_Upregulated.png", width = 2000, height = 2000, res = 300)
venn.diagram(
  x = list("CD14+ Mono" = cd14_up, "CD16+ Mono" = cd16_up),
  filename = NULL,
  fill = c("#E74C3C", "#3498DB"),
  alpha = 0.5,
  main = "Upregulated Genes in Infected Monocytes"
) |> grid::grid.draw()
dev.off()

# =============================================================================
# FIGURE 7: Violin Plots – IFN Pathway Genes
# =============================================================================

ifn_genes <- c("IFNG", "STAT1", "STAT2", "IRF7", "MX1", "ISG15")
ifn_genes <- ifn_genes[ifn_genes %in% rownames(seurat_combined)]

mono_subset <- subset(
  seurat_combined,
  subset = predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono")
)

plots <- lapply(ifn_genes, function(g) {
  VlnPlot(
    mono_subset,
    features = g,
    split.by = "Infection",
    group.by = "predicted.celltype.l2",
    pt.size = 0
  )
})

p6 <- wrap_plots(plots, ncol = 3)
ggsave("Fig6_IFN_Violin.png", p6, width = 16, height = 10, dpi = 300)

# =============================================================================
# FIGURE 8: Cell Count and Proportion Comparisons
# =============================================================================

cell_counts <- seurat_combined@meta.data %>%
  group_by(predicted.celltype.l2, Infection) %>%
  summarise(count = n(), .groups = "drop")

p7 <- ggplot(cell_counts, aes(predicted.celltype.l2, count, fill = Infection)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic()

ggsave("Fig7_CellCounts.png", p7, width = 12, height = 6, dpi = 300)

cat("\n=== Analysis Complete ===\n")

```

---
## 📚 Reference

Shahin et al. (2025).  
*Single-cell transcriptomics reveals inter-ethnic variation in immune response to Falciparum malaria.*  
**The American Journal of Human Genetics**, 112, 709–723.

---

## 👤 Author

Replication & analysis by **Opeoluwa Shodipe**  
Interests: Malaria genomics, systems immunology, molecular surveillance

*Open to collaboration.*
