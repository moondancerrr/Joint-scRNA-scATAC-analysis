library(ArchR)
library(edgeR)
library(harmony) 
library(patchwork)
library(parallel)
library(Matrix)
library(scType)
library(Seurat)
library(Signac)
library(tidyverse)
library(ggplot2)
library(wordcloud)

args <- commandArgs(TRUE)
file.h5 <- args[1]
file.frag <- args[2]
file.fig <- args[3]

multiome <- file.h5 %>% Read10X_h5()

seurat_object <- CreateSeuratObject(counts = multiome$`Gene Expression`)

p1 <- VlnPlot(seurat_object, features = c("nCount_RNA", "nCount_ATAC", "percent.mt"), ncol = 3,
    log = TRUE, pt.size = 0) + NoLegend()

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object <- subset(
	x = seurat_object, 
	subset = nFeature_RNA > 400 &
	      	 nFeature_RNA < 5000 & 
		 percent.mt < 10
)

p2 <- VlnPlot(seurat_object, features = c("nCount_RNA", "nCount_ATAC", "percent.mt"), ncol = 3,
    log = TRUE, pt.size = 0) + NoLegend()

# Save volcano plots for before and after filtering (cells with more than 400 features and fewer than 5,000 features)
file.fig %>% pdf(width = 6, height = 9)
p1 / p2
dev.off()

s.genes <- intersect(cc.genes$s.genes, rownames(seurat_object))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(seurat_object))

# To ensure that your Seurat object is appropriately normalized and that the data slot is populated:
if (!"data" %in% Assays(seurat_object)) {
    seurat_object <- NormalizeData(seurat_object)
}

# Calculate cell cycle phase scores
seurat_object <- CellCycleScoring(seurat_object, 
                                  s.features = s.genes, 
                                  g2m.features = g2m.genes, 
                                  set.ident = TRUE)

# Check the identities (cell cycle phases)
table(Idents(seurat_object))

# Scale data to standardize the expression measurements by subtracting the average expression and dividing by the standard deviation for each feature
seurat_object <- ScaleData(seurat_object, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
seurat_object <- FindVariableFeatures(seurat_object)

# Run PCA
seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = FALSE)

# Run UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# Find neighbors and clusters
seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)  # Adjust resolution for cluster granularity

# Plot UMAP with clusters
umap_plot <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5)
pdf(file.path(file.fig, "umap_clusters.pdf"))
print(umap_plot)
dev.off()

# Regress Out Mitochondrial Reads and Cell Cycle Heterogeneity
# This step will normalize the data while regressing out unwanted sources of variation (mitochondrial content and cell cycle scores)
seurat_object <- SCTransform(seurat_object, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), verbose = FALSE)

# Run PCA
seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = FALSE)

# Run UMAP based on the PCA
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# Find neighbors and perform clustering
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)  # Adjust resolution as needed

# Identify marker genes for each cluster
markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View the top markers for each cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Load necessary libraries
lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = TRUE)

# Source the scType functions from GitHub
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Prepare gene sets for your analysis
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system")

# If using SCTransform, extract the transformed data
scRNAseqData <- GetAssayData(seurat_object, layer = "data", assay = "SCT")

# Use scType to score and assign cell types
es.max <- sctype_score(
  scRNAseqData = scRNAseqData,     # Scaled scRNA-seq data
  scaled = TRUE,                   # Whether the data is already scaled
  gs = gs_list$gs_positive,        # Positive gene sets
  gs2 = gs_list$gs_negative        # Negative gene sets (optional)
)

# View the results
#View(es.max)

# Get the highest scoring cell type for each cell
seurat_object$scType_cell_type <- apply(es.max, 2, function(col) rownames(es.max)[which.max(col)])

# Graffiti Falls-inspired palette with 14 colors
graffiti_falls_14_colors <- c(
  "#F72585", "#7209B7", "#3A0CA3", "#4361EE", "#4CC9F0",  # Pinks and Purples
  "#FF006E", "#FF6B6B", "#FFD166", "#F48C06", "#FFBA08",  # Reds, Oranges, Yellows
  "#06D6A0", "#198754", "#8338EC", "#9D4EDD"              # Greens and Purples
)

# Kangaskhan-inspired palette with 14 colors
kangaskhan_palette_14_colors <- c(
  "#A0522D",  # Sienna (Brown)
  "#CD853F",  # Peru (Tan)
  "#D2B48C",  # Tan
  "#8B4513",  # SaddleBrown
  "#F4A460",  # SandyBrown
  "#FFE4C4",  # Bisque (Light Beige)
  "#8B0000",  # DarkRed
  "#000000",  # Black (Replaces White)
  "#5F9EA0",  # CadetBlue
  "#556B2F",  # DarkOliveGreen
  "#2F4F4F",  # DarkSlateGray
  "#708090",  # SlateGray
  "#696969",  # DimGray
  "#B22222"   # FireBrick (Dark Red)
)

# Apply the Graffiti Falls-inspired and Kangaskhan-inspired palettes to UMAP
combined_palette <- c(graffiti_falls_14_colors, kangaskhan_palette_14_colors)

# Visualize the UMAP with cell type labels
umap_cell_type_plot <- DimPlot(seurat_object, group.by = "scType_cell_type", reduction = "umap", label = TRUE, pt.size = 0.5)
pdf(file.path(file.fig, "umap_cell_types.pdf"))
print(umap_cell_type_plot)
dev.off()

# Identify Differentially Expressed Genes (DEGs) Across Cell Types
# Step 1: Calculate the number of cells per cluster
cell_counts <- table(Idents(seurat_object))

# Step 2: Identify clusters with 50 or more cells
valid_clusters <- names(cell_counts[cell_counts >= 50])

# Step 3: Subset the Seurat object to include only cells from valid clusters
seurat_object <- subset(seurat_object, idents = valid_clusters)

# Identify DEGs for each cell type
deg_list <- list()
for (cell_type in unique(Idents(seurat_object))) {
  deg <- FindMarkers(seurat_object, ident.1 = cell_type,
		     ident.2 = NULL,
                     test.use = "wilcox",
                     min.pct = 0.1,
                     logfc.threshold = 0.1)

  # Filter by FDR < 0.05
  deg <- deg[deg$p_val_adj < 0.05, ]

  # Store DEGs for each cell type
  deg_list[[cell_type]] <- deg
}

# Aggregate counts by cell type (for pseudobulk)
seurat_object$celltype_group <- Idents(seurat_object)

# Sum counts for pseudobulk expression
pseudobulk_counts <- AggregateExpression(seurat_object, group.by = "celltype_group")

# Create a new metadata dataframe for the pseudobulk counts
pseudobulk_metadata <- data.frame(celltype_group = colnames(pseudobulk_counts$RNA))

# Check if the pseudobulk_metadata matches the pseudobulk counts
ncol(pseudobulk_counts$RNA) == nrow(pseudobulk_metadata)  # Should be TRUE

# Create the design matrix based on the new pseudobulk metadata
design <- model.matrix(~ celltype_group, data = pseudobulk_metadata)

# Check the dimensions
nrow(design) == ncol(dge$counts)  # Should be TRUE

# Fit the model
fit <- glmFit(dge, design)

# Perform likelihood ratio test
lrt <- glmLRT(fit, coef = 2)  # Adjust coef based on your comparison

# Extract DEGs
deg_edgeR <- topTags(lrt, n = Inf)

# Filter DEGs with log2(FC) > 0.3 and FDR < 0.05
final_deg <- deg_edgeR$table[abs(deg_edgeR$table$logFC) > 0.3 & deg_edgeR$table$FDR < 0.05, ]

# Extract the data frame from the TopTags object
deg_df <- deg_edgeR$table

# Add a column to mark significant genes (FDR < 0.05 and |logFC| > 0.3)
deg_df$significant <- with(deg_df, FDR < 0.05 & abs(logFC) > 0.3)

# Create the volcano plot
volcano_plot <- ggplot(deg_df, aes(x = logFC, y = -log10(PValue), color = significant)) +
  geom_point(alpha = 0.4, size = 1.75) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot of DEGs", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal()

pdf(file.path(file.fig, "volcano_plot_degs.pdf"))
print(volcano_plot)
dev.off()

# Extract the significant DEGs
deg_df <- deg_edgeR$table
significant_deg <- deg_df[deg_df$FDR < 0.05, ]

# Set the word sizes based on absolute log fold change
words <- rownames(significant_deg)
word_sizes <- abs(significant_deg$logFC)

# Set up a color palette
palette <- brewer.pal(8, "Dark2")

# Generate the word cloud 
pdf(file.path(file.fig, "word_cloud_degs.pdf"))
wordcloud(words = words,
          freq = word_sizes,         # Use logFC to determine word size
          min.freq = 0.5,            # Minimum frequency threshold for word inclusion
          random.order = FALSE,      # Words with higher frequency appear in the center
          rot.per = 0.35,            # Proportion of words that are rotated
          colors = palette,          # Color palette for the words
          scale = c(3, 0.5))         # Scaling range for word sizes
dev.off()

pdf(file.path(file.fig, "FeaturePlot_degs.pdf"))
# Now, plot UMAP with the DEG scores for each cluster
# This assumes you want to plot the score of a specific cluster or just an average score
# For example, let's plot the DEG score for cluster 0
FeaturePlot(seurat_object, features = "DEG_score_01", reduction = "umap") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "UMAP Plot: DEG Scores for Cluster 0")
dev.off()
