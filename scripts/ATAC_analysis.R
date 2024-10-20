library(ArchR)
library(ggplot2)
library(SummarizedExperiment)
library(Matrix)
library(Seurat)

# Command-line arguments for input files and output directory
args <- commandArgs(TRUE)
file.frag <- args[1]
output_dir <- args[2]

addArchRGenome("hg38")  # Set the genome reference
addArchRThreads(threads = 1)
arrowFiles <- createArrowFiles(
    inputFiles = file.frag,
    sampleNames = "pbmc_granulocyte_sorted_10k"
)

project <- ArchRProject(
  ArrowFiles = arrowFiles,
  outputDirectory = output_dir,
  copyArrows = TRUE
)

# Add Doublet Scores to the project
project <- addDoubletScores(
    input = project,
    k = 10,                   # Number of nearest neighbors (10 is a reasonable default)
    knnMethod = "UMAP",        # The dimensionality reduction method used
    LSIMethod = 1              # Latent semantic indexing method (default is 1)
)
# Apply quality filters
project <- filterDoublets(project)  # Filter out doublets using default settings
project <- subsetArchRProject(
  ArchRProj = project,
  cells = which(
    project$TSSEnrichment > 12 &
    project$nFrags > 3000 &
    project$nFrags < 30000 &
    project$nucleosomeRatio < 2
  )
)


# Extract fragment counts per cell for Seurat object creation
cell_data <- getCellColData(project, select = c("nFrags", "TSSEnrichment"))
cell_barcodes <- rownames(cell_data)  # Store cell barcodes

# Assuming nFrags is analogous to "counts" here, use it to create the Seurat object
seurat_obj <- CreateSeuratObject(counts = as.matrix(cell_data$nFrags), project = "ATAC", meta.data = cell_data)

# Add mitochondrial content filtering
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Filter cells based on mitochondrial content
seurat_obj <- subset(seurat_obj, subset = percent.mt < 10)

# Extract the filtered cell barcodes from Seurat
filtered_barcodes <- colnames(seurat_obj)

# Use these barcodes to subset the ArchR project
project <- subsetArchRProject(
  ArchRProj = project,
  cells = filtered_barcodes
)

# Perform LSI (Linear Dimensionality Reduction)
project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  varFeatures = 20000,
  dimsToUse = 1:30
)

project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",  # Use the LSI-reduced dimensions
  name = "Clusters",             # Store the clustering result in the "Clusters" column
  resolution = 0.8               # Set an appropriate resolution for clustering
)

# Compute coverage for each group (in this case, "Clusters")
project <- addGroupCoverages(
  ArchRProj = project,
  groupBy = "Clusters"
)

# Call peaks with MACS2
project <- addReproduciblePeakSet(
  ArchRProj = project,
  groupBy = "Clusters",
  peakMethod = "MACS2",
  genomeSize = 3e9  # Set based on human genome size
)

# Compute the PeakMatrix, which contains peak accessibility information for each cell
project <- addPeakMatrix(ArchRProj = project, force = TRUE)

# Identify differentially accessible sites (DAS)
das <- getMarkerFeatures(
  ArchRProj = project,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("log10(nFrags)"),  # Correct for bias
  normBy = "ReadsInPeaks",
  testMethod = "wilcoxon",
  maxCells = 15000
)

# Extract the full matrices
FDR_matrix <- assays(das)$FDR
Log2FC_matrix <- assays(das)$Log2FC
MeanDiff_matrix <- assays(das)$MeanDiff

# Apply filtering based on at least one column meeting the criteria
# Filter DAS by FDR < 0.05, log2(FC) > 0.1, and pct > 0.1
valid_rows <- rowSums(FDR_matrix < 0.05 & abs(Log2FC_matrix) > 0.1 & MeanDiff_matrix > 0.1) > 0


# Filter the SummarizedExperiment object
das_filtered <- das[valid_rows, ]

umap_plot <- plotEmbedding(
  ArchRProj = project, 
  colorBy = "cellColData", 
  name = "Clusters", 
  embedding = "UMAP"
)

# Save the plot to a file
output_umap <- file.path(output_dir, "UMAP_ATACembedding_plot.png")
ggsave("UMAP_ATACembedding_plot.png", plot = umap_plot, width = 8, height = 6)

