library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(Signac)
library(biovizBase)
library(patchwork)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(tidyverse)

args <- commandArgs(TRUE)
file.h5 <- args[1]
file.frag <- args[2]
output_dir = arg[3]

# Read in multiome data (both ATAC and RNA)
#multiome <- Read10X_h5("path/to/multiome_data.h5")
multiome <- file.h5 %>% Read10X_h5()

# Create a Seurat object for RNA
seurat_object <- CreateSeuratObject(counts = multiome$`Gene Expression`)

# Add chromatin accessibility data (ATAC) to the Seurat object
annotations <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
genome(annotations) <- "hg38"  # Ensure that annotations are in the right genome

valid_peaks <- grepl("chr[0-9XY]+:\\d+-\\d+", rownames(multiome$Peaks))
Peaks <- multiome$Peaks[valid_peaks, ]


# Remove the `genome` parameter in `CreateChromatinAssay`
chrom_assay <- CreateChromatinAssay(
    counts = Peaks,
    sep = c(":", "-"),
    #fragments = "data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
    fragments = file.frag,
    min.cells = 10,
    min.features = 200
)

# Add chromatin accessibility data (ATAC) to the Seurat object
seurat_object[["ATAC"]] <- chrom_assay

# RNA processing (like before)
seurat_object <- NormalizeData(seurat_object)
seurat_object <- SCTransform(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object <- ScaleData(seurat_object, vars.to.regress = c("percent.mt"))
seurat_object <- RunPCA(seurat_object, npcs = 30)

# ATAC processing
seurat_object <- RunTFIDF(seurat_object, assay = "ATAC")  # Normalize the ATAC data
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 10, assay = "ATAC")
seurat_object <- RunSVD(seurat_object, assay = "ATAC")  # Perform dimensionality reduction for ATAC

# Find common cells between RNA and ATAC
common_cells <- intersect(rna_cells, atac_cells)

# Subset the Seurat object to only include common cells
seurat_object <- subset(seurat_object, cells = common_cells)

lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = TRUE)

# Source the scType functions from GitHub
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Prepare gene sets for cell type labelling
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system")

scRNAseqData <- GetAssayData(seurat_object, layer = "data", assay = "RNA")

es.max <- sctype_score(
  scRNAseqData = scRNAseqData,     # Scaled scRNA-seq data
  scaled = TRUE,                   # Whether the data is already scaled
  gs = gs_list$gs_positive,        # Positive gene sets
  gs2 = gs_list$gs_negative        # Negative gene sets (optional)
)

seurat_object$scType_cell_type <- apply(es.max, 2, function(col) rownames(es.max)[which.max(col)])

# Joint reduction using Weighted Nearest Neighbors (WNN)
seurat_object <- FindMultiModalNeighbors(
  seurat_object,
  reduction.list = list("pca", "lsi"),  # RNA and ATAC reductions
  dims.list = list(1:30, 1:30),         # Number of dimensions for each modality
  modality.weight.name = c("RNA.weight", "ATAC.weight")  # Modality weights for RNA and ATAC
)

# Get the number of cells per cluster
cell_counts <- table(Idents(seurat_object))

# Find clusters with at least 50 cells
valid_clusters <- names(cell_counts[cell_counts >= 50])

# Subset the Seurat object to include only those clusters
seurat_object <- subset(seurat_object, idents = valid_clusters)

# UMAP embedding based on WNN
seurat_object <- RunUMAP(seurat_object, nn.name = "weighted.nn")

# Visualize the integrated UMAP
pdf(file.path(output_dir, "umap_cell_types.pdf"), width = 8, height = 6)
DimPlot(seurat_object, group.by = "scType_cell_type")
dev.off()

# Peak-to-Gene Linkage Analysis
# Load the annotations again
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# Ensure that the genome version is correct
genome(annotations) <- "hg38"
Annotation(seurat_object[["ATAC"]]) <- annotations

# Manually prepend "chr" to chromosome names
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))

# Get the standard chromosomes (main chromosomes only)
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)

# Filter annotations to keep only the main chromosomes
annotations <- keepSeqlevels(annotations, main.chroms, pruning.mode = "coarse")

# Exclude chrM from the main chromosomes list
main.chroms <- setdiff(standardChromosomes(BSgenome.Hsapiens.UCSC.hg38), "chrM")
annotations <- keepSeqlevels(annotations, main.chroms, pruning.mode = "coarse")

# Calculates the GC content or other statistics for each chromatin accessibility peak
seurat_object <- RegionStats(
  object = seurat_object,
  assay = "ATAC",                   # The ATAC assay containing peak data
  genome = BSgenome.Hsapiens.UCSC.hg38  # Reference genome for human (or replace with your organism's genome)
)

# Create a link between chromatin accessibility peaks and gene expression
seurat_object <- LinkPeaks(
  object = seurat_object,
  peak.assay = "ATAC",
  expression.assay = "RNA",  # Use the RNA assay for gene expression
  gene.coords = annotations,  # Provide the gene coordinates explicitly
  genes.use = VariableFeatures(seurat_object)  # Focus on variable genes
)

# Compare RNA Expression and ATAC Accessibility

# Set the ATAC assay as the default
DefaultAssay(seurat_object) <- "ATAC"

# Plot peak-to-gene correlations
p1 <- CoveragePlot(
  object = seurat_object,
  region = "chr5-150192424-150371319",    # Plot for a specific gene
  features = "ADAM19",  # The gene whose peaks are being correlated with
  expression.assay = "RNA",       # Gene expression data
  extend.upstream = 50000,        # Extend upstream of the gene for peak visualization
  extend.downstream = 50000
)

p2 <- CoveragePlot(
  object = seurat_object,
  region = "chr1-173559805-173565987",  # Coordinates for TNFRSF18 in hg38
  features = "TNFRSF18", 
  expression.assay = "RNA", 
  extend.upstream = 50000, 
  extend.downstream = 50000
)

# Set cell type identities in the Seurat object
Idents(seurat_object) <- "scType_cell_type"

# Differential accessibility analysis
# Identify peaks that are differentially accessible in Memory B cells compared to other cell types
da_peaks <- FindMarkers(
  object = seurat_object,
  ident.1 = "Memory B cells",
  assay = "ATAC",
  slot = "data"
)

# Visualize differential peaks with gene expression correlations
# Set the default assay to RNA
DefaultAssay(seurat_object) <- "RNA"

p3 <- FeaturePlot(seurat_object, 
	    features = c("ADAM19", "chr5-150192424-150371319"),
            cols = c("blue", "red")
)

# Combine the plots using patchwork
combined_plot <- p1 / p2 / p3

# Save the combined plot to a file
pdf(file.path(output_dir, "combined_coverage_feature_plots.pdf"), width = 12, height = 10)
print(combined_plot)
dev.off()
