# Targets
all: run_analysis

# Download the filtered feature matrix (H5 file) from 10x Genomics or generate it from fastqs
#pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5:
#	cd data	
#	wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5

# Download and extract the fastqs for PBMC granulocyte sorted 10k dataset
pbmc_granulocyte_sorted_10k_fastqs:
	cd data
	wget https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_fastqs.tar
	tar -xvzf pbmc_granulocyte_sorted_10k_fastqs.tar


#Download Gencode human genome
GRCh38.p13.genome.fa.gz:
	 wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz

gencode.v38.annotation.gtf.gz:
	wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz

# Generate the reference genome using Cell Ranger ARC based on the config.yaml
GRCh38_reference: config.yaml
	gunzip GRCh38.p13.genome.fa.gz
	gunzip gencode.v38.annotation.gtf.gz
	cellranger-arc mkref --config=$^ --memgb=64 --nthreads=16

# Perform the Cell Ranger ARC count using the reference genome and library CSV
GRCh38:
	cellranger-arc count \
    	--id=sample_output \
    	--reference=GRCh38_reference \
    	--libraries=pbmc_granulocyte_sorted_10k_library.csv \
    	--localcores=8 \
    	--localmem=64

# Index the ATAC fragments file using tabix
data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi: data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
	tabix -p bed $<

# Generate figures from the R script for analysis and save them as PDFs
figures/%.pdf: data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5 data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
	@mkdir -p figures
	Rscript scripts/gex_analysis.R $^ figures/$*.pdf

run_analysis: data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz figures/
	@echo "Running ATAC analysis..."
	Rscript scripts/ATAC_analysis.R $^
	@echo "Running joint analysis..."
	Rscript scripts/joint_analysis.R $^
