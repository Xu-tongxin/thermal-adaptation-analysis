# Thermal adaptation SNP analysis and heatmap

Scripts for genome-wide SNP processing and heatmap visualization of *Phytophthora infestans* isolates.

## Files

- `run_pipeline_parallel.sh`: read quality control, alignment, duplicate removal, variant calling and SNP filtering.
- `Heatmap_Plotting.R`: heatmap visualization of a pairwise SNP-distance matrix.
- This README includes the matrix-generation commands.

## Requirements

Linux/Bash, fastp, BWA, SAMtools, GNU Parallel, bcftools, bgzip and tabix; Python with [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip); snp-dists v0.8.2; R v4.3.0 with pheatmap. The R package grDevices is included with R.

The pipeline uses 12 parallel alignment jobs and up to 48 threads in variant calling. Run it with sufficient allocated computing resources.

## Input data

- Study resequencing reads: NCBI BioProject **PRJNA1511611**.
- Reference genome: Genome Warehouse accession **GWHJXAH00000000.1**, saved as `Pinfestans_S1.fasta`.
- Accession numbers for the 50 published isolates are provided in the manuscript's Supporting Information.

Prepare `fastq/all_samples_final.txt` with one sample identifier per line. Place reads in `fastq/` as `<sample>_R1.fastq.gz` and, for paired-end samples, `<sample>_R2.fastq.gz`. Download inputs separately; sequencing data and generated matrices are not bundled here.

## 1. SNP processing

From the analysis working directory, set the input paths and build the BWA reference index before running the pipeline:

```bash
export WORKDIR="$PWD"
export REF="$WORKDIR/reference/Pinfestans_S1.fasta"
bwa index "$REF"
bash run_pipeline_parallel.sh
```

The pipeline produces `vcf/filtered_snps.vcf.gz`. Its original analysis commands and parameters are retained; the original server-specific working-directory and reference paths have been replaced with environment-variable defaults.

## 2. Pairwise SNP-distance matrix

Place `vcf2phylip.py` in the analysis directory. The recorded conversion and distance commands are:

```bash
python vcf2phylip.py -i vcf/filtered_snps.vcf.gz -f
snp-dists -c filtered_snps.min4.fasta > snp_matrix_all.csv
```

The FASTA filename above is the output filename recorded for this analysis; confirm the output name produced by your installed vcf2phylip version.

## 3. Heatmap

The plotting script expects a numeric, square pairwise SNP-distance CSV for the 70 isolates shown in Figure 2. Sample identifiers must be present in the first row and first column, with matching row and column order. Save the plotting input as `snp_matrix_no_outgroup.csv` in the R working directory.

Install pheatmap once:

```r
install.packages("pheatmap")
```

Run the plot in R:

```r
source("Heatmap_Plotting.R")
```

Rows and columns are clustered using Euclidean distance and complete linkage. The colour gradient ranges from dark red through a pale colour to dark blue.

To save a PDF explicitly:

```r
pdf("SNP_heatmap.pdf", width = 12, height = 12)
source("Heatmap_Plotting.R")
dev.off()
```
