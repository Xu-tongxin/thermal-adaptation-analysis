#!/bin/bash
set -e

WORKDIR="${WORKDIR:-$PWD}"
REF="${REF:-${WORKDIR}/reference/Pinfestans_S1.fasta}"
GENE_BED=${WORKDIR}/genes_location.bed

cd $WORKDIR
mkdir -p clean_fastq alignment_bam vcf alternative_genomes cds_extracted logs
[ ! -f "${REF}.fai" ] && samtools faidx $REF

echo "=== Step 1: QC (fastp) ==="
for sample in $(cat fastq/all_samples_final.txt); do
  if [ -f "fastq/${sample}_R2.fastq.gz" ]; then
    fastp -i fastq/${sample}_R1.fastq.gz -I fastq/${sample}_R2.fastq.gz \
         -o clean_fastq/${sample}_R1.fastq.gz -O clean_fastq/${sample}_R2.fastq.gz \
         --detect_adapter_for_pe --thread 4 --html logs/${sample}_fastp.html 2> logs/${sample}_fastp.log
  else
    fastp -i fastq/${sample}_R1.fastq.gz -o clean_fastq/${sample}_R1.fastq.gz \
         --thread 4 --html logs/${sample}_fastp.html 2> logs/${sample}_fastp.log
  fi
done

echo "=== Step 2: Alignment + MarkDuplicates (12 samples parallel) ==="
do_alignment() {
  sample=$1
  REF=$2
  RG="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA"
  if [ -f "clean_fastq/${sample}_R2.fastq.gz" ]; then
    bwa mem -R "$RG" -t 4 $REF clean_fastq/${sample}_R1.fastq.gz clean_fastq/${sample}_R2.fastq.gz | \
      samtools sort -n -@ 2 -m 4G -o - | \
      samtools fixmate -m - - | \
      samtools sort -@ 2 -m 4G -o - | \
      samtools markdup -@ 2 -r - alignment_bam/${sample}.sorted.bam
  else
    bwa mem -R "$RG" -t 4 $REF clean_fastq/${sample}_R1.fastq.gz | \
      samtools sort -@ 2 -m 4G -o - | \
      samtools markdup -@ 2 -r - alignment_bam/${sample}.sorted.bam
  fi
  samtools index alignment_bam/${sample}.sorted.bam
  echo "Finished ${sample}"
}
export -f do_alignment
export REF

cat fastq/all_samples_final.txt | parallel -j 12 "do_alignment {} $REF"

echo "=== Step 3: Variant Calling ==="
ls alignment_bam/*.sorted.bam > bam_list.txt
bcftools mpileup -f $REF -b bam_list.txt --threads 48 -a FORMAT/AD,FORMAT/DP -Q 30 | \
  bcftools call -mv -f GQ -o vcf/raw_variants.vcf

echo "=== Step 4: Strict SNP Filtering ==="
bcftools view -v snps -m2 -M2 vcf/raw_variants.vcf | \
  bcftools filter -i 'QUAL>30 && FMT/DP>4' | \
  bcftools filter -e 'F_MISSING>0.2' -o vcf/filtered_snps.vcf

bgzip -c vcf/filtered_snps.vcf > vcf/filtered_snps.vcf.gz
tabix -p vcf vcf/filtered_snps.vcf.gz

