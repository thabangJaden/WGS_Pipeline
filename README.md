# WGS_Pipeline
This repository contains a complete WGS data analysis workflow used in my MSc research.

# Whole Genome Sequencing (WGS) Variant Analysis Pipeline

## Workflow
1. **Quality Control:** `FastQC`, `MultiQC`
2. **Trimming:** `Trimmomatic`
3. **Quality Control:** `FastQC`, `MultiQC`
4. **Alignment:** `BWA-MEM`
5. **BAM Processing:** `SAMtools`, `Picard`
6. **Variant Calling:** `GATK HaplotypeCaller`
7. **Filtering:** `bcftools`
8. **Structural variant detection:** `Lumpy` 
9. **Annotation:** `VEP`, `ANNOVAR`
10. **Visualization:** IGV, R-based plots

## Example Commands
```bash
#Quality Control
fastqc *fq.gz -o FastQC_out
multiqc FastQC_out/*

#Trim if need be
trimmomatic PE -threads 24 -phred33 \
  *_1.fq.gz _2.fq.gz \
  *_R1_paired.fq.gz *_R1_unpaired.fq.gz \
  *_R2_paired.fq.gz *_R2_unpaired.fq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 HEADCROP:10 2>> trimmomatic.log


#Quality Control
fastqc -o FastQC_out_trimmed -t 4 *.fq.gz
multiqc FastQC_out_trimmed/*


# Alignment
bwa mem -M -t 24 -R "@RG\tID:PID225\tLB:lib1\tPL:ILLUMINA\tPU:unit1\tSM:PID225" references/GRCh38.primary_assembly.genome.fa \
  *_R1_paired.fq.gz *_R2_paired.fq.gz > alignment/filename_aligned.sam

# Sorting & Indexing
# Convert
samtools view -@ 8 -bS path/to/file_aligned.sam > path/to/file_aligned.bam

# Sort
samtools sort -@ 8 -o path/to/file_sorted.bam /file_rg.bam

#index
samtools index file_sorted.bam

#gatk markduplicates
gatk MarkDuplicates --REMOVE_SEQUENCING_DUPLICATES true -I file_sorted.bam -O file_dedup.bam -M file_dedup_metrics.txt
samtools index file_dedup.bam

#GATK BaseRecalibrator & ApplyBQSR
gatk BaseRecalibrator -I file_dedup.bam -R references/GRCh38.primary_assembly.genome.fa \
  --known-sites references/BQRS_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --known-sites references/BQRS_files/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O file_recal_data.table

gatk ApplyBQSR -R references/GRCh38.primary_assembly.genome.fa -I file_dedup.bam --bqsr-recal-file file_recal_data.table -O file_recal.bam

samtools index file_recal.bam


# Variant Calling
gatk HaplotypeCaller -R references/GRCh38.primary_assembly.genome.fa -I file_recal.bam --native-pair-hmm-threads 8 -O path/to/file_raw.vcf

#Variant hard-filtering **separates SNPs and INDELS
# Extract SNPs
gatk SelectVariants -V file_raw.vcf -R references/GRCh38.primary_assembly.genome.fa -select-type SNP -O file_SNPs_raw.vcf

# Filter SNPs
gatk VariantFiltration -R references/GRCh38.primary_assembly.genome.fa -V file_SNPs_raw.vcf -O file_SNPs_final.vcf \
  --filter "QD < 2.0" --filter-name "QD2" \
  --filter "QUAL < 30.0" --filter-name "QUAL30" \
  --filter "SOR > 3.0" --filter-name "SOR3" \
  --filter "FS > 60.0" --filter-name "FS60" \
  --filter "MQ < 40.0" --filter-name "MQ40" \
  --filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

#Merge filtered SNPs & INDELs
picard MergeVcfs I=file_SNPs_final.vcf I=file_INDELs_final.vcf O=file_merged.vcf

#Index and compress the final VCF
bgzip file_merged.vcf
tabix -p vcf file_merged.vcf.gz

#Lumpy SV detection
vep -i sample. lumpy.vcf -o sample.lumpy_annotated.vcf \
--vcf --cache --offline --assembly GRCh38 --fork 4 --dir_cache /path/to/vep_cache

# Annotation
vep -i file_merged.vcf.gz -o filename_annotated.vcf --vcf --cache --offline --assembly GRCh38 --everything --fork 4 --dir_cache /path/to/vep_cache

##visualization in IGV, select region of interest if necessary for your research
