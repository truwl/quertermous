# Quertermous

## Data
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113348

## WGS
- Whole-genome sequencing data were processed with the GATK best practices pipeline
- Trimmed Illumina adapters using Cutadapt v1.93
- Trimmed reads were aligned to the hg19 reference genome with Burrows-Wheeler Aligner (BWA) v0.7.124 using its bwa-mem module5 with parameters "-M -t 8 -R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tLB:lib1\tPU:unit1'".
- Duplicate reads in alignment result were marked with Picard v1.92. We performed indel
realignment and base recalibration with GATK v3.46
- The GATK HaplotypeCaller was used to generate gVCF files, which were fed into GenotypeGVCFs for joint genotype calling. We recalibrated variants using the GATK VariantRecalibrator module. Since subsequent eQTL calling (see Section 12) with RASQUAL7 required phased variants, we phased our call set with Beagle v4.18
- We first used the Beagle conform-gt module to correct any reference genotypes if
they are different from hg19. We then phased and imputed against 1000 Genomes project
phase 3 version 5a9
- Variants with imputation allelic r2 less than 0.8 and Hardy-Weinberg Equilibrium p-value less than 1x10-6 were filtered out. 

## RNA-Seq
### Prep
Demultiplexing was performed with bcl2fastq script from Illumina.
- Base quality control of demultiplexed sequences was done using FastQC v0.11.4 quality control tool
- Fastq files that correspond to the unique sample were merged using a custom script.
- Mapping of the reads was performed with STAR v2.4.0i11. In accordance to the GATK Best Practices for RNA-Seq we
used the STAR 2-pass alignment pipeline

### Alignment
- First, reads contained in the raw fastq files were mapped to GRCh37/hg19 human genome using STAR and during the first alignment pass splice junctions were discovered with high stringency
- Second pass mapping with STAR was then performed using a new index that was created with splice junction information contained in the file SJ.out.tab from the first pass STAR mapping
- Splice junctions from the first pass were used as annotation in a second pass to permit lower stringency alignment, and therefore higher sensitivity
- Prior to gene expression quantification, we used WASP12 to filter out reads that are prone to mapping bias. Read counts and RPKM were calculated with RNA-SeQC v1.1.813 using default parameters with additional flags “-n 1000 -noDoC -strictMode” using GENCODE v19 annotation14. Allele-specific read counts were generated with the createASVCF module in RASQUAL7

### Splicing
We quantified intron excision levels using LeafCutter15. In brief, we converted bam files to splice junction files using the bam2junc.sh script, and defined intron clusters using leafcutter_cluster.py with default parameters. This requires at least 30 reads supporting each cluster and at least 0.1% of reads supporting each intron within the cluster, and allows intron to have a maximum size of 100kb.

## ATAC-Seq
We used the ENCODE ATAC-seq pipeline to perform alignment and peak calling (https://github.com/kundajelab/atac_dnase_pipelines)
Steps included in the pipeline
- FASTQ files were trimmed with Cutadapt v1.93 and aligned with Bowtie2 v2.2.618 with default parameters. Duplicate reads were marked with Picard v1.126.
- The alignment were converted to ENCODE tagAlign format. Records were shifted +4 and -5 for positive-strand and minus-strand reads.
- MACS2 v2.0.819 was used to call peaks with default parameters. Each alignment was split into two pseudoreplicate (subsample of reads) and peaks were called independently.
- Irreproducible Discovery Rate (IDR) 20 analyses were performed based on pseudo-replicates with a cutoff of 0.1 to output an IDR call set, which was used for downstream analysis.
- We used WASP12 to filter out readsthat are prone to mapping bias.
