# Sponge_Epigenetics
Regulatory system of sponge Amphimedon queenslandica. 

## <div align="center">Description of the data</div>

### - ATAC-seq data:
The **dataset** comprises a comprehensive DNA accessibility profile of the Demosponge species _Amphimedon queenslandica_ through eight developmental stages, obtained using ATAC-seq. The dataset includes **40,218 consensus peaks**, with **4,751 of them showing differential accessibility** in at least one developmental stage.

### - ChiP-seq data:
Deep sequencing **(100 bp paired-end)** of the **adult libraries** – **two biological replicates** for **H3K4me3, H3K4me1, H3K36me3, H3K27me3, RNAPII, input DNA** and **no biological replicates** for **H3K27ac and total histone H3** – was performed by the Macrogen Oceania NGS Unit on **Illumina HiSeq 2000 instrument** (Illumina, San Diego, CA, United States). 

Deep sequencing **(40 bp paired-end)** of the **larva libraries** – **no biological replicates** for **H3K4me3, H3K4me1, H3K27me3, H3K27ac, RNAPII, input DNA** – was performed by the Central Analytical Research facility (CARF), Brisbane, Queensland, Australia, on **Illumina NextSeq 500 instrument** (Illumina, San Diego, CA, United States).

Paired-end 76-bp reads for **larva and adult Amphimedon queenslandica (three biological replicates for each stage)** were obtained using a **150 cycle NextSeq 500/550 high output reagent kit v2**. Raw data isn’t available, but **preprocessed data contain ATAC-seq peaks**.
### - RNA-seq data:
Raw RNA-seq [data](https://www.ncbi.nlm.nih.gov/sra/?term=SRP044247) for **adult** and **larva** _Amphimedon queenslandica_ (Illumina HiSeq 2000, PolyA selected, **paired-end**). Additionally, complementary RNA-seq [data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54364) of **_Amphimedon queenslandica_ larvae** can be used for **gene expression analysis** (Illumina HiSeq 2000, **single-end**). 

## <div align="center">Data Processing</div>
### 1) ATAC-seq Raw data processing:
There are currently 12 fastq.gz files related to larval and adult stages available on this [link](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10203/sdrf). We managed to process only 8 from each of them. We were able to do this with the following command in nextflow, and I also pointed out the components of the nextflow.contig (basically txt file which you can create by yourself for limiting cluster) and atac_seq_raw_files.csv (which was made by the instruction from [nextflow for atac-seq](https://nf-co.re/atacseq)) files for a better awareness of what's going on in the nextflow code.
**NextFlow command:**
```bash
nextflow run nf-core/atacseq -r 2.0 \
-name nf_atac_seq_raw_processing \
--input atac_seq_raw_files.csv -profile singularity\
--fasta Sponge_genome/Amphimedon_queenslandica.Aqu1.dna.toplevel.fa \
--gtf Sponge_genome/Amphimedon_queenslandica.Aqu1.56.gtf \
--outdir /home/Daniil.Sobolev/Project_directory/nf_processed_data \
--macs_gsize 880000000 -work-dir /home/Daniil.Sobolev/Project_directory/nf_intermediate/\
--email sobolevdaniil**********@gmail.com --mito_name null --aligner bowtie2
```
**atac_seq_raw_files.csv:**
```
sample,fastq_1,fastq_2,replicate
larvae1,ATAC_SEQ_RAW/Larvae1_raw/Larvae1_1.fastq.gz,ATAC_SEQ_RAW/Larvae1_raw/Larvae1_2.fastq.gz,1
larvae2,ATAC_SEQ_RAW/Larvae2_raw/Larvae2_1.fastq.gz,ATAC_SEQ_RAW/Larvae2_raw/Larvae2_2.fastq.gz,1
adult61,ATAC_SEQ_RAW/Adult61_raw/Adult61_1.fastq.gz,ATAC_SEQ_RAW/Adult61_raw/Adult61_2.fastq.gz,1
adult65,ATAC_SEQ_RAW/Adult65_raw/Adult65_1.fastq.gz,ATAC_SEQ_RAW/Adult65_raw/Adult65_2.fastq.gz,1
```
**nextflow.config:**
```
executor {
    name = 'local'
    cpus = 12
    memory = '72 GB'
}
```
**Alternative version of the command (without nextflow) requires a step-by-step script with additional file samples_atac.txt:**
```bash
#/bin/bash

while read SAMPLE; do
        echo "Running sample ${SAMPLE}"

        #running bowtie2
bowtie2 -p 4 -x genome2/amph_queen -1 reads/${SAMPLE}_1.fastq.gz -2 reads/${SAMPLE}_2.fastq.gz | samtools view -@ 4 -Sb > ${SAMPLE}.bam

        #sorting
samtools sort -o ${SAMPLE}.sorted.bam ${SAMPLE}.bam

        #mark duplicates
java -Xmx8g -jar /home/artux/.conda-envs/omics_proj/share/picard-2.18.29-0/picard.jar MarkDuplicates I=${SAMPLE}.sorted.bam O=${SAMPLE}.markdup.bam M=${SAMPLE}_dup_metrics.txt REMOVE_DUPLICATES=false



        #calculate mapping statistics
samtools stats ${SAMPLE}.markdup.bam > stats/${SAMPLE}.stats
samtools flagstat ${SAMPLE}.markdup.bam > stats/${SAMPLE}.flagstat

        #filtering of mapped reads and stat checking
samtools view -f 2 -F 1024 -q 30 -b ${SAMPLE}.markdup.bam > ${SAMPLE}.filt.bam
samtools flagstat ${SAMPLE}.filt.bam > stats/${SAMPLE}.filt.flagstat

echo "----------------------${SAMPLE}-is-done---------------------------"

done < samples_atac.txt
```
**sample_atac.txt:**
```
adult1
adult2
adult3
larvae1
larvae2
larvae3
```

**These filtered files were processed with the following command to receive .narrowPeak  and .bed files:**
```bash
macs2 callpeak -t Stage_BAM/{SAMPLE}.filt.bam \
-g 165000000 -n adult3\
--outdir Stage_peaks/ \
--format BAMPE \
--keep-dup all -q 0.05
```
To assess Fraction of reads in peaks (FRiP) (_quality metric: if FRiP > 0.3 this could be considered as a qualitive data_) we used the following code:

```bash
bedtools intersect -a Stage_BAM/{SAMPLE}.filt.bam\
-b Stage_peaks/{SAMPLE}_peaks.narrowPeak\
-bed -u | wc -l |awk -v N=`samtools view \
-c -F 4 Stage_BAM/{SAMPLE}.filt.bam` '{printf"%.4f\n", $1/N}'
```
**Command output on FRiP score for different .filt.bam files:**
<center>

| Sample | Replicate 1 | Replicate 2 | Replicate 3 |
|--------|-------------|-------------|-------------|
| Larvae | 0.2517 | 0.2181 | 0.4448 |
| Adult | 0.2243 | 0.2714 | 0.3799 |

</center>

Here we can see that only third replicate for both stages has a qualitive data, but I performed further analysis anyway.
Since most of the FRiP < 0.3 I have decided to merge replica ATAC-seq files, so first we need to assess the correlation of bam files with each other.
```bash
multiBamSummary bins -b Adult_BAM/adult1.filt.bam Adult_BAM/adult2.filt.bam Adult_BAM/adult3.filt.bam  \
Larvae_BAM/larvae1.filt.bam Larvae_BAM/larvae2.filt.bam Larvae_BAM/larvae3.filt.bam \
 -o atac-seq_multiBamSummary_output.tab \
plotCorrelation -in atac-seq_multiBamSummary_output.tab --corMethod spearman \
--skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap\
--colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_atac-seq_readCounts.png
```
![heatmap_SpearmanCorr_atac-seq_readCounts_filt_bam](https://github.com/NAGIBATOR112/Sponge_Epigenetics/assets/89070070/b1833052-f9b3-471a-8988-47fd9c314a9a)
<div align="center">Fig. 1. Between-replicate correlation heatmap. Clustering dendrogram is shown on the left.</div>

The correlation did not show values >0.9 between replicates, so the creation of a BW had to be done for each of the BAM files separately:

```bash
bamCoverage -b Stage_BAM/{SAMPLE}.filt.bam -o {SAMPLE}.bw\
--binSize 10 \
--normalizeUsing RPGC\
--smoothLength 200 --extendReads 100 \
--maxFragmentLength 140 --centerReads \
--effectiveGenomeSize 165000000 -p 8
```

### 2) ChIP-seq Raw data processing:

