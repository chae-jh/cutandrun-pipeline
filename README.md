# CUTandRUN Snakemake Pipeline 1.0.0

### Authors    
**Ayesh Awad** (@ayeshawad)
<br>
**Jeongho Chae**
<br>
*These authors are listed in alphabetical order and contributed equally to this work.*
<br>
**Date:** 2025-02-27

### Quickstart

This analysis utilises a snakemake pipeline to process CUTnRUN data. Once the pipeline has been cloned to the analysis directory (*preferably in scratch space*) using the command:

```
git clone https://sc.unc.edu/dept-fureylab/gutandrun.git
```

There's no only prerequisite for running and the command for executing the GUTandRUN pipeline is:

```
sbatch Snakemake_SLURMsubmission.sbatch
```


# Available Module Versions (from project_config.yaml)
- **fastqc**
   - Current fastqc version : 0.12.1
- **cutadapt**
   - Current cutadapt version : 4.10
- **bowtie2**
   - Current bowtie2 version : 2.4.5
- **samtools**
   - Current samtools version : 1.21
- **bedtools**
   - Current bedtools version : 2.31
- **picard**
   - Current picard version : 2.26.11 
- **deeptools**
   - Current deeptools version : 3.5.4
- **macs**
   - Current macs version : 2.2.9.1
- **python**
   - Current python version : 3.9.6


# Pipeline Rules

 1. **Rule all**
  - Defines the final expected output files.

 2. **fastQC**
  - Perform quality control on FASTQ files using FastQC
  - Generates an HTML report to assess quality before trim_adapter

 3. **trim_adapter**
  - Trim sequencing adapters from merged FASTQ files using cutadapt
  - Generates cleaned FASTQ files(adapter removed)

 4. **fastQC_trim**
  - Perform quality control on trimmed FASTQ files using FastQC
  - Generates an HTML report to assess whether quality increased after the removal of adapter

 5. **fastqScreen**
  - Check for contamination in the trimmed FASTQ files using FastqScreen
  - Creates txt, html reports to assess whether the reads are contaminated or not.

 6. **bowtie2index**
  - Generate bowtie2 index files from target, spikein genome FASTA retrieved by function get_genome_fastas

 7. **align_target**
  - Align trimmed FASTQ files to the target genome using Bowtie2.

 8. **align_spikein**
  - Align trimmed FASTQ files to the spikein genome using Bowtie2.

 9. **extract_alignment_stats_target**
  - Extract statistical information from the target-aligned sequencing data.

 10. **extract_alignment_stats_spikein**
  - Extract statistical information from the spikein-aligned sequencing data.

 11. **filter_reads_target**
  - Removes unmapped reads(-F 4) and low-quality reads(-q 30)

 12. **filter_reads_spikein**
  - Removes unmapped reads(-F 4) and low-quality reads(-q 30)

 13. **mark_duplicates_target**
  - Detect and mark PCR duplicates in BAM files(target aligned) using picard.

 14. **mark_duplicates_spikein**
  - Detect and mark PCR duplicates in BAM files(spikein aligned) using picard.

 15. **remove_duplicates**
  - Remove PCR duplicate reads from BAM files(target aligned)

 16. **compute_read_counts**
  - Count the number of reads in the target and spike-in BAM files using samtools

 17. **compute_normalization_factor**
  - Calculates a normalization factor based on the target / spike-in

 18. **sort_bam**
  - Sort and index split species BAM files(coordinate based sorting)

 19. **name_sort_bam**
  - Sort and index split species BAM files(name based sorting)

 20. **bam_to_bedgraph**
  - Converts a coordinate-sorted BAM file into a bedGraph format using bedtools

 21. **bam_to_bed**
  - Converts a name-sorted BAM file into a BEDPE format using bedtools

 22. **normalize_bedgraph**
  - Applies normalization factor to the bedGraph file

 23. **sort_bedgraph**
  - Sort normalized bedGraph file by chromosome(column 1) and start position(column 2) 

 24. **bedgraph_to_bigwig**
  - Convert sorted, normalied bedGraph file into a BigWig format

 25. **call_peaks_seacr**
  - Call peaks using SEACR to identify enriched genomic regions

 26. **call_peaks_macs2**
  - Call peaks using MACS2 to identify enriched genomic regions
 
# moveOutFiles
 Moves final results files for locally run samples to permenant space.
 set moveOutFiles to TRUE in project_config.yaml after running the pipeline, checking that everything had run correctly,
 and then rerun the pipeline using the same submission statement.
