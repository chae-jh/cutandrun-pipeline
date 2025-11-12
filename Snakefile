# Cut&Run Snakemake Pipeline 1.0.0
# Date: 2025-02-27
# Written by:
#       Ayesh Awad (ayesh_awad@med.unc.edu)
#       Jeongho Chae (chaejh@unc.edu)
#
# * These authors contributed equally to this work.



configfile: "project_config.yaml"
modules = config['module']

import pandas as pd
import os

def read_samples():
    f = pd.read_csv(config['sampleInfo'], sep="\t")
    samp_dict = {}
    for _, row in f.iterrows():
        sample_id = row['sample']
        modification = row['modification']
        if sample_id not in samp_dict:
            samp_dict[sample_id] = {}
        samp_dict[sample_id][modification] = {
            'fastq_1': row['fastq_1'],
            'fastq_2': row['fastq_2']
        }
    # IgG appears at the top
    for sample_id in samp_dict:
        samp_dict[sample_id] = {
            modification: samp_dict[sample_id][modification]
            for modification in sorted(samp_dict[sample_id], key=lambda x: x != 'IgG')
        }
    return samp_dict

# this function checks modifications and stops if it is IgG
def mods_excluding_IgG():
    return [mod for samp in SAMPLES.values() for mod in samp.keys() if mod != "IgG"]

SAMPLES = read_samples()
print(SAMPLES)

#Moves final results files for locally run samples to permenant space.
if config["moveOutFiles"]:
        projectDir = config["projectDir"]
        os.makedirs(projectDir, exist_ok=True)
        for sample in os.listdir(f"temp/"):
                sample_path = f"temp/{sample}"
                if not os.path.isdir(sample_path):
                        continue
                for modification in os.listdir(sample_path):
                        mod_path = f"{sample_path}/{modification}"
                        if not os.path.isdir(mod_path):
                                continue
                        dest_path = f"{projectDir}/{modification}/{sample}/"
                        os.makedirs(dest_path, exist_ok=True)
                        cmds = [
                            f"cp -r {mod_path}/16_sort_bam/. {dest_path}/Final_bam/",
                            f"cp -r {mod_path}/17_normalize_bam_RPGC/. {dest_path}/Bigwig/",
                            f"cp -r {mod_path}/18_normalize_bam_spikein/. {dest_path}/Bigwig/",
                            f"cp -r {mod_path}/19_call_peaks_macs2/. {dest_path}/Macs2/",
                            f"cp -r {mod_path}/1_fastQC/. {dest_path}/fastQC/",
                            f"cp -r {mod_path}/3_fastQC_trim/. {dest_path}/fastQC_trim/",
                            f"cp -r {mod_path}/4_fastqScreen/. {dest_path}/fastqScreen/"
                        ]
                        for cmd in cmds:
                                os.system(cmd)
        os._exit(0)

# create temp directory to keep intermediate files
os.makedirs("temp", exist_ok=True)

###
# BEGIN PIPELINE
###


rule all:
    input:
        expand("temp/{sample}/{modification}/1_fastQC/{sample}_{modification}_R1_fastqc.html", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]), 
        expand("temp/{sample}/{modification}/1_fastQC/{sample}_{modification}_R2_fastqc.html", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]), 
        expand("temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R1_trim.fastq.gz", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]), 
        expand("temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R2_trim.fastq.gz", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]), 
        expand("temp/{sample}/{modification}/3_fastQC_trim/{sample}_{modification}_R1_trim_fastqc.html", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]), 
        expand("temp/{sample}/{modification}/3_fastQC_trim/{sample}_{modification}_R2_trim_fastqc.html", 
               sample=SAMPLES.keys(),
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/4_fastqScreen/{sample}_{modification}_R1_trim_screen.txt", 
               sample=SAMPLES.keys(),
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/4_fastqScreen/{sample}_{modification}_R1_trim_screen.html", 
               sample=SAMPLES.keys(),
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/4_fastqScreen/{sample}_{modification}_R2_trim_screen.txt", 
               sample=SAMPLES.keys(),
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/4_fastqScreen/{sample}_{modification}_R2_trim_screen.html", 
               sample=SAMPLES.keys(),
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/5_count_reads_spikein/{sample}_{modification}_spikein_reads.txt", 
               sample=SAMPLES.keys(),
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/6_align_to_genome/{sample}_{modification}_{species}_trim.sam",
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()],
               species=["hg38", "Ecoli"]),
	    expand("temp/{sample}/{modification}/7_extract_alignment_stats/{sample}_{modification}_{species}_alignment_stats.txt",
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()],
               species=["hg38", "Ecoli"]),
        expand("temp/{sample}/{modification}/8_filter_reads/{sample}_{modification}_{species}_filtered.bam", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()],
               species=["hg38", "Ecoli"]), 
        expand("temp/{sample}/{modification}/9_count_reads_Ecoli/{sample}_{modification}_Ecoli_reads.txt", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]), 
        expand("temp/{sample}/{modification}/12_mark_duplicates/{sample}_{modification}_marked.bam", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]), 
        expand("temp/{sample}/{modification}/13_remove_duplicates/{sample}_{modification}_dedup.bam", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/16_sort_bam/{sample}_{modification}_sorted.bam", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/16_sort_bam/{sample}_{modification}_sorted.bam.bai", 
               sample=SAMPLES.keys(), 
               modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        expand("temp/{sample}/{modification}/15_compute_normalization_factor/{sample}_{modification}_norm_factor.txt", 
               sample=SAMPLES.keys(), 
               modification=mods_excluding_IgG()), 
        expand("temp/{sample}/{modification}/17_normalize_bam_RPGC/{sample}_{modification}_RPGC.bw", 
               sample=SAMPLES.keys(), 
               modification=mods_excluding_IgG()),
        expand("temp/{sample}/{modification}/18_normalize_bam_spikein/{sample}_{modification}_spikein_norm.bw", 
               sample=SAMPLES.keys(), 
               modification=mods_excluding_IgG()),
        expand("temp/{sample}/{modification}/19_call_peaks_macs2/{sample}_{modification}_peaks.broadPeak", 
               sample=SAMPLES.keys(), 
               modification=mods_excluding_IgG())
        #expand("temp/{sample}/{modification}/20_bam_to_bed/{sample}_{modification}_sorted.bed", 
        #       sample=SAMPLES.keys(), 
        #       modification=[mod for samp in SAMPLES.values() for mod in samp.keys()]),
        #expand("temp/{sample}/{modification}/21_call_peaks_sicer2", 
        #       sample=SAMPLES.keys(), 
        #       modification=mods_excluding_IgG())

#Perform quality control on FASTQ files using FastQC
#Generates an HTML report to assess quality before trim_adapter
rule fastQC:
    input:
        r1= lambda wildcards: SAMPLES[wildcards.sample][wildcards.modification]['fastq_1'],
        r2= lambda wildcards: SAMPLES[wildcards.sample][wildcards.modification]['fastq_2']
    output:
        r1_qc="temp/{sample}/{modification}/1_fastQC/{sample}_{modification}_R1_fastqc.html",
        r2_qc="temp/{sample}/{modification}/1_fastQC/{sample}_{modification}_R2_fastqc.html"
    params:
        path="temp/{sample}/{modification}/1_fastQC/"
    shell:
        """
	module load {modules[fastqcVer]}

        fastqc -o {params.path} -f fastq {input.r1} {input.r2}
        """

#Trim sequencing adapters from merged FASTQ files using cutadapt
#Generates cleaned FASTQ files(adapter removed)
rule trim_adapter:
    input:
        r1= lambda wildcards: SAMPLES[wildcards.sample][wildcards.modification]['fastq_1'],
        r2= lambda wildcards: SAMPLES[wildcards.sample][wildcards.modification]['fastq_2']
    output:
        r1="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R1_trim.fastq.gz",
        r2="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R2_trim.fastq.gz"
    params:
        a=modules["cutadapt"]["a"],
        A=modules["cutadapt"]["A"],
        m=modules["cutadapt"]["minimumLength"],
        q=modules["cutadapt"]["qualityCutoff"]
    log:
        trim_stats="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_cutadapt.log"
    shell:
        """
	    module load {modules[cutadapt][cutadaptVer]}

        cutadapt \
            -a {params.a} -A {params.A} \
            -m {params.m} -q {params.q} \
            -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.trim_stats}
        """

#Perform quality control on trimmed FASTQ files using FastQC
#Generates an HTML report to assess whether quality increased after the removal of adapter
rule fastQC_trim:
    input:
        r1="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R1_trim.fastq.gz",
        r2="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R2_trim.fastq.gz"
    output:
        r1_qc="temp/{sample}/{modification}/3_fastQC_trim/{sample}_{modification}_R1_trim_fastqc.html",
        r2_qc="temp/{sample}/{modification}/3_fastQC_trim/{sample}_{modification}_R2_trim_fastqc.html"
    params:
        outdir="temp/{sample}/{modification}/3_fastQC_trim/"
    shell:
        """
	module load {modules[fastqcVer]}

        fastqc -o {params.outdir} -f fastq {input.r1} {input.r2}
        """

#Check for contamination in the trimmed FASTQ files using FastqScreen
#Creates txt, html reports to assess whether the reads are contaminated or not.
rule fastqScreen:
    input:
        r1="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R1_trim.fastq.gz",
        r2="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R2_trim.fastq.gz"
    output:
        r1_txt="temp/{sample}/{modification}/4_fastqScreen/{sample}_{modification}_R1_trim_screen.txt",
        r1_html="temp/{sample}/{modification}/4_fastqScreen/{sample}_{modification}_R1_trim_screen.html",
        r2_txt="temp/{sample}/{modification}/4_fastqScreen/{sample}_{modification}_R2_trim_screen.txt",
        r2_html="temp/{sample}/{modification}/4_fastqScreen/{sample}_{modification}_R2_trim_screen.html"
    params:
        fqscreenPath=modules['fqscreen']['fqscreenPath'],
        fqscreenConf=modules['fqscreen']['fqscreenConf'],
        outdir="temp/{sample}/{modification}/4_fastqScreen/"
    threads: 4
    shell:
        """
        {params.fqscreenPath} --threads {threads} --force --aligner bowtie2 -conf {params.fqscreenConf} {input.r1} {input.r2} --outdir {params.outdir}
        """


rule count_reads_spikein:
    input:
        r1="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R1_trim.fastq.gz",
        r2="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R2_trim.fastq.gz"
    output:
        spikein_counts="temp/{sample}/{modification}/5_count_reads_spikein/{sample}_{modification}_spikein_reads.txt"
    params:
        barcodes_r1=[# Methylation panel
                    "TTCGCGCGTAACGACGTACCGT", "CGCGATACGACCGCGTTACGCG",
                    "CGACGTTAACGCGTTTCGTACG", "CGCGACTATCGCGCGTAACGCG",
                    "CCGTACGTCGTGTCGAACGACG", "CGATACGCGTTGGTACGCGTAA",
                    "TAGTTCGCGACACCGTTCGTCG", "TCGACGCGTAAACGGTACGTCG",
                    "TTATCGCGTCGCGACGGACGTA", "CGATCGTACGATAGCGTACCGA",
                    "CGCATATCGCGTCGTACGACCG", "ACGTTCGACCGCGGTCGTACGA",
                    "ACGATTCGACGATCGTCGACGA", "CGATAGTCGCGTCGCACGATCG",
                    "CGCCGATTACGTGTCGCGCGTA", "ATCGTACCGCGCGTATCGGTCG",
                    # Acetylation panel
                    "CGACGGATACGTACGTCCGTCG", "CGACTAACGCGTCGACGCGTAG", "ATATCGCGCGTACGACGCGATA",
                    "AGCGTCGTACGAATCGCGCGAT", "CGATTGACGCGTCGCGCGACAT", "CGCGACGCAATCCGTCGAACGC",
                    "ACGCGCGATATAACGCGTATCG", "ACGACCGTTCGACGATTCGACG", "TCGTCGATTCGATCGCGCGATT",
                    "CGCGATATGCGGCGTTCGACGG", "AATCGCGCGATTATTATGCGCG", "ACCGATACGCGTATCGCGCGAA",
                    "TTACGTCGCGACGATACGCGTT", "TTAACGTCGCGTCGGTACGCGA", "ACGACGCGATAATATCGCGCGT",
                    "TTCGCGCGATACGACGTAACGG", "TATCGCGTCGTACGATCGTCGG", "TCGCGACGAAAATATCGCGCGG",
                    "ATCGCGCGATTCGCGTAATACG", "TCGTACGACCGTTAACGCGACG", "CGTTACGCGTCCGATACGCGTT",
                    "ATCGTCGACGATACGTTACGCG", "TCGCGTATACGCTAACGCGACG"],
        barcodes_r2=[
                    # Methylation panel
                    "CGTTCGAACGTTCGTCGACGAT", "TCGCGATTACGATGTCGCGCGA",
                    "ACGCGAATCGTCGACGCGTATA", "CGCGATATCACTCGACGCGATA",
                    "CGCGAAATTCGTATACGCGTCG", "CGCGATCGGTATCGGTACGCGC",
                    "GTGATATCGCGTTAACGTCGCG", "TATCGCGCGAAACGACCGTTCG",
                    "CCGCGCGTAATGCGCGACGTTA", "CCGCGATACGACTCGTTCGTCG",
                    "GTCGCGAACTATCGTCGATTCG", "CCGCGCGTATAGTCCGAGCGTA",
                    "CGATACGCCGATCGATCGTCGG", "CCGCGCGATAAGACGCGTAACG",
                    "CGATTCGACGGTCGCGACCGTA", "TTTCGACGCGTCGATTCGGCGA",
                    # Acetylation panel
                    "CGACGTAACGCTATGCGTCGCG", "TATGCGTCGCGTACGTTACGCG", "ATATCGCGCGGATTACGTCGCG",
                    "ACGGTACGTCGAATCGCGCGAT", "TATCGCGTCGAACGATTCGCGT", "TATACGCGTCGTACGTTACGCG",
                    "ACGATCGTCGGACGATTCGACG", "ACCGTTCGTCGATCGCGCGATT", "TCGTCGACGATATTACGCGCGA",
                    "CGATTCGGCGATCGTCGTATCG", "ACGAACGTCGAACGCGCGATAT", "TTCGAGCGACGCGTTAACGTCG",
                    "CGTCGTTCGACCGTAATCGCGA", "TGTCGCGCGTATCGTACGATCG", "ACGATTCGCGTCGAATCGACGA",
                    "TTCGTACGCGACGACGTATCGG", "TCGACGTTCGTAGCGTCGTACG", "TCGCGTACGAACCGAATACGCG",
                    "GATAGTCGACGAACGCGTATCG", "TCGGTACGCTATCGGCGTATCG", "ATCGTCGACGATACGTTACGCG",
                    "TCGCGTAATCGACGTTCGAACG"]
    shell:
        """

        # here we unzip the FASTQ files for processing
        gunzip -c {input.r1} > {input.r1}.tmp
        gunzip -c {input.r2} > {input.r2}.tmp

        # here we count R1 barcodes
        total_r1=0
        for barcode in {params.barcodes_r1}; do
            count=$(grep -c "$barcode" {input.r1}.tmp)
            total_r1=$((total_r1 + count))
        done

        # here we count R2 barcodes
        total_r2=0
        for barcode in {params.barcodes_r2}; do
            count=$(grep -c "$barcode" {input.r2}.tmp)
            total_r2=$((total_r2 + count))
        done

        # here we add them together
        total_spikein_reads=$((total_r1 + total_r2))
        echo $total_spikein_reads > {output.spikein_counts}

        # here i am removing the temp files
        rm {input.r1}.tmp {input.r2}.tmp
        """
 
#Align trimmed FASTQ files to the target genome using Bowtie2.
rule align_to_genome:
    input:
        r1="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R1_trim.fastq.gz",
        r2="temp/{sample}/{modification}/2_trim_adapter/{sample}_{modification}_R2_trim.fastq.gz"
    output:
        sam="temp/{sample}/{modification}/6_align_to_genome/{sample}_{modification}_{species}_trim.sam",
        logInfo="temp/{sample}/{modification}/6_align_to_genome/{sample}_{modification}_{species}_bowtie2.txt"
    threads: 8
    params:
        refgenome=lambda wildcards: f"{modules['bowtie2']['bowtie2Index'][wildcards.species]}/{ 'hg38' if wildcards.species == 'hg38' else 'Ecoli' }"
    shell:
        """
	module load {modules[bowtie2][bowtie2Ver]}

        bowtie2 --seed 123 -p {threads} -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -x {params.refgenome} -1 {input.r1} -2 {input.r2} -S {output.sam} 2> {output.logInfo}
        """



#Extract statistical information from the aligned sequencing data.
rule extract_alignment_stats:
    input:
        sam="temp/{sample}/{modification}/6_align_to_genome/{sample}_{modification}_{species}_trim.sam"
    output:
        stats="temp/{sample}/{modification}/7_extract_alignment_stats/{sample}_{modification}_{species}_alignment_stats.txt"
    shell:
        """
	module load {modules[samtoolsVer]}

        samtools flagstat {input.sam} > {output.stats}
        """



#Removes unmapped reads(-F 4) and low-quality reads(-q 30)
rule filter_reads:
    input:
        sam="temp/{sample}/{modification}/6_align_to_genome/{sample}_{modification}_{species}_trim.sam"
    output:
        bam="temp/{sample}/{modification}/8_filter_reads/{sample}_{modification}_{species}_filtered.bam"
    shell:
        """
	module load {modules[samtoolsVer]}

        samtools view -b -F 4 -q 10 {input.sam} | samtools sort -T temp -o {output.bam}
        """

#Count the number of reads that mapped to Ecoli
rule count_reads_Ecoli:
    input:
        bam="temp/{sample}/{modification}/8_filter_reads/{sample}_{modification}_Ecoli_filtered.bam"
    output:
        Ecoli_counts="temp/{sample}/{modification}/9_count_reads_Ecoli/{sample}_{modification}_Ecoli_reads.txt"
    shell:
        """
        module load {modules[samtoolsVer]}

        samtools view -c {input.bam} > {output.Ecoli_counts}
        """

rule remove_blacklist_regions:
    input:
        bam="temp/{sample}/{modification}/8_filter_reads/{sample}_{modification}_hg38_filtered.bam",
        blacklist=config["blacklist"]
    output:
        filtered_bam="temp/{sample}/{modification}/10_remove_blacklist_regions/{sample}_{modification}_blacklist_filtered.bam"
    shell:
        """
        module load {modules[bedtoolsVer]}

        bedtools intersect -v -abam {input.bam} -b {input.blacklist} > {output.filtered_bam}
        """

rule sort_remove_blacklist_regions:
    input:
        filtered_bam="temp/{sample}/{modification}/10_remove_blacklist_regions/{sample}_{modification}_blacklist_filtered.bam"
    output:
        sort_filtered_bam="temp/{sample}/{modification}/11_sort_remove_blacklist_regions/{sample}_{modification}_sort_blacklist_filtered.bam"
    shell:
        """
        module load {modules[samtoolsVer]}

        samtools sort -T temp -o {output.sort_filtered_bam} {input.filtered_bam}
        """


#Detect and mark PCR duplicates in BAM files(target aligned) using picard.
rule mark_duplicates:
    input:
        bam="temp/{sample}/{modification}/11_sort_remove_blacklist_regions/{sample}_{modification}_sort_blacklist_filtered.bam"
    output:
        marked_bam="temp/{sample}/{modification}/12_mark_duplicates/{sample}_{modification}_marked.bam",
        metrics="temp/{sample}/{modification}/12_mark_duplicates/{sample}_{modification}_dup_metrics.txt"
    params:
        tmpDir="temp/mark_duplicates_tmp/{sample}/{modification}"
    shell:
        """
	module load {modules[picard][picardVer]}

        mkdir -p {params.tmpDir}
        export TMPDIR={params.tmpDir}
        java -Xmx16G -jar /nas/longleaf/apps/picard/2.26.11/picard-2.26.11/picard.jar MarkDuplicates \
            I={input.bam} \
            O={output.marked_bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false \
            ASSUME_SORTED=false \
            VALIDATION_STRINGENCY=LENIENT

        rm -rf {params.tmpDir}
        """


#Remove PCR duplicate reads from BAM files(target aligned)
rule remove_duplicates:
    input:
        marked_bam="temp/{sample}/{modification}/12_mark_duplicates/{sample}_{modification}_marked.bam"
    output:
        dedup_bam="temp/{sample}/{modification}/13_remove_duplicates/{sample}_{modification}_dedup.bam"
    shell:
        """
	module load {modules[samtoolsVer]}

        samtools view -b -F 1024 {input.marked_bam} > {output.dedup_bam}
	"""

#Count the number of reads in the target and spike-in BAM files using samtools
rule count_reads_hg38:
    input:
        dedup_bam="temp/{sample}/{modification}/13_remove_duplicates/{sample}_{modification}_dedup.bam"
    output:
        hg38_counts="temp/{sample}/{modification}/14_count_reads_hg38/{sample}_{modification}_hg38_counts.txt"
    shell:
        """
        module load {modules[samtoolsVer]}

        samtools view -c {input.dedup_bam} > {output.hg38_counts}
        """

rule compute_normalization_factor:
    input:
        Ecoli_counts="temp/{sample}/{modification}/9_count_reads_Ecoli/{sample}_{modification}_Ecoli_reads.txt",
        hg38_counts="temp/{sample}/{modification}/14_count_reads_hg38/{sample}_{modification}_hg38_counts.txt"
    output:
        norm_factor="temp/{sample}/{modification}/15_compute_normalization_factor/{sample}_{modification}_norm_factor.txt"
    run:
        hg38_counts=int(open(input.hg38_counts).read().strip())
        Ecoli_counts=int(open(input.Ecoli_counts).read().strip())
        if Ecoli_counts == 0:
            scale_factor = 1
            print(f"WARNING: Spike-in read count is 0 for {input.spikein_counts}. Norm factor is set to 1.")
        else:
            scale_factor = hg38_counts / (Ecoli_counts * 100)

        with open(output.norm_factor, "w") as out_file:
            out_file.write(f"{scale_factor:.6f}\n")

        print(f"Computed norm factor: {scale_factor} for normalization using Ecoli")


#Sort and index split species BAM files(coordinate based sorting)
rule sort_bam:
    input:
        bam="temp/{sample}/{modification}/13_remove_duplicates/{sample}_{modification}_dedup.bam"
    output:
        sorted_bam="temp/{sample}/{modification}/16_sort_bam/{sample}_{modification}_sorted.bam",
        sorted_bai="temp/{sample}/{modification}/16_sort_bam/{sample}_{modification}_sorted.bam.bai"
    shell:
        """
	module load {modules[samtoolsVer]}

        samtools sort -T temp -o {output.sorted_bam} {input.bam}
        samtools index {output.sorted_bam}
        """

# Normalize BAM files directly to RPGC BigWig using deepTools bamCoverage
rule normalize_bam_RPGC:
    input:
        bam="temp/{sample}/{modification}/16_sort_bam/{sample}_{modification}_sorted.bam"
    output:
        bigwig="temp/{sample}/{modification}/17_normalize_bam_RPGC/{sample}_{modification}_RPGC.bw"
    params:
        effective_genome_size=2913022398  # hg38 effective genome size
    shell:
        """
        module load {modules[deeptoolsVer]}

        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.bigwig} \
            --binSize 20 \
            --smoothLength 100 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize {params.effective_genome_size} \
            --ignoreDuplicates \
            --minMappingQuality 10
        """

# Normalize BAM files with the E coli spike ins using deepTools bamCoverage
rule normalize_bam_spikein:
    input:
        bam="temp/{sample}/{modification}/16_sort_bam/{sample}_{modification}_sorted.bam",
        scale_factor="temp/{sample}/{modification}/15_compute_normalization_factor/{sample}_{modification}_norm_factor.txt"
    output:
        bigwig="temp/{sample}/{modification}/18_normalize_bam_spikein/{sample}_{modification}_spikein_norm.bw"
    params:
        effective_genome_size=2913022398  # hg38 effective genome size
    shell:
        """
        module load {modules[deeptoolsVer]}

        scale=$(cat {input.scale_factor})

        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.bigwig} \
            --binSize 20 \
            --smoothLength 100 \
            --scaleFactor $scale \
            --normalizeUsing None \
            --effectiveGenomeSize {params.effective_genome_size} \
            --ignoreDuplicates \
            --minMappingQuality 10
        """

#Call peaks using MACS2 to identify enriched genomic regions
rule call_peaks_macs2:
    input:
        bam="temp/{sample}/{modification}/16_sort_bam/{sample}_{modification}_sorted.bam",
        control_bam="temp/{sample}/IgG/16_sort_bam/{sample}_IgG_sorted.bam"
    output:
        broadPeaks="temp/{sample}/{modification}/19_call_peaks_macs2/{sample}_{modification}_peaks.broadPeak"
    params:
        output_dir="temp/{sample}/{modification}/19_call_peaks_macs2/"
    shell:
        """
	module load {modules[macsVer]}

        macs2 callpeak -t {input.bam} -c {input.control_bam} \
        -f BAMPE -g hs --keep-dup all --outdir {params.output_dir} \
        -n {wildcards.sample}_{wildcards.modification} --broad
        """

#rule bam_to_bed:
#    input:
#        bam="temp/{sample}/{modification}/16_sort_bam/{sample}_{modification}_sorted.bam"
#    output:
#        bed="temp/{sample}/{modification}/20_bam_to_bed/{sample}_{modification}_sorted.bed"
#    shell:
#        """
#        module load {modules[bedtoolsVer]}
#
#        bamToBed -i {input.bam} > {output.bed}
#        """

#rule call_peaks_sicer2:
#    input:
#        bed="temp/{sample}/{modification}/20_bam_to_bed/{sample}_{modification}_sorted.bed",
#        control_bed="temp/{sample}/IgG/20_bam_to_bed/{sample}_IgG_sorted.bed"
#    output:
#        directory("temp/{sample}/{modification}/21_call_peaks_sicer2")
#    params:
#        sicer2_script=modules["sicer2Path"],
#        output_dir="temp/{sample}/{modification}/21_call_peaks_sicer2",
#        species="hg38",
#        window_size=200,
#        gap_size=600,
#        frag_size=150,
#        redundancy=1,
#        egf=0.74,
#        fdr=0.01,
#        cpu=8
#    shell:
#        """
#        module load {modules[pythonVer]}
#        module load {modules[bedtoolsVer]}
#
#        PYTHONPATH=/proj/fureylab/ayesh/SICER2 \
#        python3 {params.sicer2_script} \
#            --treatment_file {input.bed} \
#            --control_file {input.control_bed} \
#            --species {params.species} \
#            --window_size {params.window_size} \
#            --gap_size {params.gap_size} \
#            --fragment_size {params.frag_size} \
#            --redundancy_threshold {params.redundancy} \
#            --effective_genome_fraction {params.egf} \
#            --false_discovery_rate {params.fdr} \
#            --cpu {params.cpu} \
#            --output_directory {params.output_dir}
#        """


