# CutNRun Snakemake Pipeline 1.0.0
# Date: 2-27-2025
# Written by:
#	Ayesh Awad
#	Jeongho Chae

configfile: "project_config.yaml"
modules=config['module']


import pandas as pd

def read_samples():
	"""
	Reads the samplesheet.tsv file and returns a dictionary that looks like this:
	{
    		'Sample1': {
        		'Mod1': {'fastq_1': 'Fastq/Sample1_Mod1_R1.fastq.gz', 'fastq_2': 'Fastq/Sample1_Mod1_R2.fastq.gz'},
        		'Mod2': {'fastq_1': 'Fastq/Sample1_Mod2_R1.fastq.gz', 'fastq_2': 'Fastq/Sample1_Mod2_R2.fastq.gz'}
    		},
    		'Sample2': {
        		'Mod1': {'fastq_1': 'Fastq/Sample2_Mod1_R1.fastq.gz', 'fastq_2': 'Fastq/Sample2_Mod1_R2.fastq.gz'}
    		}
	}
	"""
	f = pd.read_csv(config['sampleInfo'], sep="\t")

	samp_dict = {}
    	for line, row in f.iterrow():
        	sample_id = row['sample']
		modification = row['modification']

		if sample_id not in samp_dict = {}

        	samp_dict[sample_id][modification] = {
            		'fastq_1': row['fastq_1'],
            		'fastq_2': row['fastq_2']
       	}
    	return samp_dict

SAMPLES = read_samples()


###
# BEGIN PIPELINE
###

rule all:
    input:
        expand(),
        expand(),
        expand()

rule fastQC:
    input:
        r1= lambda wildcards: SAMPLES[wildcards.sample][wildcards.modification]['fastq_1'],
        r2= lambda wildcards: SAMPLES[wildcards.sample][wildcards.modification]['fastq_2']
    output:
        r1_qc="FastQC/{sample}_{modification}_R1_fastqc.html",
        r2_qc="FastQC/{sample}_{modification}_R2_fastqc.html"
    envmodules:
        modules['fastqcVer']
    shell:
        """
        fastqc -o ./FastQC/ -f fastq {input.r1} {input.r2}
        """

rule trim_adapter:
    input:
        r1= lambda wildcards: SAMPLES[wildcards.sample][wildcards.modification]['fastq_1'],
        r2= lambda wildcards: SAMPLES[wildcards.sample][wildcards.modification]['fastq_2']
    output:
        r1="Fastq/{sample}_{modification}_R1_trim.fastq.gz",
        r2="Fastq/{sample}_{modification}_R2_trim.fastq.gz"
    log:
        trim_stats="Logs/{sample}_{modification}_cutadapt.log"
    envmodules:
        modules['cutadaptVer']
    shell:
        """
        cutadapt \
            -a ADAPTER_SEQ_R1 -A ADAPTER_SEQ_R2 \
            -m 20 -q 20 \
            -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.trim_stats}
        """

rule fastQC_trim:
    input:
        r1="Fastq/{sample}_{modification}_R1_trim.fastq.gz",
        r2="Fastq/{sample}_{modification}_R2_trim.fastq.gz"
    output:
        r1_qc="FastQC/{sample}_{modification}_R1_trim_fastqc.html",
        r2_qc="FastQC/{sample}_{modification}_R2_trim_fastqc.html"
    envmodules:
        modules['fastqcVer']
    shell:
        """
        fastqc -o ./FastQC/ -f fastq {input.r1} {input.r2}
        """

rule fastqScreen:
    input:
        r1="Fastq/{sample}_{modification}_R1_trim.fastq.gz",
        r2="Fastq/{sample}_{modification}_R2_trim.fastq.gz"
    output:
        r1_txt="FQscreen/{sample}_{modification}_R1_trim_screen.txt",
        r1_html="FQscreen/{sample}_{modification}_R1_trim_screen.html",
        r2_txt="FQscreen/{sample}_{modification}_R2_trim_screen.txt",
        r2_html="FQscreen/{sample}_{modification}_R2_trim_screen.html"
    params:
        fqscreenPath=modules['fqscreenPath'],
        fqscreenConf=modules['fqscreenConf']
    threads: 4
    shell:
        """
        {params.fqscreenPath} --threads {threads} --force --aligner bowtie2 -conf {params.fqscreenConf} {input.r1} --outdir ./FQscreen/
        {params.fqscreenPath} --threads {threads} --force --aligner bowtie2 -conf {params.fqscreenConf} {input.r2} --outdir ./FQscreen/
        """

def get_genome_fastas(config):
	"""
	Reads species information from the config file and returns a dictionary that looks like:
	{
        'human': '/path/human.fasta',
        'e_coli': '/path/e_coli.fasta',
        'spikein': '/path/spikein.fasta'
    	}
	"""
	fastas = {
		species_name: species['fasta']
		for species_name, species in config['genome'].items()
	}
	return fastas

rule bowtie2index:
    input:
        target_fasta=lambda wildcards: genome_fastas[wildcards.species],  # Target genome FASTA
        spikein_fasta=lambda wildcards: genome_fastas['spikein'] if 'spikein' in genome_fastas else None  # Spike-in genome FASTA (if available)
    output:
        target_combined_fa="Bowtie2Index/{species}.fa",
        target_index=expand("Bowtie2Index/{species}.{num}.bt2", species="{species}", num=["1", "2", "3", "4", "rev.1", "rev.2"]),
        spikein_combined_fa="Bowtie2Index/spikein.fa",
        spikein_index=expand("Bowtie2Index/spikein.{num}.bt2", num=["1", "2", "3", "4", "rev.1", "rev.2"]) if 'spikein' in genome_fastas else None
    params:
        module=modules['bowtie2Ver']
    shell:
        """
        # this is for the target
        module load {params.module}
        bowtie2-build {input.target_fasta} {output.target_combined_fa}

        # this is for the spike-in
        {{"bowtie2-build {input.spikein_fasta} {output.spikein_combined_fa}" if input.spikein_fasta else ""}}
        """


rule align_target:
    input:
        index=expand("Bowtie2Index/{species}.{num}.bt2", species="{species}", num=["1", "2", "3", "4", "rev.1", "rev.2"]),
        r1="Fastq/{sample}_{modification}_R1_trim.fastq.gz",
        r2="Fastq/{sample}_{modification}_R2_trim.fastq.gz"
    output:
        sam="Sam/{sample}_{modification}_{species}_trim.sam",
        logInfo="Logs/{sample}_{modification}_{species}_bowtie2.txt"
    threads: 8
    params:
        refgenome="Bowtie2Index/{species}"
    envmodules:
        modules['bowtie2Ver']
    shell:
        """
        bowtie2 --seed 123 -p {threads} -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -x {params.refgenome} -1 {input.r1} -2 {input.r2} -S {output.sam} 2> {output.logInfo}
        """

rule align_spikein:
    input:
        index=expand("Bowtie2Index/spikein.{num}.bt2", num=["1", "2", "3", "4", "rev.1", "rev.2"]),
        r1="Fastq/{sample}_{modification}_R1_trim.fastq.gz",
        r2="Fastq/{sample}_{modification}_R2_trim.fastq.gz"
    output:
        sam="Sam/{sample}_{modification}_spikein_trim.sam",
        logInfo="Logs/{sample}_{modification}_spikein_bowtie2.txt"
    threads: 8
    params:
        refgenome="Bowtie2Index/spikein"
    envmodules:
        modules['bowtie2Ver']
    shell:
        """
        bowtie2 --seed 123 -p {threads} -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -x {params.refgenome} -1 {input.r1} -2 {input.r2} -S {output.sam} 2> {output.logInfo}
        """

rule extract_alignment_stats_target:
    input:
        sam="Sam/{sample}_{modification}_{species}_trim.sam"
    output:
        stats="Logs/{sample}_{modification}_{species}_alignment_stats.txt"
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        samtools flagstat {input.sam} > {output.stats}
        """

rule extract_alignment_stats_spikein:
    input:
        sam="Sam/{sample}_{modification}_spikein_trim.sam"
    output:
        stats="Logs/{sample}_{modification}_spikein_alignment_stats.txt"
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        samtools flagstat {input.sam} > {output.stats}
        """
rule filter_reads_target:
    input:
        sam="Sam/{sample}_{modification}_{species}_trim.sam"
    output:
        bam="Bam/{sample}_{modification}_{species}_filtered.bam"
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        samtools view -b -F 4 -q 30 {input.sam} | samtools sort -o {output.bam}
        """

rule filter_reads_spikein:
    input:
        sam="Sam/{sample}_{modification}_spikein_trim.sam"
    output:
        bam="Bam/{sample}_{modification}_spikein_filtered.bam"
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        samtools view -b -F 4 -q 30 {input.sam} | samtools sort -o {output.bam}
        """
rule mark_duplicates_target:
    input:
        bam="Bam/{sample}_{modification}_{species}_filtered.bam"
    output:
        marked_bam="Bam/{sample}_{modification}_{species}_marked.bam",
        metrics="Logs/{sample}_{modification}_{species}_dup_metrics.txt"
    envmodules:
        modules['picardVer']
    shell:
        """
        picard MarkDuplicates \
            I={input.bam} \
            O={output.marked_bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT
        """

rule mark_duplicates_spikein:
    input:
        bam="Bam/{sample}_{modification}_spikein_filtered.bam"
    output:
        marked_bam="Bam/{sample}_{modification}_spikein_marked.bam",
        metrics="Logs/{sample}_{modification}_spikein_dup_metrics.txt"
    envmodules:
        modules['picardVer']
    shell:
        """
        picard MarkDuplicates \
            I={input.bam} \
            O={output.marked_bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT
        """

rule remove_duplicates_igg:
    input:
        marked_bam="Bam/{sample}_IgG_marked.bam"
    output:
        dedup_bam="Bam/{sample}_IgG_dedup.bam",
        dedup_bai="Bam/{sample}_IgG_dedup.bam.bai"
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        samtools view -b -F 1024 {input.marked_bam} | samtools sort -o {output.dedup_bam}
        samtools index {output.dedup_bam}
        """

rule keep_duplicates_histone:
    input:
        marked_bam="Bam/{sample}_{modification}_marked.bam"
    output:
        dedup_bam="Bam/{sample}_{modification}_dedup.bam",
        dedup_bai="Bam/{sample}_{modification}_dedup.bam.bai"
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        cp {input.marked_bam} {output.dedup_bam}
        samtools index {output.dedup_bam}
        """

rule compute_read_counts:
    input:
        target_bam="Bam/{sample}_{modification}_dedup.bam",
        spikein_bam="Bam/{sample}_{modification}_spikein_marked.bam"
    output:
        target_counts="Logs/{sample}_{modification}_target_counts.txt",
        spikein_counts="Logs/{sample}_{modification}_spikein_counts.txt"
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        samtools view -c {input.target_bam} > {output.target_counts}
        samtools view -c {input.spikein_bam} > {output.spikein_counts}
        """

rule compute_normalization_factor:
    input:
        target_counts="Logs/{sample}_{modification}_target_counts.txt",
        spikein_counts="Logs/{sample}_{modification}_spikein_counts.txt"
    output:
        norm_factor="Logs/{sample}_{modification}_norm_factor.txt"
    run:
        target = int(open(input.target_counts).read().strip())
        spikein = int(open(input.spikein_counts).read().strip())

        # what if there is low spike-in counts
        if spikein == 0:
            factor = 1  # this means no scaling
            print(f"WARNING: Spike-in read count is 0 for {input.spikein_counts}. Norm factor is set to 1.")
        elif spikein < 100:  # this if the spike has low counts
            factor = target / spikein
            print(f"WARNING: Very low spike-in count ({spikein}). Norm factor may not be great.")
        else:
            factor = target / spikein 

        # here we can always check the norm factor
        with open(output.norm_factor, "w") as f:
            f.write(str(factor) + "\n")

        print(f"Computed norm factor: {factor} for {input.target_counts} using spike-in {input.spikein_counts}")

rule sort_bam:
    input:
        bam="Bam/{sample}_{modification}_dedup.bam"
    output:
        sorted_bam="Bam/{sample}_{modification}_sorted.bam",
        sorted_bai="Bam/{sample}_{modification}_sorted.bam.bai"
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        samtools sort -o {output.sorted_bam} {input.bam}
        samtools index {output.sorted_bam}
        """

rule bam_to_bedgraph:
    input:
        bam="Bam/{sample}_{modification}_sorted.bam"
    output:
        bedgraph="Bedgraph/{sample}_{modification}.bedgraph"
    envmodules:
        modules['bedtoolsVer']
    shell:
        """
        bedtools genomecov -ibam {input.bam} -bg > {output.bedgraph}
        """

rule normalize_bedgraph:
    input:
        bedgraph="Bedgraph/{sample}_{modification}.bedgraph",
        norm_factor="Logs/{sample}_{modification}_norm_factor.txt"
    output:
        normalized_bedgraph="Bedgraph/{sample}_{modification}_normalized.bedgraph"
    run:
        norm_factor = float(open(input.norm_factor).read().strip())
        with open(input.bedgraph, "r") as bed, open(output.normalized_bedgraph, "w") as out:
            for line in bed:
                chrom, start, end, value = line.strip().split()
                normalized_value = float(value) * norm_factor
                out.write(f"{chrom}\t{start}\t{end}\t{normalized_value:.6f}\n")

rule sort_bedgraph:
    input:
        normalized_bedgraph="Bedgraph/{sample}_{modification}_normalized.bedgraph"
    output:
        sorted_bedgraph="Bedgraph/{sample}_{modification}_normalized_sorted.bedgraph"
    envmodules:
        modules['bedtoolsVer']
    shell:
        """
        sort -k1,1 -k2,2n {input.normalized_bedgraph} > {output.sorted_bedgraph}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph="Bedgraph/{sample}_{modification}_normalized_sorted.bedgraph",
        chrom_sizes="reference/chrom_sizes.txt"  # HERE WE NEED TO TELL IT THE CHROMOSOME SIZE
    output:
        bigwig="BigWig/{sample}_{modification}.bw"
    envmodules:
        modules['bedtoolsVer']
    shell:
        """
        bedGraphToBigWig {input.bedgraph} {input.chrom_sizes} {output.bigwig}
        """

rule call_peaks_seacr:
    input:
        bedgraph="Bedgraph/{sample}_{modification}_normalized_sorted.bedgraph",
        control_bedgraph="Bedgraph/{sample}_IgG_normalized_sorted.bedgraph"  # here i am just adding the IgG for us to compare along the modifications
    output:
        peaks="Peaks/SEACR/{sample}_{modification}_seacr.bed"
    envmodules:
        modules['seacrVer']
    shell:
        """
        SEACR_1.3.sh {input.bedgraph} {input.control_bedgraph} stringent {output.peaks}
        """

rule call_peaks_macs2:
    input:
        bam="Bam/{sample}_{modification}_sorted.bam",
        control_bam="Bam/{sample}_IgG_sorted.bam"  # again here to compare along the modifications
    output:
        narrowPeaks="Peaks/MACS2/{sample}_{modification}_peaks.narrowPeak",
        broadPeaks="Peaks/MACS2/{sample}_{modification}_peaks.broadPeak"
    envmodules:
        modules['macs2Ver']
    shell:
        """
        macs2 callpeak -t {input.bam} -c {input.control_bam} \
        -f BAMPE -g hs --keep-dup all --outdir Peaks/MACS2/ \
        -n {wildcards.sample}_{wildcards.modification} --broad
        """

