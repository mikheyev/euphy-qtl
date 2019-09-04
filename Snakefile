# Euphy QTL
from snakemake.utils import R
from scripts.split_fasta_regions import split_fasta

outDir = "/work/MikheyevU/sasha/euphy-qtl/data"
refDir = "/work/MikheyevU/sasha/euphy-qtl/ref"
SCRATCH  = "/work/scratch/sasha"
refBowtie2 = refDir + "/euphy"
ref = refDir + "/PGA_assembly.fasta"
RADsamples = set([i.split()[0] for i in open("data/radseq/reads/bcodes.txt").readlines()])

SPLITS = range(200)
REGIONS = split_fasta(ref, len(SPLITS))  # dictionary with regions to be called, with keys in SPLITS
for region in REGIONS:
	for idx,i in enumerate(REGIONS[region]):
		REGIONS[region][idx] = " -r " + str(i)
RNAseq, = glob_wildcards(outDir + "/rnaseq/reads/{sample}.R1.fastq.gz")
Poolseq, = glob_wildcards(outDir + "/poolseq/reads/{sample}_1.fastq.gz")

localrules: all, mergeVCF, stringtieMerge

rule all:
	input: expand("data/poolseq/bowtie2/{sample}.bam", sample = Poolseq)

rule trimRNAseq:
	input:
		read1 = outDir + "/rnaseq/reads/{sample}.R1.fastq.gz",
		read2 = outDir + "/rnaseq/reads/{sample}.R2.fastq.gz"
	output:
		read1 = outDir + "/rnaseq/readsTrimmed/{sample}.R1.fastq.gz",
		read2 = outDir + "/rnaseq/readsTrimmed/{sample}.R2.fastq.gz"
	resources: mem=20, time = 60*8
	shell:
		"""
		cutadapt -q 10 --trim-n  -m 25 -f fastq \
		-a PolyA=AATTGCAGTGGTATCAACGCAGAGCGGCCGC \
		-a TS=AAGCAGTGGTATCAACGCAGAGTACATGGGG \
		-g PolyArc=GCGGCCGCTCTGCGTTGATACCACTGCAATT \
		-g TSrc=CCCCATGTACTCTGCGTTGATACCACTGCTT \
		-g nextera=CTGTCTCTTATACACATCT \
		-A PolyA2=AATTGCAGTGGTATCAACGCAGAGCGGCCGC \
		-A TS2=AAGCAGTGGTATCAACGCAGAGTACATGGGG \
		-G PolyArc2=GCGGCCGCTCTGCGTTGATACCACTGCAATT \
		-G TSrc2=CCCCATGTACTCTGCGTTGATACCACTGCTT \
		-G nextera2=CTGTCTCTTATACACATCT \
		-o {output.read1} -p {output.read2} {input}
		"""

rule hisat:
	input:
		read1 = rules.trimRNAseq.output.read1,
		read2 = rules.trimRNAseq.output.read2
	output: "data/rnaseq/hisat/{sample}.bam"
	threads: 4
	resources: mem=10, time = 60*24
	shell:
		"""
		module load hisat/2.1.0		
		hisat2 -p {threads} -x {refBowtie2} -1 {input.read1} -2 {input.read2} | samtools view -Su - | samtools sort - -m 10G -T {SCRATCH}/hisat2/{wildcards.sample} -o {output} && sleep 5m
		samtools index {output}
		"""

rule stringtie:
	input: rules.hisat.output
	output: "data/rnaseq/stringtie/{sample}.gtf"
	threads: 4
	resources: mem=30, time = 60*24
	shell:
		"""
		module load stringtie/2.0.1
		stringtie -o {output} -p {threads} {input} 
		"""

rule stringtieMerge:
	input: expand("data/rnaseq/stringtie/{sample}.gtf", sample = RNAseq)
	output: "data/rnaseq/stringtie/merged.gtf"
	threads: 4
	resources: mem=10, time = 60*24
	shell:
		"""
		module load stringtie/2.0.1
		stringtie -o {output} -p {threads} {input} 
		"""

#####
# PoolSeq
#####

rule bowtiePoolseq:
	input:
		read1 = outDir + "/poolseq/reads/{sample}_1.fastq.gz",
		read2 = outDir + "/poolseq/reads/{sample}_2.fastq.gz"
	output: "data/poolseq/bowtie2/{sample}.bam"
	threads: 4
	resources: mem=10, time = 60*24
	shell:
		"""
		module load hisat/2.1.0		
		bowtie2 -p {threads} -x {refBowtie2} -1 {input.read1} -2 {input.read2} | samtools view -Su -q5 - | samtools sort - -m 10G -T {SCRATCH}/hisat2/{wildcards.sample} -o {output} && sleep 5m
		samtools index {output}
		"""

# rule tophat:
# 	input:
# 		read1 = rules.trimRNAseq.output.read1,
# 		read2 = rules.trimRNAseq.output.read2
# 	output:
# 		"data/rnaseq/tophat/{sample}/accepted_hits.bam"
# 	threads: 12
# 	resources: mem=80, time = 60*24*7
# 	shell:
# 		"""
# 		module load python/2.7.8 tophat bowtie2
# 		tophat --b2-very-sensitive --b2-N 1 --mate-inner-dist 300 --mate-std-dev 100 -o data/rnaseq/tophat/{wildcards.sample} -p {threads} {refBowtie2} {input}
# 		"""

# rule trinityGG:
# 	input: expand("data/rnaseq/tophat/{sample}/accepted_hits.bam", sample = RNAseq)
# 	output: 	outDir + "/rnaseq/trinityGG_out_dir/Trinity-GG.fasta"
# 	threads: 12
# 	resources: mem=500, time = 60*24*7
# 	params:	trinityDir = outDir + "/rnaseq/trinityGG_out_dir"

# 	shell:
# 		"""
# 		module load trinity/2.5.1 bowtie2 samtools
# 		samtools merge --threads {threads} -u {params.trinityDir}/reads.bam {input}
# 		Trinity --genome_guided_bam {params.trinityDir}/reads.bam \
# 		--max_memory {resources.mem}G --CPU {threads} --output {params.trinityDir} --genome_guided_max_intron 5000 \
# 		&& rm {params.trinityDir}/reads.bam
# 		"""

# rule trinity:
# 	input:
# 		read1 = expand(rules.trimRNAseq.output.read1, sample = RNAseq),
# 		read2 = expand(rules.trimRNAseq.output.read1, sample = RNAseq)
# 	params:
# 		read1 = lambda wildcards, input: ",".join(input.read1),
# 		read2 = lambda wildcards, input: ",".join(input.read2),
# 		trinityDir = outDir + "/rnaseq/trinity_out_dir"
# 	output: 	outDir + "/rnaseq/trinity_out_dir/Trinity.fasta"
# 	resources: mem=500, time = 60*24*7
# 	threads: 12
# 	shell:
# 		"""
# 		module load trinity/2.5.1 bowtie2
# 		Trinity --seqType fq --max_memory {resources.mem}G --output {params.trinityDir} --CPU {threads} --SS_lib_type FR --left {params.read1} --right {params.read2}
# 		"""



######
# RAD-seq
######

rule sortRAD:
	input: outDir + "/radseq/reads/Undetermined_S0_L00{lane}_R{read}_001.fastq.gz"
	output: expand(outDir + "/radseq/readsSorted/{{lane}}/{{read}}/{sample}.fq.gz", sample = RADsamples)
	resources: mem=5, time = 24*60*3
	shell: "python2 scripts/sort_barcodes.py data/radseq/reads/bcodes.txt {input} data/radseq/readsSorted/{wildcards.lane}/{wildcards.read}/"

#this script merges the results, sorts them to keep the unique sequences (including adaptor barcode)
rule mergeRAD:
	input:
		read1 = expand(outDir + "/radseq/readsSorted/{lane}/1/{{sample}}.fq.gz", lane = range(1,9)),
		read2 = expand(outDir + "/radseq/readsSorted/{lane}/2/{{sample}}.fq.gz", lane = range(1,9))
	output:
		read1 = outDir + "/radseq/readsSorted/{sample}_1.fq.gz",
		read2 = outDir + "/radseq/readsSorted/{sample}_2.fq.gz"
	resources: mem=20, time = 60*20
	shell:
		"""
		paste  \
		<(zcat {input.read1} | paste - - - -  ) \
		<(zcat {input.read2} | paste - - - -  ) \
		| awk -F"\\t" ' {{split($1, a, "_") ; seq[substr(a[2],1,4)""$2""$6] = $0}} END {{for (s in seq) print seq[s]}}' \
		| tee >(cut -f 1-4 | tr '\\t' '\\n' | gzip > {output.read1}) \
		| cut -f 5-8 | tr '\\t' '\\n' | gzip > {output.read2}
		"""

# map reads and downsample to 100 reads max per site per sample
rule mapRAD:
	input:
		read1 = outDir + "/radseq/readsSorted/{sample}_1.fq.gz",
		read2 = outDir + "/radseq/readsSorted/{sample}_2.fq.gz",
	output: outDir + "/radseq/align/{sample}.bam"
	resources: mem=20, time = 60*20
	threads: 2
	shell:
		"""
		module load bowtie2 samtools VariantBam
		bowtie2 -p {threads} --sam-rg ID:{wildcards.sample} --sam-rg LB:FASSST --sam-rg SM:{wildcards.sample} --sam-rg PL:ILLUMINA -x {refBowtie2} -1 {input.read1} -2 {input.read2} | samtools view -Su - | samtools sort - -m 10G -T {SCRATCH}/bowtie/{wildcards.sample} -o - | variant - -m 100 -b -o {output} && samtools index {output}
		"""

# rule callRad:
# 	input: expand(outDir + "/radseq/align/{sample}.bam", sample = RADsamples)
# 	output: outDir + "/radseq/var/freebayes.vcf"
# 	resources: mem=100, time = 60*24*2
# 	threads: 1
# 	shell:
# 		"""
# 		module load freebayes vcflib
# 		freebayes --min-alternate-fraction 0.2 --use-best-n-alleles 4 -m 5 -q 5 -b {outDir}/radseq/align/*.bam -f {ref} | vcffilter -f "QUAL > 40" > {output}
# 		"""

rule freeBayes:
	input: expand(outDir + "/radseq/align/{sample}.bam", sample = RADsamples)
	output: temp(outDir + "/radseq/var/split/freebayes.{region}.vcf")
	params:
		span = lambda wildcards: REGIONS[wildcards.region],
		bams = lambda wildcards, input: os.path.dirname(input[0]) + "/*.bam",
		missing = lambda wildcards, input: len(input) * 0.3
	resources: mem=10, time = 60*5
	threads: 1
	shell:
		"""
		module load freebayes/1.3.1 vcftools/0.1.12b vcflib/1.0.0-rc1
		freebayes --min-alternate-fraction 0.2 --use-best-n-alleles 3 -m 5 -q 5 -b {params.bams} {params.span} -f {ref} | vcffilter  -f "QUAL > 20 & NS > {params.missing} & AF > 0.1 & AF < 0.9" > {output}
		"""

rule mergeVCF:
	input:
		expand(outDir + "/radseq/var/split/freebayes.{region}.vcf", region = REGIONS)
	output: outDir + "/radseq/var/freebayes.vcf"
	shell:
		"""
		module load vcflib/1.0.0-rc1
		(grep "^#" {input[0]} ; cat {input} | grep -v "^#" ) | vcfuniq  > {output}
		"""
