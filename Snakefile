# Euphy QTL
from snakemake.utils import R

outDir = "data"   # symlink to /work
refDir = "ref" 	  # symlink to /work
SCRATCH  = "/work/scratch/sasha"
mcinBT2  = "ref/Mcin"


RNAseq, = glob_wildcards(outDir + "/rnaseq/reads/{sample}.R1.fastq.gz")

print(expand(outDir + "/rnaseq/readsTrimmed/{sample}.R{R}.fastq.gz", sample=RNAseq, R=[1,2]))

rule all:
	input:  expand(outDir + "/rnaseq/tophat/{sample}/accepted_hits.bam", sample = RNAseq) , outDir + "/rnaseq/trinity_out_dir/Trinity.fasta"

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

rule tophat:
	input:
		read1 = rules.trimRNAseq.output.read1,
		read2 = rules.trimRNAseq.output.read2
	output:
		"data/rnaseq/tophat/{sample}/accepted_hits.bam"
	threads: 12
	resources: mem=40, time = 60*24*7 
	shell:
		"""
		module load python/2.7.8 tophat bowtie2
		tophat --b2-very-sensitive --b2-N 1 --mate-inner-dist 300 --mate-std-dev 100 -o data/rnaseq/tophat/{wildcards.sample} -p {threads} {mcinBT2} {input}
		"""

rule trinity:
	input:
		read1 = expand(rules.trimRNAseq.output.read1, sample = RNAseq),
		read2 = expand(rules.trimRNAseq.output.read1, sample = RNAseq)
	params:
		read1 = lambda wildcards, input: ",".join(input.read1),
		read2 = lambda wildcards, input: ",".join(input.read2),
		trinityDir = outDir + "/rnaseq/trinity_out_dir" 
	output: 	outDir + "/rnaseq/trinity_out_dir/Trinity.fasta"
	resources: mem=500, time = 60*24*7 
	threads: 12
	shell: 
		"""
		module load trinity/2.5.1 bowtie2
		Trinity --seqType fq --max_memory {resources.mem}G --output {params.trinityDir} --CPU {threads} --SS_lib_type FR --left {params.read1} --right {params.read2}
		"""

# rule A:
# 	input:
# 		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
# 		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
# 	output:
# 		temp(outDir + "/meta/hosts/{sample}-{q}.txt")
# 	threads: 12
# 	shell:
# 		"""
# 		"""

