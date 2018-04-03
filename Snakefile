# nohup snakemake mapped_reads/Q-A/mapped_reads.bam --jobs 32 --cluster "sbatch -t 05:59:00" 1>log2 &
reference = "/scicore/home/neher/GROUP/data/2017_Karolinska_EV-D68/references/KX675261.fasta"
segments = ["KX675261.1"]
SVVC_dir = "./SVVC"
#run = "MIC4106"
run = "MIC3108"

rule trim:
	output:
		run+"/{sample}/trimmed_reads_1.fq.gz",
		run+"/{sample}/trimmed_reads_2.fq.gz"
	params:
		min_length = 80,
		min_length_single = 90,
		output_dir = run +"/{sample}"
	shell:
		"ml Trim_Galore &&"
		"trim_galore --length {params.min_length} --output {params.output_dir} --retain_unpaired --paired -r1 {params.min_length_single} -r2 {params.min_length_single} {params.output_dir}/*_?.fastq.gz &&"
		"cat {params.output_dir}/*val_1.fq.gz > {params.output_dir}/trimmed_reads_1.fq.gz &&"
		"cat {params.output_dir}/*val_2.fq.gz > {params.output_dir}/trimmed_reads_2.fq.gz &&"
		"rm {params.output_dir}/*val_?.fq.gz"

rule map:
	input:
		reference,
		run+"/{sample}/trimmed_reads_1.fq.gz",
		run+"/{sample}/trimmed_reads_2.fq.gz"
	output:
		"mapped_data/{sample}/mapped_reads.bam"
	shell:
		"ml BWA SAMtools &&"
		"bwa mem {input} |samtools view -Sb - > {output}"


rule pileup:
	input:
		"mapped_data/{sample}/mapped_reads.bam"
	output:
		["mapped_data/{sample}/" + "%s_allele_counts.npz"%seg for seg in segments]
	params:
		path_to_script = SVVC_dir + '/src',
		out_dir = "mapped_data/{sample}"
	shell:
		"python {params.path_to_script}/create_allele_counts.py --bam_file {input} --out_dir {params.out_dir}"

rule consensus:
	input:
		["mapped_data/{sample}/" + "%s_allele_counts.npz"%seg for seg in segments]
	output:
		"mapped_data/{sample}/consensus.fasta",
		"mapped_data/{sample}/figures/coverage.png",
		"mapped_data/{sample}/figures/diversity.png"
	params:
		path_to_script = SVVC_dir + '/src',
		out_dir = "mapped_data/{sample}"
	shell:
		"python2 {params.path_to_script}/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir}"
