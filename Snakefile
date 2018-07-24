# nohup snakemake mapped_reads/Q-A/mapped_reads.bam --jobs 32 --cluster "sbatch -t 05:59:00" 1>log2 &
reference = "/scicore/home/neher/GROUP/data/2017_Karolinska_EV-D68/references/KX675261.fasta"
segments = ["KX675261.1"]
SVVC_dir = "/scicore/home/neher/GROUP/data/2017_Karolinska_EV-D68/SVVC"
run = '.'
#gen_outdir = "/scicore/home/neher/GROUP/data/2017_Karolinska_EV-D68/mapped_data"
gen_outdir = "../../../mapped_data"


rule trim:
	output:
		run+"/trimmed_reads_1.fq.gz",
		run+"/trimmed_reads_2.fq.gz"
	params:
		min_length = 80,
		min_length_single = 90,
	shell:
		"ml Trim_Galore &&"
		"trim_galore --length {params.min_length} --output . --retain_unpaired --paired -r1 {params.min_length_single} -r2 {params.min_length_single} *_?.fastq.gz &&"
		"cat *val_1.fq.gz > trimmed_reads_1.fq.gz &&"
		"cat *val_2.fq.gz > trimmed_reads_2.fq.gz &&"
		"rm *val_?.fq.gz"

rule primers:
	input:
		rules.trim.output
	output:
		run+"/primer_trimmed_reads_1.fq.gz",
		run+"/primer_trimmed_reads_2.fq.gz"
	shell:
		"""
		ml cutadapt &&
		cutadapt -a TCGGTACCTTTGTACGCCTGTTTTA$ -a TACTAGATTGYARTCCAAARTCCCA$ -a CAACTCCAGARATGCACATNCCAGG$ -a GTTCCATARCATCAGTATCTAACCA$ -a ACCCRATTTGYTTTGAAGGYCCAGG$ -a TGTAAATGRTTTCTCCAACAGATGC$ -a CCACCTTTGTATCAATAGCNGGTGT$ -a CTAACCATTTCCGTCTAAGACTAGA$ -a TAAAACAGGCGTACAAAGGTACCGA$ -a TGGGAYTTTGGAYTRCAATCTAGTA$ -a CCTGGNATGTGCATYTCTGGAGTTG$ -a TGGTTAGATACTGATGYTATGGAAC$ -a CCTGGRCCTTCAAARCAAATYGGGT$ -a GCATCTGTTGGAGAAAYCATTTACA$ -a ACACCNGCTATTGATACAAAGGTGG$ -g ^TCTAGTCTTAGACGGAAATGGTTAG -g ^TCGGTACCTTTGTACGCCTGTTTTA -g ^TACTAGATTGYARTCCAAARTCCCA -g ^CAACTCCAGARATGCACATNCCAGG -g ^GTTCCATARCATCAGTATCTAACCA -g ^ACCCRATTTGYTTTGAAGGYCCAGG -g ^TGTAAATGRTTTCTCCAACAGATGC -g ^CCACCTTTGTATCAATAGCNGGTGT -g ^CTAACCATTTCCGTCTAAGACTAGA -g ^TAAAACAGGCGTACAAAGGTACCGA -g ^TGGGAYTTTGGAYTRCAATCTAGTA -g ^CCTGGNATGTGCATYTCTGGAGTTG -g ^TGGTTAGATACTGATGYTATGGAAC -g ^CCTGGRCCTTCAAARCAAATYGGGT -g ^GCATCTGTTGGAGAAAYCATTTACA -g ^ACACCNGCTATTGATACAAAGGTGG -g ^TCTAGTCTTAGACGGAAATGGTTAG -o {output[0]} {input[0]} &&
		cutadapt -a TCGGTACCTTTGTACGCCTGTTTTA$ -a TACTAGATTGYARTCCAAARTCCCA$ -a CAACTCCAGARATGCACATNCCAGG$ -a GTTCCATARCATCAGTATCTAACCA$ -a ACCCRATTTGYTTTGAAGGYCCAGG$ -a TGTAAATGRTTTCTCCAACAGATGC$ -a CCACCTTTGTATCAATAGCNGGTGT$ -a CTAACCATTTCCGTCTAAGACTAGA$ -a TAAAACAGGCGTACAAAGGTACCGA$ -a TGGGAYTTTGGAYTRCAATCTAGTA$ -a CCTGGNATGTGCATYTCTGGAGTTG$ -a TGGTTAGATACTGATGYTATGGAAC$ -a CCTGGRCCTTCAAARCAAATYGGGT$ -a GCATCTGTTGGAGAAAYCATTTACA$ -a ACACCNGCTATTGATACAAAGGTGG$ -g ^TCTAGTCTTAGACGGAAATGGTTAG -g ^TCGGTACCTTTGTACGCCTGTTTTA -g ^TACTAGATTGYARTCCAAARTCCCA -g ^CAACTCCAGARATGCACATNCCAGG -g ^GTTCCATARCATCAGTATCTAACCA -g ^ACCCRATTTGYTTTGAAGGYCCAGG -g ^TGTAAATGRTTTCTCCAACAGATGC -g ^CCACCTTTGTATCAATAGCNGGTGT -g ^CTAACCATTTCCGTCTAAGACTAGA -g ^TAAAACAGGCGTACAAAGGTACCGA -g ^TGGGAYTTTGGAYTRCAATCTAGTA -g ^CCTGGNATGTGCATYTCTGGAGTTG -g ^TGGTTAGATACTGATGYTATGGAAC -g ^CCTGGRCCTTCAAARCAAATYGGGT -g ^GCATCTGTTGGAGAAAYCATTTACA -g ^ACACCNGCTATTGATACAAAGGTGG -g ^TCTAGTCTTAGACGGAAATGGTTAG -o {output[1]} {input[1]}
		"""

rule map:
	input:
		reference,
		rules.primers.output
	output:
		gen_outdir + "/{sample}/mapped_reads.bam"
	shell:
		"ml BWA SAMtools &&"
		"bwa mem {input} |samtools view -Sb - > {output}"


rule pileup:
	input:
		gen_outdir + "/{sample}/mapped_reads.bam"
	output:
		[gen_outdir + "/{sample}/" + "%s_allele_counts.npz"%seg for seg in segments]
	params:
		path_to_script = SVVC_dir + '/src',
		out_dir = gen_outdir + "/{sample}"
	shell:
		"python {params.path_to_script}/create_allele_counts.py --bam_file {input} --out_dir {params.out_dir}"

rule consensus:
	input:
		[gen_outdir + "/{sample}/" + "%s_allele_counts.npz"%seg for seg in segments]
	output:
		gen_outdir + "/{sample}/consensus.fasta",
		gen_outdir + "/{sample}/figures/coverage.png",
		gen_outdir + "/{sample}/figures/diversity.png"
	params:
		path_to_script = SVVC_dir + '/src',
		out_dir = gen_outdir + "/{sample}"
	shell:
		"python2 {params.path_to_script}/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir}"


rule minor:
	input:
		[gen_outdir + "/{sample}/" + "%s_allele_counts.npz"%seg for seg in segments]
	output:
		gen_outdir + "/{sample}/minor.fasta",
	params:
		path_to_script = SVVC_dir + '/src',
		out_dir = gen_outdir + "/{sample}"
	shell:
		"""
		echo {params.path_to_script}/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir} &
		python2 {params.path_to_script}/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir}
		"""