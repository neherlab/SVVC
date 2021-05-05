# nohup snakemake mapped_reads/Q-A/mapped_reads.bam --jobs 32 --cluster "sbatch -t 05:59:00" 1>log2 &
reference = "..//reference/reference_seq.fasta"
segments = ["MN908947"]
run = ''
gen_outdir = "../results"
data_dir = "../data"

rule bwa_index:
    input:
        reference
    output:
        reference + '.bwt'
    shell:
        "bwa index {input}"

rule trim:
    input:
        r1 = data_dir+"/{sample}_1.fastq.gz",
        r2 = data_dir+"/{sample}_2.fastq.gz",
    output:
        r1 = gen_outdir+"/{sample}/trimmed_r1.fq.gz",
        r2 = gen_outdir+"/{sample}/trimmed_r2.fq.gz"
    params:
        base_name = "data/{sample}",
        outdir = gen_outdir + '/{sample}/trimming',
        min_length = 80,
        min_length_single = 90,
    shell:
        "trim_galore --length {params.min_length} --output {params.outdir} --retain_unpaired --paired -r1 {params.min_length_single} -r2 {params.min_length_single} {input.r1} {input.r2} &&"
        "mv {params.outdir}/{wildcards.sample}_1_val_1.fq.gz  {output.r1} &&"
        "mv {params.outdir}/{wildcards.sample}_2_val_2.fq.gz  {output.r2}"

rule map:
    input:
        ref = reference,
        index = rules.bwa_index.output,
        reads = rules.trim.output
    output:
        gen_outdir + "/{sample}/mapped_reads.bam"
    shell:
        "bwa mem {input.ref} {input.reads} |samtools view -Sb - > {output}"


rule pileup:
    input:
        gen_outdir + "/{sample}/mapped_reads.bam"
    output:
        gen_outdir + "/{sample}/" + "allele_counts.npz"
    params:
        path_to_script = 'src',
        out_dir = gen_outdir + "/{sample}"
    shell:
        "python3 {params.path_to_script}/create_allele_counts.py --bam_file {input} --out_dir {params.out_dir}"


rule pair_frequencies:
    input:
        gen_outdir + "/{sample}/mapped_reads.bam"
    output:
        gen_outdir + "/{sample}/" + "pair_counts.pkl.gz"
    params:
        path_to_script ='src',
        out_dir = gen_outdir + "/{sample}"
    shell:
        "python3 {params.path_to_script}/pair_statistics.py --bam_file {input} --out_dir {params.out_dir}"

rule consensus:
    input:
        gen_outdir + "/{sample}/" + "allele_counts.npz"
    output:
        gen_outdir + "/{sample}/consensus.fasta",
        gen_outdir + "/{sample}/figures/coverage.png",
        gen_outdir + "/{sample}/figures/diversity.png",
        gen_outdir + "/{sample}/minor.fasta"
    params:
        path_to_script =  'src',
        out_dir = gen_outdir + "/{sample}",
        min_freq = 0.01,
        min_cov = 1000
    shell:
        """
        echo {params.path_to_script}/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir} &
        python3 {params.path_to_script}/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir} &
        echo {params.path_to_script}/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir}  --min_freq {params.min_freq} --min_cov {params.min_cov} &
        python3 {params.path_to_script}/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir} --min_freq {params.min_freq} --min_cov {params.min_cov}
        """

