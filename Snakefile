reference = "reference/KX675261.fasta"
segments = ["KX675261.1"]
run = ''
gen_outdir = "../results"
data_dir = "../data"

rule primers:
    params:
        primerfile = "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv"
    output:
        "primers.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO, Seq
        raw_primers = pd.read_csv(params.primerfile, sep='\t')
        ref = str(SeqIO.read(reference, 'fasta').seq)
        primers = {}
        for r, row in raw_primers.iterrows():
            start = ref.find(row.seq)
            if start<0:
                start = ref.find(Seq.reverse_complement(row.seq))

            if start>0:
                primers[row.name] = {"segment":segments[0], "name":row["name"], "seq":row.seq, "start":start, "end":start+len(row.seq)}
            else:
                print(f"row {row} failed")

        pd.DataFrame(primers).T.to_csv(output[0], sep='\t', index=False)


rule concat_data:
    input:
        frags = directory(data_dir + "/{sample}")
    output:
        r1 = data_dir + "/concat/{sample}/read1.fq.gz",
        r2 = data_dir + "/concat/{sample}/read2.fq.gz"
    run:
        import glob, gzip
        from Bio import SeqIO

        for r, readfile in zip(['R1', 'R2'], [output.r1, output.r2]):
            files = sorted(glob.glob(input.frags + f'/*{r}*'))
            with gzip.open(readfile, 'wt') as out_fh:
                for in_file in files:
                    with gzip.open(in_file, 'rt') as ifh:
                        for read in SeqIO.parse(ifh, 'fastq'):
                             SeqIO.write(read, out_fh, 'fastq')


rule bwa_index:
    input:
        reference
    output:
        reference + '.bwt'
    shell:
        "bwa index {input}"

rule trim:
    input:
        r1 = data_dir+"/concat/{sample}/read1.fq.gz",
        r2 = data_dir+"/concat/{sample}/read2.fq.gz",
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
        "mv {params.outdir}/read1_val_1.fq.gz  {output.r1} &&"
        "mv {params.outdir}/read2_val_2.fq.gz  {output.r2}"

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
        reads = gen_outdir + "/{sample}/mapped_reads.bam",
        primers = "primers.tsv"
    output:
        gen_outdir + "/{sample}/" + "allele_counts.npz"
    params:
        path_to_script = 'src',
        out_dir = gen_outdir + "/{sample}"
    shell:
        "python3 {params.path_to_script}/create_allele_counts.py --bam_file {input.reads} --primers {input.primers} --out_dir {params.out_dir}"


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
        min_freq = 0.05,
        min_cov = 100
    shell:
        """
        python3 {params.path_to_script}/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir} &
        python3 {params.path_to_script}/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir} --min_freq {params.min_freq} --min_cov {params.min_cov}
        """

