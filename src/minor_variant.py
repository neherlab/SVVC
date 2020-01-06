from __future__ import print_function
from create_allele_counts import load_allele_counts
import glob, sys,os, argparse
import numpy as np
import seaborn as sns
from coverage_consensus_diversity import coverage, consensus, nuc_alpha, alpha, get_primer_mask
sns.set_style('whitegrid')

def trim_ac(ac, n_states=5):
    tmp_ac = {}
    for ref, x in ac:
        tmp_ac[ref] = x[:n_states].astype(float)/np.sum(x[:n_states], axis=0)
    return tmp_ac

# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='plot coverage, diversity and output consensus sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sample', required=True, type=str, help='the sample to analyze')
    parser.add_argument('--out_dir', required=True, type=str, help='directory to output')
    parser.add_argument('--min_cov', type=int, default=1000, help='minimal coverage to call consensus')
    parser.add_argument('--primers', type=str, help='file with primers to mask in diversity calculation')
    parser.add_argument('--min_freq', type=float, default=0.05, help='minimal frequency to accept minor variant')

    args = parser.parse_args()
    stats = {}

    ac,ins = load_allele_counts(args.sample)
    if args.primers:
        primer_masks = get_primer_mask(args.primers, ac)

    freqs = trim_ac(ac)

    sample = args.sample.split('/')[-1]

    major_freq = {ref:np.max(x, axis=0) for ref, x in freqs.items()}

    minor_seqs = {}
    any_minors = False
    seqs = []
    from Bio import SeqIO, SeqRecord, Seq
    for ref, counts in ac:
        print("ref", ref)
        consensus_seq = consensus(counts, min_cov=args.min_cov)
        cov = coverage(counts)
        div_pos = np.where((major_freq[ref]<1.0-args.min_freq)&(cov>args.min_cov))[0]
        alterations = []
        insertions_to_include = []
        for pos in div_pos:
            tmp_freqs = freqs[ref][:, pos]
            if sorted(tmp_freqs)[-2]>args.min_freq:
                ii = np.argsort(tmp_freqs)[-2]
                alterations.append([pos, nuc_alpha[ii], tmp_freqs[ii]])

        if alterations:
            print(sample, ref, 'minor variants', alterations)
            consensus_seq[[p for p,n,f in alterations]] = [n for p,n,f in alterations]
            any_minors = True

        for pos in ins[ref]:
            if cov[pos]<args.min_cov:
                continue
            total_insertion = np.sum([c.sum() for insertion, c in ins[ref][pos].items()])
            total_freq = 1.0*total_insertion/cov[pos]
            if total_freq<args.min_freq:
                continue

            insertions = [[pos, '', 1-total_freq]]
            for insertion, c in ins[ref][pos].items():
                ins_freq = 1.0*c.sum()/cov[pos]
                insertions.append([pos, insertion, ins_freq])
            insertions.sort(key=lambda x:x[2])
            if insertions[-2][2]>args.min_freq:
                insertions_to_include.append(insertions[-2])

        seq = "".join(consensus_seq)
        if insertions_to_include:
            print(sample, ref, 'minor insertions', alterations)
            complete_seq = ""
            pos = 0
            for ins_pos, ins, freq in sorted(insertions_to_include, key=lambda x:x[0]):
                complete_seq += seq[pos:ins_pos] + ins
                pos=ins_pos
                print(sample + ": inserted %s at position %d with frequency %f."%(ins, ins_pos, freq))
            complete_seq += seq[pos:]
            seq=complete_seq
            any_minors = True

        if len(ac)==1:
            seq_name = sample+'_minor'
        else:
            seq_name = sample + '_minor_' + ref
        seqs.append(SeqRecord.SeqRecord(id=seq_name, name=seq_name, description="", seq=Seq.Seq(seq)))

    if any_minors:
        SeqIO.write(seqs, args.out_dir+'/minor.fasta', 'fasta')
    else:
        os.system("touch "+args.out_dir+'/minor.fasta')
