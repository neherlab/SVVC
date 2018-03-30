from __future__ import print_function
from create_allele_counts import load_allele_counts
import glob, sys,os, argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set_style('whitegrid')

nuc_alpha = np.array(['A', 'C', 'G', 'T', '-', 'N'], dtype='S1')
alpha = nuc_alpha
gap = 200 # space between segments

def plot_coverage_concatenated(sample, ac, figure_path):
    coverage_stat = {}
    plt.figure(figsize = (18, 8))
    offset = 0
    ticks = []
    plt.title('sample %s'%sample)
    print("Sample", sample)
    for ref, counts in sorted(ac, key=lambda x:x[1].shape[-1], reverse=True):
        cov = coverage(counts)
        coverage_stat[ref] = {'cov':np.mean(cov)}
        seg = ref.split('_')[-1]
        plt.plot(offset+np.arange(counts.shape[-1]), coverage(counts), c='r')
        ticks.append([offset+counts.shape[-1]/2, seg])
        offset+=counts.shape[-1]+gap
        if (cov<100).mean()>0.20:
            print(sample, ref, "has very low coverage: %2.1f"%cov.mean(), "fraction blow 100: %1.2f"%(cov<100).mean())
        elif (cov<1000).mean()>0.20:
            print(sample, ref, "has low coverage: %2.1f"%cov.mean(), "fraction blow 100: %1.2f"%(cov<100).mean())


    #plt.xticks([t[0] for t in ticks], [t[1] for t in ticks])
    plt.yscale('log')
    plt.savefig(figure_path)
    plt.close()
    return coverage_stat

def plot_diversity(sample, ac, figure_path, primer_mask, min_cov=100, var_cutoff=0.05):
    fig, axs = plt.subplots(1,2, figsize = (20, 8))
    offset = 0
    ticks = []
    print("Sample", sample)
    diversity = {}
    for ref, counts in sorted(ac, key=lambda x:x[1].shape[-1], reverse=True):
        diversity[ref] = {}
        cov = coverage(counts)
        seg = ref.split('_')[-1]
        freq = 1.0*counts/(cov+0.001)
        div = (1-np.sum(freq**2, axis=0))*primer_mask[ref]
        div[cov<min_cov] = 0
        minor_allele = (1.0-np.max(freq, axis=0))*primer_mask[ref]
        minor_allele[cov<min_cov] = 0
        diversity[ref]['var_pos1'] = np.sum(minor_allele[0::3]>var_cutoff)
        diversity[ref]['var_pos2'] = np.sum(minor_allele[1::3]>var_cutoff)
        diversity[ref]['var_pos3'] = np.sum(minor_allele[2::3]>var_cutoff)
        diversity[ref]['mean_diversity'] = np.mean(div[cov>min_cov])
        diversity[ref]['mean_diversity>0.01'] = np.mean((div*(minor_allele>0.01))[cov>min_cov])
        diversity[ref]['mean_diversity>0.05'] = np.mean((div*(minor_allele>0.05))[cov>min_cov])
        print([np.sum(minor_allele[i::3]>var_cutoff) for i in range(3)])
        axs[0].plot(offset+np.arange(counts.shape[-1]), div, c='r')
        axs[1].plot(sorted(((1-np.max(freq, axis=0))*primer_mask[ref])[cov>min_cov]), np.linspace(1,0, np.sum(cov>min_cov)), c='r')
        offset+=counts.shape[-1]+gap
        plt.suptitle('sample %s'%sample + ' -- sites>0.05 in codon p1:%d, p2: %d, p3:%d'%(np.sum(minor_allele[0::3]>var_cutoff), np.sum(minor_allele[1::3]>0.05), np.sum(minor_allele[2::3]>0.05)))

    #plt.xticks([t[0] for t in ticks], [t[1] for t in ticks])
    axs[0].set_yscale('log')
    axs[0].set_ylabel('diversity')
    axs[1].set_xscale('log')
    axs[1].set_ylabel('fraction above')
    axs[1].set_xlabel('minor variant frequency')
    plt.savefig(figure_path)
    plt.close()
    return diversity


def coverage(ac, window = None):
    import numpy as np
    if window is None:
        return ac.sum(axis=-2)
    else:
        if np.isscalar(window):
            return np.convolve(ac.sum(axis=-2), np.ones(int(window))/float(window), 'same')
        else:
            return np.convolve(ac.sum(axis=-2), window, 'same')


def consensus(ac, min_cov=1):
    cov = coverage(ac)
    consensus_seq = alpha[np.argmax(ac, axis=-2)]
    consensus_seq[cov<min_cov]='N'
    return consensus_seq

# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='plot coverage, diversity and output consensus sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sample', required=True, type=str, help='the sample to analyze')
    parser.add_argument('--out_dir', required=True, type=str, help='directory to output')
    parser.add_argument('--primers', type=str, help='file with primers to mask in diversity calculation')
    parser.add_argument('--min_cov', type=int, default=100, help='minimal coverage to call consensus')

    args = parser.parse_args()
    stats = {}

    ac = load_allele_counts(args.sample)

    import pandas as pd
    if args.primers:
        primers = pd.read_csv(args.primers,skipinitialspace=True)

    primer_masks = {}
    for ref, counts in ac:
        primer_masks[ref] = np.ones(counts.shape[-1], dtype=int)
    if args.primers:
        for pi,p in primers.iterrows():
            if p.segment in primer_masks:
                primer_masks[p.segment][p.start:p.end]=0
            else:
                print(p.segment, "is not among the mapped segments")

    sample = args.sample.split('/')[-1]
    stats = plot_coverage_concatenated(sample, ac, args.out_dir+'/figures/%s_coverage.png'%sample)
    stats.update(plot_diversity(sample, ac, args.out_dir+"/figures/%s_diversity.png"%sample, primer_masks))
    from Bio import SeqIO, SeqRecord, Seq
    seqs=[]
    for ref, counts in ac:
        consensus_seq = consensus(counts, min_cov=args.min_cov)
        seq_name = sample+'_'+ref
        seqs.append(SeqRecord.SeqRecord(id=seq_name, name=seq_name, description="", seq=Seq.Seq("".join(consensus_seq))))
    SeqIO.write(seqs, args.out_dir+'/consensus.fasta', 'fasta')

    df = pd.DataFrame(stats)
    df.to_csv(args.out_dir+'/statistics.csv')
