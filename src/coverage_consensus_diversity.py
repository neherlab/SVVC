
from collections import defaultdict
import glob, sys,os, argparse
import numpy as np
import matplotlib
# important to use a non-interactive backend, otherwise will crash on cluster
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from create_allele_counts import load_allele_counts
import seaborn as sns
import pandas as pd
sns.set_style('whitegrid')

nuc_alpha = np.array(['A', 'C', 'G', 'T', '-', 'N'], dtype='U1')
alpha = nuc_alpha
gap = 200 # space between segments

def plot_coverage_concatenated(sample, ac, figure_path, primer_boundaries=None):
    coverage_stat = {}
    plt.figure(figsize = (18, 8))
    offset = 0
    ticks = []
    plt.title('sample %s'%sample)
    cols = ['b','r']
    colnum = 0
    print("Sample", sample)
    for ref, counts in sorted(ac, key=lambda x:x[1].shape[-1], reverse=False):#True):
        print("ref", ref)
        if ref == '':
            legKey = "New"
            colour = 'r'
        else:
            legKey = "Old"
            colour = 'b'
        #legKey = "New" if ref == '' else "Old"
        cov = coverage(counts)
        coverage_stat[ref] = {'cov':np.mean(cov)}
        seg = ref.split('_')[-1]
        plt.plot(offset+np.arange(counts.shape[-1]), coverage(counts), c=colour, label=legKey)#'r')
        colnum = colnum+1 #EBH
        ticks.append([offset+counts.shape[-1]/2, seg])
        # offset+=counts.shape[-1]+gap
        if (cov<100).mean()>0.20:
            print(sample, ref, "has very low coverage: %2.1f"%cov.mean(), "fraction blow 100: %1.2f"%(cov<100).mean())
        elif (cov<1000).mean()>0.20:
            print(sample, ref, "has low coverage: %2.1f"%cov.mean(), "fraction blow 100: %1.2f"%(cov<100).mean())
        if primer_boundaries and ref in primer_boundaries:
            for p in primer_boundaries[ref]:
                y = 2 if int(p[1])%2 else 6
                plt.plot([primer_boundaries[ref][p]['start'], primer_boundaries[ref][p]['end']],[y,y], lw=10, c=(0.7, 0.7, 0.7))

    #plt.xticks([t[0] for t in ticks], [t[1] for t in ticks])
    plt.yscale('log')
    plt.legend(loc="lower center")
    plt.savefig(figure_path)
    plt.close()
    return coverage_stat

def plot_diversity(sample, ac, figure_path, primer_mask, min_cov=100, var_cutoff=0.05, primer_boundaries=None):
    fig, axs = plt.subplots(1,2, figsize = (20, 8))
    offset = 0
    ticks = []
    print("Sample", sample)
    diversity = {}
    cols = ['r', 'b', 'g']
    for ref, counts in sorted(ac, key=lambda x:x[1].shape[-1], reverse=True):
        if ref != '':
            continue
        diversity[ref] = {}
        cov = coverage(counts)
        seg = ref.split('_')[-1]
        freq = 1.0*counts/cov
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
        max_freq = np.max(freq, axis=0)
        pos = np.arange(max_freq.shape[0])
        for ci in range(3):
            sub_set = 1-max_freq[((pos%3)==ci) & (cov>1000) & (primer_mask[ref]>0)]
            axs[1].plot(sorted(sub_set), np.linspace(1,0, len(sub_set)), c=cols[ci])

        offset+=counts.shape[-1]+gap
        plt.suptitle('sample %s'%sample + ' -- sites>0.05 in codon p1:%d, p2: %d, p3:%d'%(np.sum(minor_allele[0::3]>var_cutoff), np.sum(minor_allele[1::3]>0.05), np.sum(minor_allele[2::3]>0.05)))

        if primer_boundaries and ref in primer_boundaries:
            for p in primer_boundaries[ref]:
                y = 0.002 if int(p[1])%2 else 0.001
                axs[0].plot([primer_boundaries[ref][p]['start'], primer_boundaries[ref][p]['end']],[y,y], lw=10, c=(0.7, 0.7, 0.7))

    #plt.xticks([t[0] for t in ticks], [t[1] for t in ticks])
    axs[0].set_yscale('log')
    axs[0].set_ylabel('diversity')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_ylabel('fraction above')
    axs[1].set_xlabel('minor variant frequency')
    plt.savefig(figure_path)
    plt.close()
    return diversity

def get_primer_mask(primer_file, ac):
    import pandas as pd

    primer_masks = {}
    for ref, counts in ac:
        primer_masks[ref] = np.ones(counts.shape[-1], dtype=int)

    if primer_file:
        primers = pd.read_csv(primer_file,skipinitialspace=True)
        for pi,p in primers.iterrows():
            if p.segment in primer_masks:
                primer_masks[p.segment][p.start:p.end]=0
            else:
                print(p.segment, "is not among the mapped segments")
    return primer_masks

def get_fragment_boundaries(primer_file, ac):
    import pandas as pd

    primer_boundaries = {}
    for ref, counts in ac:
        primer_boundaries[ref] = defaultdict(dict)

    if primer_file:
        primers = pd.read_csv(primer_file,skipinitialspace=True)
        for pi,p in primers.iterrows():
            if p.segment in primer_boundaries:
                name,fr = p.loc['name'].split('_')
                if fr=='fwd':
                    primer_boundaries[p.segment][name]['start']=max(p.start, p.end)
                elif fr=='rev':
                    primer_boundaries[p.segment][name]['end']=min(p.start, p.end)
            else:
                print(p.segment, "is not among the mapped segments")
    return primer_boundaries


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
    import pandas as pd

    # Parse input args
    parser = argparse.ArgumentParser(description='plot coverage, diversity and output consensus sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sample', required=True, type=str, help='the sample to analyze')
    parser.add_argument('--out_dir', required=True, type=str, help='directory to output')
    parser.add_argument('--primers', type=str, help='file with primers to mask in diversity calculation')
    parser.add_argument('--min_cov', type=int, default=100, help='minimal coverage to call consensus')
    parser.add_argument('--all_counts', action="store_true", default=False, help="plot coverage/diversity for all count files found")

    args = parser.parse_args()
    stats = {}

    ac,ins = load_allele_counts(args.sample, allCounts=args.all_counts)
    primer_masks = get_primer_mask(args.primers, ac)
    primer_boundaries = get_fragment_boundaries(args.primers, ac)


    sample = args.sample.split('/')[-1]
    stats = plot_coverage_concatenated(sample, ac, args.out_dir+'/figures/coverage.png',
                                       primer_boundaries=primer_boundaries)
    div = plot_diversity(sample, ac, args.out_dir+"/figures/diversity.png", primer_masks,
                         primer_boundaries=primer_boundaries)
    for k, v in list(div.items()):
        stats[k].update(v)
    from Bio import SeqIO, SeqRecord, Seq
    seqs=[]
    insertions_to_include = []
    for ref, counts in ac:
        if ref != '':
            continue
        consensus_seq = consensus(counts, min_cov=args.min_cov)
        cov = coverage(counts)
        print("cov", cov)
        for pos in ins[ref]:
            if cov[pos]<args.min_cov:
                continue
            total_insertion = np.sum([c.sum() for insertion, c in list(ins[ref][pos].items())])
            total_freq = 1.0*total_insertion/cov[pos]
            max_insertion = [pos, 0,0]
            for insertion, c in list(ins[ref][pos].items()):
                ins_freq = 1.0*c.sum()/cov[pos]
                if ins_freq>max_insertion[2]:
                    max_insertion = [pos, insertion, ins_freq]
                if ins_freq>0.3:
                    print(sample + ": frequent insertion %s at position %d with frequency %f."%(insertion, pos, ins_freq))
                elif ins_freq>0.01:
                    print(sample + ": rare insertion %s at position %d with frequency %f."%(insertion, pos, ins_freq))

            # id the most frequent insertion is more common than no insertion
            if 1-total_freq<max_insertion[2]:
                insertions_to_include.append(max_insertion)

        seq = "".join(consensus_seq)
        if insertions_to_include:
            complete_seq = ""
            pos = 0
            for ins_pos, ins, freq in insertions_to_include:
                complete_seq += seq[pos:ins_pos] + ins
                pos=ins_pos
                print(sample + ": inserted %s at position %d with frequency %f."%(ins, ins_pos, freq))
            complete_seq += seq[pos:]
            seq=complete_seq

        if len(ac)==1:
            seq_name = sample
        else:
            seq_name = sample + '_' + ref
        seqs.append(SeqRecord.SeqRecord(id=seq_name, name=seq_name, description="", seq=Seq.Seq(seq)))
    SeqIO.write(seqs, args.out_dir+'/consensus.fasta', 'fasta')

    df = pd.DataFrame(stats)
    df.to_csv(args.out_dir+'/statistics.csv')
