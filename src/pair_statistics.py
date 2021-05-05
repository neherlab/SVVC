from create_allele_counts import get_primer_intervals

def pair_counts(sam_fname, paired=False, qual_min=30, max_reads=-1,
                max_isize = 700, VERBOSE = 0,
                fwd_primer_regions = None, rev_primer_regions = None):
    '''
    '''
    import numpy as np
    import pysam
    from collections import defaultdict
    from itertools import combinations

    a = {(True,True):0, (True, False):0, (False, True):0, (False, False):0}
    c = {'R':0, 'L':0, 'N':0, 'E':0}
    nuc_alpha = np.array(['A', 'C', 'G', 'T'], dtype='S1')
    # Open BAM or SAM file
    with pysam.Samfile(sam_fname) as samfile:
        ac =  []
        acc = []
        refs = {}
        read_count = 0
        for nref in range(samfile.nreferences):
            if VERBOSE: print(("allocating for:", samfile.getrname(nref), "length:", samfile.lengths[nref]))
            refs[nref]=samfile.getrname(nref)
            ac.append((samfile.getrname(nref),  np.zeros((len(nuc_alpha),samfile.lengths[nref]), dtype =int)))
            acc.append((samfile.getrname(nref), {}))

        while True:
            # find read pairs and skip secondary or supplementary alignments
            try:
                read1 = next(samfile)
                while read1.is_secondary or read1.is_supplementary:
                    read1 = next(samfile)
                read2 = next(samfile)
                while read2.is_secondary or read2.is_supplementary:
                    read2 = next(samfile)
            except:
                break

            if read1.is_unmapped or read2.is_unmapped or np.abs(read1.isize)>max_isize:
                continue
            if (read1.is_reverse==read2.is_reverse):
                continue
            if (read1.qname!=read2.qname):
                continue

            read_count+=1
            if read_count%1000==0:
                print(read_count)
                if max_reads>0 and read_count>max_reads:
                    break

            ref_name = refs[read1.rname]
            # determine which read maps to the 5p and which one the 3p end
            # pull out only positions that map, indels will be ignored in cocounts
            if read2.is_reverse:
                aln1 = np.array(read1.get_aligned_pairs(matches_only=True))
                aln2 = np.array(read2.get_aligned_pairs(matches_only=True))
                seq1 = np.fromstring(read1.seq, 'S1')[aln1[:,0]]
                qual1 = np.fromstring(read1.qual, np.int8)[aln1[:,0]] - 33
                seq2 = np.fromstring(read2.seq, 'S1')[aln2[:,0]]
                qual2 = np.fromstring(read2.qual, np.int8)[aln2[:,0]] - 33
            else:
                aln2 = np.array(read1.get_aligned_pairs(matches_only=True))
                aln1 = np.array(read2.get_aligned_pairs(matches_only=True))
                seq1 = np.fromstring(read2.seq, 'S1')[aln1[:,0]]
                qual1 = np.fromstring(read2.qual, np.int8)[aln1[:,0]] - 33
                seq2 = np.fromstring(read1.seq, 'S1')[aln2[:,0]]
                qual2 = np.fromstring(read1.qual, np.int8)[aln2[:,0]] - 33

            isize = np.abs(read1.isize)
            L1 = aln1.shape[0]
            L2 = aln2.shape[0]

            ## merge reads
            # allocate vectors
            merged_qual = np.zeros(isize, dtype=int)
            merged_seq = np.zeros(isize, dtype='S1')
            merged_pos = np.zeros((isize,2), dtype=int)

            # handle edge cases where one read in contained in the other,
            # i.e. the 5p read extends for longer than the 3p end of the 3p read
            # This can result for example from quality trimming.
            leftoverhang = aln1[0,1] - aln2[0,1]
            rightoverhang = aln1[-1,1] - aln2[-1,1]
            if leftoverhang>0: # take only the better read2
                merged_pos=aln2
                merged_qual=qual2
                merged_seq=qual2
                c['L']+=1
            elif rightoverhang>0: # take only the better read1
                merged_pos=aln1
                merged_qual=qual1
                merged_seq=qual1
                c['R']+=1
            else: # proper merging happens here
                # difference between end of aln1 and beginning of aln2 is overlap on reference
                overlap = max(0, aln1[-1,1] - aln2[0,1]+1)
                c['N']+=1

                # note that the exact coordinates might be off bc of indels
                # but what we are doing is conservate and only mapped positions
                # will be reported
                seg1 = L1 - overlap # end of non-overlap segment
                seg3 = isize - L2 + overlap # beginnning of non-overlap segment

                if seg1>0:
                    merged_pos[:seg1] = aln1[:seg1]
                    merged_qual[:seg1] = qual1[:seg1]
                    merged_seq[:seg1] = seq1[:seg1]
                else:
                    seg1=0

                merged_pos[seg3:] = aln2[overlap:]
                merged_qual[seg3:] = qual2[overlap:]
                merged_seq[seg3:] = seq2[overlap:]

                if overlap:
                    try:
                        seq_agree = (seq1[seg1:]==seq2[:overlap])&(aln1[seg1:,1]==aln2[:overlap,1])
                        better = qual1[seg1:]<qual2[:overlap]

                        from1 = np.where(seq_agree&better)[0]
                        from2 = np.where(seq_agree&(~better))[0]

                        merged_pos[seg1 + from1] = aln1[seg1 + from1]
                        merged_qual[seg1 + from1] = qual1[seg1 + from1]
                        merged_seq[seg1+from1] = seq1[seg1+from1]

                        merged_pos[seg1 + from2] = aln2[from2]
                        merged_qual[seg1 + from2] = qual2[from2]
                        merged_seq[seg1+from2] = seq2[from2]
                    except:
                        c['E']+=1
                        continue

            # mask regions in the merged read that likely derive from primer sequence
            not_primer = np.ones_like(merged_seq, 'bool')
            if rev_primer_regions:
                read_end = merged_pos[-1,1]
                for b,e in rev_primer_regions[ref_name]:
                    p_length = e-b
                    if read_end-b>0 and read_end-b<p_length:
                        not_primer[-(read_end-b):]=False
                        break

            if fwd_primer_regions:
                read_start = merged_pos[0,1]
                for b,e in fwd_primer_regions[ref_name]:
                    p_length = e-b
                    if read_start-b>0 and read_start-b<p_length:
                        not_primer[:e-read_start]=False
                        break

            counts = ac[read1.rname][1]
            cocounts = acc[read1.rname][1]
            good_ind = (merged_qual>qual_min)&not_primer
            for ni,nuc in enumerate(nuc_alpha):
                correct_state = merged_seq==nuc
                counts[ni,merged_pos[correct_state&good_ind,1]] += 1

            combo = list(zip(merged_pos[good_ind], merged_seq[good_ind]))
            for (p1, n1), (p2,n2) in combinations(combo, 2):
                posp = (p1[1], p2[1])
                p = n1+n2
                if posp not in cocounts:
                    cocounts[posp]={p:1}
                    continue
                if p not in cocounts[posp]:
                    cocounts[posp][p]=1
                else:
                    cocounts[posp][p]+=1

    return ac, acc


if __name__ == '__main__':
    import argparse, gzip
    import pickle as pickle
    parser = argparse.ArgumentParser(description='create pair counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bam_file',
                        help='bam file to pile up')

    parser.add_argument('--out_dir',
                        help='directory to save results')
    parser.add_argument('--max_reads', type=int,default=-1,
                        help='maximum number of reads to process')
    parser.add_argument('--primers', type=str, help='file with primers to mask in pile up')
    args = parser.parse_args()

    fwd_primer_intervals, rev_primer_intervals = get_primer_intervals(args.primers)

    print((fwd_primer_intervals, rev_primer_intervals))
    ac, acc = pair_counts(args.bam_file, qual_min=30, VERBOSE=3, max_isize = 600, paired=True, max_reads=args.max_reads,
                          fwd_primer_regions = fwd_primer_intervals, rev_primer_regions = rev_primer_intervals)
    acc_renamed = []
    for refname, counts in acc:
        acc_renamed.append((refname.replace('/', '_'), counts))
    acc = acc_renamed

    ac_renamed = []
    for refname, counts in ac:
        ac_renamed.append((refname.replace('/', '_'), counts))
    ac = ac_renamed

    with gzip.open(args.out_dir+'/pair_counts.pkl.gz', 'w') as fh:
        pickle.dump((ac,acc), fh)

