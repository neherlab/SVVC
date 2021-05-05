import numpy as np
nuc_alpha = np.array(['A', 'C', 'G', 'T', '-', 'N'], dtype='S1')

def insert_size_distribution(sam_fname):
    '''
    '''
    import pysam

    smax=1000
    isize_dis = np.zeros(smax)
    # Open BAM or SAM file
    with pysam.Samfile(sam_fname) as samfile:
        for i, read in enumerate(samfile):
            if np.abs(read.isize)<smax:
                isize_dis[np.abs(read.isize)]+=1
    return isize_dis


def sam_to_allele_counts(sam_fname, paired=False, qual_min=30, max_reads=-1,
                         max_isize = 700, VERBOSE = 0,
                         fwd_primer_regions = None, rev_primer_regions = None):
    '''
    calculates the allele counts for a set of mapped reads
    parameters:
    sam_fname   --   sam or bam file with mapped reads
    paired      --   differentiates between read one or two if True
                     otherwise the two reads are lumped together
    max_isize   --   maximal insert sizes to consider. this can be used to remove artifactual mappings
    qual_min    --   Ignore bases with quality less than qmin
    '''
    import pysam
    from collections import defaultdict

    alpha = nuc_alpha

    def ac_array(length, paired):
        if paired:
            return np.zeros((2,2,6,length), dtype =int)
        else:
            return np.zeros((2,6,length), dtype = int)

    # Note: the data structure for inserts is a nested dict with:
    # position --> string --> read type --> count
    #  (dict)      (dict)       (list)      (int)
    def insertion_data_structure(paired):
        if paired:
            return defaultdict(lambda: defaultdict(lambda: np.zeros((2,2), int)))
        else:
            return defaultdict(lambda: defaultdict(lambda: np.zeros(2, int)))


    # Open BAM or SAM file
    with pysam.Samfile(sam_fname) as samfile:
        ac =  []
        refs = {}
        for nref in range(samfile.nreferences):
            if VERBOSE: print(("allocating for:", samfile.getrname(nref), "length:", samfile.lengths[nref]))
            refs[nref]=samfile.getrname(nref)
            ac.append((samfile.getrname(nref), ac_array(samfile.lengths[nref], paired),
                        insertion_data_structure(paired)))

        # Iterate over single reads
        for i, read in enumerate(samfile):
            # Max number of reads
            if i == max_reads:
                if VERBOSE >= 2:
                    print(('Max reads reached:', max_reads))
                break

            if read.is_unmapped or np.abs(read.isize)>max_isize or read.is_secondary or read.is_supplementary:
                continue

            # Print output
            if (VERBOSE > 2) and (not ((i +1) % 10000)):
                print((i+1))

            # Read CIGARs (they should be clean by now)
            if paired:
                counts = ac[read.rname][1][int(read.is_read2),int(read.is_reverse)]
                insertion = ac[read.rname][2]
            else:
                counts = ac[read.rname][1][int(read.is_reverse)]
                insertion = ac[read.rname][2]

            seq = np.fromstring(read.seq, 'S1')
            qual = np.fromstring(read.qual, np.int8) - 33
            not_primer = np.ones_like(seq, 'bool')
            pos = read.pos
            # all legit reads should be FR or RF!
            if rev_primer_regions:
                if read.is_reverse or np.abs(read.isize)==seq.shape[0]:
                    read_end = pos + seq.shape[0]
                    for b,e in rev_primer_regions[refs[read.rname]]:
                        p_length = e-b
                        if read_end-b>0 and read_end-b<p_length:
                            not_primer[-(read_end-b):]=False
                            break

            if fwd_primer_regions:
                if (not read.is_reverse) or np.abs(read.isize)==seq.shape[0]:
                    for b,e in fwd_primer_regions[refs[read.rname]]:
                        p_length = e-b
                        if pos-b>0 and pos-b<p_length:
                            not_primer[:e-pos]=False
                            break

            # if pos+len(seq)>7267:
            # 	import ipdb;ipdb.set_trace()
            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(read.cigar):
                if block_type==4: # softclip
                    seq = seq[block_len:]
                    qual = qual[block_len:]
                    # not the difference here: the reported position starts after the softclip. hence the not_primer is already correct
                    not_primer = not_primer[:-block_len]
                    continue
                if block_type==5: # hard clip
                    continue

                # Check for pos: it should never exceed the length of the fragment
#                if (block_type in [0, 1, 2]) and (pos >= length):
#                    raise ValueError('Pos exceeded the length of the fragment')

                # Inline block
                if block_type == 0:
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    not_primerb = not_primer[:block_len]
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = ((seqb == a) & (qualb >= qual_min) & (not_primerb)).nonzero()[0]
                        if len(posa):
                            counts[j,pos + posa] += 1

                    # Chop off this block
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        not_primer = not_primer[block_len:]
                        pos += block_len

                # Deletion
                elif block_type == 2:
                    # Increment gap counts
                    counts[4, pos:pos + block_len] += 1
                    # Chop off pos, but not sequence
                    pos += block_len

                # Insertion
                # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                # THEN the insert, FINALLY comes seq[391:]
                elif block_type == 1:
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    not_primerb = not_primer[:block_len]
                    # Accept only high-quality inserts
                    if (qualb >= qual_min).all():
                        if paired:
                            insertion[pos][seqb.tostring().decode()][int(read.is_read2), int(read.is_reverse)] += 1
                        else:
                            insertion[pos][seqb.tostring().decode()][int(read.is_reverse)] += 1

                    # Chop off seq, but not pos
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        not_primer = not_primer[block_len:]

                # Other types of cigar?
                else:
                    if VERBOSE>2:
                        print(("unrecognized CIGAR type:", read.cigarstring))
                    #raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    return ac

def dump_allele_counts(dirname, ac, suffix=''):
    import pickle, gzip, os
    dirname = dirname.rstrip('/')+'/'
    if not os.path.isdir(dirname):
        print(("creating directory", dirname))
        try:
            os.mkdir(dirname)
        except:
            raise "creating directory failed"

    for refname, ac_array, insertions in ac:
        print(refname)
        np.savez_compressed(dirname + 'allele_counts' + suffix + '.npz', ac_array)
        with gzip.open(dirname + 'insertions' + suffix + '.pkl.gz','w') as outfile:
            pickle.dump({k:dict(v) for k,v in insertions.items()}, outfile)


def get_primer_intervals(primer_file):
    import pandas as pd
    from collections import defaultdict

    fwd_primer_intervals = defaultdict(list)
    rev_primer_intervals = defaultdict(list)

    if primer_file:
        primers = pd.read_csv(primer_file,skipinitialspace=True)
        for pi,p in primers.iterrows():
            print((pi,p.loc['name']))
            if 'fwd' in p.loc['name']:
                fwd_primer_intervals[p.segment].append(sorted((p.start, p.end)))
            else:
                rev_primer_intervals[p.segment].append(sorted((p.start, p.end)))

    return fwd_primer_intervals, rev_primer_intervals


def load_allele_counts(dirname, suffix='', allCounts=False):
    import pickle, gzip, glob
    dirname = dirname.rstrip('/')+'/'
    tmp_ac = {}
    if allCounts:
        ac_flist = glob.glob(dirname+'*allele_counts' + suffix + '.npz')
    else:
        ac_flist = glob.glob(dirname+'allele_counts' + suffix + '.npz')
    for fname in ac_flist:
        #print("reading",fname)
        tmp = '_allele_counts' + suffix + '.npz'
        refname = fname.split('/')[-1][:-len(tmp)]
        tmp_ac[refname] = list(np.load(fname).items())[0][1]

    ins_flist = glob.glob(dirname+'*insertions' + suffix + '.pkl.gz')
    tmp_ins = {}
    for fname in ins_flist:
        #print("reading",fname)
        tmp = '_insertions' + suffix + '.pkl.gz'
        refname = fname.split('/')[-1][:-len(tmp)]
        with gzip.open(fname) as fh:
            tmp_ins[refname] = pickle.load(fh)

    ac = []
    for refname in tmp_ac:
        ac.append((refname, tmp_ac[refname].sum(axis=0).sum(axis=0)))

    # ins = []
    # for refname in tmp_ins:
    #     ins.append((refname, tmp_ins[refname]))

    return ac, tmp_ins

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='create allele counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bam_file',
                        help='bam file to pile up')

    parser.add_argument('--out_dir',
                        help='directory to save results')
    parser.add_argument('--primers', type=str, help='file with primers to mask in pile up')
    args = parser.parse_args()
    fwd_primer_intervals, rev_primer_intervals = get_primer_intervals(args.primers)
    print((fwd_primer_intervals, rev_primer_intervals))
    ac = sam_to_allele_counts(args.bam_file, qual_min=30, VERBOSE=3, max_isize = 600, paired=True,
                              fwd_primer_regions = fwd_primer_intervals, rev_primer_regions = rev_primer_intervals)
    ac_renamed = []
    for refname, counts, insertions in ac:
        ac_renamed.append((refname.replace('/', '_'), counts, insertions))
    ac = ac_renamed

    dump_allele_counts(args.out_dir, ac)
