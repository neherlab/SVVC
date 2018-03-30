nuc_alpha = np.array(['A', 'C', 'G', 'T', '-', 'N'], dtype='S1')
alpha = nuc_alpha

def sam_to_allele_counts(sam_fname, paired=False, qual_min=30, max_reads=-1, max_isize = 700, VERBOSE = 0):
    '''
    calculates the allele counts for a set of mapped reads
    parameters:
    sam_fname   --   sam or bam file with mapped reads
    paired      --   differentiates between read one or two if True
                     otherwise the two reads are lumped together
    max_isize   --   maximal insert sizes to consider. this can be used to remove artifactual mappings
    qual_min    --   Ignore bases with quality less than qmin
    '''
    import numpy as np
    import pysam
    from collections import defaultdict

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
        for nref in xrange(samfile.nreferences):
            if VERBOSE: print "allocating for:", samfile.getrname(nref), "length:", samfile.lengths[nref]
            ac.append((samfile.getrname(nref), ac_array(samfile.lengths[nref], paired),
                        insertion_data_structure(paired)))

        # Iterate over single reads
        for i, read in enumerate(samfile):
            # Max number of reads
            if i == max_reads:
                if VERBOSE >= 2:
                    print 'Max reads reached:', max_reads
                break

            if read.is_unmapped or np.abs(read.isize)>max_isize:
                continue
            # Print output
            if (VERBOSE > 2) and (not ((i +1) % 10000)):
                print (i+1)

            # Read CIGARs (they should be clean by now)
            if paired:
                counts = ac[read.rname][1][int(read.is_read2),int(read.is_reverse)]
                insertion = ac[read.rname][2]
            else:
                counts = ac[read.rname][1][int(read.is_reverse)]
                insertion = ac[read.rname][2]

            seq = np.fromstring(read.seq, 'S1')
            qual = np.fromstring(read.qual, np.int8) - 33
            pos = read.pos
            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(read.cigar):
                if block_type==4: # softclip
                    seq = seq[block_len:]
                    qual = qual[block_len:]
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
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = ((seqb == a) & (qualb >= qual_min)).nonzero()[0]
                        if len(posa):
                            counts[j,pos + posa] += 1

                    # Chop off this block
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
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
                    # Accept only high-quality inserts
                    if (qualb >= qual_min).all():
                        if paired:
                            insertion[pos][seqb.tostring()][int(read.is_read2), int(read.is_reverse)] += 1
                        else:
                            insertion[pos][seqb.tostring()][int(read.is_reverse)] += 1

                    # Chop off seq, but not pos
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]

                # Other types of cigar?
                else:
                    if VERBOSE>2:
                        print "unrecognized CIGAR type:", read.cigarstring
                    #raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    return ac

def dump_allele_counts(dirname, ac, suffix=''):
    import numpy as np
    import cPickle, gzip, os
    dirname = dirname.rstrip('/')+'/'
    if not os.path.isdir(dirname):
        print "creating directory", dirname
        try:
            os.mkdir(dirname)
        except:
            raise("creating directory failed", dirname)

    for refname, ac_array, insertions in ac:
        print refname
        outname = refname.split("|")[-1]
        np.savez_compressed(dirname + outname+'_allele_counts' + suffix + '.npz', ac_array)
        with gzip.open(dirname + outname+'_insertions' + suffix + '.pkl.gz','w') as outfile:
            cPickle.dump({k:dict(v) for k,v in insertions.iteritems()}, outfile)


def load_allele_counts(dirname, suffix=''):
    import numpy as np
    import cPickle, gzip, glob
    dirname = dirname.rstrip('/')+'/'
    tmp_ac = {}
    ac_flist = glob.glob(dirname+'*_allele_counts' + suffix + '.npz')
    for fname in ac_flist:
        #print "reading",fname
        refname = fname.split('/')[-1].rstrip('_allele_counts' + suffix + '.npz')
        tmp_ac[refname] = np.load(fname).items()[0][1]

    ac = []

    for refname in tmp_ac:
        ac.append((refname, tmp_ac[refname].sum(axis=0).sum(axis=0)))

    return ac

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='create allele counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bam_file',
                        help='bam file to pile up')

    parser.add_argument('--outdir',
                        help='directory to save results')
    args = parser.parse_args()

    ac = sam_to_allele_counts(params.bam_file, qual_min=30, VERBOSE=3, max_isize = 600, paired=True)
    ac_renamed = []
    for refname, counts, insertions in ac:
        ac_renamed.append((refname.replace('/', '_'), counts, insertions))
    ac = ac_renamed

    dump_allele_counts(params.outdir, ac)
