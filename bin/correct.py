#!/fs/sz-user-supported/Linux-i686/bin/python2.5
from optparse import OptionParser
import os, random
import dna

############################################################
# correct.py
#
# Launch pipeline to correct errors in sequencing reads.
############################################################

# options
model_gc = False
qmers = True

r_dir = '/nfshomes/dakelley/research/error_correction/bin'

############################################################
# main
############################################################
def main():
    parser = OptionParser()
    parser.add_option('-r', dest='readsf', help='Fastq file of reads')
    parser.add_option('-k', dest='k', type='int', help='Size of k-mers to correct')
    parser.add_option('-p', dest='proc', type='int', default=4, help='Number of processes')
    parser.add_option('--no_count', dest='no_count', action='store_true', default=False, help='Kmers are already counted and in file expected')
    parser.add_option('--no_cut', dest='no_cut', action='store_true', default=False, help='Coverage model is optimized and cutoff is printed to file expected')
    (options, args) = parser.parse_args()

    if not options.readsf:
        parser.error('Must provide fastq file of reads with -r')
    if not options.k:
        parser.error('Must provide k-mer size with -k')

    if qmers:
        ctsf = '%s.qcts' % os.path.splitext( os.path.split(options.readsf)[1] )[0]
    else:
        ctsf = '%s.cts' % os.path.splitext( os.path.split(options.readsf)[1] )[0]

    if not options.no_count and not options.no_cut:
        #count_kmers_meryl(options.readsf, options.k, ctsf)
        count_kmers_amos(options.readsf, options.k, ctsf)

    if not options.no_cut:
        # model coverage
        if qmers:
            if model_gc:
                model_q_gc_cutoffs(ctsf, 20000)
            else:
                model_q_cutoff(ctsf, 30000)
        else:
            model_cutoff(ctsf)

    if model_gc:
        # run correct C++ code
        os.system('correct -r %s -m %s -a cutoffs.gc.txt -p %d' % (options.readsf, ctsf, options.proc))      

    else:
        cutoff = open('cutoff.txt').readline().rstrip()

        # run correct C++ code
        os.system('correct -r %s -m %s -c %s -p %d' % (options.readsf, ctsf, cutoff, options.proc))

    os.system('cat out.txt?* > out.txt')
    os.system('rm out.txt?*')


############################################################
# model_cutoff
#
# Make a histogram of kmers to give to R to learn the cutoff
############################################################
def model_cutoff(ctsf):
    # make kmer histogram
    cov_max = 0
    for line in open(ctsf):
        cov = int(line.split()[1])
        if cov > cov_max:
            cov_max = cov

    kmer_hist = [0]*cov_max
    for line in open(ctsf):
        cov = int(line.split()[1])
        kmer_hist[cov-1] += 1

    cov_out = open('kmers.hist', 'w')
    for cov in range(0,cov_max):
        if kmer_hist[cov]:
            print >> cov_out, '%d\t%d' % (cov+1,kmer_hist[cov])
    cov_out.close()

    os.system('R --slave << %s/kmer_model.z.errpar.r' % r_dir)


############################################################
# model_q_cutoff
#
# Sample kmers to give to R to learn the cutoff
############################################################
def model_q_cutoff(ctsf, sample):
    # input coverages
    covs = []
    for line in open(ctsf):
        (kmer,cov) = line.split()
        covs.append(cov)

    # sample covs
    if sample >= len(covs):
        rand_covs = covs
    else:
        rand_covs = random.sample(covs, sample)

    # print to file
    out = open('kmers.txt', 'w')
    for rc in rand_covs:
        print >> out, rc
    out.close()

    os.system('R --slave < %s/kmer_model.z.errpar.qualk.r' % r_dir)


############################################################
# model_q_gc_cutoffs
#
# Sample kmers to give to R to learn the cutoff for each
# GC value
############################################################
def model_q_gc_cutoffs(ctsf, sample): 
    print ctsf

    # input coverages
    k = 0
    for line in open(ctsf):
        (kmer,cov) = line.split()
        if k == 0:
            k = len(kmer)
            at_covs = ['']*(k+1)
        else:
            at = count_at(kmer)
            if at_covs[at]:
                at_covs[at].append(cov)
            else:
                at_covs[at] = [cov]

    for at in range(1,k):
        print '%d %d' % (at,len(at_covs[at]))

    # for each AT bin
    at_cutoffs = []
    for at in range(1,k):
        # sample covs
        if sample >= len(at_covs[at]):
            rand_covs = at_covs[at]
        else:
            rand_covs = random.sample(at_covs[at], sample)

        # print to file
        out = open('kmers.txt', 'w')
        for rc in rand_covs:
            print >> out, rc
        out.close()

        os.system('R --slave < %s/kmer_model.z.errpar.qualk.r' % r_dir)

        at_cutoffs.append( open('cutoff.txt').readline().rstrip() )
        if at in [1,k-1]:   # setting extremes to next closests
            at_cutoffs.append( open('cutoff.txt').readline().rstrip() )

        os.system('mv kmers.txt kmers.at%d.txt' % at)
        os.system('mv cutoff.txt cutoff.at%d.txt' % at)

    out = open('cutoffs.gc.txt','w')
    print >> out, '\n'.join(at_cutoffs)
    out.close()
        
    
############################################################
# count_at
#
# Count A's and T's in the given sequence
############################################################
def count_at(seq):
    return len([nt for nt in seq if nt in ['A','T']])


############################################################
# count_kmers_meryl
#
# Count kmers in the reads file using Meryl
############################################################
def count_kmers_meryl(readsf, k, ctsf, proc=1):
    # pull out sequence
    os.system('formats.py --fq2fa %s > %s.fa' % (readsf, ctsf[:-4]))
    
    # count mers
    os.system('meryl -B -C -m %d -s %s.fa -o %s -threads %d' % (k, ctsf[:-4], ctsf[:-4], proc))
    # print histogram of counts
    #os.system('meryl -Dh -s %s > %s.hist' % (kmerf,kmerf))
    # print counts
    os.system('meryl -Dt -n 1 -s %s > %s' % (ctsf[:-4],ctsf))
    

############################################################
# count_kmers_amos
#
# Count kmers in the reads file using AMOS count-kmers
############################################################
def count_kmers_amos(readsf, k, ctsf):
    kmerf = os.path.splitext( os.path.split(readsf)[1] )[0]
    os.system('cat %s | count-kmers -k %d -S > %s' % (readsf, k, ctsf))
    
            
############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
