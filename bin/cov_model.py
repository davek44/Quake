#!/fs/sz-user-supported/Linux-i686/bin/python2.5
from optparse import OptionParser
import os, random

############################################################
# cov_model.py
#
# Given a file of kmer counts, reports the cutoff to use
# to separate trusted/untrusted kmers.
############################################################

r_dir = '/nfshomes/dakelley/research/error_correction/bin'

############################################################
# main
############################################################
def main():
    usage = 'usage: %prog [options] <counts file>'
    parser = OptionParser(usage)
    parser.add_option('--int', dest='counted_kmers', action='store_true', default=False, help='Kmers were counted as integers w/o the use of quality values [default: %default]')
    parser.add_option('--gc', dest='model_gc', action='store_true', default=False, help='Model kmer coverage as a function of GC content of kmers [default: %default]')
    parser.add_option('--ratio', dest='ratio', type='int', default=1000, help='Likelihood ratio to set trusted/untrusted cutoff [default: %default]')
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide kmers counts file')
    else:
        ctsf = args[0]

    if options.counted_kmers:
        model_cutoff(ctsf, options.ratio)
        print 'Cutoff: %s' % open('cutoff.txt').readline().rstrip()
        
    else:
        if options.model_gc:
            model_q_gc_cutoffs(ctsf, 10000, options.ratio)
        else:
            model_q_cutoff(ctsf, 25000, options.ratio)
            print 'Cutoff: %s' % open('cutoff.txt').readline().rstrip()


############################################################
# model_cutoff
#
# Make a histogram of kmers to give to R to learn the cutoff
############################################################
def model_cutoff(ctsf, ratio):
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

    os.system('R --slave --args %d < %s/cov_model.r 2> r.log' % (ratio,r_dir))


############################################################
# model_q_cutoff
#
# Sample kmers to give to R to learn the cutoff
############################################################
def model_q_cutoff(ctsf, sample, ratio):
    # count number of kmer coverages
    num_covs = 0
    for line in open(ctsf):
        num_covs += 1

    # choose random kmer coverages
    if sample >= num_covs:
        rand_covs = range(num_covs)
    else:
        rand_covs = random.sample(xrange(num_covs), sample)
    rand_covs.sort()

    # print to file
    out = open('kmers.txt', 'w')
    kmer_i = 0
    rand_i = 0
    for line in open(ctsf):
        if kmer_i == rand_covs[rand_i]:
            print >> out, line.split()[1]
            rand_i += 1
            if rand_i >= sample:
                break
        kmer_i += 1
    out.close()

    os.system('R --slave --args %d < %s/cov_model_qmer.r 2> r.log' % (ratio,r_dir))


############################################################
# model_q_gc_cutoffs
#
# Sample kmers to give to R to learn the cutoff for each
# GC value
############################################################
def model_q_gc_cutoffs(ctsf, sample, ratio):
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

        os.system('R --slave --args %d < %s/cov_model_qmer.r 2> r%d.log' % (ratio,r_dir,at))

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
# __main__
############################################################
if __name__ == '__main__':
    main()
