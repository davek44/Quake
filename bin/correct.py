#!/fs/sz-user-supported/Linux-i686/bin/python2.5
from optparse import OptionParser
import os, random

############################################################
# correct.py
#
# Launch pipeline to correct errors in sequencing reads.
############################################################

r_dir = '/nfshomes/dakelley/research/error_correction/bin'

############################################################
# main
############################################################
def main():
    parser = OptionParser()
    
    parser.add_option('-r', dest='readsf', help='Fastq file of reads')
    parser.add_option('-k', dest='k', type='int', help='Size of k-mers to correct')
    parser.add_option('-p', dest='proc', type='int', help='Number of processes')
    parser.add_option('-l', dest='read_len', type='int', help='Read length')

    (options, args) = parser.parse_args()

    if not options.readsf:
        parser.error('Must provide fastq file of reads with -r')
    if not options.k:
        parser.error('Must provide k-mer size with -k')

    ctsf = '%s.cts' % os.path.splitext( os.path.split(options.readsf)[1] )[0]
    #count_kmers_meryl(options.readsf, options.k, ctsf)
    count_kmers_amos(options.readsf, options.k, ctsf)

    # sample kmer coverages
    covs = []
    for line in open(ctsf):
        covs.append(int(line.split()[1]))

    cov_out = open('kmers.txt','w')
    print >> cov_out, '\n'.join([str(x) for x in random.sample(covs,10000)])
    cov_out.close()

    # model coverage
    os.system('R --slave < %s/kmer_model.r > cutoff.txt' % r_dir)
    cutoff = int(open('cutoff.txt').readline())

    # run correct C++ code
    os.system('correct -r %s -m %s -c %d -t %d -l %d' % (options.readsf, ctsf, cutoff, options.proc, options.read_len))

    os.system('cat out.txt?* > out.txt')
    os.system('rm out.txt?*')


############################################################
# count_kmers_meryl
#
# Count kmers in the reads file using Meryl
############################################################
def count_kmers_meryl(readsf, k, ctsf):

    # pull out sequence
    os.system('formats.py --fq2fa %s > %s.fa' % (readsf, ctsf[:-4]))
    
    # count mers
    os.system('meryl -B -C -m %d -s %s.fa -o %s' % (k, ctsf[:-4], ctsf[:-4]))
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
