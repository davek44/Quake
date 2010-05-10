#!/fs/sz-user-supported/Linux-i686/bin/python2.5
from optparse import OptionParser
import os, random
import cov_model

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
    parser.add_option('-p', dest='proc', type='int', default=4, help='Number of processes')
    parser.add_option('-I', dest='illumina_qual', action='store_true', default=False, help='Interpret quality values as Illumina base 64 (as opposed to 33)')
    parser.add_option('--no_count', dest='no_count', action='store_true', default=False, help='Kmers are already counted and in file expected')
    parser.add_option('--no_cut', dest='no_cut', action='store_true', default=False, help='Coverage model is optimized and cutoff is printed to file expected')
    parser.add_option('--int', dest='counted_qmers', action='store_false', default=True, help='Kmers were counted as integers w/o the use of quality values')
    parser.add_option('--gc', dest='model_gc', action='store_true', default=False, help='Model kmer coverage as a function of GC content of kmers')

    (options, args) = parser.parse_args()

    if not options.readsf:
        parser.error('Must provide fastq file of reads with -r')
    if not options.k:
        parser.error('Must provide k-mer size with -k')

    if options.counted_qmers:
        ctsf = '%s.qcts' % os.path.splitext( os.path.split(options.readsf)[1] )[0]
    else:
        ctsf = '%s.cts' % os.path.splitext( os.path.split(options.readsf)[1] )[0]

    if options.illumina_qual:
        options.illumina_qual = '-I'
    else:
        options.illumina_qual = ''

    if not options.no_count and not options.no_cut:
        #count_kmers_meryl(options.readsf, options.k, ctsf)
        count_kmers_amos(options.readsf, options.k, ctsf, options.illumina_qual)

    if not options.no_cut:
        # model coverage
        if options.counted_qmers:
            if options.model_gc:
                cov_model.model_q_gc_cutoffs(ctsf, 10000)
            else:
                cov_model.model_q_cutoff(ctsf, 25000)
        else:
            cov_model.model_cutoff(ctsf)

    if options.model_gc:
        # run correct C++ code
        os.system('correct -r %s -m %s -a cutoffs.gc.txt -p %d %s' % (options.readsf, ctsf, options.proc, options.illumina_qual))

    else:
        cutoff = open('cutoff.txt').readline().rstrip()

        # run correct C++ code
        os.system('correct -r %s -m %s -c %s -p %d %s' % (options.readsf, ctsf, cutoff, options.proc, options.illumina_qual))


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
def count_kmers_amos(readsf, k, ctsf, iq):
    if ctsf[-4:] == 'qcts':
        os.system('cat %s | count-qmers -k %d %s > %s' % (readsf, k, iq, ctsf))
    else:
        os.system('cat %s | count-kmers -k %d %s > %s' % (readsf, k, iq, ctsf))
    
            
############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
