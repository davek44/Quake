#!/usr/bin/env python
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
from math import log, ceil
import os, random, sys, math, subprocess, pdb
import cov_model

################################################################################
# quake.py
#
# Launch full pipeline to correct errors
################################################################################

quake_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
jellyfish_dir = quake_dir

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)

    # General options
    parser.add_option('-r', dest='readsf', help='Fastq file of reads')
    parser.add_option('-f', dest='reads_listf', help='File containing fastq file names, one per line or two per line for paired end reads.')
    parser.add_option('-k', dest='k', type='int', help='Size of k-mers to correct')
    parser.add_option('-p', dest='proc', type='int', default=4, help='Number of processes [default: %default]')
    parser.add_option('-q', dest='quality_scale', type='int', default=-1, help='Quality value ascii scale, generally 64 or 33. If not specified, it will guess.')


    # Count options
    count_group = OptionGroup(parser, 'K-mer counting')
    count_group.add_option('--no_jelly', dest='no_jelly', action='store_true', default=False, help='Count k-mers using a simpler program than Jellyfish')
    count_group.add_option('--no_count', dest='no_count', action='store_true', default=False, help='Kmers are already counted and in expected file [reads file].qcts or [reads file].cts [default: %default]')    
    count_group.add_option('--int', dest='count_kmers', action='store_true', default=False, help='Count kmers as integers w/o the use of quality values [default: %default]')
    count_group.add_option('--count_only', dest='count_only', action='store_true', default=False, help=SUPPRESS_HELP)
    count_group.add_option('--hash_size', dest='hash_size', type='int', help='Jellyfish hash-size parameter. Quake will estimate using k if not given')
    parser.add_option_group(count_group)

    # Model options
    model_group = OptionGroup(parser, 'Coverage model')
    model_group.add_option('--no_cut', dest='no_cut', action='store_true', default=False, help='Coverage model is optimized and cutoff was printed to expected file cutoff.txt [default: %default]')
    model_group.add_option('--ratio', dest='ratio', type='int', default=200, help='Likelihood ratio to set trusted/untrusted cutoff.  Generally set between 10-1000 with lower numbers suggesting a lower threshold. [default: %default]')
    parser.add_option_group(model_group)

    # Correct options
    correct_group = OptionGroup(parser, 'Correction')
    correct_group.add_option('-l', dest='min_read', type='int', default=30, help='Return only reads corrected and/or trimmed to <min_read> bp')
    correct_group.add_option('-u', dest='out_err', action='store_true', default=False, help='Output error reads even if they can\'t be corrected, maintaing paired end reads')
    correct_group.add_option('-t', dest='trim_par', type='int', default=3, help='Use BWA-like trim parameter <trim_par>')
    correct_group.add_option('--headers', dest='headers', action='store_true', default=False, help='Output only the original read headers without correction messages')
    correct_group.add_option('--log', dest='log', action='store_true', default=False, help='Output a log of all corrections into *.log as "quality position new_nt old_nt"')
    parser.add_option_group(correct_group)

    # help='Model kmer coverage as a function of GC content of kmers [default: %default]'
    parser.add_option('--gc', dest='model_gc', action='store_true', default=False, help=SUPPRESS_HELP)

    (options, args) = parser.parse_args()

    if not options.readsf and not options.reads_listf:
        parser.error('Must provide fastq file of reads with -r or file with list of fastq files of reads with -f')
    if options.quality_scale == -1:
        quality_scale = guess_quality_scale(options.readsf, options.reads_listf)
    else:
        quality_scale = options.quality_scale
    if not options.k:
        parser.error('Must provide k-mer size with -k')
    elif options.k > 20:
        print >> sys.stderr, 'k=%d is a high k-mer value and should be reserved for enormous genomes (see the FAQ on the Quake website). Nevertheless, the program will proceed assuming you have chosen well.' % options.k

    if options.count_kmers:
        cts_suf = 'cts'
    else:
        cts_suf = 'qcts'
    if options.readsf:
        ctsf = '%s.%s' % (os.path.splitext( os.path.split(options.readsf)[1] )[0], cts_suf)
        reads_str = '-r %s' % options.readsf
    else:
        ctsf = '%s.%s' % (os.path.split(options.reads_listf)[1], cts_suf)
        reads_str = '-f %s' % options.reads_listf

    ############################################
    # count kmers
    ############################################
    if not options.no_count and not options.no_cut:
        if options.no_jelly:
            count_kmers(options.readsf, options.reads_listf, options.k, ctsf, quality_scale)
        else:
            jellyfish(options.readsf, options.reads_listf, options.k, ctsf, quality_scale, options.hash_size, options.proc)

        if options.count_only:
            exit(0)

    ############################################
    # model coverage
    ############################################
    if not options.no_cut:
        # clear file
        if os.path.isfile('cutoff.txt'):
            os.remove('cutoff.txt')

        if options.count_kmers:
            cov_model.model_cutoff(ctsf, options.ratio)
        else:
            if options.model_gc:
                cov_model.model_q_gc_cutoffs(ctsf, 15000, options.ratio)
            else:
                # clear file
                if os.path.isfile('kmers.txt'):
                    os.remove('kmers.txt')

                cov_model.model_q_cutoff(ctsf, 30000, options.ratio)

    ############################################
    # correct reads
    ############################################
    correct_options = make_cor_opts(options)
    if options.model_gc:        
        # run correct C++ code
        p = subprocess.Popen('%s/correct %s %s -m %s -a cutoffs.gc.txt -q %d' % (quake_dir, correct_options, reads_str, ctsf, quality_scale), shell=True)
        os.waitpid(p.pid, 0)

    else:
        cutoff = get_cutoff()            

        # run correct C++ code
        p = subprocess.Popen('%s/correct %s %s -m %s -c %s -q %d' % (quake_dir, correct_options, reads_str, ctsf, cutoff, quality_scale), shell=True)
        os.waitpid(p.pid, 0)


################################################################################
# count_kmers
#
# Count kmers in the reads file using my single-threaded
# programs
################################################################################
def count_kmers(readsf, reads_listf, k, ctsf, quality_scale):
    # find files
    fq_files = []
    if readsf:
        fq_files.append(readsf)
    else:
        for line in open(reads_listf):
            for fqf in line.split():
                fq_files.append(fqf)

    # count zipped fastq files
    fq_zipped = [fqf[-3:] == '.gz' for fqf in fq_files]

    # none zipped
    if sum(fq_zipped) == 0:
        if ctsf[-5:] == '.qcts':
            p = subprocess.Popen('cat %s | %s/count-qmers -k %d -q %d > %s' % (' '.join(fq_files), quake_dir, k, quality_scale, ctsf), shell=True)
        else:
            p = subprocess.Popen('cat %s | %s/count-kmers -k %d > %s' % (' '.join(fq_files), quake_dir, k, ctsf), shell=True)

    # all zipped
    elif sum(fq_zipped) == len(fq_zipped):
        if ctsf[-5:] == '.qcts':
            p = subprocess.Popen('gunzip -c %s | %s/count-qmers -k %d -q %d > %s' % (' '.join(fq_files), quake_dir, k, quality_scale, ctsf), shell=True)
        else:
            p = subprocess.Popen('gunzip -c %s | %s/count-kmers -k %d > %s' % (' '.join(fq_files), quake_dir, k, ctsf), shell=True)

    # mixed- boo
    else:
        print >> sys.stderr, 'Some fastq files are zipped, some are not. Please be consistent.'
        exit(1)

    os.waitpid(p.pid, 0)


################################################################################
# get_cutoff
#
# Read the chosen cutoff from disk and interpret it for the user.
################################################################################
def get_cutoff():
    if not os.path.isfile('cutoff.txt'):
        print >> sys.stderr, 'Optimization of distribution likelihood function to choose k-mer cutoff failed. Very likely you have set the value of k too high or not provided adequate coverage (>15x). Inspect the k-mer counts for a clear separation of the error and true k-mer distributions.'
        exit(1)

    cutoff = open('cutoff.txt').readline().rstrip()

    print 'Cutoff: %s' % cutoff

    if float(cutoff) < 1.0:
        print '(This is a low cutoff value indicating that you have low genome coverage. You may consider also correcting your reads by running the "correct" program directly specifying the cutoff as 1.0)'

    return cutoff


################################################################################
# guess_quality_scale
#
# Guess at ascii scale of quality values by examining
# a bunch of reads and looking for quality values < 64,
# in which case we set it to 33.
################################################################################
def guess_quality_scale(readsf, reads_listf):
    reads_to_check = 1000
    if not readsf:
        readsf = open(reads_listf).readline().split()[0]

    fqf = open(readsf)
    reads_checked = 0
    header = fqf.readline()
    while header and reads_checked < reads_to_check:
        seq = fqf.readline()
        mid = fqf.readline()
        qual = fqf.readline().rstrip()
        reads_checked += 1
        for q in qual:
            if ord(q) < 64:
                print 'Guessing quality values are on ascii 33 scale'
                return 33
        header = fqf.readline()

    print 'Guessing quality values are on ascii 64 scale'
    return 64


################################################################################
# jellyfish
#
# Run Jellyfish using integer counts only for now. I take a guess at a good
# hash table size here based on the kmer size chosen.
################################################################################
def jellyfish(readsf, reads_listf, k, ctsf, quality_scale, hash_size, proc):
    # choose parameters
    table_expansion_factor = 20
    if not hash_size:
        # guess at table size
        hash_size = int(table_expansion_factor*math.pow(4,k) / 200)

        # make power of 2
        hash_size = pow(2, ceil(log(hash_size)/log(2)))
    ct_size = 8

    # find files
    fq_files = []
    if readsf:
        output_pre = os.path.splitext(readsf)[0]
        fq_files.append(readsf)
    else:
        output_pre = os.path.splitext(reads_listf)[0]
        for line in open(reads_listf):
            for fqf in line.split():
                fq_files.append(fqf)

    # count zipped fastq files
    fq_zipped = [fqf[-3:] == '.gz' for fqf in fq_files]

    # none zipped
    if sum(fq_zipped) == 0:
        if ctsf[-4:] == 'qcts':
            p = subprocess.Popen('%s/jellyfish count -q --quality-start %d -c %d -o %s.db -m %d -t %d -s %d --both-strands %s' % (jellyfish_dir, quality_scale, ct_size, output_pre, k, proc, hash_size, ' '.join(fq_files)), shell=True)
        else:
            p = subprocess.Popen('%s/jellyfish count -c %d -o %s.db -m %d -t %d -s %d --both-strands %s' % (jellyfish_dir, ct_size, output_pre, k, proc, hash_size, ' '.join(fq_files)), shell=True)

    # all zipped
    elif sum(fq_zipped) == len(fq_zipped):
        if ctsf[-4:] == 'qcts':
            p = subprocess.Popen('gunzip -c %s | %s/jellyfish count -q --quality-start %d -c %d -o %s.db -m %d -t %d -s %d --both-strands /dev/fd/0' % (' '.join(fq_files), jellyfish_dir, quality_scale, ct_size, output_pre, k, proc, hash_size), shell=True)
        else:
            p = subprocess.Popen('gunzip -c %s | %s/jellyfish count -c %d -o %s.db -m %d -t %d -s %d --both-strands /dev/fd/0' % (' '.join(fq_files), jellyfish_dir, ct_size, output_pre, k, proc, hash_size), shell=True)

    # mixed- boo
    else:
        print >> sys.stderr, 'Some fastq files are zipped, some are not. Please be consistent.'
        exit(1)

    os.waitpid(p.pid, 0)

    # merge
    max_db = 0
    while os.path.isfile('%s.db_%d' % (output_pre,max_db)):
        max_db += 1
    if max_db == 1:
        os.rename('%s.db_0' % output_pre, '%s.dbm' % output_pre)
    else:
        if ctsf[-4:] == 'qcts':
            # expand hash table
            hash_size *= max_db

            # merge db
            p = subprocess.Popen('%s/jellyfish qmerge -s %d -m %d -o %s.dbm %s' % (jellyfish_dir, hash_size, k, output_pre, ' '.join(['%s.db_%d' % (output_pre,i) for i in range(max_db)])), shell=True)
            os.waitpid(p.pid, 0)
            
            # rename file
            os.rename('%s.dbm_0' % output_pre, '%s.dbm' % output_pre)

        else:
            # merge db
            p = subprocess.Popen('%s/jellyfish merge -o %s.dbm %s' % (jellyfish_dir, output_pre, ' '.join(['%s.db_%d' % (output_pre,i) for i in range(max_db)])), shell=True)
            os.waitpid(p.pid, 0)

    # produce actual counts
    if ctsf[-4:] == 'qcts':
        p = subprocess.Popen('%s/jellyfish qdump -c %s.dbm > %s' % (jellyfish_dir, output_pre, ctsf), shell=True)
    else:
        p = subprocess.Popen('%s/jellyfish dump -c %s.dbm > %s' % (jellyfish_dir, output_pre, ctsf), shell=True)
    os.waitpid(p.pid, 0)


################################################################################
# make_cor_opts
#
# Interpret the python command line options to pass along to the correct program
################################################################################
def make_cor_opts(options):
    correct_options = '-k %d -p %d -l %d -t %d' % (options.k, options.proc, options.min_read, options.trim_par)
    if options.out_err:
        correct_options += ' -u'
    if options.log:
        correct_options += ' --log'
    if options.headers:
        correct_options += ' --headers'
    return correct_options

############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
