#!/fs/sz-user-supported/Linux-i686/bin/python2.5
from optparse import OptionParser
import os

############################################################
# shrec.py
#
# Functions to run Shrec's error correction tool
# as well as parse the results.
############################################################

############################################################
# main
############################################################
def main():
    parser = OptionParser()
    parser.add_option('-r', dest='readsf', default='test_reads.fq', help='Reads fastq file')
    parser.add_option('-c', dest='cutoff', default=4, type='int', help='Kmer cutoff')

    (options, args) = parser.parse_args()

############################################################
# run
#
# Run Shrec.  Reads file should be fasta.
############################################################
def run(readsf, genome_size):
    os.system('java -Xmx3072m Shrec %s %d correct.fa discard.fa strictness=5.0' % (readsf,genome_size))


############################################################
# accuracy
#
# Count the following:
# 1. Error reads properly corrected.
# 2. Error reads improperly corrected.
# 3. Error reads thrown away.
# 4. Error reads ignored and kept.
# 5. OK reads improperly corrected.
# 6. OK reads thrown away.
############################################################
def accuracy(infof):
    # get actual error reads
    error_reads = {}
    for line in open(infof):
        a = line.split()
        error_reads[a[0]] = a[1]

    # stats
    err_proper = 0
    err_improper = 0
    err_toss = 0
    err_kept = 0
    ok_improper = 0
    ok_toss = 0

    # parse discarded reads
    for line in open('discard.fa'):
        if line[0] == '>':
            read_id = line[1:].rstrip()
            if error_reads.has_key(read_id):
                err_toss += 1
            else:
                ok_toss += 1

    # parse corrected reads
    for line in open('correct.fa'):
        if line[0] == '>':
            if line.find('correct') == -1:
                read_id = line[1:].rstrip()
                correction = False
            else:
                read_id = line[1:].split()[0]
                correction = True

        else:
            seq = line.rstrip()

            if correction:
                # read was corrected
                if error_reads.has_key(read_id):
                    # read had an error
                    if seq == error_reads[read_id]:                        
                        err_proper += 1
                    else:
                        err_improper += 1
                else:
                    # but had no error
                    ok_improper += 1

            else:
                # read was kept
                if error_reads.has_key(read_id):
                    # but had an error
                    err_kept += 1
                
    print 'Error reads properly corrected\t%d' % err_proper
    print 'Error reads improperly corrected\t%d' % err_improper
    print 'Error reads thrown away\t%d' % err_toss
    print 'Error reads ignored and kept\t%d' % err_kept
    print 'OK reads improperly corrected\t%d' % ok_improper
    print 'OK reads thrown away\t%d' % ok_toss

    print 'Bad reads: %d' % (err_improper+err_kept+ok_improper)
        
