#!/fs/sz-user-supported/Linux-i686/bin/python2.5
from optparse import OptionParser
import os

############################################################
# soap.py
#
# Functions to run SOAP de novo's error correction tool
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
# Run SOAP de novo error correction
############################################################
def run(readsf, kmer, cutoff=''):
    open('soap_in.txt','w').write('%s\n' % readsf)
    os.system('time KmerFreq -i soap_in.txt -o soap_out -s %d' % kmer)
    if cutoff:
        os.system('time Corrector -i soap_in.txt -r soap_out.freq -s %d -k %d -e %d' % (kmer,cutoff,cutoff))
    else:
        os.system('time Corrector -i soap_in.txt -r soap_out.freq -s %d' % kmer)


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
def accuracy(infof, corf):
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

    # get corrected reads
    corf = open(corf)
    header = corf.readline()[1:]
    exp_read_id = 0
    while header:
        read_id = '@'+header.split()[0]
        seq = corf.readline().rstrip()

        skip_read_id = '@read_%d' % exp_read_id
        while read_id != skip_read_id:
            # skip read was tossed
            if error_reads.has_key(skip_read_id):
                err_toss += 1
            else:
                ok_toss += 1

            exp_read_id += 1
            skip_read_id = '@read_%d' % exp_read_id
            

        if header.find('no error found') != -1:
            # read was allowed
            if error_reads.has_key(read_id):
                err_kept += 1
        else:
            # read was corrected
            if not error_reads.has_key(read_id):
                ok_improper += 1
            else:
                if seq == error_reads[read_id]:
                    err_proper += 1
                else:
                    err_improper += 1                

        exp_read_id += 1
        header = corf.readline()[1:]

    print 'Error reads properly corrected\t%d' % err_proper
    print 'Error reads improperly corrected\t%d' % err_improper
    print 'Error reads thrown away\t%d' % err_toss
    print 'Error reads ignored and kept\t%d' % err_kept
    print 'OK reads improperly corrected\t%d' % ok_improper
    print 'OK reads thrown away\t%d' % ok_toss

    print 'Bad reads: %d' % (err_improper+err_kept+ok_improper)
        
