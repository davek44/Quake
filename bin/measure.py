#!/fs/sz-user-supported/Linux-i686/bin/python2.5
from optparse import OptionParser
import random, math, os
import dna, formats
import soap, shrec

############################################################
# measure.py
#
#
############################################################

refgenome_file = '/fs/szdata/ncbi/genomes/Bacteria/Helicobacter_pylori_G27/NC_011333.fna' # 1.65 Mb
#errorprofile_file = '/nfshomes/dakelley/research/hpylori/data/reads/HPKX_1039_AG0C1.1.fq'
errorprofile_file = '/nfshomes/dakelley/research/error_correction/data/turkey.fq'

proc = 4
num_reads = 900000 # ~40x
kmer = 15
rand_ntnt = False

############################################################
# main
############################################################
def main():
    parser = OptionParser()
    parser.add_option('-a', dest='accuracy', default=False, action='store_true', help='Just compute accuracy of existing read corrections')
    parser.add_option('--me', dest='me', default=False, action='store_true', help='Do my corrector')
    parser.add_option('--soap', dest='soap', default=False, action='store_true', help='Do SOAP')
    parser.add_option('--shrec', dest='shrec', default=False, action='store_true', help='Do SHREC')
    parser.add_option('--all', dest='all', default=False, action='store_true', help='Do all')
    (options, args) = parser.parse_args()

    if not (options.me or options.soap or options.shrec or options.all):
        print 'Not going to correct any reads- no program provided'

    ########################################
    # generate test data
    ########################################
    if not options.accuracy:
        # load error profiles
        err_profs = get_error_profiles(num_reads)

        # load genome
        genome_dict = dna.fasta2dict(refgenome_file)

        # make nt->nt transition matrix
        ntnt_matrix = make_ntntmatrix(rand_ntnt)

        # simulate reads
        simulate_error_reads(genome_dict.values()[0], err_profs, ntnt_matrix)
    
    ########################################
    # my correct
    ########################################
    if options.all or options.me:
        if not options.accuracy:
            # correct reads
            os.system('time correct.py -r test_reads.fq -k %d -p %d' % (kmer,proc))

        # compute accuracy
        my_accuracy()

    ########################################
    # SOAP de novo
    ########################################
    if options.all or options.soap:
        if not options.accuracy:
            # steal my cutoff
            cutoff = int(open('cutoff.txt').readline())

            # correct reads
            soap.run('test_reads.fq', kmer, cutoff)

        # compute accuracy
        soap.accuracy('test_reads.info', 'test_reads.fq.corr')

    ########################################
    # EULER
    ########################################    
    #formats.fq2fa('test_reads.fq', 'test_reads.fa')
    #euler.run('test_reads.fa', kmer)

    ########################################
    # SHREC
    ########################################
    if options.all or options.shrec:
        if not options.accuracy:
            # convert to fasta
            formats.fq2fa('test_reads.fq', 'test_reads.fa')

            # correct reads
            shrec.run('test_reads.fa', 1650000)
        
        # compute accuracy
        shrec.accuracy('test_reads.info')


############################################################
# get_error_profiles
#
# Read in actual error profiles from a fastq file.  If rand
# is true, read in all error profiles and randomly sample
# the number desired.  (May have memory issues for a huge
# file.)
############################################################
def get_error_profiles(num, rand=False):
    err_profs = []
    fqf = open(errorprofile_file)
    line = fqf.readline()
    while line:
        line = fqf.readline()
        line = fqf.readline()
        line = fqf.readline()
        err_profs.append(line.rstrip())
        if not rand:
            num = num - 1
            if num <= 0:
                break
        line = fqf.readline()

    if rand:
        return random.sample(err_profs, num)
    else:
        return err_profs

############################################################
# make_ntntmatrix
#
# Make the nt->nt transition matrix for errors
############################################################
def make_ntntmatrix(rand_ntnt):
    nts = ['A','C','G','T']

    # choose randomly, or set to 1/3 if !=
    ntnt_matrix = []
    for nt in range(4):
        ntnt_matrix.append([0.0]*4)
        for nt2 in range(4):
            if nt != nt2:
                if rand_ntnt:
                    ntnt_matrix[nt][nt2] = random.random()
                else:
                    ntnt_matrix[nt][nt2] = 1.0/3.0

        # normalize
        if rand_ntnt:
            ntsum = sum(ntnt_matrix[nt])
            for nt2 in range(4):
                ntnt_matrix[nt][nt2] /= ntsum

    # print to file
    ntnt_out = open('ntnt.out','w')
    print >> ntnt_out, '\t'+'\t'.join(nts)
    for nt in range(4):
        print >> ntnt_out, '%s\t%s' % (nts[nt],'\t'.join([str(x) for x in ntnt_matrix[nt]]))
    ntnt_out.close()

    return ntnt_matrix

############################################################
# simulate_error_reads
#
# Randomly simulate a bunch of reads using the error
# profiles given.  Only bother keeping those with errors
# and save information about the errors in order to check.
############################################################
def simulate_error_reads(genome, err_profs, ntnt_matrix):
    readlen = len(err_profs[0])
    nts = {'A':0, 'C':1, 'G':2, 'T':3}

    reads_out = open('test_reads.fq','w')
    info_out = open('test_reads.info','w')

    for i in range(len(err_profs)):
        # simulate read
        start = random.randint(0, len(genome)-readlen)
        seq = genome[start:start+readlen]

        # reverse complement
        if random.choice([True,False]):
            seq = dna.rc(seq)

        # mutate
        mseq = list(seq)
        for j in range(readlen):
            err_prob = math.pow(10.0,-(ord(err_profs[i][j])-33.0)/10.0)
            if err_prob == 1.0:
                mseq[j] = 'N'
            elif random.random() < err_prob:
                mseq[j] = myrand_choice(['A','C','G','T'],ntnt_matrix[nts[mseq[j]]])                    
        mseq = ''.join(mseq)        

        # print to fastq
        header = '@read_'+str(i)
        print >> reads_out, '%s\n%s\n+\n%s' % (header,mseq,err_profs[i])

        if seq != mseq:
            print >> info_out, '%s\t%s' % (header,seq)

    reads_out.close()
    info_out.close()


############################################################
# myrand_choice
#
# Modeled after random.choice to choose one of the values in
# the list 'vals', but now the choice is definited by the
# discrete distribution given in 'p'.
############################################################
def myrand_choice(vals, p):
    r = random.random()
    psum = 0.0
    for i in range(len(p)):
        psum += p[i]
        if psum > r:
            return vals[i]

    print 'Probabilities did not sum to 1'
    print p
    return ''


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
def my_accuracy():
    # get actual error reads
    error_reads = {}
    for line in open('test_reads.info'):
        a = line.split()
        error_reads[a[0]] = a[1]

    # stats
    err_proper = 0
    err_improper = 0
    err_improperf = open('err_improper.txt','w')
    err_toss = 0
    err_kept = 0
    err_keptf = open('err_kept.txt','w')
    ok_improper = 0
    ok_improperf = open('ok_improper.txt','w')
    ok_toss = 0
    ok_tossf = open('ok_toss.txt','w')

    # get corrected reads
    for line in open('out.txt'):
        (header,seq,cseq) = line.split()

        if len(cseq) > 1:
            # read was corrected
            if not error_reads.has_key(header):
                ok_improper += 1
                print >> ok_improperf, header
            else:
                if cseq == error_reads[header]:
                    err_proper += 1
                else:
                    err_improper += 1
                    print >> err_improperf, header
                    
        else:
            # read was thrown away
            if not error_reads.has_key(header):
                ok_toss += 1
                print >> ok_tossf, header
            else:
                err_toss += 1

        if error_reads.has_key(header):
            del error_reads[header]

    err_kept += len(error_reads)
    for header in error_reads:
        print >> err_keptf, header

    print 'Error reads properly corrected\t%d' % err_proper
    print 'Error reads improperly corrected\t%d' % err_improper
    print 'Error reads thrown away\t%d' % err_toss
    print 'Error reads ignored and kept\t%d' % err_kept
    print 'OK reads improperly corrected\t%d' % ok_improper
    print 'OK reads thrown away\t%d' % ok_toss
    #print 'Bad reads: %d' % (err_improper+err_kept+ok_improper)


############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
