#!/fs/sz-user-supported/Linux-i686/bin/python2.5
import math, random, os, dna

############################################################
# test.py
#
# Test my error correction program with simulated reads
#
# Notes:
# -"Correction should have been made" can appear wrongly
#   if a read has two errors and one of them lies in a 
#   repetitive region such that some of its surrounding
#   kmers are trusted by chance.
############################################################

#errorprofile_file = '/nfshomes/dakelley/research/hpylori/data/reads/HPKX_1039_AG0C1.1.fq'
#errorprofile_file = '/nfshomes/dakelley/research/error_correction/data/ecoli/SRR001665_1.fastq'
errorprofile_file = '/fs/szattic-asmg4/Bees/Bombus_impatiens/s_3_1_sequence.txt'
read_len = 124
illumina_qual = True

like_t = .00001
like_spread_t = .2
trimq = 3
k = 19

############################################################
# main
############################################################
def main():
    print 'REMEBER YOU HAVE TO CHANGE THE OUTPUT FROM [-,.]\'S TO PRINTING AMBIGUOUS CORRECTED READS'

    # get error profiles
    num_reads = 50000
    err_profs = get_error_profiles(num_reads, False)

    # make genome
    genome = ''
    genome_size = 100000
    for i in range(genome_size):
        genome += random.choice(['A','C','G','T'])
    gf = open('genome.fa','w')
    gf.write('>genome\n%s\n' % genome)
    gf.close()
    #genome = dna.fasta2dict('genome.fa').values()[0]

    # simulate reads
    err_reads = simulate_error_reads(genome, err_profs)

    # find trusted kmers
    trusted_kmers(genome)

    # run ./correct
    #os.system('time ./build_bithash -m genome.cts -c 99 -o genome.bh.out')
    #os.system('time ./correct -r err_reads.fq -m genome.bh.out -p 1 2> out.txt')

    if illumina_qual:
        os.system('time ./correct -r err_reads.fq -m genome.cts -c 99 -I -p 1 2> out.txt')
    else:
        os.system('time ./correct -r err_reads.fq -m genome.cts -c 99 -p 1 2> out.txt')

    # compare corrected reads
    compare_corrections(err_reads, genome)


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
# simulate_error_reads
#
# Randomly simulate a bunch of reads using the error
# profiles given.  Only bother keeping those with errors
# and same information about the errors in order to check.
############################################################
def simulate_error_reads(genome, err_profs):
    reads_out = open('err_reads.fq','w')
    err_reads = {}
    for i in range(len(err_profs)):
        # simulate read
        start = random.randint(0, len(genome)-read_len)
        seq = genome[start:start+read_len]

        # mutate
        mseq = list(seq)
        for j in range(read_len):
            if illumina_qual:
                err_prob = math.pow(10.0,-(ord(err_profs[i][j])-64.0)/10.0)
            else:
                err_prob = math.pow(10.0,-(ord(err_profs[i][j])-33.0)/10.0)
            if err_prob > .75:
                mseq[j] = 'N'
            elif random.random() < err_prob:
                if mseq[j] == 'A':
                    mseq[j] = random.choice(['C','G','T'])
                elif mseq[j] == 'C':
                    mseq[j] = random.choice(['A','G','T'])
                elif mseq[j] == 'G':
                    mseq[j] = random.choice(['A','C','T'])
                elif mseq[j] == 'T':
                    mseq[j] = random.choice(['A','C','G'])
                    
        mseq = ''.join(mseq)

        if mseq != seq or True:
            header = '@'+str(i)
            err_reads[header] = {'orig':seq, 'err':mseq, 'qual':err_profs[i]}
            print >> reads_out, '%s\n%s\n+\n%s' % (header,mseq,err_profs[i])

    reads_out.close()
    return err_reads
        

############################################################
# trusted_kmers
#
# Print out a meryl-like out .cts file with all kmers
# in the genome.
############################################################
def trusted_kmers(genome):
    tkout = open('genome.cts','w')
    for i in range(len(genome)-k+1):
        print >> tkout, '%s\t100.0' % genome[i:i+k]
    tkout.close()

############################################################
# compare_corrections
#
# Open the corrected read output file and compare to
# the actual errors
############################################################
def compare_corrections(err_reads, genome):
    right_correction = 0
    right_wrong_correction = 0
    right_ambiguous = 0
    right_threshold = 0
    right_trim = 0
    correctable_abstained = 0
    
    for line in open('out.txt'):
        a = line.split()

        # corrected read
        if not err_reads.has_key(a[0]):
            # already seen as ambiguous or badnt forfeit
            continue
        er = err_reads[a[0]]

        
        if len(a) == 3:
            # check trimming
            if len(a[2]) > 1 and len(a[2]) != len(er['err']):
               if check_trim(a[2], er):
                   right_trim += 1
               else:
                   print 'Bad trim: %s %s %s %s %s' % (a[0],er['orig'],er['err'],er['qual'],a[2])

            if er['orig'] == a[2]:
                # proper correction
                rlike = likelihood(er['orig'], er['err'], er['qual'])
                if rlike > like_t:
                    right_correction += 1
                else:
                    print 'correction should have failed due to low likelihood: %s %s %s %s %f' % (a[0], er['orig'], er['err'], er['qual'], rlike)

            elif a[2] == '-':
                # no correction found
                rlike = likelihood(er['orig'], er['err'], er['qual'])
                if rlike > like_t and er['err'].find('N') == -1:
                    print 'correction should have been made: %s %s %s %s %f' % (a[0],er['orig'],er['err'],er['qual'],rlike)
                else:
                    right_threshold += 1

            elif a[2] == '.':
                # read too crappy to bother
                rlike = likelihood(er['orig'], er['err'], er['qual'])
                if rlike > like_t:
                    correctable_abstained += 1

            else:
                if len(a[2]) == len(er['orig']):
                    # no trim so corrected improperly
                    if(not check_trust(a[2], genome)):
                        print 'improper correction should have failed due to lack of trust %s %s %s' % (a[0], a[1], a[2])

                    clike = likelihood(a[2], a[1], er['qual'])
                    rlike = likelihood(er['orig'], er['err'], er['qual'])

                    if clike*like_spread_t > rlike:
                        # proper action
                        right_wrong_correction += 1

                    elif clike > rlike and clike*like_spread_t < rlike:
                        # action is best, but ambiguous
                        print 'improper correction should have failed due to ambiguity: %s %s %s %s %f %f' % (a[0], er['orig'], er['err'], er['qual'], rlike, clike)

                    elif clike < like_t:
                        print 'improper correction should have failed due to low likelihood: %s %s %s %s %f' % (a[0], er['orig'], er['err'], er['qual'], clike)

                else:
                    # trimmed
                    if a[2] == er['orig'][:len(a[2])]:
                        # proper correction
                        rlike = likelihood(er['orig'][:len(a[2])], er['err'][:len(a[2])], er['qual'][:len(a[2])])
                        if rlike > like_t:
                            right_correction += 1
                        else:
                            print 'trimmed correction should have failed due to low likelihood: %s %s %s %s %f' % (a[0], er['orig'], er['err'], er['qual'], rlike)

                    else:
                        # corrected improperly
                        if(not check_trust(a[2], genome)):
                            print 'trimmed correction should have failed due to lack of trust %s %s %s' % (a[0], a[1], a[2])

                        clike = likelihood(a[2], a[1][:len(a[2])], er['qual'][:len(a[2])])
                        rlike = likelihood(er['orig'][:len(a[2])], er['err'][:len(a[2])], er['qual'][:len(a[2])])

                        if clike*like_spread_t > rlike:
                            # proper action
                            right_wrong_correction += 1

                        elif clike > rlike and clike*like_spread_t < rlike:
                            # action is best, but ambiguous
                            print 'trimmed correction should have failed due to ambiguity: %s %s %s %s %f %f' % (a[0], er['orig'], er['err'], er['qual'], rlike, clike)

                        elif clike < like_t:
                            print 'trimmed correction should have failed due to low likelihood: %s %s %s %s %f' % (a[0], er['orig'], er['err'], er['qual'], clike)


        elif len(a) == 4:
            # ambiguous read
            r1like = likelihood(a[2], a[1], er['qual'])
            r2like = likelihood(a[3], a[1], er['qual'])
            
            if r1like*like_spread_t > r2like:
                print 'correction should not have failed due to ambiguity: %s %s %s %s %f %f' % (a[0], a[2], a[3], er['qual'], r1like, r2like)
            else:
                right_ambiguous += 1

        del err_reads[a[0]]


    print 'Right correction: %d' % right_correction
    print 'Optimal but incorrect correction: %d' % right_wrong_correction
    print 'Optimal ambiguity restraint in not correcting: %d' % right_ambiguous
    print 'Optimal threshold restraint in not correcting: %d' % right_threshold
    print 'Gave up, but read was correctable over threshold: %d' % correctable_abstained
    print 'Right trim: %d' % right_trim


############################################################
# likelihood
#
# Return the likelihood of observing 'obs_read' given 
# the actual read and quality values.
# 
# Note that if the obs_q = .25, the likelihood ratio is 1.0
############################################################
def likelihood(actual_read, obs_read, quals):
    like = 1.0
    for i in range(len(actual_read)):
        if actual_read[i] != obs_read[i]:
            if illumina_qual:
                obs_q = max(.25, 1.0 - math.pow(10.0,-(ord(quals[i])-64.0)/10.0))
            else:
                obs_q = max(.25, 1.0 - math.pow(10.0,-(ord(quals[i])-33.0)/10.0))
            change_q = (1.0 - obs_q) / 3.0
            
            like *= change_q / obs_q
    return like


############################################################
# check_trust
#
# Check for the existence of the kmers in seq in the 
# genome.
############################################################
def check_trust(seq, genome):
    for i in range(len(seq)-k+1):
        if genome.find(seq[i:i+k]) == -1 and genome.find(dna.rc(seq[i:i+k])) == -1:
            return False
    return True

############################################################
# check_trim
#
# Check that trimming was properly done.  The true trim
# is the BWA trim point minus any correctable errors.  So
# vaid trim points are the BWA point as well as any number
# of correctable errors
############################################################
def check_trim(seq, error_read):
    # determine BWA trim point
    BWA_trim = len(error_read['qual'])
    max_score = 0
    for i in range(len(error_read['qual'])):
        score = 0
        for q in range(i,len(error_read['qual'])):
            if illumina_qual:
                score += (trimq - (ord(error_read['qual'][q])-64))
            else:
                score += (trimq - (ord(error_read['qual'][q])-33))
        #print '%d: %d' % (i,score)
        if score >= max_score:
            max_score = score
            BWA_trim = i    

    #print len(seq)
    #print BWA_trim
    trim_pt = len(seq)
    if trim_pt == BWA_trim:
        return True

    elif error_read['orig'][trim_pt+k-1] != error_read['err'][trim_pt+k-1]:
        return True

    else:
        print 'Quaec trim: %d' % trim_pt
        print 'BWA trim: %d' % BWA_trim
        print 'Error trims: ',
        for i in range(BWA_trim):
            if error_read['orig'][i] != error_read['err'][i]:
                print '%d,' % (i-k+1),
        print ''
        
        return False


############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
