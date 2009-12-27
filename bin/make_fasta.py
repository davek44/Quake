#!/fs/sz-user-supported/Linux-i686/bin/python2.5
import pdb

############################################################
# make_fasta.py
#
#
############################################################


############################################################
# main
############################################################
def main():
    corrected = {}
    for line in open('out.txt'):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()
        corrected[a[0]] = a[2]
        
    fqf = open('/nfshomes/dakelley/research/hpylori/data/reads/HPKX_1039_AG0C2.2.fq')
    header = fqf.readline().rstrip()
    while header:
        seq = fqf.readline().rstrip()
        skip = fqf.readline().rstrip()
        qual = fqf.readline().rstrip()

        if corrected.has_key(header):
            if corrected[header] not in ['-','.']:
                print header
                print corrected[header]
                print '+'
                print qual

        else:
            print header
            print seq
            print '+'
            print qual

        header = fqf.readline().rstrip()
        

############################################################
# __main__
############################################################
if __name__ == '__main__':
    #pdb.runcall(main)
    main()
