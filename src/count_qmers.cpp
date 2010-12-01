#include  <string>
#include  <string.h>
#include  <vector>
#include  <iostream>
#include  <fstream>
#include  <math.h>
#include  "count.h"
#include "qmer_hash.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// count_qmers
//
// Count q-mers using my qmer_hash space efficient hash table, which is
// especially helpful for Hadoop clusters.
//
// To determine the size of the table, the user can enter the size directly
// or it will guess based on k.  However, the size will be reduced if a memory
// limit is established.
////////////////////////////////////////////////////////////////////////////////

static void CountMers (const string & s, const string & q, qmer_hash & mer_table);

////////////////////////////////////////////////////////////////////////////////
// Additional options
////////////////////////////////////////////////////////////////////////////////
int quality_scale = -1;
unsigned long long table_size = 0;
unsigned int count_max = 500;

////////////////////////////////////////////////////////////////////////////////
// Usage
////////////////////////////////////////////////////////////////////////////////
static void Usage (char * command)
{
  fprintf(stderr, 
	  "\n"
	  ".USAGE.\n"
	  "  count_qmers [-f fastq] [-k len] [options]\n"
	  "\n"
	  ".DESCRIPTION.\n"
	  "  Count kmers in a fastq file. Output is to stdout in simple nmer"
	  "  count format: mer count\n"
	  "\n.OPTIONS.\n"
	  "  -f <fastq> fastq file to count\n"
	  "  -k <len>   Length of kmer \n"
	  "  -l <limit> Gigabyte limit on RAM. If limited, the output will\n"
	  "             contain redundancies\n"
	  "  -t <size>  Define hash table size explicitly. [Default: chosen via k]\n"
	  "  -m <max>   Maximum k-mer count. [Default: 500]\n"
	  "  -q <num>   Quality value ascii scale, generally 64 or 33.  If\n"
	  "             not specified, it will guess.\n"
	  "\n");
  return;
}

////////////////////////////////////////////////////////////////////////////////
// parse_command_line
////////////////////////////////////////////////////////////////////////////////
static void parse_command_line(int argc, char **argv) {
  bool errflg = false;
  int ch;
  optarg = NULL;
  char* p;
  
  ////////////////////////////////////////
  // parse args
  ////////////////////////////////////////
  while(!errflg && ((ch = getopt(argc, argv, myopts)) != EOF)) {
    switch(ch) {

    case 'f':
      fastqfile = strdup(optarg);
      break;
      
    case 'k':
      Kmer_Len = strtol(optarg, &p, 10);
      if(p == optarg || Kmer_Len < 1) {
	fprintf(stderr, "Bad kmer length value \%s\"\n", optarg);
	errflg = true;
      }
      break;

    case 't':
      table_size = strtol(optarg, &p, 15);
      if(p == optarg || table_size <= 0) {
	fprintf(stderr, "Bad table size \%s\"\n", optarg);
	errflg = true;
      }
      break;

    case 'm':
      count_max = strtol(optarg, &p, 15);
      if(p == optarg || count_max <= 0) {
	fprintf(stderr, "Bad max count \%s\"\n", optarg);
	errflg = true;
      }
      break;

    case 'l':
      gb_limit = strtod(optarg, &p);
      break;

    case 'q': 
      quality_scale = int(strtol(optarg, &p, 10));
      if(p == optarg || quality_scale < -1) {
	fprintf(stderr, "Bad quality value scale \"%s\"\n",optarg);
	errflg = true;
      }
      break;
      
    case 'h':
      Usage(argv[0]);
      exit(EXIT_FAILURE);  

    case  '?' :
      fprintf (stderr, "Unrecognized option -%c\n", optopt);      
    default:
      errflg = true; 
    }
  }

  // return errors
  if(errflg || optind != argc) {
    Usage(argv[0]);
    exit(EXIT_FAILURE);
  }

  ////////////////////////////////////////
  // user input errors
  ////////////////////////////////////////
  if (Kmer_Len > 26 || Kmer_Len < 1) {
       cerr << "Kmer length must be <= 26" << endl;
       exit(1);
  }
  
  Forward_Mask = ((long long unsigned) 1 << (2 * Kmer_Len)) - 1;
}


////////////////////////////////////////////////////////////////////////////////
// guess_quality_scale
//
// Guess at ascii scale of quality values by examining
// a bunch of reads and looking for quality values < 64,
// in which case we set it to 33.
////////////////////////////////////////////////////////////////////////////////
static void guess_quality_scale(char* fqf) {
  string header, seq, mid, strqual;
  int reads_to_check = 1000;
  int reads_checked = 0;
  ifstream reads_in(fqf);
  while(getline(reads_in, header)) {
    getline(reads_in, seq);
    getline(reads_in, mid);
    getline(reads_in, strqual);
    
    for(int i = 0; i < strqual.size(); i++) {
      if(strqual[i] < 64) {
	cerr << "Guessing quality values are on ascii 33 scale" << endl;
	quality_scale = 33;
	reads_in.close();
	return;
      }
    }

    if(++reads_checked >= reads_to_check)
      break;
  }
  reads_in.close();
  cerr << "Guessing quality values are on ascii 64 scale" << endl;
  quality_scale = 64;
}

//////////////////////////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////////////////////////
int  main (int argc, char * argv [])
{
  parse_command_line(argc, argv);

  if(table_size == 0) {
       // guess based on k
       table_size = 30 * ((unsigned long long)pow(4, Kmer_Len) / 200);
  }

  // limit usurps table size
  if(gb_limit != 0) {
       // use user size only if it's smaller
       unsigned long long limit_size = (unsigned long long)((gb_limit*1073741824)/8);
       if(limit_size < table_size)
	    table_size = limit_size;
  }

  // but make size a power of 2
  unsigned long long size2 = 1024;
  while(size2 < table_size)
       size2 <<= 1;
  size2 >>= 1;

  cerr << "Table size: " << size2 << endl;

  qmer_hash mer_table(size2, Kmer_Len, count_max);

  FILE * fp;
  if(strcmp(fastqfile,"-") == 0) {
    fp = stdin;
    if(quality_scale == -1) {
      cerr << "Cannot guess at quality scale on reads from stdin- assuming 64." << endl;
      quality_scale = 64;
    }
  } else {
    cerr << fastqfile << endl;
    fp = fopen(fastqfile, "r");
    if (!fp)
      {
        cerr << "Couldn't open " << fastqfile << endl;
        exit(1);
      }
    if(quality_scale == -1)
      guess_quality_scale(fastqfile);
  }

  cerr << "Processing sequences..." << endl;

  string s, q, tag;
  unsigned int proc_seq = 0;
  
  while(Fastq_Read(fp, s, tag, q)) {
    CountMers(s, q, mer_table);
    if(mer_table.load() > 0.9) {
      // print table
      cerr << COUNT << " sequences processed, " << LEN << " bp scanned" << endl;
      mer_table.print();
      // clear table
      mer_table.clear();
    }
    if(++proc_seq == 1000000) {
      cerr << ".";
      proc_seq = 0;
    }
  }

  cerr << COUNT << " sequences processed, " << LEN << " bp scanned" << endl;
  
  if (BAD_CHAR)
     cerr << "WARNING: Input had " << BAD_CHAR << " non-DNA (ACGT) characters whose kmers were not counted" << endl;

  // final print
  mer_table.print();

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// compare_mers
//
// Return true if the first mer is lexicographically before the second
////////////////////////////////////////////////////////////////////////////////
bool compare_mers(Mer_t mer1, Mer_t mer2) {
     unsigned long long nt_mask = Forward_Mask - ((1ULL << (2*Kmer_Len-2)) - 1);
     unsigned long long nt1, nt2;
     for(int i = 0; i < Kmer_Len; i++) {
	  nt1 = nt_mask & mer1;
	  nt2 = nt_mask & mer2;

	  if(nt1 < nt2)
	       return true;
	  else if(nt2 < nt1)
	       return false;

	  mer1 <<= 2;
	  mer2 <<= 2;
     }     
     return false;
}


////////////////////////////////////////////////////////////////////////////////
// CountMers
//
// I edited this function to detect non ACGT's and ignore
// the Kmer_Len affected kmers.
////////////////////////////////////////////////////////////////////////////////
static void  CountMers (const string & s, const string & q, qmer_hash & mer_table)
{
   Mer_t  fwd_mer, rev_mer;
   int  i, j, n;
   int non_acgt_buffer = 0;

   // convert quality values
   vector<double> quals;
   double quality = 1.0;
   for(i = 0; i < q.size(); i++) {
      quals.push_back(max(.25, 1.0-pow(10.0,-(q[i]-quality_scale)/10.0)));
   }

   InitMer(fwd_mer);
   InitMer(rev_mer);

   n = s . length ();

   COUNT++;
   LEN+=n;

   if  (n < Kmer_Len) { return; }

   for  (i = 0;  i < Kmer_Len - 1;  i ++)
   {
     Forward_Add_Ch (fwd_mer, s [i]);
     Reverse_Add_Ch (rev_mer, s [i]);

     quality *= quals[i];

     if(strchr(bintoascii, s[i]) == NULL)
       non_acgt_buffer = Kmer_Len;
     else if(non_acgt_buffer > 0)
       non_acgt_buffer--;
   }

   while (i < n)
   {
     Forward_Add_Ch (fwd_mer, s [i]);
     Reverse_Add_Ch (rev_mer, s [i]);

     if(i == Kmer_Len-1)
       quality *= quals[i];
     else
       quality *= (quals[i] / quals[i - Kmer_Len]);

     if(strchr(bintoascii, s[i]) == NULL)
       non_acgt_buffer = Kmer_Len;
     else if(non_acgt_buffer > 0)
       non_acgt_buffer--;

     if(non_acgt_buffer == 0 && quality > .005) {
       string f,r;

       if(compare_mers(fwd_mer, rev_mer))
	    mer_table.add(fwd_mer, quality);
       else
	    mer_table.add(rev_mer, quality);
     }
     
     i++;
   }

   return;
}
