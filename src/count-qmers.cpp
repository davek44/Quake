#include  <string>
#include  <cstring>
#include  <vector>
#include  <iostream>
#include  <fstream>
#include  <math.h>
#include  "count.h"

using namespace std;
using namespace HASHMAP;

typedef hash_map<Mer_t, double, hash<unsigned long> > MerTable_t;

static void CountMers (const string & s, const string & q, MerTable_t & mer_table);
static void PrintMers(const MerTable_t & mer_table, int min_count);

////////////////////////////////////////////////////////////
// Additional options
////////////////////////////////////////////////////////////
int quality_scale = -1;

//////////////////////////////////////////////////////////////////////
// Usage
//////////////////////////////////////////////////////////////////////
static void Usage (char * command)
{
  fprintf(stderr, 
	  "\n"
	  ".USAGE.\n"
	  "  count-qmers [-f fastq] [-k len] [options]\n"
	  "\n"
	  ".DESCRIPTION.\n"
	  "  Count kmers in a fastq file. Output is to stdout in simple nmer"
	  "  count format: mer count\n"
	  "\n.OPTIONS.\n"
	  "  -f <fastq> fastq file to count\n"
	  "  -k <len>   Length of kmer \n"
	  "  -m <min>   Minimum count to report (default: >0)\n"
	  "  -l <limit> Gigabyte limit on RAM. If limited, the output will\n"
	  "             contain redundancies\n"
	  "  -q <num>   Quality value ascii scale, generally 64 or 33.  If\n"
	  "             not specified, it will guess.\n"
	  "\n");
  return;
}


//////////////////////////////////////////////////////////////////////
// parse_command_line
//////////////////////////////////////////////////////////////////////
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

    case 'l':
      gb_limit = strtod(optarg, &p);
      break;

    case 'm':
      min_count = strtol(optarg, &p, 10);
      if(p == optarg || min_count <= 0) {
	fprintf(stderr, "Bad min count value \%s\"\n", optarg);
	errflg = true;
      }
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
  if (Kmer_Len > 31 || Kmer_Len < 1)
    {
      cerr << "Kmer length must be <= 31" << endl;
      exit(1);
    }
  Forward_Mask = ((long long unsigned) 1 << (2 * Kmer_Len)) - 1;
}


////////////////////////////////////////////////////////////
// guess_quality_scale
//
// Guess at ascii scale of quality values by examining
// a bunch of reads and looking for quality values < 64,
// in which case we set it to 33.
////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////
int  main (int argc, char * argv [])
{
  parse_command_line(argc, argv);

  MerTable_t mer_table;

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
  unsigned long mb_limit = (unsigned long)(1024.0*gb_limit);
  unsigned long kmer_limit = mb_limit * 1048576UL / (unsigned long)bytes_per_kmer;
  unsigned int proc_seq = 0;
  
  while(Fastq_Read(fp, s, tag, q)) {
    CountMers(s, q, mer_table);
    if(gb_limit > 0 && mer_table.size() > kmer_limit) {
      // print table
      cerr << COUNT << " sequences processed, " << LEN << " bp scanned" << endl;
      PrintMers(mer_table, min_count);
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
    {
      cerr << "WARNING: Input had " << BAD_CHAR << " non-DNA (ACGT) characters whose kmers were not counted" << endl;
    }

  PrintMers(mer_table, min_count);

  return 0;
}


////////////////////////////////////////////////////////////
// CountMers
//
// I edited this function to detect non ACGT's and ignore
// the Kmer_Len affected kmers.
////////////////////////////////////////////////////////////
static void  CountMers (const string & s, const string & q, MerTable_t & mer_table)
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

   MerTable_t::iterator fi;

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

     if(non_acgt_buffer == 0 && quality > .0001) {
       string f,r;

       MerToAscii(fwd_mer, f);
       MerToAscii(rev_mer, r);
       
       if (f < r)
       {
	 fi = mer_table.find(fwd_mer);
	 if (fi == mer_table.end()) { fi=mer_table.insert(make_pair(fwd_mer,0)).first; }
       }
       else
       {
	 fi = mer_table.find(rev_mer);
	 if (fi == mer_table.end()) { fi=mer_table.insert(make_pair(rev_mer,0)).first; }
       }

       fi->second += quality;
     }
     
     i++;
   }

   return;
}


static void PrintMers(const MerTable_t & mer_table, int min_count)
{
  cerr << mer_table.size() << " total distinct mers" << endl;
  string mer;
  int printed = 0;

  MerTable_t::const_iterator fi;
  for (fi = mer_table.begin(); fi != mer_table.end(); fi++)
  {
    if (fi->second >= min_count)
    {
      MerToAscii(fi->first, mer);
      if (PRINT_SIMPLE)
      {
        printf("%s\t%f\n", mer.c_str(), fi->second);
      }
      else
      {
        printf(">%d\n%f\n", fi->second, mer.c_str());
      }
      printed++;
    }
  }

  cerr << printed << " mers occur at least " << min_count << " times" << endl;
}
