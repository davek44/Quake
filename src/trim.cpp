#include "Read.h"
#include "bithash.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <getopt.h>
#include <omp.h>
#include <cstdlib>
#include <iomanip>

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////
const static char* myopts = "r:o:t:q:p:I";
// -r, fastq file of reads
static char* fastqf = NULL;
// -o, output file
static char* outf = "out.txt";
// -q
static int trimq = 3;
// -I
bool Read::illumina_qual = false;
// -t
static int trim_t = 30;
// -p, number of threads
static int threads = 4;

// Note: to not trim, set trimq=0 and trim_t>read_length-k

// constants
static const char* nts = "ACGTN";

static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  build-icm [options] output_file < input-file\n"
           "\n"
           "Read sequences from standard input and output to  output-file\n"
           "the interpolated context model built from them.\n"
           "Input also can be piped into the program, e.g.,\n"
           "  cat abc.in | build-icm xyz.icm\n"
           "If <output-file> is \"-\", then output goes to standard output\n"
           "\n"
           "Options:\n"
           " -d <num>\n"
           "    Set depth of model to <num>\n"
           " -F\n"
           "    Ignore input strings with in-frame stop codons\n"
           " -h\n"
           "    Print this message\n"
           " -p <num>\n"
           "    Set period of model to <num>\n"
           " -r\n"
           "    Use the reverse of input strings to build the model\n"
           " -t\n"
           "    Output model as text (for debugging only)\n"
           " -v <num>\n"
           "    Set verbose level; higher is more diagnostic printouts\n"
           " -w <num>\n"
           "    Set length of model window to <num>\n"
           "\n");

   return;
  }

////////////////////////////////////////////////////////////
// parse_command_line
////////////////////////////////////////////////////////////
static void parse_command_line(int argc, char **argv) {
  bool errflg = false;
  int ch;
  optarg = NULL;
  char* p;
  
  // parse args
  while(!errflg && ((ch = getopt(argc, argv, myopts)) != EOF)) {
    switch(ch) {
    case 'r':
      fastqf = strdup(optarg);
      break;

    case 'o':
      outf = strdup(optarg);
      break;

    case 'q':
      trimq = int(strtol(optarg, &p, 10));
      if(p == optarg || trimq < 0) {
	fprintf(stderr, "Bad trim quality value \"%s\"\n",optarg);
	errflg = true;
      }
      break;

    case 't':
      trim_t = int(strtol(optarg, &p, 10));
      if(p == optarg || trim_t < 1) {
	fprintf(stderr, "Bad trim threshold \"%s\"\n",optarg);
	errflg = true;
      }
      break;

    case 'I':
      Read::illumina_qual = true;
      break;

    case 'p':
      threads = int(strtol(optarg, &p, 10));
      if(p == optarg || threads <= 0) {
	fprintf(stderr, "Bad number of threads \"%s\"\n",optarg);
	errflg = true;
      }
      break;

    case  '?' :
      fprintf (stderr, "Unrecognized option -%c\n", optopt);

    default:
      errflg = true;
    }
  }

  // for some reason, optind is not advancing properly so this
  // always returns an error

  // return errors
  /*
  if(errflg || optind != argc-1) {
    Usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  */
}

////////////////////////////////////////////////////////////
// pa_params
////////////////////////////////////////////////////////////
void pa_params(vector<streampos> & starts, vector<streampos> & counts) {
  // count number of sequences
  int N = 0;
  ifstream reads_in(fastqf);
  string toss;
  while(getline(reads_in, toss))
    N++;
  reads_in.close();
  N /= 4;

  // determine counts per thread
  counts.push_back(N - (threads-1)*(N/threads));
  for(int i = 1; i < threads; i++) {
    counts.push_back(N / threads);
  }

  // find start points
  reads_in.open(fastqf);
  starts.push_back(reads_in.tellg());
  int s = 0;
  int t = 0;
  while(getline(reads_in,toss)) {
    // sequence
    getline(reads_in, toss);
    // +
    getline(reads_in, toss);
    // quality
    getline(reads_in, toss);

    if(++s == counts[t] && s != N) {
      starts.push_back(reads_in.tellg());
      s = 0;
      t++;
    }
  }
}

////////////////////////////////////////////////////////////
// trim_reads
////////////////////////////////////////////////////////////
static void trim_reads(vector<streampos> & starts, vector<streampos> & counts) {
#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();

    // output
    char* toutf = strdup(outf);
    char strtid[10];
    sprintf(strtid,"%d",tid);
    strcat(toutf, strtid);
    ofstream reads_out(toutf);

    // input
    ifstream reads_in(fastqf);
    reads_in.seekg(starts[tid]);
    
    string header,ntseq,strqual,mid;
    char* nti;
    Read *r;
    vector<int> untrusted;  // dummy
    vector<correction> cor; // dummy

    int tcount = 0;
    while(tcount++ < counts[tid] && getline(reads_in, header)) {
      //cout << header << endl;

      // get sequence
      getline(reads_in, ntseq);
      //cout << ntseq << endl;

      // convert ntseq to iseq
      vector<unsigned int> iseq;
      for(int i = 0; i < ntseq.size(); i++) {
	nti = strchr(nts, ntseq[i]);	
	iseq.push_back(nti - nts);
      }     
    
      // get quality values
      getline(reads_in,mid);
      //cout << mid << endl;
      getline(reads_in,strqual);
      //cout << strqual << endl;
      
      // trim
      r = new Read(header, &iseq[0], strqual, untrusted, iseq.size());      
      r->trim(trimq);
      ntseq = r->print_corrected(cor);

      // print if large enough
      if(ntseq.size() >= trim_t) {
	reads_out << header << endl << ntseq << endl << mid << endl << strqual.substr(0, ntseq.size()) << endl;
      }
      
      delete r;
    }
    reads_in.close();
  }
}

////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  // set up parallelism
  omp_set_num_threads(threads);
  vector<streampos> starts;
  vector<streampos> counts;
  pa_params(starts, counts);
  
  for(int i = 0; i < counts.size(); i++) {
    cout << starts[i] << " " << counts[i] << endl;
  }

  trim_reads(starts, counts);

  return 0;
}
