#include "Read.h"
#include "bithash.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <string.h>
#include <getopt.h>
#include <omp.h>
#include <cstdlib>
#include <iomanip>
#include <sys/stat.h>

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////
const static char* myopts = "r:t:q:p:Ih";
// -r, fastq file of reads
static char* fastqf = NULL;
// -q
static int trimq = 3;
// -I
bool Read::illumina_qual = false;
// -t
static int trim_t = 30;
// -p, number of threads
static int threads = 4;

// constants
static const char* nts = "ACGTN";

////////////////////////////////////////////////////////////
// Usage
//
//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.
////////////////////////////////////////////////////////////
static void  Usage(char * command)
{
  fprintf (stderr,
           "USAGE:  trim [options]\n"
           "\n"
           "Trims reads in a fastq file.\n"
           "\n"
	   "Options:\n"
	   " -r <file>\n"
	   "    Fastq file of reads to trim\n"
	   " -p <num>\n"
	   "    Use <num> openMP threads\n"
	   " -t <num>=30\n"
	   "    Return only reads trimmed to >= <num> bp\n"
	   " -q <num>=3\n"
	   "    Use BWA trim parameter <num>\n"
	   " -I\n"
	   "    Use 64 scale Illumina quality values (else base 33)\n"
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

    case 'h':
      Usage(argv[0]);
      exit(EXIT_FAILURE);

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

  ////////////////////////////////////////
  // correct user input errors
  ////////////////////////////////////////
  if(fastqf == NULL) {
    cerr << "Must provide a fastq file of reads with -r" << endl;
    exit(EXIT_FAILURE);
  }
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
// combine_output
//
// Combine output files into one.
////////////////////////////////////////////////////////////
static void combine_output(const char* outf) {
  ofstream combine_out(outf);
  char strt[10];
  char toutf[50];
  string line;
  struct stat st_file_info;
  
  for(int t = 0; t < threads; t++) {
    strcpy(toutf, outf);
    sprintf(strt,"%d",t);
    strcat(toutf, strt);

    // if file exists, add to single output
    if(stat(toutf, &st_file_info) == 0) {
      ifstream thread_out(toutf);
      while(getline(thread_out, line))
	combine_out << line << endl;
      thread_out.close();
      
      remove(toutf);
    }
  }
}

////////////////////////////////////////////////////////////
// trim_reads
////////////////////////////////////////////////////////////
static void trim_reads(vector<streampos> & starts, vector<streampos> & counts) {
  //format output file
  string fqf_str(fastqf);
  int suffix_index = fqf_str.rfind(".");
  string prefix = fqf_str.substr(0, suffix_index);
  string suffix = fqf_str.substr(suffix_index, fqf_str.size()-suffix_index);
  string outf = prefix + string(".trim") + suffix;

#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();

    // output
    string toutf(outf);
    stringstream tconvert;
    tconvert << tid;
    toutf += tconvert.str();
    ofstream reads_out(toutf.c_str());

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
    reads_out.close();
  }

  combine_output(outf.c_str());
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
