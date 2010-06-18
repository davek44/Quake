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
// -t
static int trimq = 3;
// -q
int Read::quality_scale = -1;
// -l
static int trim_t = 30;
// -p, number of threads
static int threads = 4;

// constants
static const char* nts = "ACGTN";

static unsigned int chunks_per_thread = 200;

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
           "Trims reads in a fastq file.  This program will shuffle the\n"
	   "sequences.  If there is demand for a version that maintains\n"
	   "the order, I can fix that - dakelley@umiacs.umd.edu\n"
           "\n"
	   "Options:\n"
	   " -r <file>\n"
	   "    Fastq file of reads to trim\n"
	   " -p <num>\n"
	   "    Use <num> openMP threads\n"
	   " -l <num>=30\n"
	   "    Return only reads corrected and/or trimmed to >= <num>\n"
	   "    bp\n"
	   " -q <num>\n"
	   "    Quality value ascii scale.  Older reads may be 33.  If\n"
	   "    not specified, it will guess.\n"
	   " -t <num>=3\n"
	   "    Use BWA trim parameter <num>\n"
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

    case 't':
      trimq = int(strtol(optarg, &p, 10));
      if(p == optarg || trimq < 0) {
	fprintf(stderr, "Bad trim quality value \"%s\"\n",optarg);
	errflg = true;
      }
      break;

    case 'l':
      trim_t = int(strtol(optarg, &p, 10));
      if(p == optarg || trim_t < 1) {
	fprintf(stderr, "Bad trim threshold \"%s\"\n",optarg);
	errflg = true;
      }
      break;

    case 'q': 
      Read::quality_scale = int(strtol(optarg, &p, 10));
      if(p == optarg || Read::quality_scale < 0) {
	fprintf(stderr, "Bad quality value scale \"%s\"\n",optarg);
	errflg = true;
      }
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
void pa_params(vector<streampos> & starts, vector<unsigned long long> & counts) {
  // count number of sequences
  unsigned long long N = 0;
  ifstream reads_in(fastqf);
  string toss;
  while(getline(reads_in, toss))
    N++;
  reads_in.close();
  N /= 4ULL;

  if(threads*chunks_per_thread > N) {
    // use 1 thread for everything
    counts.push_back(N);
    starts.push_back(0);   
    omp_set_num_threads(1);

  } else {
    // determine counts per thread
    unsigned long long sum = 0;
    for(int i = 0; i < threads*chunks_per_thread-1; i++) {
      counts.push_back(N / (threads*chunks_per_thread));
      sum += counts.back();
    }
    counts.push_back(N - sum);

    // find start points
    reads_in.open(fastqf);
    starts.push_back(reads_in.tellg());
    unsigned long long s = 0;
    unsigned int t = 0;
    while(getline(reads_in,toss)) {
      // sequence
      getline(reads_in, toss);
      // +
      getline(reads_in, toss);
      // quality
      getline(reads_in, toss);
      
      if(++s == counts[t] && t < counts.size()-1) {
	starts.push_back(reads_in.tellg());
	s = 0;
	t++;
      }
      
      // set up parallelism
      omp_set_num_threads(threads);
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
  
  for(int t = 0; t < threads*chunks_per_thread; t++) {
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

  combine_out.close();
}


////////////////////////////////////////////////////////////
// guess_quality_scale
//
// Guess at ascii scale of quality values by examining
// a bunch of reads and looking for quality values < 64,
// in which case we set it to 33.
////////////////////////////////////////////////////////////
static void guess_quality_scale() {
  string header, seq, mid, strqual;
  int reads_to_check = 1000;
  int reads_checked = 0;
  ifstream reads_in(fastqf);
  while(getline(reads_in, header)) {
    getline(reads_in, seq);
    getline(reads_in, mid);
    getline(reads_in, strqual);
    
    for(int i = 0; i < strqual.size(); i++) {
      if(strqual[i] < 64) {
	cerr << "Guessing quality values are on ascii 33 scale" << endl;
	Read::quality_scale = 33;
	reads_in.close();
	return;
      }
    }

    if(++reads_checked >= reads_to_check)
      break;
  }
  reads_in.close();
  cerr << "Guessing quality values are on ascii 64 scale" << endl;
  Read::quality_scale = 64;
}


////////////////////////////////////////////////////////////
// trim_reads
////////////////////////////////////////////////////////////
static void trim_reads(vector<streampos> & starts, vector<unsigned long long> & counts) {
  //format output file
  string fqf_str(fastqf);
  int suffix_index = fqf_str.rfind(".");
  string prefix = fqf_str.substr(0, suffix_index);
  string suffix = fqf_str.substr(suffix_index, fqf_str.size()-suffix_index);
  string outf = prefix + string(".trim") + suffix;

  unsigned int chunk = 0;
#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();

    // input
    ifstream reads_in(fastqf);
    
     unsigned int tchunk;
    string header,ntseq,strqual,mid;
    char* nti;
    Read *r;
    vector<int> untrusted;  // dummy
    vector<correction> cor; // dummy
    
    #pragma omp critical
    tchunk = chunk++;

    while(tchunk < starts.size()) {
      reads_in.seekg(starts[tchunk]);

      // output
      string toutf(outf);
      stringstream tconvert;
      tconvert << tchunk;
      toutf += tconvert.str();
      ofstream reads_out(toutf.c_str());
      
      unsigned long long tcount = 0;
      while(tcount++ < counts[tchunk] && getline(reads_in, header)) {
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
      reads_out.close();

      #pragma omp critical
      tchunk = chunk++;
    }
    reads_in.close();
  }
  combine_output(outf.c_str());
}


////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  guess_quality_scale();

  // set up parallelism
  vector<streampos> starts;
  vector<unsigned long long> counts;
  pa_params(starts, counts);
    
  /*
  for(int i = 0; i < counts.size(); i++) {
    cout << i << " " << starts[i] << " " << counts[i] << endl;
  }
  */
  
  trim_reads(starts, counts);

  return 0;
}
