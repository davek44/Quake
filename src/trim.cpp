#include "Read.h"
#include "bithash.h"
#include "edit.h"
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
const static char* myopts = "r:f:t:q:l:p:h";
// -r, fastq file of reads
//char* fastqf;
// -f, file of fastq files of reads
//char* file_of_fastqf;
// -t
//static int trimq = 3;
// -q
//int Read::quality_scale;
// -l
static int trim_t = 30;
// -p, number of threads
//int threads;

// constants
static const char* nts = "ACGTN";

//unsigned int chunks_per_thread;

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
	   "    Fastq file of reads\n"
	   " -f <file>\n"
	   "    File containing fastq file names, one per line or\n"
	   "    two per line for paired end reads.\n"
	   " -p <num>\n"
	   "    Use <num> openMP threads\n"
	   " -l <num>=30\n"
	   "    Return only reads corrected and/or trimmed to >= <num>\n"
	   "    bp\n"
	   " -q <num>\n"
	   "    Quality value ascii scale, generally 64 or 33. If not\n"
	   "    specified, it will guess.\n"
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

    case 'f':
      file_of_fastqf = strdup(optarg);
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
      if(p == optarg || Read::quality_scale < -1) {
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
  if(fastqf == NULL && file_of_fastqf == NULL) {
    cerr << "Must provide a fastq file of reads (-r) or a file containing a list of fastq files of reads (-f)" << endl;
    exit(EXIT_FAILURE);
  }
}


////////////////////////////////////////////////////////////
// trim_reads
////////////////////////////////////////////////////////////
static void trim_reads(string fqf, int pe_code, vector<streampos> & starts, vector<unsigned long long> & counts) {
  //format output file
  string path_suffix = split(fqf,'/').back();
  string out_dir("."+path_suffix);
  mkdir(out_dir.c_str(), S_IRWXU);

  unsigned int chunk = 0;
#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();

    // input
    ifstream reads_in(fqf.c_str());
    
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
      string toutf(out_dir+"/");
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
	} else if(pe_code > 0) {
	  reads_out << header << " error" << endl << ntseq << endl << mid << endl << strqual.substr(0, ntseq.size()) << endl;
	}
	
	delete r;
      }
      reads_out.close();

      #pragma omp critical
      tchunk = chunk++;
    }
    reads_in.close();
  }

  if(pe_code == 0)
    combine_output(fqf, string("trim"));
}


////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  // make list of files
  vector<string> fastqfs;
  vector<int> pairedend_codes;
  parse_fastq(fastqfs, pairedend_codes);

  // determine quality value scale
  if(Read::quality_scale == -1)
    guess_quality_scale(fastqfs[0]);

  // process each file
  string fqf;
  for(int f = 0; f < fastqfs.size(); f++) {
    fqf = fastqfs[f];
    cout << fqf << endl; 

    // split up file
    vector<streampos> starts;
    vector<unsigned long long> counts;
    chunkify_fastq(fqf, starts, counts);
    
    /*
    for(int i = 0; i < counts.size(); i++)
      cout << i << " " << starts[i] << " " << counts[i] << endl;
    */
    
    trim_reads(fqf, pairedend_codes[f], starts, counts);

    // combine paired end
    if(pairedend_codes[f] == 2)
      combine_output_paired(fastqfs[f-1], fqf, string("trim"));
  }

  return 0;
}
