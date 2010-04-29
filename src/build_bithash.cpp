#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <getopt.h>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include "bithash.h"

using namespace::std;

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////
const static char* myopts = "m:c:o:";
// -m, kmer count file
static char* merf = NULL;
// -c, kmer count trusted cutoff
static double cutoff = NULL;
// -a, AT cutoff
static char* ATcutf = NULL;
// -o, bithash output file
static char* outf = "bithash.out";

static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  build_bithash [options]\n"
           "\n"
	    "Load kmer counts and build bithash data structure from\n"
	    "trusted kmers.\n"
           "\n"
           "Options:\n"
	   " -m <file>\n"
	   "    File containg kmer counts in format `seq\tcount`.\n"
	   "    Can also be piped in with '-'\n"
	   " -c <num>\n"
	   "    Separate trusted/untrusted kmers at cutoff <num>\n"
	   " -a <file>\n"
	   "    Separate trusted/untrusted kmers as a function of\n"
	   "    AT content, with cutoffs found in <file>, one per line\n"
	   " -o <file>\n"
	   "    Bithash will be dumped as binary to <file>\n"
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
    case 'm':
      merf = strdup(optarg);
      break;

    case 'c':
      cutoff = double(strtod(optarg, &p));
      if(p == optarg || cutoff < 0) {
	fprintf(stderr, "Bad mer cutoff value \"%s\"\n",optarg);
	errflg = true;
      }
      break;

    case 'a':
      ATcutf = strdup(optarg);
      break;

    case 'o':
      outf = strdup(optarg);
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

  ////////////////////////////////////////
  // correct user input errors
  ////////////////////////////////////////  
  if(cutoff == 0 && ATcutf == NULL) {
    cerr << "Must provide a trusted/untrusted kmer cutoff (-c) or a file containing the cutoff as a function of the AT content (-a)" << endl;
    exit(EXIT_FAILURE);
  }
  if(merf == NULL) {
    cerr << "Must provide a file of kmer counts (-m)" << endl;
    exit(EXIT_FAILURE);
  }
}

////////////////////////////////////////////////////////////
// load_AT_cutoffs
//
// Load AT cutoffs from file
////////////////////////////////////////////////////////////
vector<double> load_AT_cutoffs() {
  vector<double> cutoffs;
  ifstream cut_in(ATcutf);
  string line;
  double cut;
  
  while(getline(cut_in, line)) {
    stringstream ss(stringstream::in | stringstream::out);
    ss << line;
    ss >> cut;
    cutoffs.push_back(cut);
  }

  if(cutoffs.size() != (k+1)) {
    cerr << "Must specify " << (k+1) << " AT cutoffs in " << ATcutf << endl;
    exit(EXIT_FAILURE);
  }

  return cutoffs;
}


////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);
  
  // make trusted kmer data structure
  bithash *trusted = new bithash();
  if(ATcutf != NULL) {
    if(strcmp(merf,"-") == 0)
      trusted->tab_file_load(cin, load_AT_cutoffs(), NULL);
    else {
      ifstream mer_in(merf);
      trusted->tab_file_load(mer_in, load_AT_cutoffs(), NULL);
    }
  } else {
    if(strcmp(merf,"-") == 0) {
      trusted->tab_file_load(cin, cutoff, NULL);
    } else {
      ifstream mer_in(merf);
      trusted->tab_file_load(mer_in, cutoff, NULL);
    }
  }
  cout << trusted->num_kmers() << " trusted kmers" << endl;
  
  // write to file  
  trusted->binary_file_output(outf);
    
  return 0;
}
