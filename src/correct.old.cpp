#include "prefix_tree.h"
#include "Read.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////
const static char* myopts = "r:m:oc:l";
// -r, fastq file of reads
static char* fastqf = NULL;
// -m, mer counts
static char* merf = NULL;
// -o, output file
static char* outf = "out.txt";
// -c, cutoff between trusted and untrusted mers
static int cutoff;
// -l, read length
static int read_len = 36;

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

    case 'm':
      merf = strdup(optarg);
      break;

    case 'o':
      outf = strdup(optarg);
      break;

    case 'c':
      cutoff = int(strtol(optarg, &p, 10));
      if(p == optarg || cutoff < 0) {
	fprintf(stderr, "Bad mer cutoff value \"%s\"\n",optarg);
	errflg = true;
      }
      break;

    case 'l':
      read_len = int(strtol(optarg, &p, 10));
      if(p == optarg || read_len <= 0) {
	fprintf(stderr, "Bad read length value \"%s\"\n",optarg);
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
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  prefix_tree trusted;
  trusted.file_load(merf, cutoff);

  ifstream reads_in(fastqf);
  ofstream reads_out(outf);
  string header;
  string ntseq;
  string strqual;
  int iseq[read_len];
  char* nti;
  int error_reads = 0;
  int fixed_reads = 0;
  const char* nts = "ACGTN";

  int sample = 0;
  while(getline(reads_in, header)) {
    if(sample > 500000)
      break;
    sample++;

    // get sequence
    getline(reads_in, ntseq);
    //cout << ntseq << endl;

    // convert ntseq to iseq
    for(int i = 0; i < read_len; i++) {
      nti = strchr(nts, ntseq[i]);
      iseq[i] = nti - nts;
      /*
      if(nti)
	iseq[i] = nti - nts;
      else
	iseq[i] = 0;   // change all non ACGT's to A
      */
    }

    // find untrusted kmers
    vector<int> untrusted;
    for(int i = 0; i < read_len-k+1; i++) {
      if(!trusted.check(&iseq[i])) {
	untrusted.push_back(i);
      }
    }
    
    // get quality values
    getline(reads_in,strqual);
    getline(reads_in,strqual);

    // save error reads
    if(untrusted.size() > 0) {
      Read* r = new Read(header, &iseq[0], strqual, untrusted, read_len);
      error_reads++;
      if(r->correct(&trusted, reads_out))
	fixed_reads++;
      delete r;
    }
  }

  cout << "Error reads: " << error_reads << endl;
  cout << "Fixed reads: " << fixed_reads << endl;

  return 0;
}
