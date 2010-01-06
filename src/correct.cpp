//#include "prefix_tree.h"
#include "bithash.h"
#include "Read.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <omp.h>

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
// -p, number of threads
static int threads = 1;

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
// pa_params
////////////////////////////////////////////////////////////
void pa_params(vector<int> & starts, vector<int> & counts) {
  // count number of sequences
  int N = 0;
  ifstream reads_in(fastqf);
  string toss;
  while(getline(reads_in, toss))
    N++;
  reads_in.close();
  N /= 4;

  // determine counts per thread
  int sum = 0;
  for(int i = 0; i < threads-1; i++) {
    counts.push_back(N / threads);
    sum += N/threads;
  }
  counts.push_back(N - sum);

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
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  vector<int> starts;
  vector<int> counts;
  pa_params(starts, counts);
  for(int i = 0; i < threads; i++)
    cout << counts[i] << " " << starts[i] << endl;

  int error_reads = 0;
  int fixed_reads = 0;

  //prefix_tree trusted;
  bithash trusted;
  trusted.file_load(merf, cutoff);


  omp_set_num_threads(threads);
  //cout << omp_get_max_threads() << " threads" << endl;

#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();
    //cout << "Initializing thread " << tid << endl;

    char* toutf = strdup(outf);
    char strtid[10];
    sprintf(strtid,"%d",tid);
    strcat(toutf, strtid);
    ofstream reads_out(toutf);

    ifstream reads_in(fastqf);
    reads_in.seekg(starts[tid]);
    
    string header,ntseq,strqual;
    unsigned int iseq[read_len];
    char* nti;
    const char* nts = "ACGTN";
    Read *r;

    int tcount = 0;
    while(getline(reads_in, header)) {
      //cout << header << endl;

      // get sequence
      getline(reads_in, ntseq);
      //cout << ntseq << endl;

      // convert ntseq to iseq
      for(int i = 0; i < read_len; i++) {
	nti = strchr(nts, ntseq[i]);
	iseq[i] = nti - nts;
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
      //cout << strqual << endl;
      getline(reads_in,strqual);
      //cout << strqual << endl;
      
      // fix error reads
      if(untrusted.size() > 0) {
	r = new Read(header, &iseq[0], strqual, untrusted, read_len);
	error_reads++;
	if(r->correct(&trusted, reads_out))
	  fixed_reads++;
	delete r;
      }

      if(++tcount == counts[tid])
	break;
    }
    reads_in.close();
  }

  cout << "Error reads: " << error_reads << endl;
  cout << "Fixed reads: " << fixed_reads << endl;

  return 0;
}
