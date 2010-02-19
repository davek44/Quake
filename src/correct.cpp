//#include "prefix_tree.h"
#include "bithash.h"
//#include "dawg.h"
#include "Read.h"
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
const static char* myopts = "r:m:o:c:t:l:";
// -r, fastq file of reads
static char* fastqf = NULL;
// -m, mer counts
static char* merf = NULL;
// -o, output file
static char* outf = "out.txt";
// -c, cutoff between trusted and untrusted mers
static int cutoff;
// -t
static int trimq = 3;
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

    case 't':
      trimq = int(strtol(optarg, &p, 10));
      if(p == optarg || trimq < 0) {
	fprintf(stderr, "Bad trim quality value \"%s\"\n",optarg);
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
// correct_reads
////////////////////////////////////////////////////////////
static void correct_reads(bithash * trusted, vector<int> & starts, vector<int> & counts, double (&ntnt_prob)[4][4]) {
  int error_reads = 0;
  int fixed_reads = 0;

#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();
    char* toutf = strdup(outf);
    char strtid[10];
    sprintf(strtid,"%d",tid);
    strcat(toutf, strtid);
    ofstream reads_out(toutf);

    ifstream reads_in(fastqf);
    reads_in.seekg(starts[tid]);
    
    string header,ntseq,strqual;    
    char* nti;
    Read *r;

    int tcount = 0;
    while(getline(reads_in, header)) {
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
      
      vector<int> untrusted;
      for(int i = 0; i < iseq.size()-k+1; i++) {
	if(!trusted->check(&iseq[i])) {
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
	r = new Read(header, &iseq[0], strqual, untrusted, iseq.size());
	error_reads++;

	if(r->trim(trimq, reads_out) || r->correct(trusted, reads_out, ntnt_prob))
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
}

////////////////////////////////////////////////////////////
// learn_errors
//
// Correct reads using a much stricter filter in order
// to count the nt->nt errors and learn the errors
// probabilities
////////////////////////////////////////////////////////////
static void learn_errors(bithash * trusted, vector<int> & starts, vector<int> & counts, double (&ntnt_prob)[4][4]) {
  int ntnt_counts[4][4] = {0};
  unsigned int samples = 0;
#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();
    char* toutf = strdup(outf);
    char strtid[10];
    sprintf(strtid,"%d",tid);
    strcat(toutf, strtid);
    ofstream reads_out(toutf);

    ifstream reads_in(fastqf);
    reads_in.seekg(starts[tid]);
    
    string header,ntseq,strqual;
    char* nti;
    Read *r;
    //correction* cor;

    int tcount = 0;
    while(getline(reads_in, header)) {
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
      
      vector<int> untrusted;
      for(int i = 0; i < iseq.size()-k+1; i++) {
	if(!trusted->check(&iseq[i])) {
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
	//cout << "Processing " << header;
	//for(int i = 0; i < untrusted.size(); i++) {
	//  cout << " " << untrusted[i];	  
	//}
	//cout << endl;

	r = new Read(header, &iseq[0], strqual, untrusted, iseq.size());

	// if trimmed to no errors, trusted_read doesn't exist
	if(!r->trim(trimq, reads_out) && r->correct(trusted, reads_out, ntnt_prob, false)) {
	  for(int c = 0; c < r->trusted_read->corrections.size(); c++) {
	    correction cor = r->trusted_read->corrections[c];
	    if(iseq[cor.index] < 4) {
	      // real -> observed, (real approximated by correction)
	      ntnt_counts[cor.to][iseq[cor.index]]++;
	      //ntnt_counts[iseq[cor->index]][cor->to]++;
	      samples++;
	    }
	  }
	}
	delete r;
      }

      if(++tcount == counts[tid] || samples > 35000)
	break;
    }
    reads_in.close();
  }

  cout << "Count matrix:" << endl;
  cout << "\tA\tC\tG\tT" << endl;
  for(int i = 0; i < 4; i++) {
    cout << nts[i];
    for(int j = 0; j < 4; j++) {
      if(i == j)
	cout << "\t-";
      else
	cout << "\t" << ntnt_counts[i][j];
    }
    cout << endl;
  }
  cout << endl;

  // make counts into probabilities
  int ntsum;
  for(int i = 0; i < 4; i++) {
    // sum all corrections from this nt
    ntsum = 0;
    for(int j = 0; j < 4; j++) {
      ntsum += ntnt_counts[i][j];
    }

    // normalize counts by sum
    for(int j = 0; j < 4; j++) {
      ntnt_prob[i][j] = (double)ntnt_counts[i][j] / double(ntsum);
    }
  }
}


////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  // set up parallelism
  omp_set_num_threads(threads);
  vector<int> starts;
  vector<int> counts;
  pa_params(starts, counts);

  // make trusted kmer data structure
  //prefix_tree *trusted = new prefix_tree;
  bithash *trusted = new bithash();
  //dawg *trusted = new dawg();
  trusted->tab_file_load(merf, cutoff);

  // learn nt->nt transitions
  double ntnt_prob[4][4] = {0};
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      if(i != j)
	ntnt_prob[i][j] = 1.0/3.0;
    }
  }
  learn_errors(trusted, starts, counts, ntnt_prob);

  cout << "New error matrix:" << endl;
  cout << "\tA\tC\tG\tT" << endl;
  for(int i = 0; i < 4; i++) {
    cout << nts[i];
    for(int j = 0; j < 4; j++) {
      if(i == j)
	cout << "\t-";
      else
	cout << "\t" << setprecision(4) << ntnt_prob[i][j];
    }
    cout << endl;
  }

  correct_reads(trusted, starts, counts, ntnt_prob);  

  return 0;
}
