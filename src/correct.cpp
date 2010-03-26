#include "bithash.h"
#include "Read.h"
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

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////
const static char* myopts = "r:f:m:c:a:t:q:p:z:ICh";
// -r, fastq file of reads
static char* fastqf = NULL;
// -f, file of fastq files of reads
static char* file_of_fastqf = NULL;
// -m, mer counts
static char* merf = NULL;
// -c, cutoff between trusted and untrusted mers
static double cutoff = 0;
// -a, AT cutoff
static char* ATcutf = NULL;
// -q
static int trimq = 3;
// -I
bool Read::illumina_qual = false;
// -t
static int trim_t = 30;
// -p, number of threads
static int threads = 4;
// -z, zip mode directory
static char* zipd = NULL;
// -C, Contrail output
static bool contrail_out = false;

// Note: to not trim, set trimq=0 and trim_t>read_length-k

// constants
#define TESTING true
static const char* nts = "ACGTN";

static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

{
  fprintf (stderr,
           "USAGE:  correct [options]\n"
           "\n"
	   "Correct sequencing errors in fastq file provided with -r\n"
	   "and output trusted and corrected reads to\n"
	   "<fastq-prefix>.cor.fastq.\n"
           "\n"
           "Options:\n"
           " -r <file>\n"
	   "    Fastq file of reads\n"
	   " -f <file>\n"
	   "    File containing fastq file names, one per line\n"
	   " -m <file>\n"
	   "    File containg kmer counts in format `seq\tcount`.\n"
	   "    Can also be piped in with '-'\n"
	   " -c <num>\n"
	   "    Separate trusted/untrusted kmers at cutoff <num>\n"
	   " -a <file>\n"
	   "    Separate trusted/untrusted kmers as a function of\n"
	   "    AT content, with cutoffs found in <file>, one per line\n"
	   " -p <num>\n"
	   "    Use <num> openMP threads\n"
	   " -t <num>=30\n"
	   "    Return only reads trimmed to >= <num> bp\n"
	   " -q <num>=3\n"
	   "    Use BWA trim parameter <num>\n"
	   " -I\n"
	   "    Use 64 scale Illumina quality values\n"
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

    case 'm':
      merf = strdup(optarg);
      break;

    case 'z':
      zipd = strdup(optarg);
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

    case 'C':
      contrail_out = true;
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

  if(merf == NULL) {
    cerr << "Must provide a file of kmer counts (-m)" << endl;
    exit(EXIT_FAILURE);
  }  
}

////////////////////////////////////////////////////////////
// pa_params
////////////////////////////////////////////////////////////
void pa_params(string fqf, vector<streampos> & starts, vector<unsigned long long> & counts) {
  // count number of sequences
  unsigned long long N = 0;
  ifstream reads_in(fqf.c_str());
  string toss;
  while(getline(reads_in, toss))
    N++;
  reads_in.close();
  N /= 4ULL;

  // determine counts per thread
  unsigned long long sum = 0;
  for(int i = 0; i < threads-1; i++) {
    counts.push_back(N / threads);
    sum += N/threads;
  }
  counts.push_back(N - sum);

  // find start points
  reads_in.open(fqf.c_str());
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
  char mycmd[100];

  // clean file
  ofstream out_all(outf);
  out_all.close();
  
  char strt[10];
  for(int t = 0; t < threads; t++) {
    sprintf(strt,"%d",t);
    
    // cat thread output
    strcpy(mycmd, "cat ");
    strcat(mycmd, outf);
    strcat(mycmd, strt);
    strcat(mycmd, " >> ");
    strcat(mycmd, outf);
    system(mycmd);
    
    // rm thread output
    strcpy(mycmd, "rm ");
    strcat(mycmd, outf);
    strcat(mycmd, strt);
    system(mycmd);
  }
}

////////////////////////////////////////////////////////////
// correct_reads
////////////////////////////////////////////////////////////
static void correct_reads(string fqf, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double (&ntnt_prob)[4][4]) {
  // format output file
  int suffix_index = fqf.rfind(".");
  string prefix = fqf.substr(0,suffix_index);
  string suffix = fqf.substr(suffix_index, fqf.size()-suffix_index);
  string outf = prefix + string(".cor") + suffix;
  //cout << outf << endl;
  
#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();
    
    string toutf(outf);
    stringstream tconvert;
    tconvert << tid;
    toutf += tconvert.str();
    ofstream reads_out(toutf.c_str());

    /*
    char* toutf = strdup(outf.c_str());
    char strtid[10];
    sprintf(strtid,"%d",tid);
    strcat(toutf, strtid);
    cout << toutf << endl;
    ofstream reads_out(toutf);
    */

    ifstream reads_in(fqf.c_str());
    reads_in.seekg(starts[tid]);
    
    string header,ntseq,mid,strqual,corseq;
    char* nti;
    Read *r;

    unsigned long long tcount = 0;
    while(getline(reads_in, header)) {
      //cout << tid << " " << header << endl;

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
      getline(reads_in,mid);
      //cout << mid << endl;
      getline(reads_in,strqual);
      //cout << strqual << endl;
      
      // fix error reads
      if(untrusted.size() > 0) {
	r = new Read(header, &iseq[0], strqual, untrusted, iseq.size());

	// try to trim
	corseq = r->trim(trimq);	

	if(r->untrusted.empty()) {
	  // long enough?
	  if(corseq.size() >= trim_t) {
	    if(contrail_out)
	      reads_out << header << "\t" << corseq << endl;
	    else
	      reads_out << header << endl << corseq << endl << mid << endl << strqual.substr(0,corseq.size()) << endl;

	    if(TESTING)
	      cerr << header << "\t" << ntseq << "\t" << corseq << endl;
	  
	  // else throw away
	  } else {
	    if(TESTING)
	      cerr << header << "\t" << ntseq << "\t." << endl;
	  }

	} else {
	  // if still untrusted, correct
	  corseq = r->correct(trusted, ntnt_prob);

	  // if trimmed to long enough
	  if(corseq.size() >= trim_t) {
	    if(contrail_out)
	      reads_out << header << "\t" << corseq << endl;
	    else
	      reads_out << header << endl << corseq << endl << mid << endl << strqual.substr(0,corseq.size()) << endl;

	    if(TESTING)
	      cerr << header << "\t" << ntseq << "\t" << corseq << endl;

	  // else throw away
	  } else {
	    if(TESTING)
	      cerr << header << "\t" << ntseq << "\t-" << endl;
	  }
	}

	delete r;
      } else {
	// output read as is
	if(contrail_out)
	  reads_out << header << "\t" << ntseq << endl;
	else
	  reads_out << header << endl << ntseq << endl << mid << endl << strqual << endl;
      }

      if(++tcount == counts[tid])
	break;
    }
    reads_in.close();
  }
  
  //combine_output(outf.c_str());
}

////////////////////////////////////////////////////////////
// learn_errors
//
// Correct reads using a much stricter filter in order
// to count the nt->nt errors and learn the errors
// probabilities
////////////////////////////////////////////////////////////
static void learn_errors(string fqf, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double (&ntnt_prob)[4][4]) {
  int ntnt_counts[4][4] = {0};
  unsigned int samples = 0;

#pragma omp parallel //shared(trusted)
  {
    int tid = omp_get_thread_num();

    ifstream reads_in(fqf.c_str());
    reads_in.seekg(starts[tid]);
    
    string header,ntseq,strqual,corseq;
    char* nti;
    Read *r;    
    //correction* cor;

    unsigned long long tcount = 0;
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

	// try to trim
	corseq = r->trim(trimq);
	if(!r->untrusted.empty()) {
	  // if still untrusted, correct
	  corseq = r->correct(trusted, ntnt_prob, true);

	  // if trimmed to long enough
	  if(corseq.size() >= trim_t) {
	    if(r->trusted_read != 0) { // else no guarantee there was a correction
	      for(int c = 0; c < r->trusted_read->corrections.size(); c++) {
		correction cor = r->trusted_read->corrections[c];
		if(iseq[cor.index] < 4) {
		  // real -> observed, (real approximated by correction)
		  //ntnt_counts[cor.to][iseq[cor.index]]++;
		  
		  // hmm well its best to compute P(real|observed)
		  ntnt_counts[iseq[cor.index]][cor.to]++;
		  samples++;
		}
	      }
	    }
	  }
	}
	delete r;
      }

      if(++tcount == counts[tid] || samples > 50000)
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
// unzip_fastq
//
// Copy the fastq file to the current directory and unzip it
////////////////////////////////////////////////////////////
void unzip_fastq(const char* fqf) {
  char mycmd[100];

  // cp to cwd
  strcpy(mycmd, "cp ");
  strcat(mycmd, zipd);
  strcat(mycmd, "/");
  strcat(mycmd, fqf);
  strcat(mycmd, " .");
  system(mycmd);

  // gunzip
  strcpy(mycmd, "gunzip ");
  strcat(mycmd, fqf);
  system(mycmd);
}

////////////////////////////////////////////////////////////
// zip_fastq
//
// Remove the fastq file, zip the corrected fastq file, and
// move it to zipd
////////////////////////////////////////////////////////////
void zip_fastq(const char* fqf) {
  char mycmd[100];

  // rm fqf
  strcpy(mycmd, "rm ");
  strcat(mycmd, fqf);
  system(mycmd);

  // determine output file
  string fqf_str(fqf);
  int suffix_index = fqf_str.rfind(".");
  string prefix = fqf_str.substr(0,suffix_index);
  string suffix = fqf_str.substr(suffix_index, fqf_str.size()-suffix_index);
  string outf = prefix + string(".cor") + suffix;

  // zip
  strcpy(mycmd, "gzip ");
  strcat(mycmd, outf.c_str());
  system(mycmd);

  // move
  strcpy(mycmd, "mv ");
  strcat(mycmd, outf.c_str());
  strcat(mycmd, ".gz ");
  strcat(mycmd, zipd);
  system(mycmd);
}

////////////////////////////////////////////////////////////
// is_kmer_file
//
// Return true if the file is in the format "kmer\tcount",
// else its probably a binary bithash dump file.
////////////////////////////////////////////////////////////
bool is_kmer_file(char* f) {
  ifstream f_in(f);
  string line;
  getline(f_in, line);

  char nt;
  for(int i = 0; i < 5; i++) {
    nt = line[0];
    if(!(nt == 'A' || nt == 'C' || nt == 'G' || nt == 'T')) {
      cout << "Reading -m file as bithash binary" << endl;
      return false;
    }
  }
  return true;
}

////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  // set up parallelism
  omp_set_num_threads(threads);

  // make trusted kmer data structure
  bithash *trusted = new bithash();
  if(is_kmer_file(merf)) {
    if(cutoff == 0 && ATcutf == NULL) {
      cerr << "Must provide a trusted/untrusted kmer cutoff (-c) or a file containing the cutoff as a function of the AT content (-a)" << endl;
      exit(EXIT_FAILURE);
    }

    if(ATcutf != NULL) {
      if(strcmp(merf,"-") == 0)
	trusted->tab_file_load(cin, load_AT_cutoffs());
      else {
	ifstream mer_in(merf);
	trusted->tab_file_load(mer_in, load_AT_cutoffs());
      }
    } else {
      if(strcmp(merf,"-") == 0) {
	trusted->tab_file_load(cin, cutoff);
      } else {
	ifstream mer_in(merf);
	trusted->tab_file_load(mer_in, cutoff);
      }
    }
  } else {
    // from from file
    trusted->binary_file_input_lowmem(merf);
    cout << trusted->num_kmers() << " trusted kmers" << endl;
  }

  // make list of files
  vector<string> fastqfs;
  if(file_of_fastqf != NULL) {
    ifstream ff(file_of_fastqf);
    string next_fastqf;

    while(getline(ff, next_fastqf))
      fastqfs.push_back(next_fastqf);
  } else
    fastqfs.push_back(string(fastqf));

  // process each
  string fqf;
  for(int f = 0; f < fastqfs.size(); f++) {
    fqf = fastqfs[f];
    cout << fqf << endl;

    if(zipd != NULL) {
      unzip_fastq(fqf.c_str());
      fqf = fqf.substr(0, fqf.size()-3);
    }

    // split file
    vector<streampos> starts;
    vector<unsigned long long> counts;
    pa_params(fqf, starts, counts);

    // learn nt->nt transitions
    double ntnt_prob[4][4] = {0};
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
	if(i != j)
	  ntnt_prob[i][j] = 1.0/3.0;
      }
    }
    learn_errors(fqf, trusted, starts, counts, ntnt_prob);

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

    correct_reads(fqf, trusted, starts, counts, ntnt_prob);

    if(zipd != NULL)
      zip_fastq(fqf.c_str());
  }

  return 0;
}
