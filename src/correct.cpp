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
#include <sys/stat.h>

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////
const static char* myopts = "r:f:m:b:c:a:t:q:p:z:ICuh";
static struct option  long_options [] = {
  {"headers", 0, 0, 1000},
  {0, 0, 0, 0}
};

// -r, fastq file of reads
static char* fastqf = NULL;
// -f, file of fastq files of reads
static char* file_of_fastqf = NULL;

// -m, mer counts
static char* merf = NULL;
// -b, bithash
static char* bithashf = NULL;

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

// --headers, Print only normal headers
static bool orig_headers = false;

// -C, Contrail output
static bool contrail_out = false;
// -u, output uncorrected reads
static bool uncorrected_out = false;

static unsigned int chunks_per_thread = 200;

// Note: to not trim, set trimq=0 and trim_t>read_length-k

// constants
#define TESTING false
static const char* nts = "ACGTN";
static const unsigned int max_qual = 50;

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
	   "    File containing kmer counts in format `seq\tcount`.\n"
	   "    Can also be piped in with '-'\n"
	   " -b <file>\n"
	   "    File containing saved bithash.\n"
	   " -c <num>\n"
	   "    Separate trusted/untrusted kmers at cutoff <num>\n"
	   " -a <file>\n"
	   "    Separate trusted/untrusted kmers as a function of AT\n"
	   "    content, with cutoffs found in <file>, one per line\n"
	   " -p <num>\n"
	   "    Use <num> openMP threads\n"
	   " -t <num>=30\n"
	   "    Return only reads trimmed to >= <num> bp\n"
	   " -q <num>=3\n"
	   "    Use BWA trim parameter <num>\n"
	   " -I\n"
	   "    Use 64 scale Illumina quality values (else base 33)\n"
	   " -u\n"
	   "    Output errors reads even if they can't be corrected\n"
	   " --headers\n"
	   "    Output only the original read headers without\n"
	   "    correction messages\n"
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
  int option_index = 0;
  char* p;
  
  // parse args
  while(!errflg && ((ch = getopt_long(argc, argv, myopts, long_options, &option_index)) != EOF)) {
  //while(!errflg && ((ch = getopt(argc, argv, myopts)) != EOF)) {
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

    case 'b':
      bithashf = strdup(optarg);
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

    case 'u':
      uncorrected_out = true;
      break;  

    case 'p':
      threads = int(strtol(optarg, &p, 10));
      if(p == optarg || threads <= 0) {
	fprintf(stderr, "Bad number of threads \"%s\"\n",optarg);
	errflg = true;
      }
      break;

    case 1000:
      orig_headers = true;
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

  if(merf != NULL) {
    if(cutoff == 0 && ATcutf == NULL) {
      cerr << "Must provide a trusted/untrusted kmer cutoff (-c) or a file containing the cutoff as a function of the AT content (-a)" << endl;
      exit(EXIT_FAILURE);
    }
  } else if(bithashf == NULL) {
    cerr << "Must provide a file of kmer counts (-m) or a saved bithash (-b)" << endl;
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

  /*
  char mycmd[30000];

  // cat thread output
  strcpy(mycmd, "cat ");  
  char strt[10];
  for(int t = 0; t < threads; t++) {
    strcat(mycmd, outf);
    sprintf(strt,"%d",t);
    strcat(mycmd, strt);
    strcat(mycmd, " ");
  }
  // into single file
  strcat(mycmd, "> ");
  strcat(mycmd, outf);
    
  // rm thread output
  strcat(mycmd, "; rm ");
  for(int t = 0; t < threads; t++) {
    strcat(mycmd, outf);
    sprintf(strt,"%d",t);
    strcat(mycmd, strt);
    strcat(mycmd, " ");
  }

  if(!final)
    // run in background
    strcat(mycmd, "&");

  // execute
  cout << mycmd << endl;
  system(mycmd);
  */
}


////////////////////////////////////////////////////////////
// regress_probs
//
// Use ntnt_counts to perform nonparametric regression
// on ntnt_prob across quality values.
////////////////////////////////////////////////////////////
void regress_probs(double ntnt_prob[max_qual][4][4], unsigned int ntnt_counts[max_qual][4][4]) {
  double sigma = 3.0;
  double sigma2 = pow(sigma, 2);
  
  // count # occurrences for each (quality=q,actual=a) tuple
  unsigned int actual_counts[max_qual][4] = {0};
  for(int q = 1; q < max_qual; q++)
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
	actual_counts[q][i] += ntnt_counts[q][i][j];

  // regress
  double ntdsum;
  for(int q = 1; q < max_qual; q++) {
    for(int i = 0; i < 4; i++) {
      //ntdsum = 0;
      for(int j = 0; j < 4; j++) {
	double pnum = 0;
	double pden = 0;
	for(int qr = 1; qr < max_qual; qr++) {
	  pnum += ntnt_counts[qr][i][j] * exp(-pow((double)(qr - q), 2)/(2*sigma2));
	  pden += actual_counts[qr][i] * exp(-pow((double)(qr - q), 2)/(2*sigma2));
	}
	ntnt_prob[q][i][j] = pnum / pden;
	//ntdsum += ntnt_prob[q][i][j];
      }

      // re-normalize to sum to 1
      //for(int j = 0; j < 4; j++)
      //ntnt_prob[q][i][j] /= ntdsum;
    }
  }
}


////////////////////////////////////////////////////////////
// output_model
//
// Print the error model to the file error_model.txt
////////////////////////////////////////////////////////////
void output_model(double ntnt_prob[max_qual][4][4], unsigned int ntnt_counts[max_qual][4][4]) {
  ofstream mod_out("error_model.txt");

  unsigned int ntsum;
  for(int q = 1; q < max_qual; q++) {
    mod_out << "Quality = " << q << endl;

    // counts
    mod_out << "\tA\tC\tG\tT" << endl;
    for(int i = 0; i < 4; i++) {
      mod_out << nts[i];

      ntsum = 0;
      for(int j = 0; j < 4; j++)
	ntsum += ntnt_counts[q][i][j];

      for(int j = 0; j < 4; j++) {
	if(i == j)
	  mod_out << "\t-";
	else if(ntsum > 0)
	  mod_out << "\t" << ((double)ntnt_counts[q][i][j] / (double)ntsum) << "(" << ntnt_counts[q][i][j] << ")";
	else
	  mod_out << "\t0";
      }
      mod_out << endl;
    }

    // probs
    mod_out << "\tA\tC\tG\tT" << endl;
    for(int i = 0; i < 4; i++) {
      mod_out << nts[i];
      for(int j = 0; j < 4; j++) {
	if(i == j)
	  mod_out << "\t-";
	else
	  mod_out << "\t" << ntnt_prob[q][i][j];
      }
      mod_out << endl;
    }
    mod_out << endl;    
  }
}


////////////////////////////////////////////////////////////
// output_read
//
// Output the given possibly corrected and/or trimmed
// read according to the given options.
////////////////////////////////////////////////////////////
static void output_read(ofstream & reads_out, string header, string ntseq, string mid, string strqual, string corseq) {
  if(corseq.size() >= trim_t) {
    // check for changes
    bool corrected = false;
    for(int i = 0; i < corseq.size(); i++) {
      if(corseq[i] != ntseq[i]) {
	corrected = true;
	// set qual to crap
	if(Read::illumina_qual)
	  strqual[i] = 'B';
	else
	  strqual[i] = '#';
      }
    }
    // update header
    if(!orig_headers) {
      if(corrected)
	header += " correct";
      unsigned int trimlen = ntseq.size()-corseq.size();
      if(trimlen > 0) {
	stringstream trim_inter;
	trim_inter << trimlen;
	header += " trim=" + trim_inter.str();
      }
    }
    // print
    if(contrail_out)
      reads_out << header << "\t" << corseq << endl;
    else
      reads_out << header << endl << corseq << endl << mid << endl << strqual.substr(0,corseq.size()) << endl;
    if(TESTING)
      cerr << header << "\t" << ntseq << "\t" << corseq << endl;
  } else {
    if(uncorrected_out) {
      // update header
      if(!orig_headers)
	header += " error";
      //print
      if(contrail_out)
	reads_out << header << "\t" << ntseq << endl;
      else
	reads_out << header << endl << ntseq << endl << mid << endl << strqual << endl;	  
    }
    if(TESTING)
      cerr << header << "\t" << ntseq << "\t-" << endl; // or . if it's only trimmed?
  }
}


////////////////////////////////////////////////////////////
// correct_reads
////////////////////////////////////////////////////////////
//static void correct_reads(string fqf, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double (&ntnt_prob)[4][4], double prior_prob[4]) {
static void correct_reads(string fqf, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double ntnt_prob[max_qual][4][4], double prior_prob[4]) {
  // format output file
  int suffix_index = fqf.rfind(".");
  string prefix = fqf.substr(0,suffix_index);
  string suffix = fqf.substr(suffix_index, fqf.size()-suffix_index);
  string outf = prefix + string(".cor") + suffix;
  //cout << outf << endl;

  unsigned int chunk = 0;
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
    ifstream reads_in(fqf.c_str());
    
    unsigned int tchunk;
    string header,ntseq,mid,strqual,corseq;
    char* nti;
    Read *r;

    while(chunk < starts.size()) {
      #pragma omp critical
      tchunk = chunk++;

      reads_in.seekg(starts[tchunk]);

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

	  // if still untrusted, correct
	  if(!r->untrusted.empty())
	    corseq = r->correct(trusted, ntnt_prob, prior_prob);

	  // output read w/ trim and corrections
	  output_read(reads_out, header, ntseq, mid, strqual, corseq);
	  
	  delete r;
	} else {
	  // output read as is
	  if(contrail_out)
	    reads_out << header << "\t" << ntseq << endl;
	  else
	    reads_out << header << endl << ntseq << endl << mid << endl << strqual << endl;
	}
	
	if(++tcount == counts[tchunk])
	  break;
      }
    }
    reads_in.close();
    reads_out.close();
  }
  
  combine_output(outf.c_str());
}


////////////////////////////////////////////////////////////
// learn_errors
//
// Correct reads using a much stricter filter in order
// to count the nt->nt errors and learn the errors
// probabilities
////////////////////////////////////////////////////////////
//static void learn_errors(string fqf, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double (&ntnt_prob)[4][4], double prior_prob[4]) {
static void learn_errors(string fqf, bithash * trusted, vector<streampos> & starts, vector<unsigned long long> & counts, double ntnt_prob[max_qual][4][4], double prior_prob[4]) {
  unsigned int ntnt_counts[max_qual][4][4] = {0};
  unsigned int samples = 0;

  unsigned int chunk = 0;
#pragma omp parallel //shared(trusted)
  {    
    unsigned int tchunk;
    string header,ntseq,strqual,corseq;
    char* nti;
    Read *r;    
    ifstream reads_in(fqf.c_str());
    
    while(chunk < threads*chunks_per_thread) {
#pragma omp critical
      tchunk = chunk++;     
      
      reads_in.seekg(starts[tchunk]);
      
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

	vector<int> iqual;
	for(int i = 0; i < strqual.size(); i++) {
	  if(Read::illumina_qual)
	    iqual.push_back(strqual[i]-64);
	  else
	    iqual.push_back(strqual[i]-33);
	}
	
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
	    corseq = r->correct(trusted, ntnt_prob, prior_prob, true);
	    
	    // if trimmed to long enough
	    if(corseq.size() >= trim_t) {
	      if(r->trusted_read != 0) { // else no guarantee there was a correction
		for(int c = 0; c < r->trusted_read->corrections.size(); c++) {
		  correction cor = r->trusted_read->corrections[c];
		  if(iseq[cor.index] < 4) {
		    // P(obs=o|actual=a,a!=o) for Bayes
		    ntnt_counts[iqual[cor.index]][cor.to][iseq[cor.index]]++;
		    
		    // P(actual=a|obs=o,a!=o)
		    //ntnt_counts[iseq[cor.index]][cor.to]++;
		    samples++;
		  }
		}
	      }
	    }
	  }
	  delete r;
	}
	
	if(++tcount == counts[tchunk] || samples > 200000)
	  break;
      }
    }
    reads_in.close();
  }

  regress_probs(ntnt_prob, ntnt_counts);

  output_model(ntnt_prob, ntnt_counts);
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

  // zcat to cwd
  strcpy(mycmd, "zcat ");
  strcat(mycmd, zipd);
  strcat(mycmd, "/");
  strcat(mycmd, fqf);
  strcat(mycmd, " > ");
  strncat(mycmd, fqf, strstr(fqf, ".gz") - fqf);
  system(mycmd);

  /*
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
  */
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
// main
////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  // count AT's and GC's
  unsigned long long atgc[2] = {0};

  // make trusted kmer data structure
  bithash *trusted = new bithash();

  // kmer counts
  if(merf != NULL) {
    if(ATcutf != NULL) {
      if(strcmp(merf,"-") == 0)
	trusted->tab_file_load(cin, load_AT_cutoffs(), atgc);
      else {
	ifstream mer_in(merf);
	trusted->tab_file_load(mer_in, load_AT_cutoffs(), atgc);
      }
    } else {
      if(strcmp(merf,"-") == 0) {
	trusted->tab_file_load(cin, cutoff, atgc);
      } else {
	ifstream mer_in(merf);
	trusted->tab_file_load(mer_in, cutoff, atgc);
      }
    }

  // saved bithash
  } else if(bithashf != NULL) {
    if(strcmp(bithashf,"-") == 0) {
      cerr << "Saved bithash cannot be piped in.  Please specify file." << endl;
      exit(EXIT_FAILURE);
    } else
      trusted->binary_file_input(bithashf, atgc);
  }
  
  cout << trusted->num_kmers() << " trusted kmers" << endl;

  double prior_prob[4];
  prior_prob[0] = (double)atgc[0] / (double)(atgc[0]+atgc[1]) / 2.0;
  prior_prob[1] = .5 - prior_prob[0];
  prior_prob[2] = prior_prob[1];
  prior_prob[3] = prior_prob[0];
  
  cout << "AT: " << atgc[0] << " GC: " << atgc[1] << endl;
  cout << "AT% = " << (2*prior_prob[0]) << endl;

  // make list of files
  vector<string> fastqfs;
  if(file_of_fastqf != NULL) {
    ifstream ff(file_of_fastqf);
    string next_fastqf;

    while(getline(ff, next_fastqf) && next_fastqf.size() > 0)
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
    double ntnt_prob[max_qual][4][4] = {0};
    for(int q = 0; q < max_qual; q++)
      for(int i = 0; i < 4; i++)
	for(int j = 0; j < 4; j++)
	  if(i != j)
	    ntnt_prob[q][i][j] = 1.0/3.0;

    learn_errors(fqf, trusted, starts, counts, ntnt_prob, prior_prob);

    correct_reads(fqf, trusted, starts, counts, ntnt_prob, prior_prob);

    if(zipd != NULL)
      zip_fastq(fqf.c_str());
  }

  return 0;
}
