#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include <ext/hash_map>
namespace Sgi = ::__gnu_cxx;         // GCC 4.0 and later
#define HASHMAP __gnu_cxx
using namespace::std;

////////////////////////////////////////////////////////////
// hash_map
////////////////////////////////////////////////////////////
struct eqstr {
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1,s2) == 0;
  }
};

namespace HASHMAP                                                                                 
{                                                                                             
  template<> struct hash< std::string >                                                       
  {                                                                                           
    size_t operator()( const std::string& x ) const                                           
    {                                                                                         
      return HASHMAP::hash< const char* >()( x.c_str() );                                              
    }                                                                                         
  };
}          

//typedef HASHMAP::hash_map<const char*, const char*, HASHMAP::hash<const char*>, eqstr> seq_hash; 
//typedef HASHMAP::hash_map<const char*, const char*, HASHMAP::hash<const char*>, equal_to<string> > seq_hash; 
typedef HASHMAP::hash_map<string, string> seq_hash; 

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////
const static char* myopts = "r:c:IC";
// -r, fastq file of reads
static char* fastqf = NULL;
// -c, fastq file of corrected reads
static char* corf = NULL;
// -I
static bool illumina_qual = false;
// -C, Contrail output
static bool contrail_out = false;


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

    case 'c':
      corf = strdup(optarg);
      break;

    case 'I':
      illumina_qual = true;
      break;

    case 'C':
      contrail_out = true;
      break;

    case  '?' :
      fprintf (stderr, "Unrecognized option -%c\n", optopt);

    default:
      errflg = true;
    }
  }

  if(fastqf == NULL) {
    cerr << "Must provide original read fastq file with -r" << endl;
    exit(EXIT_FAILURE);

  } else if(corf == NULL) {
    // infer correction file
    string fqf_str(fastqf);
    unsigned int suffix_index = fqf_str.rfind(".");
    corf = new char[fqf_str.size()+5];
    strncpy(corf, fastqf, suffix_index);
    strcat(corf, ".cor");
    strncat(corf, &fastqf[suffix_index], fqf_str.size()-suffix_index);
    cout << "Correction file assumed to be " << corf << endl;
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
  string header, seq, mid, qual;
  const char * nts = "ACGT";

  unsigned long oread_count = 0;
  unsigned long cread_count = 0;
  unsigned long correct_bp = 0;
  unsigned long trim_bp = 0;
  unsigned long ntnt_counts[4][4] = {0};
  vector<unsigned long> pos_corrections;
  vector<unsigned long> pos_trims;  // to do this right, must hash on trim tag

  parse_command_line(argc, argv);

  ////////////////////////////////////////
  // hash corrections
  ////////////////////////////////////////
  seq_hash corrected_reads;
  ifstream correctionsf(corf);
  unsigned int ci;
  const char* crhead;
  while(getline(correctionsf, header)) {
    if(contrail_out) {
      unsigned int line_tab = header.find('\t');
      seq = header.substr(line_tab+1, header.size()-1-line_tab);
      //qual = garbage...
      header = header.substr(0, line_tab);
    } else {
      getline(correctionsf, seq);
      getline(correctionsf, mid);
      getline(correctionsf, qual);
    }

    cread_count++;

    ci = header.find("correct");
    //cout << header << " " << ci << endl;
    if(ci != -1)
      corrected_reads.insert(make_pair(header.substr(0,ci-1), seq));
  }
  correctionsf.close();

  ////////////////////////////////////////
  // parse original file
  ////////////////////////////////////////
  ifstream originalf(fastqf);
  seq_hash::iterator fi;
  string cseq;
  while(getline(originalf, header)) {
    getline(originalf, seq);
    getline(originalf, mid);
    getline(originalf, qual);

    oread_count++;

    // set up positional counts
    if(pos_corrections.empty()) {
      for(int p = 0; p < seq.size(); p++)
	pos_corrections.push_back(0);
    }
    if(pos_trims.empty()) {
      for(int p = 0; p < seq.size(); p++)
	pos_trims.push_back(0);
    }

    // if corrected
    fi = corrected_reads.find(header);
    if(fi != corrected_reads.end()) {
      cseq = fi->second;      

      int lastbp = correct_bp;

      int i;
      for(i = 0; i < cseq.size(); i++) {
	// if correction
	if(seq[i] != cseq[i]) {
	  correct_bp++;

	  // count nt->nt
	  if(seq[i] != 'N')
	    ntnt_counts[strchr(nts, cseq[i]) - nts][strchr(nts, seq[i]) - nts]++;

	  // count position - corrected and trimmed
	  pos_corrections[i]++;
	}
      }
      // count trims
      trim_bp += seq.size()-i;
      for(; i < seq.size(); i++)
	pos_trims[i]++;

      if(correct_bp == lastbp)
	cout << "No corrections for " << header << endl;
    }
  }

  ////////////////////////////////////////
  // print stats
  ////////////////////////////////////////
  cout << "Total reads: " << oread_count << endl << endl;
  //cout << "Validated reads: " << (cread_count - ... hash trims to get this right  
  cout << "Tossed reads: " << (oread_count-cread_count) << endl;
  cout << "Corrected reads: " << corrected_reads.size() << endl;
  cout << "Corrected bp: " << correct_bp << endl;
  cout << "Trimmed bp: " << trim_bp << " (incomplete)" << endl;

  cout << "nt->nt error rate:" << endl;
  cout << "\tA\tC\tG\tT" << endl;
  for(int i = 0; i < 4; i++) {
    cout << nts[i];
    for(int j = 0; j < 4; j++)
      cout << "\t" << ntnt_counts[i][j];
    cout << endl;
  }

  cout << "errors by position:" << endl;
  for(int i = 0; i < pos_corrections.size(); i++)
    cout << (i+1) << " " << pos_corrections[i] << endl;

  /* hash trims to get this right
  cout << "trims by position:" << endl;
  for(int i = 0; i < pos_trims.size(); i++)
    cout << (i+1) << " " << pos_trims[i] << endl;
  */

  return 0;
}
