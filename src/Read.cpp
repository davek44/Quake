#include "Read.h"
#include "bithash.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <set>
#include <queue>

#define TESTING false
#define DEBUG false

////////////////////////////////////////////////////////////
// corrections_compare
//
// Simple class to compare to corrected_read's in the
// priority queue
////////////////////////////////////////////////////////////
class corrections_compare {
public:
  //  corrections_compare() {};
  bool operator() (const corrected_read* lhs, const corrected_read* rhs) const {
    //return lhs->likelihood < rhs->likelihood;
    if(lhs->likelihood < rhs->likelihood)
      return true;
    else if(lhs->likelihood > rhs->likelihood)
      return false;
    else
      return lhs->region_edits > rhs->region_edits;
  }
};

////////////////////////////////////////////////////////////
// Read (constructor)
//
// Make shallow copies of sequence and untrusted, and
// convert quality value string to array of probabilities
////////////////////////////////////////////////////////////
Read::Read(const string & h, const unsigned int* s, const string & q, vector<int> & u, const int rl)
  :untrusted(u) {

  header = h;
  read_length = rl;
  trim_length = rl;
  seq = new unsigned int[read_length];
  prob = new float[read_length];
  for(int i = 0; i < read_length; i++) {
    seq[i] = s[i];    
    // quality values of 0,1 lead to p < .25
    if(illumina_qual)
      prob[i] = max(.25, 1.0-pow(10.0,-(q[i]-64)/10.0));
    else
      prob[i] = max(.25, 1.0-pow(10.0,-(q[i]-33)/10.0));
  }
  trusted_read = 0;
  global_like = 1.0;
}

Read::~Read() {
  delete[] seq;
  delete[] prob;
  if(trusted_read != 0)
    delete trusted_read;
}

////////////////////////////////////////////////////////////
// trim
//
// Trim the end of the read the way BWA does it.
// Returns the trimmed read as a string
////////////////////////////////////////////////////////////
string Read::trim(int t) {
  // find trim index
  int phredq;
  int current_trimfunc = 0;
  int max_trimfunc = 0;
  trim_length = read_length; // already set in constructor but ok
  for(int i = read_length-1; i >= 0; i--) {
    phredq = floor(.5-10*log(1.0 - prob[i])/log(10));
    current_trimfunc += (t - phredq);
    if(current_trimfunc > max_trimfunc) {
      max_trimfunc = current_trimfunc;
      trim_length = i;
    }
  }

  // update untrusted
  for(int i = untrusted.size()-1; i >= 0; i--) {
    if(untrusted[i] >= trim_length)
      untrusted.pop_back();
  }

  vector<correction> no_cors;
  return print_corrected(no_cors);
}


////////////////////////////////////////////////////////////
// single_correct
//
// Find the set of corrections with maximum likelihood
// that result in all trusted kmers.
//
// Assumes a short read so obsolete.
////////////////////////////////////////////////////////////
/*
bool Read::single_correct(bithash *trusted, ofstream & out, double (&ntnt_prob)[4][4], bool learning) {
  if(correct_subset(untrusted, trusted, out, ntnt_prob, learning)) {
    out << header << "\t" << print_seq() << "\t" << print_corrected(trusted_read->corrections) << endl;
    return true;
  } else
    return false;
}
*/

////////////////////////////////////////////////////////////
// correct_subset
//
// Find the set of corrections with maximum likelihood that
// result in all trusted kmers in the region defined by
// the given untrusted kmer indices.
//
// Will print output if the read will not be corrected,
// but otherwise abstains.
//
// Corrections can be accessed through 'trusted_read'
////////////////////////////////////////////////////////////
bool Read::correct_subset(vector<int> untrusted_subset, bithash *trusted, double (&ntnt_prob)[4][4], bool learning) {  
  /*
  if(header == "@19732") {
    cout << "Untrusted: " << untrusted_subset.size() << endl;
    for(int i = 0; i < untrusted_subset.size(); i++)
      cout << untrusted_subset[i] << " ";
    cout << endl;
  }
  */

  // determine region to consider
  vector<short> region = error_region(untrusted_subset);   

  // filter really bad reads  
  float avgprob = 0;
  for(int i = 0; i < region.size(); i++)
    avgprob += prob[region[i]];
  avgprob /= region.size();
  
  int badnt = 0;
  for(int i = 0; i < region.size(); i++)
    if(prob[region[i]] < .95)
      badnt += 1;

  bool forfeit_easily = false;
  if(learning) {
    if(badnt >= 4) {
      //out << header << "\t" << print_seq() << "\t." << endl;
      return false;
    }
  }
  else if(badnt >= 8) {
    forfeit_easily = true;
  }

  // sort by quality
  quality_quicksort(region, 0, region.size()-1);

  // data structure for corrected_reads sorted by likelihood
  //priority_queue< corrected_read*, vector<corrected_read*>, corrections_compare > cpq;
  vector<corrected_read*> cpq;
  corrections_compare cpq_comp;

  int max_pq = 0;

  // add initial reads
  corrected_read *cr, *next_cr;
  short edit_i;
  float like;
  bitset<bitsize> bituntrusted;
  for(int i = 0; i < untrusted_subset.size(); i++)
    bituntrusted.set(untrusted_subset[i]);

  bool cr_added = true;  // once an iteration passes w/ no corrected reads added, we can stop
  for(short region_edit = 0; region_edit < region.size() && cr_added; region_edit++) {
    edit_i = region[region_edit];
    cr_added = false;

    for(short nt = 0; nt < 4; nt++) {
      if(seq[edit_i] != nt) {
	//like = (1.0-prob[edit_i]) * ntnt_prob[nt][seq[edit_i]] / prob[edit_i];
	like = (1.0-prob[edit_i]) * ntnt_prob[seq[edit_i]][nt] / prob[edit_i];
	
	next_cr = new corrected_read(bituntrusted, like, region_edit+1);
	next_cr->corrections.push_back(correction(edit_i, nt));
      
	// add to priority queue
	//cpq.push(next_cr);
	cpq.push_back(next_cr);
	push_heap(cpq.begin(), cpq.end(), cpq_comp);
	cr_added = true;
      }
    }
  }

  // initialize likelihood parameters
  trusted_read = 0;
  float trusted_likelihood;
  int untrusted_count;

  ////////////////////////////////////////
  // process corrected reads
  ////////////////////////////////////////
  while(cpq.size() > 0) {
    if(cpq.size() > max_pq)
      max_pq = cpq.size();
    
    // give up on read if it was marked and no trusted reads have been found by 30k
    if(forfeit_easily && trusted_read == 0 && cpq.size() > 30000) {
      if(TESTING)
	cerr << header << "\t" << print_seq() << "\t." << endl;
      break;
    }
    
    // give up on read if pq is too big
    if(cpq.size() > 400000) {
      if(DEBUG)
	cout << "queue is too large for " << header << endl;
      if(trusted_read != 0) {
	delete trusted_read;
	trusted_read = 0;
      }
      break;
    }

    //cr = cpq.top();
    cr = cpq[0];
    //cpq.pop();
    pop_heap(cpq.begin(), cpq.end(), cpq_comp);
    cpq.pop_back();

    // save for later comparison
    untrusted_count = (signed int)cr->untrusted.count();
    
    if(trusted_read != 0) {
      // if a corrected read exists, compare likelihoods and if likelihood is too low, break loop return true
      if(cr->likelihood < trusted_likelihood*trust_spread_t) {
	delete cr;
	break;
      }
    } else {
      // if no corrected read exists and likelihood is too low, break loop return false
      if(global_like*cr->likelihood < correct_min_t) {
	delete cr;
	break;
      }
    }
    
    // check for all kmers trusted
    if(check_trust(cr, trusted)) {
      if(trusted_read == 0) {
	// if yes, and first trusted read, save
	trusted_read = cr;
	trusted_likelihood = cr->likelihood;
      } else {
	// output ambiguous corrections for testing
	if(TESTING)
	  cerr << header << "\t" << print_seq() << "\t" << print_corrected(trusted_read->corrections) << "\t" << print_corrected(cr->corrections) << endl;
	
	// if yes, and if trusted read exists delete trusted_read, break loop
	delete trusted_read;
	delete cr;
	trusted_read = 0;
	break;
      }
    }
    
    /*
    if(header == "@19732") {
      cout << cr->likelihood << "\t";
      for(int c = 0; c < cr->corrections.size(); c++) {
	cout << " (" << cr->corrections[c].index << "," << cr->corrections[c].to << ")";
      }
      cout << "\t";
      for(int c = 0; c < read_length; c++) {
	if(cr->untrusted[c])
	  cout << 1;
	else
	  cout << 0;
      }
      cout << endl;
    }
    */
    
    // if untrusted sharply increases, just bail
    if(((signed int)cr->untrusted.count() - untrusted_count) < k/3) {    

      // determine next nt flips
      bool cr_added = true;  // once an iteration passes w/ no corrected reads added, we can stop
      for(short region_edit = cr->region_edits; region_edit < region.size() && cr_added; region_edit++) {
	edit_i = region[region_edit];
	cr_added = false;
	
	// add relatives
	for(short nt = 0; nt < 4; nt++) {
	  // if actual edit, 
	  if(seq[edit_i] != nt) {
	    // calculate new likelihood
	    // error is real->observed, (with real approximated by correction)
	    //like = cr->likelihood * (1.0-prob[edit_i]) * ntnt_prob[nt][seq[edit_i]] / prob[edit_i];
	    like = cr->likelihood * (1.0-prob[edit_i]) * ntnt_prob[seq[edit_i]][nt] / prob[edit_i];
	    
	    // if thresholds ok, add new correction
	    if(trusted_read != 0) {
	      if(like < trusted_likelihood*trust_spread_t)
		continue;
	    } else {
	      // must consider spread or risk missing a case of ambiguity
	      if(global_like*like < correct_min_t*trust_spread_t)
		continue;
	    }
	    
	    next_cr = new corrected_read(cr->corrections, cr->untrusted, like, region_edit+1);
	    next_cr->corrections.push_back(correction(edit_i, nt));
	  
	    // add to priority queue
	    //cpq.push(next_cr);
	    cpq.push_back(next_cr);
	    push_heap(cpq.begin(), cpq.end(), cpq_comp);
	    cr_added = true;
	  }
	}
      }
    }

    // if not the saved max trusted, delete
    if(trusted_read != cr) {
      delete cr;
    }
  }

  // delete memory from cpq
  /*
  while(cpq.size() > 0) {
    cr = cpq.top();
    cpq.pop();
    delete cr;
  }
  */

  for(int i = 0; i < cpq.size(); i++)
    delete cpq[i];

  if(trusted_read != 0)
    return true;
  else
    return false;
}

////////////////////////////////////////////////////////////
// print_seq
////////////////////////////////////////////////////////////
string Read::print_seq() {
  char nts[5] = {'A','C','G','T','N'};
  string sseq;
  for(int i = 0; i < read_length; i++)
    sseq.push_back(nts[seq[i]]);
  return sseq;
}

////////////////////////////////////////////////////////////
// print_corrected
//
// Print read with corrections and trimming.
////////////////////////////////////////////////////////////
string Read::print_corrected(vector<correction> & cor) {
  return print_corrected(cor, trim_length);
}
string Read::print_corrected(vector<correction> & cor, int print_nt) {
  char nts[5] = {'A','C','G','T','N'};
  string sseq;
  int correct_i;
  for(int i = 0; i < print_nt; i++) {
    correct_i = -1;
    for(int c = 0; c < cor.size(); c++) {
      if(cor[c].index == i)
	correct_i = c;
    }
    if(correct_i != -1)
      sseq.push_back(nts[cor[correct_i].to]);
    else
      sseq.push_back(nts[seq[i]]);
  }
  return sseq;
}
  

////////////////////////////////////////////////////////////
// correct
//
// Perform correction by breaking up untrusted kmers
// into connected components and correcting them
// independently.
////////////////////////////////////////////////////////////
string Read::correct(bithash *trusted, double (&ntnt_prob)[4][4], bool learning) {
  // current connected component
  vector<int> cc_untrusted;
  cc_untrusted.push_back(untrusted[0]);
  
  vector<correction> multi_cors;
  for(int i = 1; i < untrusted.size(); i++) {

    // if kmer from last untrusted doesn't reach next
    if(untrusted[i-1]+k-1 < untrusted[i]) {

      // correct current component
      if(correct_subset(cc_untrusted, trusted, ntnt_prob, learning)) {
	global_like *= trusted_read->likelihood;

	// store corrections
	for(int c = 0; c < trusted_read->corrections.size(); c++)
	  multi_cors.push_back(trusted_read->corrections[c]);

	// start new component
	cc_untrusted.clear();
      } else {
	// cannot correct
	if(cc_untrusted[0]-k+1 >= 1) {
	  // but can trim
	  return print_corrected(multi_cors, cc_untrusted[0]-k+1);
	} else
	  return string("");
      }
    }
    // add to current component
    cc_untrusted.push_back(untrusted[i]);
  }

  // finish last component
  if(correct_subset(cc_untrusted, trusted, ntnt_prob, learning)) {
    global_like *= trusted_read->likelihood;

    // store corrections
    for(int c = 0; c < trusted_read->corrections.size(); c++)
      multi_cors.push_back(trusted_read->corrections[c]);
    
  } else {
    // cannot correct
    if(cc_untrusted[0]-k+1 >= 1) {
      // but can trim
      return print_corrected(multi_cors, cc_untrusted[0]-k+1);
    } else
      return string("");
  }

  // create new trusted read (mostly for learn_errors)
  corrected_read * tmp = trusted_read;
  trusted_read = new corrected_read(multi_cors, tmp->untrusted, global_like, 0);
  delete tmp;

  // print read with all corrections
  return print_corrected(multi_cors);
}


////////////////////////////////////////////////////////////
// error_region
//
// Find region of the read to consider for errors based
// on the pattern of untrusted kmers
////////////////////////////////////////////////////////////
vector<short> Read::error_region(vector<int> untrusted_subset) {
  // find read indexes to consider
  vector<short> region;
  if(untrusted_intersect(untrusted_subset, region)) {
    // if front kmers can reach region, there may be more
    // errors in the front

    short f = region.front();
    short b = region.back();

    if(k-1 >= f) {
      // extend to front
      for(short i = f-1; i >= 0; i--)
	region.push_back(i);
    }
    if(read_length-k <= b) {
      // extend to back
      for(short i = b+1; i < read_length; i++)
	region.push_back(i);
    }
  } else
    untrusted_union(untrusted_subset, region);

  return region;
}
    
////////////////////////////////////////////////////////////
// untrusted_intersect
//
// Compute the intersection of the untrusted kmers as
// start,end return true if it's non-empty or false
// otherwise
////////////////////////////////////////////////////////////
bool Read::untrusted_intersect(vector<int> untrusted_subset, vector<short> & region) {
  int start = 0;
  int end = read_length-1;

  int u;
  for(int i = 0; i < untrusted_subset.size(); i++) {
    u = untrusted_subset[i];

    // if overlap
    if(start <= u+k-1 && u <= end) {
      // take intersection
      start = max(start, u);
      end = min(end, u+k-1);
    } else {
      // intersection is empty
      return false;   
    }
  }
    
  // intersection is non-empty
  for(short i = start; i <= end; i++)
    region.push_back(i);
  return true;
}

////////////////////////////////////////////////////////////
// untrusted_union
//
// Compute the union of the untrusted kmers, though not
////////////////////////////////////////////////////////////
void Read::untrusted_union(vector<int> untrusted_subset, vector<short> & region) {
  short u;
  set<short> region_set;
  for(int i = 0; i < untrusted_subset.size(); i++) {
    u = untrusted_subset[i];

    for(short ui = u; ui < u+k; ui++)
      region_set.insert(ui);
  }

  set<short>::iterator it;
  for(it = region_set.begin(); it != region_set.end(); it++)
    region.push_back(*it);      
}

////////////////////////////////////////////////////////////
// quality_quicksort
//
// Sort the indexes from lowest probability of an accurate
// basecall to highest
////////////////////////////////////////////////////////////
void Read::quality_quicksort(vector<short> & indexes, int left, int right) {
  int i = left, j = right;  
  short tmp;
  float pivot = prob[indexes[(left + right) / 2]];
 
  /* partition */
  while (i <= j) {
    while (prob[indexes[i]] < pivot)
      i++;    
    while (prob[indexes[j]] > pivot)      
      j--;
    if (i <= j) {
      tmp = indexes[i];
      indexes[i] = indexes[j];
      indexes[j] = tmp;
      i++;
      j--;
    }
  }

  /* recursion */
  if (left < j)
    quality_quicksort(indexes, left, j);
  if (i < right)
    quality_quicksort(indexes, i, right);  
}

////////////////////////////////////////////////////////////
// check_trust
//
// Given a corrected read and data structure holding
// trusted kmers, update the corrected_reads's vector
// of untrusted kmers and return true if it's now empty
////////////////////////////////////////////////////////////
bool Read::check_trust(corrected_read *cr, bithash *trusted) {
  // original read HAS errors
  if(cr->corrections.empty())
    return false;

  // make corrections to sequence, saving nt's to fix later
  vector<int> seqsave;
  int i;
  for(i = 0; i < cr->corrections.size(); i++) {
    seqsave.push_back(seq[cr->corrections[i].index]);
    seq[cr->corrections[i].index] = cr->corrections[i].to;
  }

  int edit = cr->corrections.back().index;
  int kmer_start = max(0, edit-k+1);
  int kmer_end = min(edit, read_length-k);

  bool non_acgt = false;
  for(i = kmer_start; i < kmer_end+k; i++)
    if(seq[i] >=4)
      non_acgt = true;

  non_acgt = true;
  if(non_acgt) {
    // easier to just check kmers one by one
    for(i = kmer_start; i <= kmer_end; i++)
      // check kmer
      cr->untrusted.set(i, !trusted->check(&seq[i]));

  } else {
    // check affected kmers
    long long unsigned kmermap;
    // check first kmer and save map value  
    cr->untrusted.set(kmer_start, !trusted->check(&seq[kmer_start], kmermap));
    for(i = kmer_start+1; i <= kmer_end; i++) {
      // check kmer using map value
      cr->untrusted.set(i, !trusted->check(kmermap, seq[i-1], seq[i+k-1]));
    }
  }

  // fix sequence
  for(i = 0; i < cr->corrections.size(); i++)
    seq[cr->corrections[i].index] = seqsave[i];

  return(cr->untrusted.none());
}
