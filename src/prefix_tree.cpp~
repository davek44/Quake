#include "prefix_tree.h"
#include <cstring>
#include <iostream>
#include <fstream>

using namespace::std;

prefix_tree::prefix_tree() {
  root = new pt_node();
}

prefix_tree::~prefix_tree() {
  delete root;
}

////////////////////////////////////////////////////////////
// add
//
// Add a single sequence to the tree
////////////////////////////////////////////////////////////
void prefix_tree::add(int kmer[k]) {
  pt_node* current = root;
  for(int i = 0; i < k; i++) {
    if(current->child[kmer[i]] == NULL) {
      current->child[kmer[i]] = new pt_node();
    }
    current = current->child[kmer[i]];
  }
}

////////////////////////////////////////////////////////////
// check
//
// Check for the presence of a sequence in the tree
////////////////////////////////////////////////////////////
bool prefix_tree::check(int kmer[k]) {
  pt_node* current = root;
  //cout << "check: ";
  for(int i = 0; i < k; i++) {
    //cout << kmer[i] << " ";
    if(kmer[i] >= 4 || current->child[kmer[i]] == NULL) {
      //cout << endl;
      return false;
    }
    current = current->child[kmer[i]];
  }
  //cout << endl;
  return true;
} 

////////////////////////////////////////////////////////////
// file_load
//
// Make a prefix_tree from kmers in the file given that
// occur >= "boundary" times
////////////////////////////////////////////////////////////
void prefix_tree::file_load(const char* merf, const int boundary) {
  const char* nts = "ACGT";
  ifstream mer_in(merf);
  string line;
  int count;
  int seq[k];
  int tmpseq[k];
  bool add_kmer = false;

  while(getline(mer_in, line)) {
    if(line[0] == '>') {
      // get count
      count = atoi(line.substr(1).c_str());
      //cout << count << endl;
      
      // compare to boundary
      if(count >= boundary) {\
	add_kmer = true;
      } else {
	add_kmer = false;
      }

    } else if(add_kmer) {
      // convert sequence
      for(int i = 0; i < k; i++) {
	seq[i] = strchr(nts, line[i]) - nts;
	//cout << seq[i] << " ";
      }
      //cout << endl;

      // add to tree
      add(seq);

      // reverse
      for(int i = 0; i < k; i++)
	tmpseq[i] = seq[k-1-i];
      // complement
      for(int i = 0; i < k; i++)
	seq[i] = 3 - tmpseq[i];

      // add to tree
      add(seq);
    }
  }
}

int prefix_tree::count_nodes() {
  return count_nodes_recurse(root);
}

int prefix_tree::count_nodes_recurse(pt_node* current) {
  int count = 1;
  for(int nt = 0; nt < 4; nt++) {
    if(current->child[nt] != NULL) {
      count += count_nodes_recurse(current->child[nt]);
    }
  }
  return count;
}
