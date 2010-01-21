#ifndef PREFIX_TREE_H
#define PREFIX_TREE_H

#include "node.h"

const int k = 15;

class prefix_tree {
 public:
  prefix_tree();
  ~prefix_tree();
  void add(unsigned int kmer[k]);
  bool check(unsigned int kmer[k]);
  void file_load(const char* merf, const int boundary);
  int count_nodes();
 private:
  int count_nodes_recurse(node* current);
  node* root;
};

#endif
