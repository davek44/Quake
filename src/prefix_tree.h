#ifndef PREFIX_TREE_H
#define PREFIX_TREE_H

#include "pt_node.h"

const int k = 15;

class prefix_tree {
 public:
  prefix_tree();
  ~prefix_tree();
  void add(int kmer[k]);
  bool check(int kmer[k]);
  void file_load(const char* merf, const int boundary);
  int count_nodes();
 private:
  int count_nodes_recurse(pt_node* current);
  pt_node* root;
};

#endif
