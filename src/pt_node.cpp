#include "pt_node.h"
#include <cstring>

pt_node::pt_node() {
  for(int i = 0; i < 4; i++)
    child[i] = NULL;
}

pt_node::~pt_node() {
  for(int i = 0; i < 4; i++) {
    if(child[i] != NULL)
      delete child[i];
  }
}
