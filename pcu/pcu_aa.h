/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_AA_H
#define PCU_AA_H

#include <stdbool.h>

typedef struct pcu_aa_node_struct* pcu_aa_tree;

struct pcu_aa_node_struct
{
  pcu_aa_tree left;
  pcu_aa_tree right;
  int level;
};
typedef struct pcu_aa_node_struct pcu_aa_node;

typedef bool pcu_aa_less(pcu_aa_node* a, pcu_aa_node* b);

void pcu_make_aa(pcu_aa_tree* t);
bool pcu_aa_empty(pcu_aa_tree t);
pcu_aa_node* pcu_aa_insert(pcu_aa_node* x, pcu_aa_tree* t, pcu_aa_less* less);
pcu_aa_node* pcu_aa_remove(pcu_aa_node* x, pcu_aa_tree* t, pcu_aa_less* less);
pcu_aa_node* pcu_aa_find(pcu_aa_node* x, pcu_aa_tree t, pcu_aa_less* less);
int pcu_aa_count(pcu_aa_tree t);

#endif
