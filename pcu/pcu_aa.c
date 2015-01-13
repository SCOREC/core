/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
/*
   These functions and data structures are from:
    A. Andersson, "Balanced Search Trees Made Simple," pp. 60-71, 1993.
   The "aa" prefix is for Arne Andersson.
   He deserves full credit for their functionality and simplicity.
   I have modified them such that nodes are intrusive rather than
   containing a key which is the reason for the "less" function pointers
   and the more complicated remove function.
*/

#include <stdlib.h>
#include <assert.h>
#include "pcu_aa.h"

static pcu_aa_node pcu_aa_bottom = 
{
  .left = &pcu_aa_bottom,
  .right = &pcu_aa_bottom,
  .level = 0
};

void pcu_make_aa(pcu_aa_tree* t)
{
  *t = &pcu_aa_bottom;
}

bool pcu_aa_empty(pcu_aa_tree t)
{
  return t == &pcu_aa_bottom;
}

static void skew(pcu_aa_tree* t)
{
  pcu_aa_tree temp;
  if ((*t)->left->level == (*t)->level)
  { /* rotate right */
    temp = *t;
    *t = (*t)->left;
    temp->left = (*t)->right;
    (*t)->right = temp;
  }
}

static void split(pcu_aa_tree* t)
{
  pcu_aa_tree temp;
  if ((*t)->right->right->level == (*t)->level)
  { /* rotate left */
    temp = *t;
    *t = (*t)->right;
    temp->right = (*t)->left;
    (*t)->left = temp;
    ++((*t)->level);
  }
}

pcu_aa_node* pcu_aa_insert(pcu_aa_node* x, pcu_aa_tree* t, pcu_aa_less* less)
{
  pcu_aa_node* result;
  if (*t == &pcu_aa_bottom)
  {
    result = x;
    *t = x;
    (*t)->left = &pcu_aa_bottom;
    (*t)->right = &pcu_aa_bottom;
    (*t)->level = 1;
  }
  else
  {
    if (less(x,*t))
      result = pcu_aa_insert(x,&((*t)->left),less);
    else if (less(*t,x))
      result = pcu_aa_insert(x,&((*t)->right),less);
    else
      result = *t;
    skew(t);
    split(t);
  }
  return result;
}

struct remove_vars
{
  pcu_aa_node* x;
  pcu_aa_less* less;
  pcu_aa_node* deleted;
  pcu_aa_node* last;
  pcu_aa_node* successor;
};

static void remove_successor(pcu_aa_tree* t, struct remove_vars* v)
{
  assert(*t);
  if (*t != &pcu_aa_bottom)
  { /* search down the tree and set pointers last and deleted */
    v->last = *t;
    if (v->less(v->x,*t))
      remove_successor(&((*t)->left),v);
    else
    {
      v->deleted = *t;
      remove_successor(&((*t)->right),v);
    }
  }
  /* at the bottom of the tree we remove the element (if it is present) */
  if ((*t == v->last)&&(v->deleted != &pcu_aa_bottom)
    &&(!v->less(v->x,v->deleted))&&(!v->less(v->deleted,v->x)))
  {
    v->successor = *t;
    *t = (*t)->right;
  }
  /* on the way back, we rebalance */
  else if (((*t)->left->level < (*t)->level-1)
       ||((*t)->right->level < (*t)->level-1))
  {
    --((*t)->level);
    if ((*t)->right->level > (*t)->level)
      (*t)->right->level = (*t)->level;
    skew(t);
    skew(&((*t)->right));
    skew(&((*t)->right->right));
    split(t);
    split(&((*t)->right));
  }
}

static void swap_successor(pcu_aa_tree* t, struct remove_vars* v)
{
  if (v->less(v->deleted,*t))
    swap_successor(&((*t)->left),v);
  else if (v->less(*t,v->deleted))
    swap_successor(&((*t)->right),v);
  else
  {
    if (v->deleted != *t)
      abort();
    *t = v->successor;
    (*t)->left = v->deleted->left;
    (*t)->right = v->deleted->right;
    (*t)->level = v->deleted->level;
  }
}

pcu_aa_node* pcu_aa_remove(pcu_aa_node* x, pcu_aa_tree* t, pcu_aa_less* less)
{
  struct remove_vars v;
  v.x = x;
  v.less = less;
  v.deleted = &pcu_aa_bottom;
  v.last = 0;
  v.successor = 0;
  remove_successor(t,&v);
  if (!v.successor)
    return 0;
  if (v.successor != v.deleted)
    swap_successor(t,&v);
  return v.deleted;
}

pcu_aa_node* pcu_aa_find(pcu_aa_node* x, pcu_aa_tree t, pcu_aa_less* less)
{
  if (t == &pcu_aa_bottom)
    return 0;
  if (less(x,t))
    return pcu_aa_find(x,t->left,less);
  else if (less(t,x))
    return pcu_aa_find(x,t->right,less);
  else
    return t;
}

int pcu_aa_count(pcu_aa_tree t)
{
  if (pcu_aa_empty(t))
    return 0;
  return pcu_aa_count(t->left) + pcu_aa_count(t->right) + 1;
}

