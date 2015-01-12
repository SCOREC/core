/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_order.h"
#include "pcu_aa.h"
#include "pcu_msg.h"
#include <assert.h>

struct message {
  pcu_aa_node node;
  int from;
  pcu_buffer buf;
};

struct pcu_order_struct {
  pcu_aa_tree tree;
  struct message** array;
  int count;
  bool ready;
  int at;
};

static void init_order(pcu_order o)
{
  pcu_make_aa(&o->tree);
  o->array = 0;
  o->count = 0;
  o->ready = false;
  o->at = 0;
}

pcu_order pcu_order_new(void)
{
  pcu_order o;
  PCU_MALLOC(o,1);
  init_order(o);
  return o;
}

static void free_message(struct message* m)
{
  pcu_free_buffer(&m->buf);
  pcu_free(m);
}

static void free_messages(pcu_aa_tree* t)
{
  if (pcu_aa_empty(*t))
    return;
  free_messages(&((*t)->left));
  free_messages(&((*t)->right));
  struct message* m;
  m = (struct message*) *t;
  free_message(m);
  pcu_make_aa(t);
}

static void dtor_order(pcu_order o)
{
  free_messages(&o->tree);
  pcu_free(o->array);
}

void pcu_order_free(pcu_order o)
{
  dtor_order(o);
  pcu_free(o);
}

static struct message* take_message(pcu_msg* t)
{
  struct message* m;
  PCU_MALLOC(m,1);
  m->from = t->received.peer;
  m->buf = t->received.buffer; /* steal the buffer */
  pcu_make_buffer(&t->received.buffer);
  return m;
}

static bool message_less(pcu_aa_node* a, pcu_aa_node* b)
{
  struct message* ma;
  struct message* mb;
  ma = (struct message*) a;
  mb = (struct message*) b;
  return ma->from < mb->from;
}

static void fill(pcu_order o, pcu_aa_tree t)
{
  if (pcu_aa_empty(t))
    return;
  fill(o, t->left);
  o->array[o->at++] = (struct message*) t;
  fill(o, t->right);
}

static void prepare(pcu_order o, pcu_msg* t)
{
  struct message* m;
  while (pcu_msg_receive(t)) {
    m = take_message(t);
    pcu_aa_insert(&m->node, &o->tree, message_less);
  }
  o->count = pcu_aa_count(o->tree);
  PCU_MALLOC(o->array, o->count);
  o->at = 0;
  fill(o, o->tree);
  o->at = -1;
  o->ready = true;
}

bool pcu_order_receive(pcu_order o, pcu_msg* m)
{
  if (!o->ready)
    prepare(o, m);
  o->at++;
  if (o->at == o->count) {
    dtor_order(o);
    init_order(o);
    return false;
  }
  pcu_begin_buffer(&o->array[o->at]->buf);
  return true;
}

void* pcu_order_unpack(pcu_order o, size_t size)
{
  return pcu_walk_buffer(&o->array[o->at]->buf, size);
}

bool pcu_order_unpacked(pcu_order o)
{
/* compatibility with pcu_msg_unpacked before pcu_msg_receive */
  if (!o->ready)
    return true;
  return pcu_buffer_walked(&o->array[o->at]->buf);
}

int pcu_order_received_from(pcu_order o)
{
  return o->array[o->at]->from;
}

size_t pcu_order_received_size(pcu_order o)
{
  return o->array[o->at]->buf.capacity;
}

