/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "mds_net.h"
#include <PCU_C.h>
#include <string.h>
#include <stdlib.h>
#include <pcu_util.h>

void mds_create_net(struct mds_net* net)
{
  memset(net, 0, sizeof(*net));
}

void mds_destroy_net(struct mds_net* net, struct mds* m)
{
  int t;
  mds_id i;
  for (t = 0; t < MDS_TYPES; ++t) {
    if (net->data[t])
      for (i = 0; i < m->cap[t]; ++i)
        free(net->data[t][i]);
    free(net->data[t]);
  }
}

struct mds_copies* mds_make_copies(int n)
{
  struct mds_copies* c;
  c = malloc(sizeof(struct mds_copies) + (n - 1) * sizeof(struct mds_copy));
  c->n = n;
  return c;
}

void mds_set_copies(struct mds_net* net, struct mds* m, mds_id e,
    struct mds_copies* c)
{
  struct mds_copies** p;
  int t;
  mds_id i;
  t = mds_type(e);
  i = mds_index(e);
  if (!net->data[t]) {
    if (c)
      net->data[t] = calloc(m->cap[t], sizeof(*(net->data[t])));
    else
      return;
  }
  p = &net->data[t][i];
  PCU_ALWAYS_ASSERT(p);
  if (!*p && c)
    ++net->n[t];
  else if (*p && !c)
    --net->n[t];
  free(*p);
  *p = c;
  if (!net->n[t]) {
    free(net->data[t]);
    net->data[t] = NULL;
  }
}

struct mds_copies* mds_get_copies(struct mds_net* net, mds_id e)
{
  int t = mds_type(e);
  if (!net->data[t])
    return NULL;
  return net->data[t][mds_index(e)];
}

void mds_grow_net(
    struct mds_net* net,
    struct mds* m,
    mds_id old_cap[MDS_TYPES])
{
  int t;
  mds_id i;
  for (t = 0; t < MDS_TYPES; ++t)
    if (net->data[t]) {
      net->data[t] = realloc(net->data[t],
          m->cap[t] * sizeof(struct mds_copies*));
      for (i = old_cap[t]; i < m->cap[t]; ++i)
        net->data[t][i] = NULL;
    }
}

static int find_place(struct mds_copies* cs, int p)
{
  int i;
  for (i = 0; i < cs->n; ++i)
    if (cs->c[i].p > p)
      return i;
  return cs->n;
}

void mds_add_copy(struct mds_net* net, struct mds* m, mds_id e,
    struct mds_copy c)
{
  struct mds_copies* cs;
  int t;
  int p;
  mds_id i;
  t = mds_type(e);
  i = mds_index(e);
  cs = mds_get_copies(net, e);
  if (cs) {
    p = find_place(cs, c.p);
    cs = realloc(cs, sizeof(struct mds_copies) +
        (cs->n) * sizeof(struct mds_copy));
/* insert sorted by moving greater items up by one */
    memmove(&cs->c[p + 1], &cs->c[p], (cs->n - p) * sizeof(struct mds_copy));
    cs->c[p] = c;
    ++cs->n;
    net->data[t][i] = cs;
  } else {
    cs = mds_make_copies(1);
    cs->c[0] = c;
    mds_set_copies(net, m, e, cs);
  }
}

static int find_peer(struct mds_links* ln, unsigned p)
{
  unsigned i;
  for (i = 0; i < ln->np; ++i)
    if (ln->p[i] == p)
      return i;
  return -1;
}

static int insert_peer(struct mds_links* ln, int p)
{
  int i;
  i = find_peer(ln, p);
  if (i != -1)
    return i;
  i = ln->np;
  ++(ln->np);
  ln->p = realloc(ln->p, ln->np * sizeof(unsigned));
  ln->n = realloc(ln->n, ln->np * sizeof(unsigned));
  ln->l = realloc(ln->l, ln->np * sizeof(unsigned*));
  ln->p[i] = p;
  ln->n[i] = 0;
  ln->l[i] = NULL;
  return i;
}

static void note_peer(struct mds_links* ln, int p)
{
  int i;
  i = insert_peer(ln, p);
  ++(ln->n[i]);
}

static void for_type_net(PCU_t h, struct mds_net* net, struct mds* m,
    int t, void (*f)(PCU_t h, mds_id i, struct mds_copy c, void* u), void* u)
{
  mds_id i;
  int j;
  struct mds_copies* cs;
  for (i = 0; i < m->end[t]; ++i) {
    cs = mds_get_copies(net, mds_identify(t, i));
    if (!cs)
      continue;
    for (j = 0; j < cs->n; ++j)
      f(h, i, cs->c[j], u);
  }
}

static void note_remote_link(PCU_t h, mds_id i, struct mds_copy c, void* u)
{
  (void)i;
  if (c.p != PCU_Comm_Self(h))
    note_peer(u, c.p);
}

/* allocate the link arrays based on np and n */
static void alloc_links(struct mds_links* ln)
{
  unsigned i;
  for (i = 0; i < ln->np; ++i)
    ln->l[i] = malloc(ln->n[i] * sizeof(unsigned));
}

static void take_remote_link(PCU_t h, mds_id i, struct mds_copy c, void* u)
{
  struct mds_links* ln = u;
  int pi;
  if (PCU_Comm_Self(h) < c.p) {
    pi = find_peer(ln, c.p);
    ln->l[pi][ln->n[pi]] = i;
    /* use n as position keeper */
    ++ln->n[pi];
  }
}

static void send_remote_link(PCU_t h, mds_id i, struct mds_copy c, void* u)
{
  (void)i;
  (void)u;
  if (PCU_Comm_Self(h) < c.p)
    PCU_COMM_PACK(h, c.p, c.e);
}

static void recv_links(PCU_t h, struct mds_links* ln)
{
  int from;
  mds_id* tmp;
  int pi;
  unsigned i;
  from = PCU_Comm_Sender(h);
  pi = find_peer(ln, from);
  tmp = PCU_Comm_Extract(h, ln->n[pi] * sizeof(mds_id));
  for (i = 0; i < ln->n[pi]; ++i)
    ln->l[pi][i] = mds_index(tmp[i]);
}

void mds_get_type_links(PCU_t h, struct mds_net* net, struct mds* m,
    int t, struct mds_links* ln)
{
  unsigned i;
  /* count remote links in np, p and n */
  for_type_net(h, net, m, t, note_remote_link, ln);
  alloc_links(ln);
  for (i = 0; i < ln->np; ++i)
    if (((unsigned)PCU_Comm_Self(h)) < ln->p[i])
      /* zero n's for behavior in take_copy */
      ln->n[i] = 0;
  /* record indices in local order for owned boundaries */
  for_type_net(h, net, m, t, take_remote_link, ln);
  PCU_Comm_Begin(h);
  /* send indices in local order for owned boundaries */
  for_type_net(h, net, m, t, send_remote_link, NULL);
  PCU_Comm_Send(h);
  while (PCU_Comm_Listen(h))
  /* recv indices in remote order for non-owned boundaries */
    recv_links(h, ln);
}

void mds_set_type_links(PCU_t h, struct mds_net* net, struct mds* m,
    int t, struct mds_links* ln)
{
  unsigned i;
  unsigned j;
  unsigned* in;
  struct mds_copy c;
  PCU_Comm_Begin(h);
  for (i = 0; i < ln->np; ++i) {
    PCU_ALWAYS_ASSERT(ln->l);
    for (j = 0; j < ln->n[i]; ++j)
      PCU_COMM_PACK(h, ln->p[i], ln->l[i][j]);
  }
  PCU_Comm_Send(h);
  while (PCU_Comm_Listen(h)) {
    c.p = PCU_Comm_Sender(h);
    PCU_ALWAYS_ASSERT(c.p != PCU_Comm_Self(h));
    i = find_peer(ln, c.p);
    in = PCU_Comm_Extract(h, ln->n[i] * sizeof(unsigned));
    for (j = 0; j < ln->n[i]; ++j) {
      c.e = mds_identify(t, in[j]);
      mds_add_copy(net, m, mds_identify(t, ln->l[i][j]), c);
    }
  }
}

void mds_free_links(struct mds_links* ln)
{
  unsigned i;
  free(ln->n);
  free(ln->p);
  for (i = 0; i < ln->np; ++i)
    free(ln->l[i]);
  free(ln->l);
}

int mds_net_empty(struct mds_net* net)
{
  int t;
  for (t = 0; t < MDS_TYPES; ++t)
    if (net->data[t])
      return 0;
  return 1;
}

static void note_local_link(PCU_t h, mds_id i, struct mds_copy c, void* u)
{
  if (c.p == PCU_Comm_Self(h)) {
    if (i < mds_index(c.e))
      note_peer(u, PCU_Comm_Self(h));
    else /* hack id to store self-receivers */
      note_peer(u, PCU_Comm_Peers(h));
  }
}

static void take_local_link(PCU_t h, mds_id i, struct mds_copy c, void* u)
{
  struct mds_links* ln = u;
  int self = find_peer(ln, PCU_Comm_Self(h));
  int other = find_peer(ln, PCU_Comm_Peers(h));
  mds_id j = mds_index(c.e);
  if ((PCU_Comm_Self(h) == c.p) && (i < j)) {
    ln->l[self][ln->n[self]] = i;
    ln->l[other][ln->n[other]] = j;
    /* use ns as (redundant) position keepers */
    ++ln->n[self];
    ++ln->n[other];
  }
}

void mds_get_local_matches(PCU_t h, struct mds_net* net, struct mds* m,
                         int t, struct mds_links* ln)
{
  int self, other;
  for_type_net(h, net, m, t, note_local_link, ln);
  self = find_peer(ln, PCU_Comm_Self(h));
  if (self == -1)
    return;
  other = find_peer(ln, PCU_Comm_Peers(h));
  PCU_ALWAYS_ASSERT(ln->n[self] == ln->n[other]);
  ln->l[self] = malloc(ln->n[self] * sizeof(unsigned));
  ln->l[other] = malloc(ln->n[other] * sizeof(unsigned));
  ln->n[self] = 0;
  ln->n[other] = 0;
  for_type_net(h, net, m, t, take_local_link, ln);
}

void mds_set_local_matches(PCU_t h, struct mds_net* net, struct mds* m,
                         int t, struct mds_links* ln)
{
  int self, other;
  unsigned i;
  mds_id a, b;
  struct mds_copy c;
  c.p = PCU_Comm_Self(h);
  self = find_peer(ln, PCU_Comm_Self(h));
  if (self == -1)
    return;
  other = find_peer(ln, PCU_Comm_Peers(h));
  PCU_ALWAYS_ASSERT(ln->n != 0);
  PCU_ALWAYS_ASSERT(ln->n[self] == ln->n[other]);
  for (i = 0; i < ln->n[self]; ++i) {
    PCU_ALWAYS_ASSERT(ln->l != 0);
    a = mds_identify(t, ln->l[self][i]);
    b = mds_identify(t, ln->l[other][i]);
    c.e = b;
    mds_add_copy(net, m, a, c);
    c.e = a;
    mds_add_copy(net, m, b, c);
  }
}

void mds_free_local_links(PCU_t h, struct mds_links* ln)
{
  int self, other;
  self = find_peer(ln, PCU_Comm_Self(h));
  if (self == -1)
    return;
  other = find_peer(ln, PCU_Comm_Peers(h));
  PCU_ALWAYS_ASSERT(ln->n != 0);
  ln->n[self] = ln->n[other] = 0;
  PCU_ALWAYS_ASSERT(ln->l != 0);
  free(ln->l[self]);
  free(ln->l[other]);
  ln->l[self] = ln->l[other] = NULL;
}
