/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "mds_tag.h"
#include <stdlib.h>
#include <string.h>

void mds_create_tags(struct mds_tags* ts)
{
  ts->first = NULL;
}

void mds_destroy_tags(struct mds_tags* ts)
{
  while (ts->first)
    mds_destroy_tag(ts,ts->first);
}

static void grow_tag(
    struct mds_tag* tag,
    struct mds* m,
    mds_id old_cap[MDS_TYPES])
{
  int t;
  mds_id i;
  mds_id has[2];
  for (t = 0; t < MDS_TYPES; ++t) {
    if ( ! tag->has[t])
      continue;
    has[0] = (old_cap[t] / 8) + 1;
    has[1] = (m->cap[t] / 8) + 1;
    tag->has[t] = realloc(tag->has[t], has[1]);
    for (i = has[0]; i < has[1]; ++i)
      tag->has[t][i] = 0;
    tag->data[t] = realloc(tag->data[t],
        tag->bytes * m->cap[t]);
  }
}

void mds_grow_tags(
    struct mds_tags* ts,
    struct mds* m,
    mds_id old_cap[MDS_TYPES])
{
  struct mds_tag* t;
  for (t = ts->first; t; t = t->next)
    grow_tag(t,m,old_cap);
}

struct mds_tag* mds_create_tag(
    struct mds_tags* ts,
    const char* name,
    int bytes,
    int user_type)
{
  int l;
  struct mds_tag* t;
  t = calloc(1,sizeof(*t));
  t->next = ts->first;
  ts->first = t;
  t->bytes = bytes;
  t->user_type = user_type;
  l = strlen(name);
  t->name = malloc(l + 1);
  strcpy(t->name,name);
  return t;
}

void mds_destroy_tag(struct mds_tags* ts, struct mds_tag* t)
{
  struct mds_tag** p;
  int i;
  for (p = &(ts->first); *p != t; p = &((*p)->next));
  *p = (*p)->next;
  for (i = 0; i < MDS_TYPES; ++i)
    free(t->data[i]);
  for (i = 0; i < MDS_TYPES; ++i)
    free(t->has[i]);
  free(t->name);
  free(t);
}

void* mds_get_tag(struct mds_tag* tag, mds_id e)
{
  return tag->data[mds_type(e)] + tag->bytes * mds_index(e);
}

struct mds_tag* mds_find_tag(struct mds_tags* ts, const char* name)
{
  struct mds_tag* p;
  for (p = ts->first; p; p = p->next)
    if ( ! strcmp(p->name,name))
      return p;
  return 0;
}

int mds_has_tag(struct mds_tag* tag, mds_id e)
{
  int t;
  mds_id i;
  mds_id c;
  int b;
  unsigned char v;
  t = mds_type(e);
  if ( ! tag->has[t])
    return 0;
  i = mds_index(e);
  c = i / 8;
  b = i % 8;
  v = tag->has[t][c] & (1 << b);
  return v != 0;
}

void mds_give_tag(struct mds_tag* tag, struct mds* m, mds_id e)
{
  int t;
  mds_id i;
  mds_id c;
  int b;
  unsigned char* has;
  t = mds_type(e);
  if ( ! tag->has[t]) {
    tag->has[t] = calloc((m->cap[t] / 8) + 1, 1);
    tag->data[t] = malloc(tag->bytes * m->cap[t]);
  }
  i = mds_index(e);
  c = i / 8;
  b = i % 8;
  has = tag->has[t] + c;
  *has |= (1<<b);
}

void mds_take_tag(struct mds_tag* tag, mds_id e)
{
  int t;
  mds_id i;
  mds_id c;
  int b;
  unsigned char* has;
  t = mds_type(e);
  i = mds_index(e);
  c = i / 8;
  b = i % 8;
  if (!tag->has[t])
    return;
  has = tag->has[t] + c;
  *has &= ~(1 << b);
}

static struct mds_tag** find_prev(struct mds_tags* ts, struct mds_tag* t)
{
  struct mds_tag* p;
  if (ts->first == t)
    return &ts->first;
  for (p = ts->first; p; p = p->next)
    if (p->next == t)
      return &p->next;
  return NULL;
}

/* users that rely on struct mds_tag* (such as apf::MeshTag*)
   not changing during the tag lifetime are
   disappointed when those pointers are invalidated by mds_reorder.
   Thus, this HACK swaps mds_tag structs so that the new mesh
   has the same struct mds_tag pointers as the old one */

void mds_swap_tag_structs(struct mds_tags* as, struct mds_tag** a,
    struct mds_tags* bs, struct mds_tag** b)
{
  struct mds_tag** pa;
  struct mds_tag** pb;
  struct mds_tag tmp;
  struct mds_tag* tmp_p;
  pa = find_prev(as, *a);
  pb = find_prev(bs, *b);
  tmp = **a;
  **a = **b;
  **b = tmp;
  *pa = *b;
  *pb = *a;
  tmp_p = *a;
  *a = *b;
  *b = tmp_p;
}
