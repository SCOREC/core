/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "agm.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

struct ents {
  int n[AGM_ENT_TYPES];
  int cap[AGM_ENT_TYPES];
  int* first_use[AGM_ENT_TYPES];
  int* first_bdry[AGM_ENT_TYPES];
  int* last_bdry[AGM_ENT_TYPES];
};

struct uses {
  int n[AGM_USE_TYPES];
  int cap[AGM_USE_TYPES];
  int* used[AGM_USE_TYPES];
  int* next_use_of[AGM_USE_TYPES];
  int* next_use_by[AGM_USE_TYPES];
  int* user[AGM_USE_TYPES];
};

struct bdrys {
  int n[AGM_BDRY_TYPES];
  int cap[AGM_BDRY_TYPES];
  int* next_bdry[AGM_BDRY_TYPES];
  int* first_use[AGM_BDRY_TYPES];
  int* last_use[AGM_BDRY_TYPES];
  int* bounds[AGM_BDRY_TYPES];
};

enum { MAX_SUBTYPES = 4 };

struct agm_tag {
  struct agm_tag* next;
  struct agm* m;
  int bytes;
  char* data[AGM_OBJ_TYPES][MAX_SUBTYPES];
};

struct agm_tags {
  struct agm_tag* first;
};

struct agm {
  struct ents ents;
  struct uses uses;
  struct bdrys bdrys;
  struct agm_tags tags;
};

struct agm* agm_new(void)
{
  struct agm* m = calloc(1,sizeof(*m));
  return m;
}

static void free_tag(struct agm_tag* t)
{
  int i,j;
  for (i = 0; i < AGM_OBJ_TYPES; ++i)
    for (j = 0; j < MAX_SUBTYPES; ++j)
      free(t->data[i][j]);
  free(t);
}

static void free_tags(struct agm_tags* t)
{
  struct agm_tag* old;
  while (t->first) {
    old = t->first;
    t->first = old->next;
    free_tag(old);
  }
}

static void free_ents(struct ents* e)
{
  for (enum agm_ent_type t = 0; t < AGM_ENT_TYPES; ++t) {
    free(e->first_use[t]);
    free(e->first_bdry[t]);
    free(e->last_bdry[t]);
  }
}

static void free_uses(struct uses* u)
{
  for (enum agm_use_type t = 0; t < AGM_USE_TYPES; ++t) {
    free(u->used[t]);
    free(u->next_use_of[t]);
    free(u->next_use_by[t]);
    free(u->user[t]);
  }
}

static void free_bdrys(struct bdrys* b)
{
  for (enum agm_bdry_type t = 0; t < AGM_BDRY_TYPES; ++t) {
    free(b->next_bdry[t]);
    free(b->first_use[t]);
    free(b->last_use[t]);
    free(b->bounds[t]);
  }
}

void agm_free(struct agm* m)
{
  free_ents(&m->ents);
  free_uses(&m->uses);
  free_bdrys(&m->bdrys);
  free_tags(&m->tags);
  free(m);
}

static void grow(int* n)
{
  *n = ((*n + 1) * 3) / 2;
}

static void resize(int** a, int cap)
{
  *a = realloc(*a, cap * sizeof(int));
}

static int get_cap(struct agm* m, enum agm_obj_type o, int subtype)
{
  switch (o) {
    case AGM_ENTITY:
      return m->ents.cap[subtype];
    case AGM_USE:
      return m->uses.cap[subtype];
    case AGM_BOUNDARY:
      return m->bdrys.cap[subtype];
    default:
      assert(!"bad obj type");
  }

  return 0; // Should never be reached (but stops compilers from complaining)
}

static void grow_tag(struct agm_tag* t, enum agm_obj_type o, int subtype)
{
  char** data;
  data = &t->data[o][subtype];
  if (*data)
    *data = realloc(*data, t->bytes * get_cap(t->m, o, subtype));
}

static void grow_tags(struct agm_tags* ts, enum agm_obj_type o, int subtype)
{
  struct agm_tag* t;
  for (t = ts->first; t; t = t->next)
    grow_tag(t, o, subtype);
}

void agm_reserve(struct agm* m, enum agm_ent_type t, int n)
{
  struct ents* e;
  e = &m->ents;
  assert(e->n[t] <= n);
  resize(&e->first_use[t], n);
  resize(&e->first_bdry[t], n);
  resize(&e->last_bdry[t], n);
  grow_tags(&m->tags, AGM_ENTITY, t);
  e->cap[t] = n;
}

struct agm_ent agm_add_ent(struct agm* m, enum agm_ent_type t)
{
  struct ents* e;
  struct agm_ent ent;
  e = &m->ents;
  ent.type = t;
  if (e->n[t] == e->cap[t]) {
    grow(&e->cap[t]);
    agm_reserve(m, t, e->cap[t]);
  }
  ent.id = e->n[t]++;
  e->first_use[t][ent.id] = -1;
  e->first_bdry[t][ent.id] = -1;
  e->last_bdry[t][ent.id] = -1;
  return ent;
}

static enum agm_bdry_type const ent_bdry_types[AGM_ENT_TYPES] = {
  AGM_BDRY_TYPES,
  AGM_ENDPOINTS,
  AGM_LOOP,
  AGM_SHELL
};

struct agm_bdry agm_add_bdry(struct agm* m, struct agm_ent e)
{
  struct bdrys* b;
  struct agm_bdry bdry;
  int prev;
  enum agm_bdry_type t;
  b = &m->bdrys;
  bdry.type = t = ent_bdry_types[e.type];
  if (b->n[t] == b->cap[t]) {
    grow(&b->cap[t]);
    resize(&b->next_bdry[t], b->cap[t]);
    resize(&b->first_use[t], b->cap[t]);
    resize(&b->last_use[t], b->cap[t]);
    resize(&b->bounds[t], b->cap[t]);
    grow_tags(&m->tags, AGM_BOUNDARY, t);
  }
  bdry.id = b->n[t]++;
  b->next_bdry[t][bdry.id] = -1;
  if (m->ents.first_bdry[e.type][e.id] == -1) {
    m->ents.first_bdry[e.type][e.id] = bdry.id;
    m->ents.last_bdry[e.type][e.id] = bdry.id;
  } else {
    prev = m->ents.last_bdry[e.type][e.id];
    b->next_bdry[t][prev] = bdry.id;
    m->ents.last_bdry[e.type][e.id] = bdry.id;
  }
  b->first_use[t][bdry.id] = -1;
  b->last_use[t][bdry.id] = -1;
  b->bounds[t][bdry.id] = e.id;
  return bdry;
}

static enum agm_use_type const bdry_use_types[AGM_BDRY_TYPES] = {
  AGM_VERTEX_USE,
  AGM_EDGE_USE,
  AGM_FACE_USE
};

static enum agm_ent_type const use_ent_types[AGM_USE_TYPES] = {
  AGM_VERTEX,
  AGM_EDGE,
  AGM_FACE
};

static void check_ent(struct agm* m, struct agm_ent e)
{
  assert(e.type >= 0);
  assert(e.type < AGM_ENT_TYPES);
  assert(e.id >= 0);
  assert(e.id < m->ents.n[e.type]);
}

struct agm_use agm_add_use(struct agm* m, struct agm_bdry b, struct agm_ent of)
{
  struct uses* u;
  struct agm_use use;
  int prev;
  enum agm_use_type t;
  check_ent(m, of);
  u = &m->uses;
  use.type = t = bdry_use_types[b.type];
  assert(of.type == use_ent_types[use.type]);
  if (u->n[t] == u->cap[t]) {
    grow(&u->cap[t]);
    resize(&u->used[t], u->cap[t]);
    resize(&u->next_use_of[t], u->cap[t]);
    resize(&u->next_use_by[t], u->cap[t]);
    resize(&u->user[t], u->cap[t]);
    grow_tags(&m->tags, AGM_USE, t);
  }
  use.id = u->n[t]++;
  u->used[t][use.id] = of.id;
  u->next_use_of[t][use.id] = m->ents.first_use[of.type][of.id];
  m->ents.first_use[of.type][of.id] = use.id;
  u->next_use_by[t][use.id] = -1;
  if (m->bdrys.first_use[b.type][b.id] == -1) {
    m->bdrys.first_use[b.type][b.id] = use.id;
    m->bdrys.last_use[b.type][b.id] = use.id;
  } else {
    prev = m->bdrys.last_use[b.type][b.id];
    u->next_use_by[t][prev] = use.id;
    m->bdrys.last_use[b.type][b.id] = use.id;
  }
  u->user[t][use.id] = b.id;
  return use;

}

int agm_ent_count(struct agm* m, enum agm_ent_type t)
{
  return m->ents.n[t];
}

int agm_use_count(struct agm* m, enum agm_use_type t)
{
  return m->uses.n[t];
}

int agm_bdry_count(struct agm* m, enum agm_bdry_type t)
{
  return m->bdrys.n[t];
}

static void limit_ent(struct agm* m, struct agm_ent* e)
{
  if (e->id >= m->ents.n[e->type])
    e->id = -1;
}

struct agm_ent agm_first_ent(struct agm* m, enum agm_ent_type t)
{
  struct agm_ent e;
  e.type = t;
  e.id = 0;
  limit_ent(m, &e);
  return e;
}

struct agm_ent agm_next_ent(struct agm* m, struct agm_ent e)
{
  ++(e.id);
  limit_ent(m, &e);
  return e;
}

int agm_ent_null(struct agm_ent e)
{
  return e.id == -1;
}

int agm_ent_eq(struct agm_ent a, struct agm_ent b)
{
  return a.type == b.type && a.id == b.id;
}

int agm_use_null(struct agm_use u)
{
  return u.id == -1;
}

int agm_use_eq(struct agm_use a, struct agm_use b)
{
  return a.type == b.type && a.id == b.id;
}

int agm_bdry_null(struct agm_bdry b)
{
  return b.id == -1;
}

int agm_bdry_eq(struct agm_bdry a, struct agm_bdry b)
{
  return a.type == b.type && a.id == b.id;
}

static enum agm_use_type const ent_use_types[AGM_ENT_TYPES] = {
  AGM_VERTEX_USE,
  AGM_EDGE_USE,
  AGM_FACE_USE,
  AGM_USE_TYPES
};

struct agm_use agm_first_use_of(struct agm* m, struct agm_ent e)
{
  struct agm_use u;
  u.type = ent_use_types[e.type];
  u.id = m->ents.first_use[e.type][e.id];
  return u;
}

struct agm_use agm_next_use_of(struct agm* m, struct agm_use u)
{
  u.id = m->uses.next_use_of[u.type][u.id];
  return u;
}

struct agm_use agm_first_use_by(struct agm* m, struct agm_bdry b)
{
  struct agm_use u;
  u.type = bdry_use_types[b.type];
  u.id = m->bdrys.first_use[b.type][b.id];
  return u;
}

struct agm_use agm_next_use_by(struct agm* m, struct agm_use u)
{
  u.id = m->uses.next_use_by[u.type][u.id];
  return u;
}

struct agm_ent agm_used(struct agm* m, struct agm_use u)
{
  struct agm_ent e;
  e.type = use_ent_types[u.type];
  e.id = m->uses.used[u.type][u.id];
  return e;
}

static enum agm_bdry_type const use_bdry_types[AGM_USE_TYPES] = {
  AGM_ENDPOINTS,
  AGM_LOOP,
  AGM_SHELL
};

struct agm_bdry agm_user(struct agm* m, struct agm_use u)
{
  struct agm_bdry b;
  b.type = use_bdry_types[u.type];
  b.id = m->uses.user[u.type][u.id];
  return b;
}

static enum agm_ent_type const bdry_ent_types[AGM_BDRY_TYPES] = {
  AGM_EDGE,
  AGM_FACE,
  AGM_REGION
};

struct agm_ent agm_bounds(struct agm* m, struct agm_bdry b)
{
  struct agm_ent e;
  e.type = bdry_ent_types[b.type];
  e.id = m->bdrys.bounds[b.type][b.id];
  return e;
}

struct agm_bdry agm_first_bdry_of(struct agm* m, struct agm_ent e)
{
  struct agm_bdry b;
  b.type = ent_bdry_types[e.type];
  b.id = m->ents.first_bdry[e.type][e.id];
  return b;
}

struct agm_bdry agm_next_bdry_of(struct agm* m, struct agm_bdry b)
{
  b.id = m->bdrys.next_bdry[b.type][b.id];
  return b;
}

int agm_use_count_of(struct agm* m, struct agm_ent e)
{
  struct agm_use u;
  int i = 0;
  for (u = agm_first_use_of(m, e); !agm_use_null(u);
       u = agm_next_use_of(m, u))
    ++i;
  return i;
}

int agm_use_count_by(struct agm* m, struct agm_bdry b)
{
  struct agm_use u;
  int i = 0;
  for (u = agm_first_use_by(m, b); !agm_use_null(u);
       u = agm_next_use_by(m, u))
    ++i;
  return i;
}

int agm_bdry_count_of(struct agm* m, struct agm_ent e)
{
  struct agm_bdry b;
  int i = 0;
  for (b = agm_first_bdry_of(m, e); !agm_bdry_null(b);
       b = agm_next_bdry_of(m, b))
    ++i;
  return i;
}

int agm_down_count(struct agm* m, struct agm_ent e)
{
  struct agm_bdry b;
  int i = 0;
  for (b = agm_first_bdry_of(m, e); !agm_bdry_null(b);
       b = agm_next_bdry_of(m, b))
    i += agm_use_count_by(m, b);
  return i;
}

struct agm_tag* agm_new_tag(struct agm* m, int bytes)
{
  struct agm_tag* t;
  t = calloc(1, sizeof(*t));
  t->next = m->tags.first;
  m->tags.first = t;
  t->bytes = bytes;
  t->m = m;
  return t;
}

void* agm_tag_at(struct agm_tag* t, enum agm_obj_type o,
    int subtype, int index)
{
  char** data;
  data = &t->data[o][subtype];
  if (!(*data))
    *data = malloc(t->bytes * get_cap(t->m, o, subtype));
  assert(*data);
  return *data + t->bytes * index;
}

/* yes, these tables are the identity map for now,
   but we will cause trouble later if we just
   cast back and forth everywhere */

enum agm_ent_type agm_type_from_dim(int dim)
{
  static enum agm_ent_type const tab[4] = {
    AGM_VERTEX,
    AGM_EDGE,
    AGM_FACE,
    AGM_REGION
  };
  return tab[dim];
}

int agm_dim_from_type(enum agm_ent_type t)
{
  static int const tab[AGM_ENT_TYPES] = {0,1,2,3};
  return tab[t];
}

struct agm_use agm_find_use_by_bdry(struct agm* m, struct agm_ent of,
    struct agm_bdry by)
{
  struct agm_use u;
  for (u = agm_first_use_of(m, of);
       !agm_use_null(u);
       u = agm_next_use_of(m, u))
    if (agm_bdry_eq(agm_user(m, u), by))
      return u;
  u.id = -1;
  return u;
}

struct agm_use agm_find_use_by_ent(struct agm* m, struct agm_ent of,
    struct agm_ent by)
{
  struct agm_use u;
  for (u = agm_first_use_of(m, of);
       !agm_use_null(u);
       u = agm_next_use_of(m, u))
    if (agm_ent_eq(agm_bounds(m, agm_user(m, u)), by))
      return u;
  u.id = -1;
  return u;
}

static int find_path(struct agm* m, struct agm_ent from, struct agm_ent to,
    struct agm_use* path, int depth)
{
  if (!depth)
    return agm_ent_eq(from, to);
  for (*path = agm_first_use_of(m, from);
       !agm_use_null(*path);
       *path = agm_next_use_of(m, *path))
    if (find_path(m, agm_bounds(m, agm_user(m, *path)), to,
                  path + 1, depth - 1))
      return 1;
  return 0;
}

int agm_find_path(struct agm* m, struct agm_ent from, struct agm_ent to,
    struct agm_use path[4])
{
  int depth = to.type - from.type;
  assert(depth >= 0);
  if (!find_path(m, from, to, path, depth))
    return -1;
  return depth;
}
