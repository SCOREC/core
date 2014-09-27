/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "mds.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static void* mds_realloc(void* p, size_t n)
{
  if ((!p)&&(!n))
    return NULL;
  if (n)
    p = realloc(p,n);
  else {
    free(p);
    p = NULL;
  }
  assert((p) || (!n));
  return p;
}

#define REALLOC(p,n) ((p)=mds_realloc(p,(n)*sizeof(*(p))))
#define ZERO(o) memset(&(o),0,sizeof(o))

#define MDS_LIVE -2

int const mds_dim[MDS_TYPES] =
{0 /* MDS_VERTEX */
,1 /* MDS_EDGE */
,2 /* MDS_TRIANGLE */
,2 /* MDS_QUADRILATERAL */
,3 /* MDS_WEDGE */
,3 /* MDS_PYRAMID */
,3 /* MDS_TETRAHEDRON */
,3 /* MDS_HEXAHEDRON */
};

int const mds_degree[MDS_TYPES][4] =
{{1, 0,0,0} /* MDS_VERTEX */
,{2, 1,0,0} /* MDS_EDGE */
,{3, 3,1,0} /* MDS_TRIANGLE */
,{4, 4,1,0} /* MDS_QUADRILATERAL */
,{6, 9,5,1} /* MDS_WEDGE */
,{5, 8,5,1} /* MDS_PYRAMID */
,{4, 6,4,1} /* MDS_TETRAHEDRON */
,{8,12,6,1} /* MDS_HEXAHEDRON */
};

static int const e0[] = {MDS_VERTEX,MDS_VERTEX};

static int const t0[] = {MDS_VERTEX,MDS_VERTEX,MDS_VERTEX};
static int const t1[] = {MDS_EDGE,MDS_EDGE,MDS_EDGE};

static int const t01[] =     {0,1,1,2,2,0};
static int const t10[] = {2,0,0,1,1,2};

static int const q0[] = {MDS_VERTEX,MDS_VERTEX,MDS_VERTEX,MDS_VERTEX};
static int const q1[] = {MDS_EDGE,MDS_EDGE,MDS_EDGE,MDS_EDGE};

static int const q01[] =     {0,1,1,2,2,3,3,0};
static int const q10[] = {3,0,0,1,1,2,2,3};

static int const T0[] = {MDS_VERTEX,MDS_VERTEX,MDS_VERTEX,MDS_VERTEX};
static int const T1[] = {MDS_EDGE,MDS_EDGE,MDS_EDGE
                        ,MDS_EDGE,MDS_EDGE,MDS_EDGE};
static int const T2[] = {MDS_TRIANGLE,MDS_TRIANGLE,MDS_TRIANGLE,MDS_TRIANGLE};

static int const T01[] = {0,1,1,2,2,0
                         ,0,3,1,3,2,3};
static int const T10[] = {2,0,0,1,1,2,3,4};
static int const T12[] = {0,1,2
                         ,0,4,3
                         ,1,5,4
                         ,2,3,5};
static int const T21[] = {0,1
                         ,0,2
                         ,0,3
                         ,1,3
                         ,1,2
                         ,2,3};

static int const W0[] = {MDS_VERTEX,MDS_VERTEX,MDS_VERTEX
                        ,MDS_VERTEX,MDS_VERTEX,MDS_VERTEX};
static int const W1[] = {MDS_EDGE,MDS_EDGE,MDS_EDGE
                        ,MDS_EDGE,MDS_EDGE,MDS_EDGE
                        ,MDS_EDGE,MDS_EDGE,MDS_EDGE};
static int const W2[] = {MDS_TRIANGLE
                        ,MDS_QUADRILATERAL,MDS_QUADRILATERAL,MDS_QUADRILATERAL
                        ,MDS_TRIANGLE};

static int const W01[] = {0,1,1,2,2,0
                         ,0,3,1,4,2,5
                         ,3,4,4,5,5,3};
static int const W10[] = {0,2,0,1,1,2
                         ,6,8,6,7,7,8};
static int const W12[] = {0,1,2
                         ,0,4,6,3
                         ,1,4,7,5
                         ,2,3,8,5
                         ,6,7,8};
static int const W21[] = {0,1,0,2,0,3
                         ,1,3,1,2,2,3
                         ,1,4,2,4,3,4};

static int const P0[] = {MDS_VERTEX,MDS_VERTEX,MDS_VERTEX,MDS_VERTEX
                        ,MDS_VERTEX};
static int const P1[] = {MDS_EDGE,MDS_EDGE,MDS_EDGE,MDS_EDGE
                        ,MDS_EDGE,MDS_EDGE,MDS_EDGE,MDS_EDGE};
static int const P2[] = {MDS_QUADRILATERAL
                        ,MDS_TRIANGLE,MDS_TRIANGLE,MDS_TRIANGLE,MDS_TRIANGLE};

static int const P01[] = {0,1,1,2
                         ,2,3,3,0
                         ,0,4,1,4
                         ,2,4,3,4};
static int const P10[] = {0,3,0,1
                         ,1,2,2,3
                         ,4,5};
static int const P12[] = {0,1,2,3
                         ,0,5,4
                         ,1,6,5
                         ,2,7,6
                         ,3,4,7};
static int const P21[] = {0,1,0,2,0,3
                         ,0,4,1,4,1,2
                         ,2,3,3,4};

static int const H0[] = {MDS_VERTEX,MDS_VERTEX,MDS_VERTEX,MDS_VERTEX
                        ,MDS_VERTEX,MDS_VERTEX,MDS_VERTEX,MDS_VERTEX};
static int const H1[] = {MDS_EDGE,MDS_EDGE,MDS_EDGE,MDS_EDGE
                        ,MDS_EDGE,MDS_EDGE,MDS_EDGE,MDS_EDGE
                        ,MDS_EDGE,MDS_EDGE,MDS_EDGE,MDS_EDGE};
static int const H2[] = {MDS_QUADRILATERAL,MDS_QUADRILATERAL,MDS_QUADRILATERAL
                        ,MDS_QUADRILATERAL,MDS_QUADRILATERAL,MDS_QUADRILATERAL};

static int const H01[] = {0,1,1,2,2,3,3,0
                         ,0,4,1,5,2,6,3,7
                         ,4,5,5,6,6,7,7,4};
static int const H10[] = {0,3,0,1,1, 2,2, 3
                         ,4,8,5,9,6,10,7,11};
static int const H12[] = {3,2, 1, 0
                         ,0,5, 8, 4
                         ,1,6, 9, 5
                         ,2,7,10, 6
                         ,3,4,11, 7
                         ,8,9,10,11};
static int const H21[] = {0,1,0,2,0,3,0,4
                         ,1,4,1,2,2,3,3,4
                         ,1,5,2,5,3,5,4,5};

int const* mds_types[MDS_TYPES][4] =
{{0 ,0 ,0 ,0}
,{e0,0 ,0 ,0}
,{t0,t1,0 ,0}
,{q0,q1,0 ,0}
,{W0,W1,W2,0}
,{P0,P1,P2,0}
,{T0,T1,T2,0}
,{H0,H1,H2,0}
};
static int const* convs[MDS_TYPES][4][4] =
{{{0,0  ,0,0},{0  ,0,0  ,0},{0,0  ,0,0},{0,0,0,0}}
,{{0,0  ,0,0},{0  ,0,0  ,0},{0,0  ,0,0},{0,0,0,0}}
,{{0,t01,0,0},{t10,0,0  ,0},{0,0  ,0,0},{0,0,0,0}}
,{{0,q01,0,0},{q10,0,0  ,0},{0,0  ,0,0},{0,0,0,0}}
,{{0,W01,0,0},{W10,0,W12,0},{0,W21,0,0},{0,0,0,0}}
,{{0,P01,0,0},{P10,0,P12,0},{0,P21,0,0},{0,0,0,0}}
,{{0,T01,0,0},{T10,0,T12,0},{0,T21,0,0},{0,0,0,0}}
,{{0,H01,0,0},{H10,0,H12,0},{0,H21,0,0},{0,0,0,0}}
};

static void resize_down(struct mds* m, int from, int to,
    mds_id new_cap[MDS_TYPES])
{
  int t;
  int deg;
  for (t = 0; t < MDS_TYPES; ++t)
    if (mds_dim[t] == from) {
      deg = mds_degree[t][to];
      REALLOC(m->down[to][t],new_cap[t] * deg);
    }
}

static void resize_up(struct mds* m, int from, int to,
    mds_id old_cap[MDS_TYPES],
    mds_id new_cap[MDS_TYPES])
{
  int t;
  int deg;
  mds_id i;
  for (t = 0; t < MDS_TYPES; ++t) {
    if (mds_dim[t] == to) {
      deg = mds_degree[t][from];
      REALLOC(m->up[from][t],new_cap[t] * deg);
    } else if (mds_dim[t] == from) {
      REALLOC(m->first_up[to][t],new_cap[t]);
      for (i = old_cap[t]; i < new_cap[t]; ++i) {
        assert(m->first_up[to][t]);
        m->first_up[to][t][i] = MDS_NONE;
      }
    }
  }
}

static void resize_adjacency(struct mds* m, int from, int to,
    mds_id old_cap[MDS_TYPES],
    mds_id new_cap[MDS_TYPES])
{
  if (from < to)
    resize_up(m,from,to,old_cap,new_cap);
  else
    resize_down(m,from,to,new_cap);
}

static void alloc_adjacency(struct mds* m, int from_dim, int to_dim)
{
  mds_id zero_cap[MDS_TYPES] = {0};
  resize_adjacency(m,from_dim,to_dim,zero_cap,m->cap);
}

void mds_remove_adjacency(struct mds* m, int from_dim, int to_dim)
{
  mds_id zero_cap[MDS_TYPES] = {0};
  resize_adjacency(m,from_dim,to_dim,m->cap,zero_cap);
  m->mrm[from_dim][to_dim] = 0;
}

static void resize_adjacencies(struct mds* m, mds_id old_cap[MDS_TYPES])
{
  int i,j;
  for (i = 0; i <= 3; ++i)
  for (j = 0; j <= 3; ++j)
    if (m->mrm[i][j])
      resize_adjacency(m,i,j,old_cap,m->cap);
}

static void resize_free(struct mds* m)
{
  int t;
  for (t = 0; t < MDS_TYPES; ++t)
    REALLOC(m->free[t],m->cap[t]);
}

static void resize(struct mds* m, mds_id old_cap[MDS_TYPES])
{
  resize_adjacencies(m,old_cap);
  resize_free(m);
}

void mds_create(struct mds* m, int d, mds_id cap[MDS_TYPES])
{
  int i,j;
  mds_id zero_cap[MDS_TYPES] = {0};
  ZERO(*m);
  m->d = d;
  for (i = 0; i < MDS_TYPES; ++i)
    m->cap[i] = cap[i];
  for (i = 0; i <= d; ++i)
  for (j = 0; j <= d; ++j)
    if (abs(i-j)==1)
      m->mrm[i][j] = 1;
  resize(m,zero_cap);
  for (i = 0; i < MDS_TYPES; ++i)
    m->first_free[i] = MDS_NONE;
}

void mds_destroy(struct mds* m)
{
  int i;
  mds_id old_cap[MDS_TYPES];
  for (i = 0; i < MDS_TYPES; ++i)
    old_cap[i] = m->cap[i];
  ZERO(m->cap);
  resize(m,old_cap);
}

#define ID(t,i) ((i)*MDS_TYPES + (t))
#define TYPE(id) ((id) % MDS_TYPES)
#define INDEX(id) ((id) / MDS_TYPES)

int mds_type(mds_id e)
{
  return TYPE(e);
}

mds_id mds_index(mds_id e)
{
  return INDEX(e);
}

mds_id mds_identify(int type, mds_id idx)
{
  return ID(type,idx);
}

static mds_id* at_id(mds_id* a[MDS_TYPES], mds_id x)
{
  return &(a[TYPE(x)][INDEX(x)]);
}

static void relate_down(struct mds* m, mds_id from, mds_id* to)
{
  int t;
  mds_id i;
  int j;
  int deg;
  int to_dim;
  t = TYPE(from);
  i = INDEX(from);
  to_dim = mds_dim[TYPE(to[0])];
  deg = mds_degree[t][to_dim];
  for (j = 0; j < deg; ++j)
    *at_id(m->down[to_dim],ID(t,i * deg + j)) = to[j];
}

static void relate_up(struct mds* m, mds_id e, mds_id node)
{
  int node_dim;
  int e_dim;
  mds_id* head;
  mds_id* nodep;
  e_dim = mds_dim[TYPE(e)];
  node_dim = mds_dim[TYPE(node)];
  head = at_id(m->first_up[node_dim],e);
  nodep = at_id(m->up[e_dim],node);
  *nodep = *head;
  *head = node;
}

static void unrelate_up(struct mds* m, mds_id e, mds_id node)
{
  int node_dim;
  int e_dim;
  mds_id* prev;
  e_dim = mds_dim[TYPE(e)];
  node_dim = mds_dim[TYPE(node)];
  prev = at_id(m->first_up[node_dim],e);
  while (*prev != node)
    prev = at_id(m->up[e_dim],*prev);
  *prev = *at_id(m->up[e_dim],node);
}

static void relate_back_up(struct mds* m, mds_id* from, mds_id to)
{
  int t;
  mds_id i;
  int j;
  int deg;
  int from_dim;
  mds_id x;
  t = TYPE(to);
  i = INDEX(to);
  from_dim = mds_dim[TYPE(from[0])];
  deg = mds_degree[t][from_dim];
  for (j = 0; j < deg; ++j) {
    x = ID(t, i * deg + j);
    relate_up(m,from[j],x);
  }
}

static void unrelate_back_up(struct mds* m, mds_id* from, mds_id to)
{
  int t;
  mds_id i;
  int j;
  int deg;
  int from_dim;
  mds_id x;
  t = TYPE(to);
  i = INDEX(to);
  from_dim = mds_dim[TYPE(from[0])];
  deg = mds_degree[t][from_dim];
  for (j = 0; j < deg; ++j) {
    x = ID(t, i * deg + j);
    unrelate_up(m,from[j],x);
  }
}

static void relate_both(struct mds* m, mds_id* down, mds_id up)
{
  relate_down(m,up,down);
  relate_back_up(m,down,up);
}

static void look_up(struct mds* m, mds_id const e, int d, struct mds_set* s)
{
  mds_id* n;
  mds_id nv;
  int t;
  mds_id i;
  int deg;
  mds_id* es = s->e;
  mds_id** p;
  t = TYPE(e);
  i = INDEX(e);
  p = m->first_up[d];
  n = p[t] + i;
  nv = *n;
  d = mds_dim[t];
  p = m->up[d];
  while (nv != MDS_NONE) {
    t = TYPE(nv);
    i = INDEX(nv);
    deg = mds_degree[t][d];
    *es = ID(t, i / deg);
    ++es;
    n = p[t] + i;
    nv = *n;
  }
  s->n = es - s->e;
}

struct down {
  int n;
  mds_id* e;
};

static struct down reach_down(struct mds* m, mds_id e, int d)
{
  struct down out;
  int t;
  mds_id i;
  t = TYPE(e);
  i = INDEX(e);
  out.n = mds_degree[t][d];
  out.e = m->down[d][t] + i * out.n;
  return out;
}

static void look_down(struct mds* m, mds_id e, int d, struct mds_set* s)
{
  struct down dn;
  int j;
  dn = reach_down(m,e,d);
  s->n = dn.n;
  for (j = 0; j < dn.n; ++j)
    s->e[j] = dn.e[j];
}

static void look(struct mds* m, mds_id e, int d, struct mds_set* s)
{
  int from_dim;
  from_dim = mds_dim[TYPE(e)];
  if (d < from_dim) {
    look_down(m,e,d,s);
  } else if (from_dim == d) {
    s->n = 1; s->e[0] = e;
  } else {
    look_up(m,e,d,s);
  }
}

static void copy_set(struct mds_set* to, struct mds_set* from)
{
  int i;
  to->n = from->n;
  for (i = 0; i < from->n; ++i)
    to->e[i] = from->e[i];
}

static int contains(struct mds_set* s, mds_id id)
{
  int i;
  for (i = 0; i < s->n; ++i)
    if (s->e[i] == id)
      return 1;
  return 0;
}

static void intersect(struct mds_set* s, struct mds_set* with)
{
  int i;
  int j = 0;
  for (i = 0; i < s->n; ++i)
    if (contains(with,s->e[i]))
      s->e[j++] = s->e[i];
  s->n = j;
}

static void unite(struct mds_set* s, struct mds_set* with)
{
  int i;
  int j = s->n;
  for (i = 0; i < with->n; ++i)
    if ( ! contains(s,with->e[i])) {
      assert(j < MDS_SET_MAX);
      s->e[j++] = with->e[i];
    }
  s->n = j;
}

static mds_id common_down(struct mds* m, mds_id a, mds_id b, int d)
{
  struct down da;
  struct down db;
  int i,j;
  da = reach_down(m,a,d);
  db = reach_down(m,b,d);
  for (i = 0; i < da.n; ++i)
  for (j = 0; j < db.n; ++j)
    if (da.e[i] == db.e[j])
      return da.e[i];
  return MDS_NONE;
}

static mds_id common_up(struct mds* m, struct mds_set* s, int d)
{
  struct mds_set found;
  struct mds_set adjacent;
  int i;
  assert(0 < s->n);
  for (i = 0; i < s->n; ++i)
    if (s->e[i] == MDS_NONE)
      return MDS_NONE;
  look_up(m,s->e[0],d,&found);
  for (i = 1; i < s->n; ++i) {
    look_up(m,s->e[i],d,&adjacent);
    intersect(&found,&adjacent);
  }
  assert(found.n <= 1);
  if (found.n)
    return found.e[0];
  return MDS_NONE;
}

static void grow(struct mds* m, int t)
{
  int i;
  mds_id old_cap[MDS_TYPES];
  for (i = 0; i < MDS_TYPES; ++i)
    old_cap[i] = m->cap[i];
  m->cap[t] = ((old_cap[t] + 2) * 3) / 2;
  resize(m,old_cap);
}

static mds_id fill_hole(struct mds* m, int t)
{
  mds_id *head;
  mds_id *node;
  mds_id i;
  head = &(m->first_free[t]);
  i = *head;
  node = &(m->free[t][i]);
  *head = *node;
  *node = MDS_LIVE;
  return ID(t,i);
}

static mds_id alloc_ent(struct mds* m, int t)
{
  mds_id id;
  if (m->n[t] == m->cap[t])
    grow(m,t);
  ++(m->n[t]);
  if (m->first_free[t] == MDS_NONE)
    id = ID(t,m->n[t] - 1);
  else
    id = fill_hole(m,t);
  if (INDEX(id) == m->end[t]) {
    m->free[t][m->end[t]] = MDS_LIVE;
    ++(m->end[t]);
  }
  assert(m->end[t] >= m->n[t]);
  return id;
}

static void free_ent(struct mds* m, mds_id e)
{
  mds_id *head;
  mds_id *node;
  int t;
  mds_id i;
  t = TYPE(e);
  i = INDEX(e);
  head = &(m->first_free[t]);
  node = &(m->free[t][i]);
  *node = *head;
  *head = i;
  --(m->n[t]);
}

static mds_id add_ent(struct mds* m, int t, mds_id* from)
{
  mds_id id;
  id = alloc_ent(m,t);
  relate_both(m,from,id);
  return id;
}

static mds_id add_or_find_ent(struct mds* m, int t, struct mds_set* from,
    int can_add)
{
  mds_id id;
  id = common_up(m,from,mds_dim[t]);
  if (id != MDS_NONE)
    return id;
  if (can_add)
    return add_ent(m,t,from->e);
  return MDS_NONE;
}

static void check_ent(struct mds* m, mds_id e)
{
  int t;
  mds_id i;
  assert(e >= 0);
  t = TYPE(e);
  assert(t < MDS_TYPES);
  i = INDEX(e);
  assert(i < m->end[t]);
  assert(m->free[t][i] == MDS_LIVE);
}

static void unrelate_ent(struct mds* m, mds_id e)
{
  struct mds_set down;
  look_down(m,e,mds_dim[TYPE(e)] - 1,&down);
  unrelate_back_up(m,down.e,e);
}

void mds_destroy_entity(struct mds* m, mds_id e)
{
  check_ent(m,e);
  if (TYPE(e) != MDS_VERTEX)
    unrelate_ent(m,e);
  free_ent(m,e);
}

static void step_up(struct mds* m,
    struct mds_set* from_s, int from_dim,
    struct mds_set* to_s, int to_dim,
    int t, int make)
{
  struct mds_set s;
  int i;
  int j;
  int const* ct;
  int const* ci;
  int tt;
  ct = mds_types[t][to_dim];
  ci = convs[t][from_dim][to_dim];
  to_s->n = mds_degree[t][to_dim];
  for (i = 0; i < to_s->n; ++i) {
    tt = *ct;
    s.n = mds_degree[tt][from_dim];
    for (j = 0; j < s.n; ++j) {
      s.e[j] = from_s->e[*ci];
      ++ci;
    }
    to_s->e[i] = add_or_find_ent(m,tt,&s,make);
    ++ct;
  }
}

static void step_down(struct mds* m,
    struct mds_set* from_s, int from_dim,
    struct mds_set* to_s, int to_dim,
    int t)
{
  int i;
  int const* ci;
  mds_id a;
  mds_id b;
  ci = convs[t][from_dim][to_dim];
  to_s->n = mds_degree[t][to_dim];
  for (i = 0; i < to_s->n; ++i) {
    a = from_s->e[ci[2*i]];
    b = from_s->e[ci[2*i + 1]];
    to_s->e[i] = common_down(m,a,b,to_dim);
  }
}

static void convert_up(struct mds* m,
    struct mds_set* from_s, int from_dim,
    struct mds_set* to_s, int to_dim,
    int t, int make)
{
  struct mds_set sets[2];
  struct mds_set* s[2];
  struct mds_set* tmp;
  s[0] = sets;
  s[1] = sets + 1;
  copy_set(s[0],from_s);
  for (; from_dim != to_dim; ++from_dim) {
    step_up(m,s[0],from_dim,s[1],from_dim + 1,t,make);
    tmp = s[0];
    s[0] = s[1];
    s[1] = tmp;
  }
  copy_set(to_s,s[0]);
}

static void convert_down(struct mds* m,
    struct mds_set* from_s, int from_dim,
    struct mds_set* to_s, int to_dim,
    int t)
{
  struct mds_set sets[2];
  struct mds_set* s[2];
  struct mds_set* tmp;
  s[0] = sets;
  s[1] = sets + 1;
  copy_set(s[0],from_s);
  for (; from_dim != to_dim; --from_dim) {
    step_down(m,s[0],from_dim,s[1],from_dim - 1,t);
    tmp = s[0];
    s[0] = s[1];
    s[1] = tmp;
  }
  copy_set(to_s,s[0]);
}

static void check_set(struct mds* m, struct mds_set* s)
{
  int i;
  int dim;
  assert(0 < s->n);
  for (i = 0; i < s->n; ++i)
    check_ent(m,s->e[i]);
  dim = mds_dim[TYPE(s->e[0])];
  for (i = 1; i < s->n; ++i)
    assert(dim == mds_dim[TYPE(s->e[i])]);
}

static mds_id get_ent_far(struct mds* m, int t, struct mds_set* in, int make)
{
  int from_dim;
  int to_dim;
  struct mds_set down;
  from_dim = mds_dim[TYPE(in->e[0])];
  to_dim = mds_dim[t];
  convert_up(m,in,from_dim,&down,to_dim - 1,t,make);
  return add_or_find_ent(m,t,&down,make);
}

static mds_id get_ent(struct mds* m, int t, mds_id* from, int make)
{
  int from_dim;
  int to_dim;
  struct mds_set in;
  int i;
  check_ent(m,from[0]);
  from_dim = mds_dim[TYPE(from[0])];
  in.n = mds_degree[t][from_dim];
  for (i = 0; i < in.n; ++i)
    in.e[i] = from[i];
  check_set(m,&in);
  to_dim = mds_dim[t];
  assert(from_dim < to_dim);
  if (from_dim + 1 == to_dim)
    return add_or_find_ent(m,t,&in,make);
  return get_ent_far(m,t,&in,make);
}

mds_id mds_create_entity(struct mds* m, int t, mds_id* from)
{
  if (t == MDS_VERTEX)
    return alloc_ent(m,t);
  return get_ent(m,t,from,1);
}

mds_id mds_find_entity(struct mds* m, int t, mds_id* from)
{
  return get_ent(m,t,from,0);
}

static void expand_once(struct mds* m, struct mds_set* from, struct mds_set* to)
{
  int i;
  int dim;
  struct mds_set up;
  to->n = 0;
  assert(from->n);
  dim = mds_dim[TYPE(from->e[0])];
  for (i = 0; i < from->n; ++i) {
    look_up(m,from->e[i],dim + 1,&up);
    unite(to,&up);
  }
}

static void get_up(struct mds* m, mds_id e, int d, struct mds_set* out)
{
  struct mds_set sets[2];
  struct mds_set* s[2];
  struct mds_set* tmp;
  int dim;
  s[0] = sets;
  s[1] = sets + 1;
  s[0]->n = 1;
  s[0]->e[0] = e;
  dim = mds_dim[TYPE(e)];
  for (; dim < d; ++dim) {
    expand_once(m,s[0],s[1]);
    tmp = s[0];
    s[0] = s[1];
    s[1] = tmp;
  }
  copy_set(out,s[0]);
}

static void get_down(struct mds* m, mds_id e, int d, struct mds_set* out)
{
  struct mds_set in;
  int t;
  int from_dim;
  t = TYPE(e);
  from_dim = mds_dim[t];
  look_down(m,e,from_dim - 1,&in);
  convert_down(m,&in,from_dim - 1,out,d,t);
}

void mds_get_adjacent(struct mds* m, mds_id e, int d, struct mds_set* s)
{
  int e_dim;
  if (d > m->d) {
    s->n = 0;
    return;
  }
  check_ent(m,e);
  e_dim = mds_dim[TYPE(e)];
  if ((e_dim == d) || m->mrm[e_dim][d]) {
    look(m,e,d,s);
    return;
  } if (d < e_dim) {
    get_down(m,e,d,s);
    return;
  }
  get_up(m,e,d,s);
}

static mds_id skip(struct mds* m, mds_id e)
{
  int t;
  int d;
  mds_id i;
  t = TYPE(e);
  d = mds_dim[t];
  i = INDEX(e);
  for (; t < MDS_TYPES; ++t) {
    if (mds_dim[t] == d) {
      for (; i < m->end[t]; ++i) {
        if (m->free[t][i] == MDS_LIVE)
          return ID(t,i);
      }
      i = 0;
    }
  }
  return MDS_NONE;
}

mds_id mds_begin(struct mds* m, int d)
{
  int t;
  for (t = 0; t < MDS_TYPES; ++t)
    if (mds_dim[t] == d)
      return skip(m,ID(t,0));
  return MDS_NONE;
}

mds_id mds_next(struct mds* m, mds_id e)
{
  return skip(m,ID(TYPE(e),INDEX(e) + 1));
}

void mds_add_adjacency(struct mds* m, int from_dim, int to_dim)
{
  mds_id e;
  struct mds_set adj;
  alloc_adjacency(m,from_dim,to_dim);
  if (from_dim < to_dim)
    for (e = mds_begin(m,to_dim);
         e != MDS_NONE;
         e = mds_next(m,e)) {
      mds_get_adjacent(m,e,from_dim,&adj);
      relate_back_up(m,adj.e,e);
    }
  else
    for (e = mds_begin(m,from_dim);
         e != MDS_NONE;
         e = mds_next(m,e)) {
      mds_get_adjacent(m,e,to_dim,&adj);
      relate_down(m,e,adj.e);
    }
  m->mrm[from_dim][to_dim] = 1;
}

int mds_has_up(struct mds* m, mds_id e)
{
  int d;
  d = mds_dim[TYPE(e)];
  if (d == m->d)
    return 0;
  return *at_id(m->first_up[d + 1],e) != MDS_NONE;
}

static void increase_dimension(struct mds* m)
{
  int old_d;
  old_d = m->d;
  ++(m->d);
  assert(old_d < m->d);
  mds_add_adjacency(m, old_d, m->d);
  mds_add_adjacency(m, m->d, old_d);
}

static void decrease_dimension(struct mds* m)
{
  int old_d;
  old_d = m->d;
  --(m->d);
  mds_remove_adjacency(m, old_d, m->d);
  mds_remove_adjacency(m, m->d, old_d);
}

void mds_change_dimension(struct mds* m, int d)
{
  while (m->d < d)
    increase_dimension(m);
  while (m->d > d)
    decrease_dimension(m);
}
