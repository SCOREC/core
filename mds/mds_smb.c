/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "mds_apf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <PCU.h>
#include <pcu_io.h>
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/

enum { SMB_VERSION = 4 };

enum {
  SMB_VERT,
  SMB_EDGE,
  SMB_TRI,
  SMB_QUAD,
  SMB_HEX,
  SMB_PRIS,
  SMB_PYR,
  SMB_TET,
  SMB_TYPES
};

enum {
  SMB_INT,
  SMB_DBL
};

static int smb2mds(int smb_type)
{
  int const table[SMB_TYPES] =
  {MDS_VERTEX
  ,MDS_EDGE
  ,MDS_TRIANGLE
  ,MDS_QUADRILATERAL
  ,MDS_HEXAHEDRON
  ,MDS_WEDGE
  ,MDS_PYRAMID
  ,MDS_TETRAHEDRON};
  return table[smb_type];
}

static int mds2smb(int mds_type)
{
  int const table[MDS_TYPES] =
  {SMB_VERT
  ,SMB_EDGE
  ,SMB_TRI
  ,SMB_QUAD
  ,SMB_PRIS
  ,SMB_PYR
  ,SMB_TET
  ,SMB_HEX};
  return table[mds_type];
}

static int down_degree(int t)
{
  return mds_degree[t][mds_dim[t] - 1];
}

static void read_links(struct pcu_file* f, struct mds_links* l)
{
  unsigned i;
  PCU_READ_UNSIGNED(f, l->np);
  if (!l->np)
    return;
  l->p = malloc(l->np * sizeof(unsigned));
  pcu_read_unsigneds(f, l->p, l->np);
  l->n = malloc(l->np * sizeof(unsigned));
  l->l = malloc(l->np * sizeof(unsigned*));
  pcu_read_unsigneds(f, l->n, l->np);
  for (i = 0; i < l->np; ++i) {
    l->l[i] = malloc(l->n[i] * sizeof(unsigned));
    pcu_read_unsigneds(f, l->l[i], l->n[i]);
  }
}

static void write_links(struct pcu_file* f, struct mds_links* l)
{
  unsigned i;
  PCU_WRITE_UNSIGNED(f, l->np);
  if (!l->np)
    return;
  pcu_write_unsigneds(f, l->p, l->np);
  pcu_write_unsigneds(f, l->n, l->np);
  for (i = 0; i < l->np; ++i)
    pcu_write_unsigneds(f, l->l[i], l->n[i]);
}

static void read_header(struct pcu_file* f, unsigned* version, unsigned* dim)
{
  unsigned magic, np;
  PCU_READ_UNSIGNED(f, magic);
  PCU_READ_UNSIGNED(f, *version);
  assert(*version <= SMB_VERSION);
  PCU_READ_UNSIGNED(f, *dim);
  PCU_READ_UNSIGNED(f, np);
  if (*version >= 1)
    assert(np == (unsigned)PCU_Comm_Peers());
}

static void write_header(struct pcu_file* f, unsigned dim)
{
  unsigned magic = 0;
  unsigned version = SMB_VERSION;
  unsigned np;
  PCU_WRITE_UNSIGNED(f, magic);
  PCU_WRITE_UNSIGNED(f, version);
  PCU_WRITE_UNSIGNED(f, dim);
  np = (unsigned)PCU_Comm_Peers();
  PCU_WRITE_UNSIGNED(f, np);
}

static void make_verts(struct mds_apf* m)
{
  mds_id i;
  for (i = 0; i < m->mds.cap[MDS_VERTEX]; ++i)
    mds_create_entity(&m->mds, MDS_VERTEX, NULL);
}

static void read_conn(struct pcu_file* f, struct mds_apf* m)
{
  unsigned* conn;
  struct mds_set down;
  int const* dt;
  mds_id cap;
  size_t size;
  int type_mds;
  int i;
  mds_id j;
  int k;
  for (i = 1; i < SMB_TYPES; ++i) {
    type_mds = smb2mds(i);
    down.n = down_degree(type_mds);
    cap = m->mds.cap[type_mds];
    dt = mds_types[type_mds][mds_dim[type_mds] - 1];
    size = down.n * cap;
    conn = malloc(size * sizeof(*conn));
    pcu_read_unsigneds(f, conn, size);
    for (j = 0; j < cap; ++j) {
      for (k = 0; k < down.n; ++k)
        down.e[k] = mds_identify(dt[k], conn[j * down.n + k]);
      mds_create_entity(&m->mds, type_mds, down.e);
    }
    free(conn);
    assert(m->mds.n[type_mds] == m->mds.cap[type_mds]);
  }
}

static void write_conn(struct pcu_file* f, struct mds_apf* m)
{
  unsigned* conn;
  struct mds_set down;
  mds_id end;
  size_t size;
  int type_mds;
  int i;
  mds_id j;
  int k;
  for (i = 1; i < SMB_TYPES; ++i) {
    type_mds = smb2mds(i);
    down.n = down_degree(type_mds);
    end = m->mds.end[type_mds];
    size = down.n * end;
    conn = malloc(size * sizeof(*conn));
    for (j = 0; j < end; ++j) {
      mds_get_adjacent(&m->mds, mds_identify(type_mds, j),
          mds_dim[type_mds] - 1, &down);
      for (k = 0; k < down.n; ++k)
        conn[j * down.n + k] = mds_index(down.e[k]);
    }
    pcu_write_unsigneds(f, conn, size);
    free(conn);
  }
}

static void read_remotes(struct pcu_file* f, struct mds_apf* m)
{
  struct mds_links ln = MDS_LINKS_INIT;
  read_links(f, &ln);
  mds_set_type_links(&m->remotes, &m->mds, MDS_VERTEX, &ln);
  mds_free_links(&ln);
}

static void write_remotes(struct pcu_file* f, struct mds_apf* m)
{
  struct mds_links ln = MDS_LINKS_INIT;
  mds_get_type_links(&m->remotes, &m->mds, MDS_VERTEX, &ln);
  write_links(f, &ln);
  mds_free_links(&ln);
}

static void read_class(struct pcu_file* f, struct mds_apf* m)
{
  mds_id cap;
  size_t size;
  int type_mds;
  unsigned* class;
  int i,j;
  for (i = 0; i < SMB_TYPES; ++i) {
    type_mds = smb2mds(i);
    cap = m->mds.cap[type_mds];
    size = 2 * cap;
    class = malloc(size * sizeof(*class));
    pcu_read_unsigneds(f, class, size);
    for (j = 0; j < cap; ++j) {
      m->model[type_mds][j] =
        mds_find_model(m, class[2 * j + 1], class[2 * j]);
      assert(m->model[type_mds][j]);
    }
    free(class);
  }
}

static void write_class(struct pcu_file* f, struct mds_apf* m)
{
  mds_id end;
  size_t size;
  int type_mds;
  unsigned* class;
  struct gmi_ent* model;
  int i,j;
  for (i = 0; i < SMB_TYPES; ++i) {
    type_mds = smb2mds(i);
    end = m->mds.end[type_mds];
    size = 2 * end;
    class = malloc(size * sizeof(*class));
    for (j = 0; j < end; ++j) {
      model = m->model[type_mds][j];
      class[2 * j + 1] = mds_model_dim(m, model);
      class[2 * j] = mds_model_id(m, model);
    }
    pcu_write_unsigneds(f, class, size);
    free(class);
  }
}

static struct mds_tag* read_tag_header(struct pcu_file* f, struct mds_apf* m)
{
  unsigned type, count;
  char* name;
  struct mds_tag* t;
  int type_apf[2];
  size_t bytes[2];
  type_apf[SMB_INT] = mds_apf_int;
  type_apf[SMB_DBL] = mds_apf_double;
  bytes[SMB_INT] = sizeof(int);
  bytes[SMB_DBL] = sizeof(double);
  PCU_READ_UNSIGNED(f, type);
  PCU_READ_UNSIGNED(f, count);
  pcu_read_string(f, &name);
  t = mds_create_tag(&m->tags, name,
      count * bytes[type], type_apf[type]);
  free(name);
  return t;
}

static void write_tag_header(struct pcu_file* f, struct mds_tag* t)
{
  unsigned type, count;
  int type_smb[2];
  size_t bytes[2];
  type_smb[mds_apf_int] = SMB_INT;
  type_smb[mds_apf_double] = SMB_DBL;
  bytes[mds_apf_int] = sizeof(int);
  bytes[mds_apf_double] = sizeof(double);
  type = type_smb[t->user_type];
  count = t->bytes / bytes[t->user_type];
  PCU_WRITE_UNSIGNED(f, type);
  PCU_WRITE_UNSIGNED(f, count);
  pcu_write_string(f, t->name);
}

static void read_int_tag(struct pcu_file* f, struct mds_apf* m,
                         struct mds_tag* tag, unsigned count, int t)
{
  unsigned* ids;
  unsigned* tmp;
  int size;
  unsigned i;
  int j;
  mds_id e;
  int* p;
  unsigned* q;
  ids = malloc(count * sizeof(*ids));
  size = tag->bytes / sizeof(int);
  tmp = malloc(size * count * sizeof(*tmp));
  pcu_read_unsigneds(f, ids, count);
  pcu_read_unsigneds(f, tmp, size * count);
  for (i = 0; i < count; ++i) {
    e = mds_identify(t, ids[i]);
    mds_give_tag(tag, &m->mds, e);
    p = mds_get_tag(tag, e);
    q = tmp + i * size;
    for (j = 0; j < size; ++j)
      p[j] = q[j];
  }
  free(tmp);
  free(ids);
}

static mds_id count_tagged(struct mds_apf* m, struct mds_tag* tag, int t)
{
  mds_id i;
  mds_id count = 0;
  for (i = 0; i < m->mds.end[t]; ++i)
    if (mds_has_tag(tag, mds_identify(t, i)))
      ++count;
  return count;
}

static void write_int_tag(struct pcu_file* f, struct mds_apf* m,
                          struct mds_tag* tag, unsigned count, int t)
{
  unsigned* ids;
  unsigned* tmp;
  int size;
  unsigned i;
  int j;
  unsigned k;
  mds_id e;
  int* p;
  unsigned* q;
  ids = malloc(count * sizeof(*ids));
  size = tag->bytes / sizeof(int);
  tmp = malloc(size * count * sizeof(*tmp));
  k = 0;
  for (i = 0; i < (unsigned)(m->mds.end[t]); ++i) {
    e = mds_identify(t, i);
    if (!mds_has_tag(tag, e))
      continue;
    ids[k] = i;
    p = mds_get_tag(tag, e);
    q = tmp + k * size;
    for (j = 0; j < size; ++j)
      q[j] = p[j];
    ++k;
  }
  assert(k == count);
  pcu_write_unsigneds(f, ids, count);
  pcu_write_unsigneds(f, tmp, size * count);
  free(tmp);
  free(ids);
}

static void read_dbl_tag(struct pcu_file* f, struct mds_apf* m,
    struct mds_tag* tag, unsigned count, int t)
{
  unsigned* ids;
  double* tmp;
  int size;
  unsigned i;
  mds_id e;
  double* p;
  double* q;
  ids = malloc(count * sizeof(*ids));
  size = tag->bytes / sizeof(double);
  tmp = malloc(size * count * sizeof(*tmp));
  pcu_read_unsigneds(f, ids, count);
  pcu_read_doubles(f, tmp, size * count);
  for (i = 0; i < count; ++i) {
    e = mds_identify(t, ids[i]);
    mds_give_tag(tag, &m->mds, e);
    p = mds_get_tag(tag, e);
    q = tmp + i * size;
    memcpy(p, q, size * sizeof(double));
  }
  free(tmp);
  free(ids);
}

static void write_dbl_tag(struct pcu_file* f, struct mds_apf* m,
                          struct mds_tag* tag, unsigned count, int t)
{
  unsigned* ids;
  double* tmp;
  int size;
  unsigned i;
  int j;
  unsigned k;
  mds_id e;
  double* p;
  double* q;
  ids = malloc(count * sizeof(*ids));
  size = tag->bytes / sizeof(double);
  tmp = malloc(size * count * sizeof(*tmp));
  k = 0;
  for (i = 0; i < (unsigned)(m->mds.end[t]); ++i) {
    e = mds_identify(t, i);
    if (!mds_has_tag(tag, e))
      continue;
    ids[k] = i;
    p = mds_get_tag(tag, e);
    q = tmp + k * size;
    for (j = 0; j < size; ++j)
      q[j] = p[j];
    ++k;
  }
  assert(k == count);
  pcu_write_unsigneds(f, ids, count);
  pcu_write_doubles(f, tmp, size * count);
  free(tmp);
  free(ids);
}

static void read_tags(struct pcu_file* f, struct mds_apf* m)
{
  unsigned n;
  unsigned* sizes;
  struct mds_tag** tags;
  unsigned i,j;
  int type_mds;
  PCU_READ_UNSIGNED(f,n);
  tags = malloc(n * sizeof(*tags));
  sizes = malloc(n * sizeof(*sizes));
  for (i = 0; i < n; ++i)
    tags[i] = read_tag_header(f, m);
  for (i = 0; i < SMB_TYPES; ++i) {
    pcu_read_unsigneds(f, sizes, n);
    type_mds = smb2mds(i);
    for (j = 0; j < n; ++j) {
      if (tags[j]->user_type == mds_apf_int)
        read_int_tag(f, m, tags[j], sizes[j], type_mds);
      else
        read_dbl_tag(f, m, tags[j], sizes[j], type_mds);
    }
  }
  free(tags);
  free(sizes);
}

static void write_tags(struct pcu_file* f, struct mds_apf* m)
{
  unsigned n;
  unsigned* sizes;
  struct mds_tag* t;
  int i,j;
  int type_mds;
  n = 0;
  for (t = m->tags.first; t; t = t->next)
    if (t->user_type != mds_apf_long)
      ++n;
  PCU_WRITE_UNSIGNED(f,n);
  sizes = malloc(n * sizeof(*sizes));
  for (t = m->tags.first; t; t = t->next)
    if (t->user_type != mds_apf_long)
      write_tag_header(f, t);
  for (i = 0; i < SMB_TYPES; ++i) {
    type_mds = smb2mds(i);
    j = 0;
    for (t = m->tags.first; t; t = t->next) {
      if (t->user_type == mds_apf_long)
        continue;
      sizes[j++] = count_tagged(m, t, type_mds);
    }
    pcu_write_unsigneds(f, sizes, n);
    j = 0;
    for (t = m->tags.first; t; t = t->next) {
      if (t->user_type == mds_apf_int)
        write_int_tag(f, m, t, sizes[j++], type_mds);
      else if (t->user_type == mds_apf_double)
        write_dbl_tag(f, m, t, sizes[j++], type_mds);
    }
  }
  free(sizes);
}

static void read_type_matches(struct pcu_file* f, struct mds_apf* m, int t)
{
  struct mds_links ln = MDS_LINKS_INIT;
  read_links(f, &ln);
  mds_set_local_matches(&m->matches, &m->mds, t, &ln);
  mds_free_local_links(&ln);
  mds_set_type_links(&m->matches, &m->mds, t, &ln);
  mds_free_links(&ln);
}

static void write_type_matches(struct pcu_file* f, struct mds_apf* m, int t)
{
  struct mds_links ln = MDS_LINKS_INIT;
  mds_get_type_links(&m->matches, &m->mds, t, &ln);
  mds_get_local_matches(&m->matches, &m->mds, t, &ln);
  write_links(f, &ln);
  mds_free_links(&ln);
}

static void read_matches_old(struct pcu_file* f, struct mds_apf* m)
{
  int t;
  for (t = 0; t < MDS_HEXAHEDRON; ++t)
    read_type_matches(f, m, t);
}

static void read_matches_new(struct pcu_file* f, struct mds_apf* m)
{
  int t;
  for (t = 0; t < SMB_TYPES; ++t)
    read_type_matches(f, m, smb2mds(t));
}

static void write_matches(struct pcu_file* f, struct mds_apf* m)
{
  int t;
  for (t = 0; t < SMB_TYPES; ++t)
    write_type_matches(f, m, smb2mds(t));
}

static struct mds_apf* read_smb(struct gmi_model* model, const char* filename, int zip)
{
  struct mds_apf* m;
  struct pcu_file* f;
  unsigned version;
  unsigned dim;
  unsigned n[SMB_TYPES];
  mds_id cap[MDS_TYPES];
  int i;
  f = pcu_fopen(filename, 0, zip);
  assert(f);
  read_header(f, &version, &dim);
  pcu_read_unsigneds(f, n, SMB_TYPES);
  for (i = 0; i < MDS_TYPES; ++i)
    cap[i] = n[mds2smb(i)];
  m = mds_apf_create(model, dim, cap);
  make_verts(m);
  read_conn(f, m);
  pcu_read_doubles(f, &m->point[0][0], 3 * n[SMB_VERT]);
  if (version >= 2)
    pcu_read_doubles(f, &m->param[0][0], 2 * n[SMB_VERT]);
  read_remotes(f, m);
  read_class(f, m);
  read_tags(f, m);
  if (version >= 4)
    read_matches_new(f, m);
  else if (version >= 3)
    read_matches_old(f, m);
  pcu_fclose(f);
  return m;
}

static void write_coords(struct pcu_file* f, struct mds_apf* m)
{
  size_t count;
  count = m->mds.end[MDS_VERTEX] * 3;
  pcu_write_doubles(f, &m->point[0][0], count);
  count = m->mds.end[MDS_VERTEX] * 2;
  pcu_write_doubles(f, &m->param[0][0], count);
}

static void write_smb(struct mds_apf* m, const char* filename,
                      int zip)
{
  struct pcu_file* f;
  unsigned n[SMB_TYPES] = {0};
  int i;
  f = pcu_fopen(filename, 1, zip);
  assert(f);
  write_header(f, m->mds.d);
  for (i = 0; i < MDS_TYPES; ++i)
    n[mds2smb(i)] = m->mds.end[i];
  pcu_write_unsigneds(f, n, SMB_TYPES);
  write_conn(f, m);
  write_coords(f, m);
  write_remotes(f, m);
  write_class(f, m);
  write_tags(f, m);
  write_matches(f, m);
  pcu_fclose(f);
}

static int ends_with(const char* s, const char* w)
{
  int ls = strlen(s);
  int lw = strlen(w);
  if (ls < lw)
    return 0;
  return !strcmp(s + ls - lw, w);
}

static int starts_with(const char* s, const char* w)
{
  int ls = strlen(s);
  int lw = strlen(w);
  if (ls < lw)
    return 0;
  return !strncmp(s, w, lw);
}

#define SMB_FANOUT 2048

static char* handle_path(const char* in, int is_write, int* zip)
{
  size_t n;
  char* tmp;
  char* out;
  int li = strlen(in);
  mode_t const dir_perm = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int subdir;
  int self = PCU_Comm_Self();
  n = li + 256;
  tmp = malloc(n);
  out = malloc(n);
  strcpy(tmp,in);
  if (starts_with(tmp, "bz2:")) {
    *zip = 1;
    memmove(tmp, tmp + 4, li - 3);
  } else {
    *zip = 0;
  }
  if (ends_with(tmp, "/")) {
    if (is_write && (!self))
      mkdir(tmp, dir_perm);
    PCU_Barrier();
    if (PCU_Comm_Peers() > SMB_FANOUT) {
      subdir = self / SMB_FANOUT;
      snprintf(out, n, "%s%d/", tmp, subdir);
      if (is_write && (self % SMB_FANOUT == 0))
        mkdir(out, dir_perm);
      PCU_Barrier();
      strcpy(tmp, out);
    }
  } else if (ends_with(tmp, ".smb")) {
    tmp[li - 4] = '\0';
  } else {
    fprintf(stderr,"invalid smb path %s\n",tmp);
    abort();
  }
  snprintf(out, n, "%s%d.smb", tmp, self);
  free(tmp);
  return out;
}

struct mds_apf* mds_read_smb(struct gmi_model* model, const char* pathname)
{
  char* filename;
  int zip;
  struct mds_apf* m;
  filename = handle_path(pathname, 0, &zip);
  m = read_smb(model, filename, zip);
  free(filename);
  return m;
}

static int is_compact(struct mds_apf* m)
{
  int t;
  for (t = 0; t < MDS_TYPES; ++t)
    if (m->mds.n[t] != m->mds.end[t])
      return 0;
  return 1;
}

struct mds_apf* mds_write_smb(struct mds_apf* m, const char* pathname)
{
  char* filename;
  int zip;
  if (PCU_Or(!is_compact(m)))
    m = mds_reorder(m);
  filename = handle_path(pathname, 1, &zip);
  write_smb(m, filename, zip);
  free(filename);
  return m;
}

