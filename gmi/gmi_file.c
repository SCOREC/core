#include "gmi_base.h"
#include "gmi_lookup.h"
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

void gmi_fscanf(FILE* f, int n, const char* format, ...)
{
  va_list ap;
  int r;
  va_start(ap, format);
  r = vfscanf(f, format, ap);
  va_end(ap);
  assert(r == n);
}

int gmi_getline(char** line, size_t* cap, FILE* f)
{
#ifdef __SUNPRO_C
  /* Oracle Solaris Studio doesn't come with getline, so
     we can instead use this completely non-standard-conformant
     version that works well enough for our purposes */
  int c;
  int i = 0;
  while (1) {
    c = fgetc(f);
    if (c == EOF)
      return -1;
    if (i == *cap) {
      *cap = ((*cap + 1) * 3) / 2;
      *line = realloc(*line, *cap);
    }
    (*line)[i++] = c;
    if (c == '\n')
      break;
  }
  if (i == *cap) {
    *cap = ((*cap + 1) * 3) / 2;
    *line = realloc(*line, *cap);
  }
  (*line)[i++] = '\0';
  return i;
#else
  return (int) getline(line, cap, f);
#endif
}

void gmi_base_read_dmg(struct gmi_base* m, FILE* f)
{
  int n[4];
  int i,j,k;
  int loops, shells;
  int faces, edges;
  struct agm_ent e;
  struct agm_bdry b;
  int used[2];
  struct agm_ent d;
  int tag;
  /* read entity counts */
  gmi_fscanf(f, 4, "%d %d %d %d", &n[3], &n[2], &n[1], &n[0]);
  gmi_base_init(m);
  /* bounding box */
  gmi_fscanf(f, 0, "%*f %*f %*f");
  gmi_fscanf(f, 0, "%*f %*f %*f");
  /* vertices */
  gmi_base_reserve(m, AGM_VERTEX, n[0]);
  for (i = 0; i < n[0]; ++i) {
    gmi_fscanf(f, 1, "%d %*f %*f %*f", &tag);
    e = agm_add_ent(m->topo, AGM_VERTEX);
    gmi_set_lookup(m->lookup, e, tag);
  }
  gmi_freeze_lookup(m->lookup, 0);
  /* edges */
  gmi_base_reserve(m, AGM_EDGE, n[1]);
  for (i = 0; i < n[1]; ++i) {
    gmi_fscanf(f, 3, "%d %d %d", &tag, &used[0], &used[1]);
    e = agm_add_ent(m->topo, AGM_EDGE);
    gmi_set_lookup(m->lookup, e, tag);
    b = agm_add_bdry(m->topo, e);
    for (j = 0; j < 2; ++j) {
      d = gmi_look_up(m->lookup, AGM_VERTEX, used[j]);
      if (!agm_ent_null(d))
        agm_add_use(m->topo, b, d);
    }
  }
  gmi_freeze_lookup(m->lookup, 1);
  /* faces */
  gmi_base_reserve(m, AGM_FACE, n[2]);
  for (i = 0; i < n[2]; ++i) {
    gmi_fscanf(f, 2, "%d %d", &tag, &loops);
    e = agm_add_ent(m->topo, AGM_FACE);
    gmi_set_lookup(m->lookup, e, tag);
    for (j = 0; j < loops; ++j) {
      gmi_fscanf(f, 1, "%d", &edges);
      b = agm_add_bdry(m->topo, e);
      for (k = 0; k < edges; ++k) {
        /* tag, direction */
        gmi_fscanf(f, 1, "%d %*d", &tag);
        d = gmi_look_up(m->lookup, AGM_EDGE, tag);
        agm_add_use(m->topo, b, d);
      }
    }
  }
  gmi_freeze_lookup(m->lookup, 2);
  /* regions */
  gmi_base_reserve(m, AGM_REGION, n[3]);
  for (i = 0; i < n[3]; ++i) {
    gmi_fscanf(f, 2, "%d %d", &tag, &shells);
    e = agm_add_ent(m->topo, AGM_REGION);
    gmi_set_lookup(m->lookup, e, tag);
    for (j = 0; j < shells; ++j) {
      gmi_fscanf(f, 1, "%d", &faces);
      b = agm_add_bdry(m->topo, e);
      for (k = 0; k < faces; ++k) {
        /* tag, direction */
        gmi_fscanf(f, 1, "%d %*d", &tag);
        d = gmi_look_up(m->lookup, AGM_FACE, tag);
        agm_add_use(m->topo, b, d);
      }
    }
  }
  gmi_freeze_lookup(m->lookup, 3);
}

static int starts_with(char const* with, char const* s)
{
  int lw;
  int ls;
  lw = strlen(with);
  ls = strlen(s);
  if (ls < lw)
    return 0;
  return strncmp(with, s, lw) == 0;
}

static void seek_marker(FILE* f, char const* marker)
{
  char* line = 0;
  size_t linecap = 0;
  while (-1 != gmi_getline(&line, &linecap, f))
    if (starts_with(marker, line))
      return;
}

void gmi_base_read_tess(struct gmi_base* m, FILE* f)
{
  int n;
  int tag;
  int i;
  struct agm_ent e;
  struct agm_bdry b;
  int j;
  struct agm_ent d;
  int ndown;
  gmi_base_init(m);
  seek_marker(f, " **vertex");
  gmi_fscanf(f, 1, "%d", &n);
  gmi_base_reserve(m, AGM_VERTEX, n);
  for (i = 0; i < n; ++i) {
    gmi_fscanf(f, 1, "%d %*f %*f %*f %*d", &tag);
    e = agm_add_ent(m->topo, AGM_VERTEX);
    gmi_set_lookup(m->lookup, e, tag);
  }
  gmi_freeze_lookup(m->lookup, 0);
  seek_marker(f, " **edge");
  gmi_fscanf(f, 1, "%d", &n);
  gmi_base_reserve(m, AGM_EDGE, n);
  for (i = 0; i < n; ++i) {
    gmi_fscanf(f, 1, "%d", &tag);
    e = agm_add_ent(m->topo, AGM_EDGE);
    gmi_set_lookup(m->lookup, e, tag);
    b = agm_add_bdry(m->topo, e);
    for (j = 0; j < 2; ++j) {
      gmi_fscanf(f, 1, "%d", &tag);
      d = gmi_look_up(m->lookup, AGM_VERTEX, tag);
      agm_add_use(m->topo, b, d);
    }
    gmi_fscanf(f, 0, "%*d");
  }
  gmi_freeze_lookup(m->lookup, 1);
  seek_marker(f, " **face");
  gmi_fscanf(f, 1, "%d", &n);
  gmi_base_reserve(m, AGM_FACE, n);
  for (i = 0; i < n; ++i) {
    gmi_fscanf(f, 1, "%d", &tag);
    e = agm_add_ent(m->topo, AGM_FACE);
    gmi_set_lookup(m->lookup, e, tag);
    b = agm_add_bdry(m->topo, e);
    gmi_fscanf(f, 1, "%d", &ndown);
    for (j = 0; j < ndown; ++j)
      gmi_fscanf(f, 0, "%*d"); /* skip face->vert data */
    gmi_fscanf(f, 1, "%d", &ndown);
    for (j = 0; j < ndown; ++j) {
      gmi_fscanf(f, 1, "%d", &tag);
      d = gmi_look_up(m->lookup, AGM_EDGE, abs(tag));
      agm_add_use(m->topo, b, d);
    }
    gmi_fscanf(f, 0, "%*f %*f %*f %*f"); /* face plane equation */
    gmi_fscanf(f, 0, "%*d %*d %*f %*f %*f"); /* regularization */
  }
  gmi_freeze_lookup(m->lookup, 2);
  seek_marker(f, " **polyhedron");
  gmi_fscanf(f, 1, "%d", &n);
  gmi_base_reserve(m, AGM_REGION, n);
  for (i = 0; i < n; ++i) {
    gmi_fscanf(f, 1, "%d", &tag);
    e = agm_add_ent(m->topo, AGM_REGION);
    gmi_set_lookup(m->lookup, e, tag);
    b = agm_add_bdry(m->topo, e);
    gmi_fscanf(f, 1, "%d", &ndown);
    for (j = 0; j < ndown; ++j) {
      gmi_fscanf(f, 1, "%d", &tag);
      d = gmi_look_up(m->lookup, AGM_FACE, abs(tag));
      agm_add_use(m->topo, b, d);
    }
  }
  gmi_freeze_lookup(m->lookup, 3);
}
