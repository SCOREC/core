/*****************************************************************************

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

******************************************************************************/
#include "gmi_analytic.h"
#include "gmi_null.h"
#include "gmi_mesh.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef uint8_t (*ptr_to_2_u8s)[2];
typedef double (*ptr_to_2x2_doubles)[2][2];

struct gmi_analytic
{
  struct gmi_mesh mesh;
  gmi_analytic_fun* f[4];
  ptr_to_2_u8s periodic[4];
  ptr_to_2x2_doubles ranges[4];
};

#define REALLOC(a,n) ((a)=realloc((a),(n)*sizeof(*(a))))

void gmi_add_analytic(struct gmi_model* m, int dim, int tag,
    gmi_analytic_fun f, int* periodic, double (*ranges)[2])
{
  struct gmi_analytic* m2;
  struct gmi_ent* e;
  int n;
  int i;
  int j;
  m2 = (struct gmi_analytic*)m;
  n = m->n[dim];
  e = gmi_null_find(m, dim, tag);
  if (n != m->n[dim]) {
    REALLOC(m2->f[dim], m->n[dim]);
    REALLOC(m2->periodic[dim], m->n[dim]);
    REALLOC(m2->ranges[dim], m->n[dim]);
  }
  i = gmi_index(e);
  m2->f[dim][i] = f;
  for (j = 0; j < dim; ++j) {
    m2->periodic[dim][i][j] = periodic[j];
    m2->ranges[dim][i][j][0] = ranges[j][0];
    m2->ranges[dim][i][j][1] = ranges[j][1];
  }
  for (; j < 2; ++j) {
    m2->periodic[dim][i][j] = 0;
    m2->ranges[dim][i][j][0] = 0;
    m2->ranges[dim][i][j][1] = 0;
  }
}

static struct gmi_ent* find(struct gmi_model* m, int dim, int tag)
{
  int i;
  struct gmi_mesh* m2 = (struct gmi_mesh*)m;
  for (i = 0; i < m->n[dim]; ++i)
    if (m2->tags[dim][i] == tag)
      return gmi_identify(dim, i);
  return 0;
}

static void eval(struct gmi_model* m, struct gmi_ent* e,
      double const p[2], double x[3])
{
  struct gmi_analytic* m2;
  m2 = (struct gmi_analytic*)m;
  m2->f[gmi_dim(m, e)][gmi_index(e)](p, x);
}

static void reparam(struct gmi_model* m, struct gmi_ent* from,
      double const from_p[2], struct gmi_ent* to, double to_p[2])
{
  to_p[0] = 0;
  to_p[1] = 0;
}

static int periodic(struct gmi_model* m, struct gmi_ent* e, int dim)
{
  struct gmi_analytic* m2 = (struct gmi_analytic*)m;
  return m2->periodic[gmi_dim(m, e)][gmi_index(e)][dim];
}

static void range(struct gmi_model* m, struct gmi_ent* e, int dim, double r[2])
{
  struct gmi_analytic* m2 = (struct gmi_analytic*)m;
  memcpy(r, m2->ranges[gmi_dim(m, e)][gmi_index(e)][dim], 2 * sizeof(double));
}

static void destroy(struct gmi_model* m)
{
  struct gmi_analytic* m2;
  int i;
  m2 = (struct gmi_analytic*)m;
  for (i = 0; i < 4; ++i) {
    free(m2->f[i]);
    free(m2->periodic[i]);
    free(m2->ranges[i]);
  }
  gmi_mesh_destroy(m);
}

static struct gmi_model_ops ops = {
  .begin   = gmi_mesh_begin,
  .next    = gmi_mesh_next,
  .end     = gmi_mesh_end,
  .dim     = gmi_mesh_dim,
  .tag     = gmi_mesh_tag,
  .find    = find,
  .eval    = eval,
  .reparam = reparam,
  .periodic = periodic,
  .range = range,
  .destroy = destroy
};

struct gmi_model* gmi_make_analytic(void)
{
  struct gmi_analytic* m;
  m = calloc(1, sizeof(*m));
  m->mesh.model.ops = &ops;
  return &m->mesh.model;
}
