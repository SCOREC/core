#include "gmi_lookup.h"
#include "agm.h"
#include <stdlib.h>
#include <assert.h>

struct entry {
  int index;
  int tag;
};

struct gmi_lookup {
  struct entry* sorted[AGM_ENT_TYPES];
  struct agm_tag* tag;
  struct agm* topo;
};

static int comp_entries(const void* va, const void* vb)
{
  const struct entry* a = va;
  const struct entry* b = vb;
  return a->tag - b->tag;
}

static int* get_tag(struct gmi_lookup* l, struct agm_ent e)
{
  return agm_tag_at(l->tag, AGM_ENTITY, e.type, e.id);
}

void gmi_freeze_lookup(struct gmi_lookup* l, enum agm_ent_type t)
{
  struct entry* es;
  int n;
  struct agm_ent e;
  struct entry* y;
  assert(!(l->sorted[t]));
  n = agm_ent_count(l->topo, t);
  es = malloc(n * sizeof(*es));
  y = es;
  for (e = agm_first_ent(l->topo, t);
       !agm_ent_null(e);
       e = agm_next_ent(l->topo, e)) {
    y->index = e.id;
    y->tag = *(get_tag(l, e));
    ++y;
  }
  assert(y - es == n);
  qsort(es, n, sizeof(*es), comp_entries);
  l->sorted[t] = es;
}

struct gmi_lookup* gmi_new_lookup(struct agm* topo)
{
  struct gmi_lookup* l;
  l = calloc(1, sizeof(*l));
  l->tag = agm_new_tag(topo, sizeof(int));
  l->topo = topo;
  return l;
}

void gmi_set_lookup(struct gmi_lookup* l, struct agm_ent e, int tag)
{
  *(get_tag(l, e)) = tag;
}

int gmi_get_lookup(struct gmi_lookup* l, struct agm_ent e)
{
  return *(get_tag(l, e));
}

struct agm_ent gmi_look_up(struct gmi_lookup* l, enum agm_ent_type t, int tag)
{
  struct entry key;
  struct entry* found;
  struct agm_ent e;
  if (l->sorted[t]) {
    e.type = t;
    e.id = -1;
    key.tag = tag;
    found = bsearch(&key, l->sorted[t], agm_ent_count(l->topo, t), sizeof(key),
                    comp_entries);
    if (found)
      e.id = found->index;
  } else {
    for (e = agm_first_ent(l->topo, t);
         !agm_ent_null(e);
         e = agm_next_ent(l->topo, e))
      if (*(get_tag(l, e)) == tag)
        return e;
  }
  return e;
}

void gmi_unfreeze_lookups(struct gmi_lookup* l)
{
  enum agm_ent_type t;
  for (t = 0; t < AGM_ENT_TYPES; ++t) {
    free(l->sorted[t]);
    l->sorted[t] = 0;
  }
}

void gmi_free_lookup(struct gmi_lookup* l)
{
  gmi_unfreeze_lookups(l);
  free(l);
}
