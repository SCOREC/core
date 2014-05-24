/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_H
#define GMI_H

struct gmi_ent;

struct gmi_iter;

struct gmi_model;

struct gmi_model_ops {
  struct gmi_iter* (*begin)(struct gmi_model* m, int dim);
  struct gmi_ent* (*next)(struct gmi_model* m, struct gmi_iter* i);
  void (*end)(struct gmi_model* m, struct gmi_iter* i);
  int (*dim)(struct gmi_model* m, struct gmi_ent* e);
  int (*tag)(struct gmi_model* m, struct gmi_ent* e);
  struct gmi_ent* (*find)(struct gmi_model* m, int dim, int tag);
  void (*eval)(struct gmi_model* m, struct gmi_ent* e,
      double const p[2], double x[3]);
  void (*eval_du)(struct gmi_model* m, struct gmi_ent* e,
      double const p[2], double t[3]);
  void (*eval_dv)(struct gmi_model* m, struct gmi_ent* e,
      double const p[2], double t[3]);
  void (*reparam)(struct gmi_model* m, struct gmi_ent* from,
      double const from_p[2], struct gmi_ent* to, double to_p[2]);
  int (*periodic)(struct gmi_model* m, struct gmi_ent* e, int dim);
  void (*range)(struct gmi_model* m, struct gmi_ent* e, int dim,
      double r[2]);
  void (*destroy)(struct gmi_model* m);
};

struct gmi_model {
  struct gmi_model_ops const* ops;
  int n[4];
};

typedef struct gmi_model* (*gmi_creator)(const char* filename);

#ifdef __cplusplus
extern "C" {
#endif

struct gmi_iter* gmi_begin(struct gmi_model* m, int dim);
struct gmi_ent* gmi_next(struct gmi_model* m, struct gmi_iter* i);
void gmi_end(struct gmi_model* m, struct gmi_iter* i);
int gmi_dim(struct gmi_model* m, struct gmi_ent* e);
int gmi_tag(struct gmi_model* m, struct gmi_ent* e);
struct gmi_ent* gmi_find(struct gmi_model* m, int dim, int tag);
int gmi_can_eval(struct gmi_model* m);
void gmi_eval(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double x[3]);
void gmi_eval_du(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double t[3]);
void gmi_eval_dv(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double t[3]);
void gmi_reparam(struct gmi_model* m, struct gmi_ent* from,
    double const from_p[2], struct gmi_ent* to, double to_p[2]);
int gmi_periodic(struct gmi_model* m, struct gmi_ent* e, int dim);
void gmi_range(struct gmi_model* m, struct gmi_ent* e, int dim,
    double r[2]);
void gmi_destroy(struct gmi_model* m);

void gmi_register(gmi_creator f, const char* ext);
struct gmi_model* gmi_load(const char* filename);

void gmi_fail(const char* why);

#ifdef __cplusplus
}
#endif

#endif
