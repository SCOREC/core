/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_H
#define GMI_H

/** \page gmi GMI
  The Geometric Model Interface (GMI) is designed to allow SCOREC codes
  to query many different kinds of geometric models from one unified
  set of APIs.

  - The main GMI interface is in gmi.h
  - The built-in meshmodel system is in gmi_mesh.h
  - The built-in analytic model is in gmi_analytic.h
  - The don't-use null model is in gmi_null.h
  */

/** \file gmi.h
  \brief abstract Geometric Model Interface */

#include <stdio.h>

struct gmi_ent;

struct gmi_iter;

struct gmi_model;

/** \brief a set of model entities
 \details users should create these only with gmi_make_set
          and always call gmi_free_set after receiving one */
struct gmi_set {
  /** \brief number of model entities */
  int n;
  /** \brief array of model entity pointers
   \details vlas at the end of structs are ok by C99, but C++ still hates them.
   trick it with a 1, like the good old days of the linux kernel */
  struct gmi_ent* e[1];
};

/** \brief model interface definition
  \details this is the C equivalent of a C++ class with
  pure virtual methods. It contains pointers to implementations
  for all the different geometric model operations supported
  by GMI. Each model points to one of these definitions
  so that GMI knows which code to use to interact with said model.
  Implementations are allowed to omit any of these pointers,
  in which case either that functionality is prohibited
  (and the program crashes upon attempted use)
  or a comment may describe the default behavior if the pointer is omitted.
  Some functionality is known to be rare, so specific APIs are in
  place to check whether a particular model supports those functions. */
struct gmi_model_ops {
  /** \brief implement gmi_begin */
  struct gmi_iter* (*begin)(struct gmi_model* m, int dim);
  /** \brief implement gmi_begin */
  struct gmi_ent* (*next)(struct gmi_model* m, struct gmi_iter* i);
  /** \brief implement gmi_end */
  void (*end)(struct gmi_model* m, struct gmi_iter* i);
  /** \brief implement gmi_dim */
  int (*dim)(struct gmi_model* m, struct gmi_ent* e);
  /** \brief implement gmi_tag */
  int (*tag)(struct gmi_model* m, struct gmi_ent* e);
  /** \brief implement gmi_find */
  struct gmi_ent* (*find)(struct gmi_model* m, int dim, int tag);
  /** \brief implement gmi_adjacent
   \details if omitted then gmi_adjacent returns gmi_make_set(0) */
  struct gmi_set* (*adjacent)(struct gmi_model* m, struct gmi_ent* e, int dim);
  /** \brief implement gmi_eval
   \details if omitted then gmi_can_eval returns false */
  void (*eval)(struct gmi_model* m, struct gmi_ent* e,
      double const p[2], double x[3]);
  /** \brief implement gmi_reparam */
  void (*reparam)(struct gmi_model* m, struct gmi_ent* from,
      double const from_p[2], struct gmi_ent* to, double to_p[2]);
  /** \brief implement gmi_periodic */
  int (*periodic)(struct gmi_model* m, struct gmi_ent* e, int dim);
  /** \brief implement gmi_range */
  void (*range)(struct gmi_model* m, struct gmi_ent* e, int dim,
      double r[2]);
  /** \brief implement gmi_destroy */
  void (*destroy)(struct gmi_model* m);
};

/** \brief the basic structure for all GMI models */
struct gmi_model {
  /** \brief pointer to interface definition */
  struct gmi_model_ops const* ops;
  /** \brief number of model entities per dimension */
  int n[4];
};

/** \brief model from file constructor, give to gmi_register */
typedef struct gmi_model* (*gmi_creator)(const char* filename);

#ifdef __cplusplus
extern "C" {
#endif

/** \brief allocate a gmi_set with (n) elements */
struct gmi_set* gmi_make_set(int n);
/** \brief free a gmi_set */
void gmi_free_set(struct gmi_set* s);

/** \brief begin an iterator over model entities of one dimension
 \details call gmi_end on this iterator afterwards */
struct gmi_iter* gmi_begin(struct gmi_model* m, int dim);
/** \brief dereference and then increment an interator
 \returns 0 if past the end, otherwise a valid entity */
struct gmi_ent* gmi_next(struct gmi_model* m, struct gmi_iter* i);
/** \brief free an iterator */
void gmi_end(struct gmi_model* m, struct gmi_iter* i);
/** \brief get the dimension of a model entity */
int gmi_dim(struct gmi_model* m, struct gmi_ent* e);
/** \brief get the tag of a model entity */
int gmi_tag(struct gmi_model* m, struct gmi_ent* e);
/** \brief lookup a model entity by dimension and tag */
struct gmi_ent* gmi_find(struct gmi_model* m, int dim, int tag);
/** \brief query model entity adjacencies
 \details currently only one-level adjacencies are supported by most
          implementations */
struct gmi_set* gmi_adjacent(struct gmi_model* m, struct gmi_ent* e, int dim);
/** \brief check whether the model implements gmi_eval */
int gmi_can_eval(struct gmi_model* m);
/** \brief evaluate the parametric definition of a model boundary entity
  \param p ignored for vertices. for edges, p[0] should be the edge
           parametric coordinate. for faces, p should contain the
           parametric u,v face coordinates
  \param x the resulting point in space */
void gmi_eval(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double x[3]);
/** \brief re-parameterize from one model entity to another
  \param from the model entity to start from
  \param from_p the parametric coordinates on entity (from),
                see gmi_eval
  \param to the model entity to reparameterize onto
  \param to_p the resulting parametric coordinates, again
              in the form described by gmi_eval */
void gmi_reparam(struct gmi_model* m, struct gmi_ent* from,
    double const from_p[2], struct gmi_ent* to, double to_p[2]);
/** \brief return true iff the model entity is periodic around this dimension */
int gmi_periodic(struct gmi_model* m, struct gmi_ent* e, int dim);
/** \brief return the range of parametric coordinates along this dimension */
void gmi_range(struct gmi_model* m, struct gmi_ent* e, int dim,
    double r[2]);
/** \brief destroy a geometric model */
void gmi_destroy(struct gmi_model* m);

/** \brief register a new geometric modeler
  \details this function registers a geometric modeler in
           GMI's global list. the modeler declares that it
           is responsible for files with extension (ext).
           Subsequent calls to gmi_load will check for extension
           (ext), and if it matches then the gmi_creator (f)
           is called with the given filename and the resulting
           model is returned.
  \param ext the model file extension, without the dot */
void gmi_register(gmi_creator f, const char* ext);
/** \brief load a geometric model file
  \details see gmi_register. the filename is checked against
           all registered extensions. The first match
           triggers a gmi_creator to load the file.
           This function will fail and abort if there are
           no matches */
struct gmi_model* gmi_load(const char* filename);

/** \brief write a dmg (meshmodel) file
  \details the .dmg format is SCOREC's meshmodel format
  which contains only model entities and topology, and
  can be loaded by the gmi_mesh.h system as its own
  structure.
  GMI can write this file format from any gmi_model
  object that implements basic iteration and dim/tag queries. */
void gmi_write_dmg(struct gmi_model* m, const char* filename);

/** \brief print the message as a gmi failure and abort
  \details this is for GMI internal use, not public users */
void gmi_fail(const char* why) __attribute__((noreturn));

/** \brief fscanf wrapper that checks return values
  \details programmers often fail to check the return
  value of fscanf, and some compiler configurations will
  complain about that. this function calls fscanf(f,format,...)
  and then requires that the return value is equal to n. */
void gmi_fscanf(FILE* f, int n, const char* format, ...);

int gmi_getline(char** line, size_t* cap, FILE* f);

#ifdef __cplusplus
}
#endif

#endif
