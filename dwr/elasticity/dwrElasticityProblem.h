/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_ELASTICITYPROBLEM_H
#define DWR_ELASTICITYPROBLEM_H

/** \file dwrElasticityProblem.h
  * \brief Provides an interface to define a dual variational problem for
  * linear elasticity */

#include <apfDynamicArray.h>

class Epetra_Map;

namespace apf {
class Mesh;
class Field;
class ModelEntity;
template <class T> class NumberingOf;
typedef NumberingOf<long> GlobalNumbering;
}

namespace dwr {

class LinearSystem;

/** \brief structure for defining the boundary conditions of the dual problem.
  * \details these correspond to model entities where dirichlet boundary
  * condititions have been applied to the primal problem */
struct ElasticityDBC {
  /** \brief dimension of geometric entity */
  int dim;
  /** \brief dimension-unique tag associated with the geometric entity */
  int tag;
  /** \brief the component of the displacement vector that is constrained
    * \details 0 for x, 1 for y, 2 for z */
  int component;
};

/** \brief interface for solving elasticity duals
  * \details solves a dual variational problem for isotropic, linear
  * elasticity. currently only single material problems are supported */
class ElasticityProblem
{
  public:
    ElasticityProblem();
    ~ElasticityProblem();
    /** \brief elastic modulus */
    double E;
    /** \brief Poisson's ration */
    double nu;
    /** \brief quadrature degree used for dual solve */
    int quadratureDegree;
    /** \brief primal displacement field */
    apf::Field* primal;
    /** \brief dual displacement field
      * \details the elasticity problem will create this */
    apf::Field* dual;
    /** \brief boundary condition information */
    apf::DynamicArray<ElasticityDBC> dbc;
    /** \brief solve the dual problem */
    apf::Field* computeDual();

  private:
    apf::Mesh* mesh_;
    long nge_;
    apf::GlobalNumbering* gn_;
    Epetra_Map* owned_;
    Epetra_Map* overlap_;
    LinearSystem* ls_;
    void validate();
    void setup();
    void assemble();
    void solve();
};

}

#endif
