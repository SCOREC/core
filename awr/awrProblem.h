/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWR_PROBLEM_H
#define AWR_PROBLEM_H

/* Teuchos forward declarations */
namespace Teuchos { 
class ParameterList;
}

/* apf forward declarations */
namespace apf {
class Mesh;
class Field;
class MeshEntity;
class VectorElement;
typedef VectorElement MeshElement;
template <class T> class NumberingOf;
typedef NumberingOf<int> Numbering;
typedef NumberingOf<long> GlobalNumbering;
class DynamicMatrix;
}

namespace awr {

/* save some typing */
using Teuchos::ParameterList;

/* awr forward declarations */
class LinearSystem;

/* main problem interface */
class Problem
{
  public:

    Problem(ParameterList& p, apf::Mesh* m);

    virtual ~Problem() = 0;

    void setup();

    void assemble();

    void solve();

  protected:

    /* specific parameter sublists */
    ParameterList& problemList_;
    ParameterList& qoiList_;
    ParameterList& bcList_;

    apf::Mesh* mesh_;
    int integrationOrder_;

    /* linear algebra system */
    long numGlobalEqs_;
    LinearSystem* ls_;

    /* primal and adjoint fields */
    apf::Field* primal_;
    apf::Field* adjoint_;
    int numComponents_;

    /* numbering for dof holders */
    apf::Numbering* numbering_;
    apf::GlobalNumbering* globalNumbering_;

    /* setup methods */
    virtual void validateProblemList() = 0;
    virtual void setPrimalField() = 0;
    void createAdjointField();
    void createNumbering();
    void computeNumGlobalEqs();
    void globalizeNumbering();

    /* assemble methods */
    void processBC();
    virtual void createIntegrator() = 0;
    virtual void processKe(apf::MeshElement* me,
                           int& numNodes,
                           apf::DynamicMatrix& Ke) = 0;
};

Problem* createProblem(ParameterList& p, apf::Mesh* m);

}

#endif
