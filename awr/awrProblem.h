/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWR_PROBLEM_H
#define AWR_PROBLEM_H

namespace Teuchos { 
class ParameterList;
}

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

using Teuchos::ParameterList;

class LinearSystem;
class QoI;

void print(const char* format, ...) __attribute__((format(printf,1,2)));

class Problem
{
  public:
    Problem(ParameterList& p, apf::Mesh* m);
    virtual ~Problem() = 0;
    void setup();
    void assemble();
    void solve();
  protected:
    ParameterList& problemList_;
    ParameterList& qoiList_;
    ParameterList& bcList_;
    apf::Mesh* mesh_;
    int integrationOrder_;
    long numGlobalEqs_;
    LinearSystem* ls_;
    apf::Field* primal_;
    apf::Field* adjoint_;
    int numComponents_;
    apf::Numbering* numbering_;
    apf::GlobalNumbering* globalNumbering_;
    QoI* qoi_;
    virtual void validateProblemList() = 0;
    virtual void setPrimalField() = 0;
    virtual void createIntegrator() = 0;
    virtual void processKe(apf::MeshEntity* e, apf::DynamicMatrix& Ke) = 0;
};

Problem* createProblem(ParameterList& p, apf::Mesh* m);

}

#endif
