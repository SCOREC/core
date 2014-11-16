/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWR_QOI_H
#define AWR_QOI_H

namespace Teuchos { 
class ParameterList;
}

namespace apf {
class Mesh;
class Field;
class MeshEntity;
class DynamicVector;
}

namespace awr {

using Teuchos::ParameterList;

class QoI
{
  public:
    QoI(ParameterList& p, apf::Mesh* m, apf::Field* f);
    virtual ~QoI() = 0;
    void setup();
  protected:
    ParameterList& qoiList_;
    apf::Mesh* mesh_;
    int integrationOrder_;
    apf::Field* primal_;
    virtual void validateQoIList() = 0;
    virtual void createIntegrator() = 0;
    virtual void processFe(apf::MeshEntity* e, apf::DynamicVector& Fe) = 0;
};

QoI* createQoI(ParameterList& p, apf::Mesh* m, apf::Field* f);

}

#endif
