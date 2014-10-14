/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRLHS_H
#define AWRLHS_H

#include "apf.h"
#include <apfDynamicMatrix.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

class LHS
{
  public:
    LHS(apf::Mesh* m, const Teuchos::ParameterList& p);
    virtual ~LHS() {};
    virtual void 
    evaluateElementLHS(apf::MeshEntity* e,
                       apf::DynamicMatrix& k) = 0;
    apf::Mesh* getMesh() { return mesh_; };
  protected:
    apf::Mesh* mesh_;
    Teuchos::ParameterList params_;
  private:
    LHS(const LHS&);
    LHS& operator=(const LHS&);
};

}

#endif
