/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRQOI_H
#define AWRQOI_H

#include "apf.h"
#include <apfDynamicVector.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

class QOI
{
  public:
    QOI(apf::Mesh* m, const Teuchos::ParameterList& p);
    virtual ~QOI() {};
    virtual void 
    evaluateElementQOI(apf::MeshEntity* e,
                       apf::DynamicVector& f) = 0;
  protected:
    apf::Mesh* mesh_;
    Teuchos::ParameterList params_;
  private:
    QOI(const QOI&);
    QOI& operator=(const QOI&);
};

}

#endif
