/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRQOIFACTORY_H
#define AWRQOIFACTORY_H

#include "awrQOI.h"
#include <Teuchos_RCP.hpp>

namespace awr {

class QOIFactory
{
  public:
    QOIFactory(apf::Mesh* m, const Teuchos::ParameterList& p);
    virtual ~QOIFactory() {};
    virtual Teuchos::RCP<QOI> create();
  protected:
    apf::Mesh* mesh_;
    Teuchos::ParameterList params_;
  private:
    QOIFactory(const QOIFactory&);
    QOIFactory& operator=(const QOIFactory&);
};

}

#endif
