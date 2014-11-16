/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */


#include "awrQoI.h"
#include "awrDomainIntegral.h"
#include <Teuchos_ParameterList.hpp>

namespace awr {

QoI::QoI(ParameterList& p, apf::Mesh* m, apf::Field* f) :
  qoiList_(p),
  mesh_(m),
  primal_(f)
{
  integrationOrder_ = qoiList_.get("Integration Order",2);
}

QoI::~QoI()
{
}

void rejectQoIName(const char* name)
{
  fprintf(stderr,"AWR problem error\n");
  fprintf(stderr,"Unknown quantity of interest name *%s*\n",name);
  abort();
}

QoI* createQoI(ParameterList& p, apf::Mesh* m, apf::Field* f)
{
  std::string name = p.get<std::string>("Name","");
  std::cout << name << std::endl;
  QoI* qoi;
  if (name == "Domain Integral")
    qoi = new DomainIntegral(p,m,f);
  else
    rejectQoIName(name.c_str());
  return qoi;
}

}
