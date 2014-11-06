/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrPoisson.h"
#include "awrLinearSystem.h"
#include <apfNumbering.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

Problem::Problem(ParameterList& p, apf::Mesh* m) :
  problemList_(p.sublist("Adjoint Problem")),
  qoiList_(p.sublist("Quantity of Interest")),
  bcList_(p.sublist("Boundary Conditions")),
  mesh_(m)
{
  integrationOrder_ = problemList_.get("Integration Order",2);
}

Problem::~Problem()
{
  delete ls_;
  apf::destroyGlobalNumbering(globalNumbering_);
}

void rejectProblemName(const char* name)
{
  fprintf(stderr,"AWR problem error\n");
  fprintf(stderr,"Unknown problem name *%s*\n",name);
  abort();
}

void rejectSublists(const char* msg)
{
  fprintf(stderr,"AWR sublist error\n");
  fprintf(stderr,"sublist *%s* is not defined\n",msg);
  abort();
}

void validateSublists(ParameterList& p)
{
  if (! p.isSublist("Adjoint Problem"))
    rejectSublists("Adjoint Problem");
  if (! p.isSublist("Quantity of Interest"))
    rejectSublists("Quantity of Interest");
  if (! p.isSublist("Boundary Conditions"))
    rejectSublists("Boundary Conditions");
}

Problem* createProblem(ParameterList& p, apf::Mesh* m)
{
  validateSublists(p);
  ParameterList& pp = p.sublist("Adjoint Problem");
  std::string name = pp.get<std::string>("Name","");
  Problem* problem;
  if (name == "Poisson")
    problem = new PoissonProblem(p,m); 
  else
    rejectProblemName(name.c_str());
  return problem;
}

}
