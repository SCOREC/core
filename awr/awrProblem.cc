/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrPoisson.h"
#include "awrNonlinearPoisson.h"
#include "awrQoI.h"
#include "awrLinearSystem.h"
#include <PCU.h>
#include <apfNumbering.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

void print(const char* format, ...)
{
  if (PCU_Comm_Self())
    return;
  printf("\nAWR: ");
  va_list ap;
  va_start(ap,format);
  vfprintf(stdout,format,ap);
  va_end(ap);
  printf("\n");
}

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
  delete qoi_;
  apf::destroyGlobalNumbering(globalNumbering_);
}

void rejectProblemName(const char* name)
{
  print("Unknown problem name *%s*\n",name);
  abort();
}

void rejectSublists(const char* msg)
{
  print("sublist *%s* is not defined\n",msg);
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
  print("solving %s adjoint problem",name.c_str());
  Problem* problem;
  if (name == "Poisson")
    problem = new PoissonProblem(p,m); 
  else if (name == "Nonlinear Poisson")
    problem = new NonlinearPoissonProblem(p,m);
  else
    rejectProblemName(name.c_str());
  return problem;
}

}
