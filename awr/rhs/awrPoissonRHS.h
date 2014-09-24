/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRPOISSONRHS_H
#define AWRPOISSONRHS_H

#include "awrRHS.h"
#include "apfMatrix.h"
#include "apfDynamicVector.h"

namespace awr {

class PoissonRHS : public RHS
{
  public:
    /** constructor */
    PoissonRHS(const Teuchos::ParameterList& p);
    /** destructor */
    virtual ~PoissonRHS() {};
    /** evaluate element level stiffness matrix */
    virtual void
    evaluateElementRHS(apf::MeshEntity* element,
                       apf::Field* primal_solution,
                       int integration_order,
                       apf::DynamicMatrix& k);
  protected:
    /** number of spatial dimensions */
    int num_dims_;
    /** number of element nodes */
    int num_nodes_;
    /** number of element quadrature points */
    int num_qp_;
    /** shape functions */
    apf::DynamicVector N;
    /** shape function derivatives */
    apf::DynamicMatrix dN;
    /** element Jacobian transformation */
    apf::Matrix3x3 jac;
  private:
    PoissonRHS(const PoissonRHS&);
    PoissonRHS& operator=(const PoissonRHS&);
};

}

#endif
