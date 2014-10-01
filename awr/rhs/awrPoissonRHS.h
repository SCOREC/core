/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRPOISSONRHS_H
#define AWRPOISSONRHS_H

#include "awrRHS.h"
#include "../awrBasisUtils.h"
#include "apfMatrix.h"
#include "apfDynamicVector.h"

namespace awr {

class PoissonRHS : public RHS
{
  public:
    /** constructor */
    PoissonRHS(apf::Mesh* m, const Teuchos::ParameterList& p);
    /** destructor */
    virtual ~PoissonRHS() {};
    /** evaluate element level stiffness matrix */
    virtual void
    evaluateElementRHS(apf::MeshEntity* element,
                       apf::DynamicMatrix& k);
  protected:
    /** number of spatial dimensions */
    int num_dims_;
    /** number of element nodes */
    int num_nodes_;
    /** number of element quadrature points */
    int num_qp_;
    /** integration order */
    int integration_order_;
    /** primal solution field */
    apf::Field* sol_;
    /** gradient of basis functions */
    NodeQPVector grad_bf_;
    /** weighted gradient of basis functions */
    NodeQPVector w_grad_bf_;
    /** initialize basis stuff */
    void init();
    /** validate parameters */
    void validateParameters();
  private:
    PoissonRHS(const PoissonRHS&);
    PoissonRHS& operator=(const PoissonRHS&);
};

}

#endif
