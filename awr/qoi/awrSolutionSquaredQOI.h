/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRSOLUTIONSQUAREDQOI_H
#define AWRSOLUTIONSQUAREDQOI_H

#include "awrQOI.h"
#include "../awrBasisUtils.h"

namespace awr {

class SolutionSquaredQOI : public QOI
{
  public:
    /** constructor */
    SolutionSquaredQOI(apf::Mesh* m, const Teuchos::ParameterList& p);
    /** destructor */
    virtual ~SolutionSquaredQOI() {};
    /** evaluate element level force vector */
    virtual void
    evaluateElementQOI(apf::MeshEntity* element,
                       apf::DynamicVector& f);
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
    /** initialize basis stuff */
    void init();
    /** validate parameters */
    void validateParameters();
  private:
    SolutionSquaredQOI(const SolutionSquaredQOI&);
    SolutionSquaredQOI& operator=(const SolutionSquaredQOI&);
};

}

#endif
