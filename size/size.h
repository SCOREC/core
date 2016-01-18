/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef SIZE_H
#define SIZE_H

/** \file size.h
  * \brief Size field specification interface */

namespace apf {
class Field;
}

/** \namespace size
  * \brief All size field specification functions */
namespace size {

/** \namespace mech
  * \details size field specification for mechanics residuals */
namespace mech
{

/** \brief input for mechanics error estimation */
struct Input
{
  /** \brief the relevant stress measure
    * \details Cauchy stress tensor for small strain problems and the first
    * Piola-Kirchhoff stress tensor for finite strain problems */
  apf::Field* stress;
  /** \brief the solution to an auxiliary adjoint problem */
  apf::Field* adjoint;
};

/** \brief estimate the error for a mechanics residual
  * \param in user-constructed input
  * \returns updated field value for the adjoint-weighted residual
  * \details uses an adjoint weighted residual method. problems with
  * Neumann boudnary conditions are not yet supported */
apf::Field* estimateError(Input const& in);

}

}

#endif
