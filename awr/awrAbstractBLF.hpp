/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWR_ABSTRACTBLF_H
#define AWR_ABSTRACTBLF_H

#include "apf.h"

namespace awr {

/** \brief Abstract interface for bilinear forms */
class AbstractBLF
{
  /** default constructor */
  AbstractBLF() {};

  /** destructor */
  virtual ~AbstractBLF() {};

  /** evaluate discretized element-level bilinear form */
  virtual void
  evaluateBLF(const apf::MeshEntity* element,
              const apf::MeshEntity* enriched_solution,
              apf::DynamicVector<double> blf) = 0;
}

}

#endif
