/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFVTK_H
#define APFVTK_H

#include "apfField.h"


namespace apf {

class HasAll : public FieldOp
{
  public:
    virtual bool inEntity(MeshEntity* e)
    {
      if (!f->getData()->hasEntity(e))
        ok = false;
      return false;
    }
    bool run(FieldBase* f_)
    {
      f = f_;
      ok = true;
      this->apply(f);
      return ok;
    }
  private:
    bool ok;
    FieldBase* f;
};



bool isPrintable(FieldBase* f);


} // namespace apf

#endif
