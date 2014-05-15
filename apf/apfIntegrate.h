/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFINTEGRATE_H
#define APFINTEGRATE_H

#include "apfVector.h"

namespace apf {

struct IntegrationPoint
{
  IntegrationPoint(Vector3 const& p, double w):
    param(p),weight(w)
  {}
  Vector3 param;
  double weight;
};

class Integration
{
  public:
    virtual ~Integration() {}
    virtual int countPoints() const = 0;
    virtual IntegrationPoint const* getPoint(int i) const = 0;
    virtual int getAccuracy() const = 0;
};

class EntityIntegration
{
  public:
    virtual ~EntityIntegration() {}
    Integration const* getAccurate(int minimumAccuracy) const;
    virtual int countIntegrations() const = 0;
    virtual Integration const* getIntegration(int i) const = 0;
};

EntityIntegration const* getIntegration(int meshEntityType);

}//namespace apf

#endif
