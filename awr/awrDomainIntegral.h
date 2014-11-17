/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWR_DOMAININTEGRAL_H
#define AWR_DOMAININTEGRAL_H

namespace awr {

class QoI;
class DomainIntegrator;

class DomainIntegral : public QoI
{
  public:
    DomainIntegral(ParameterList& p, apf::Mesh* m, apf::Field* f);
    ~DomainIntegral();
  private:
    DomainIntegrator* integrator_;
    void validateQoIList();
    void createIntegrator();
    void processFe(apf::MeshEntity* e, apf::DynamicVector& Fe);
};

}

#endif
