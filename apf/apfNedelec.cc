/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfShape.h"
#include "apfMesh.h"
#include "apfFieldOf.h"
#include "apfElement.h"
#include "apfVectorElement.h"
#include <pcu_util.h>

#include <iostream>

namespace apf {

class Nedelec: public FieldShape {
  public:
    const char* getName() const { return "Nedelec"; }
    Nedelec(int order) : P(order) {
      // TODO
      ;
    }
    class Vertex : public apf::EntityShape
    {
    public:
      void getValues(apf::Mesh*, apf::MeshEntity*,
	  apf::Vector3 const&, apf::NewArray<double>& values) const
      {
	(void)values;
	// TODO inform the user that this is not implemented and abort()
      }
      void getLocalGradients(apf::Mesh*, apf::MeshEntity*,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      }
      int countNodes() const {return 0;}
      void alignSharedNodes(apf::Mesh*,
	  apf::MeshEntity*, apf::MeshEntity*, int order[])
      {
	(void)order;
      }
    };
    class Edge : public apf::EntityShape
    {
    public:
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<double>&) const
      {
	// TODO inform the user that this is not implemented and abort()
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
	// TODO inform the user that this is not implemented and abort()
      }
      int countNodes() const {return 0; // TODO update this}
      void alignSharedNodes(apf::Mesh*,
	  apf::MeshEntity*, apf::MeshEntity*, int order[])
      {
	(void)order;
      }
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	// TODO: to be completed
      }
      void getLocalVectorGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Matrix3x3>&) const
      {
      	// TODO: to be completed
      }
    };
    int getOrder() {return P;}
  private:
    int P;
};

}
