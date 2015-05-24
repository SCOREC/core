#ifndef DSP_SMOOTHERS_H
#define DSP_SMOOTHERS_H

#include <apfMesh.h>
#include <set>

namespace dsp {
  
  typedef std::set<apf::ModelEntity*> Boundary;
  
  class Smoother {
  public:
    virtual ~Smoother();
    virtual void preprocess(apf::Mesh* m, Boundary& fixed, Boundary& moving, vector < apf::MeshEntity* >& V_total, vector < apf::Vector3 >& D_total, int& in_0, int& fb_0);
    virtual void smooth(apf::Field* df, vector < apf::MeshEntity* >& V_total, vector < apf::Vector3 >& D_total, int& in_0, int& fb_0) = 0;
    virtual void cleanup();
    static Smoother* makeLaplacian();
    static Smoother* makeSemiSpring();
    static Smoother* makeEmpty();
  };
  
}

#endif