#ifndef DSP_SMOOTHERS_H
#define DSP_SMOOTHERS_H

#include <apfMesh.h>
#include <set>
#include <vector>

using namespace std;

namespace dsp {
  
  typedef std::set<apf::ModelEntity*> Boundary;
  
  class Smoother {
  public:
    virtual ~Smoother();
    virtual void preprocess(apf::Mesh* m, Boundary& fixed, Boundary& moving, vector < apf::MeshEntity* >& V_total,
        int& in_0, int& fb_0);
    virtual void smooth(apf::Field* df, vector < apf::MeshEntity* >& V_total, int in_0, int fb_0) = 0;
    virtual void cleanup(apf::Mesh* m);
    static Smoother* makeLaplacian();
    static Smoother* makeSemiSpring();
    static Smoother* makeElastic();
    static Smoother* makeEmpty();
  };
  
}

#endif
