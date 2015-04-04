#ifndef DSP_SMOOTHERS_H
#define DSP_SMOOTHERS_H

#include <apfMesh.h>
#include <set>

namespace dsp {

typedef std::set<apf::ModelEntity*> Boundary;

class Smoother {
  public:
    virtual ~Smoother();
    virtual void smooth(apf::Field* df, Boundary& fixed, Boundary& moving) = 0;
    static Smoother* makeLaplacian();
};

}

#endif
