#ifndef DSP_ADAPTERS_H
#define DSP_ADAPTERS_H

#include <apfMesh2.h>

namespace dsp {

class Adapter {
  public:
    virtual ~Adapter();
    virtual void adapt(apf::Mesh2* m) = 0;
    static Adapter* makeUniform(double size);
};

}

#endif
