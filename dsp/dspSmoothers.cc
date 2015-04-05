#include "dspSmoothers.h"
#include <apf.h>

namespace dsp {

Smoother::~Smoother()
{
}

class LaplacianSmoother : public Smoother {
  public:
    void smooth(apf::Field* df, Boundary& fixed, Boundary& moving)
    {
      apf::Mesh* m = apf::getMesh(df);
      (void)m;
      (void)df;
      (void)fixed;
      (void)moving;
    }
};

Smoother* Smoother::makeLaplacian()
{
  return new LaplacianSmoother();
}

}
