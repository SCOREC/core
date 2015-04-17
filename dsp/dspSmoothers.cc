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
      /* this needs to be filled in with the actual
         laplacian nodal averaging algorithm */
      (void)m;
      (void)df;
      (void)fixed;
      (void)moving;
    }
};

class EmptySmoother : public Smoother {
  public:
    void smooth(apf::Field* df, Boundary& fixed, Boundary& moving)
    {
      (void)df;
      (void)fixed;
      (void)moving;
    }
};

Smoother* Smoother::makeLaplacian()
{
  return new LaplacianSmoother();
}

Smoother* Smoother::makeEmpty()
{
  return new EmptySmoother();
}

}
