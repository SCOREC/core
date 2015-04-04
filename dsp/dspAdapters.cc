#include "dspAdapters.h"
#include <ma.h>

namespace dsp {

Adapter::~Adapter()
{
}

class UniformAdapter : public Adapter {
  public:
  UniformAdapter(double s):myFunction(s) {}
  class MyFunction : public ma::IsotropicFunction {
    public:
    MyFunction(double s):size(s) {}
    virtual double getValue(ma::Entity* v)
    {
      (void)v;
      return size;
    }
    double size;
  };
  virtual void adapt(apf::Mesh2* m)
  {
    ma::adapt(m, &myFunction);
  }
  MyFunction myFunction;
};

Adapter* Adapter::makeUniform(double size)
{
  return new UniformAdapter(size);
}

}
