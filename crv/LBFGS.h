
#ifndef LBFGS_H
#define LBFGS_H

#include <iostream>
#include <pcu_util.h>
#include <vector>

namespace crv{

class ObjFunction
{
  public:
  ObjFunction(){};
  
  public:  
  virtual double getValue(std::vector<double> &x) = 0;
  virtual std::vector<double> getGrad(std::vector<double> &x) = 0;
  //virtual int getVectorDim() = 0;
  virtual int getSpaceDim() = 0;
};

class LBFGS
{
public:
  LBFGS(double inTol, int inIter, const std::vector<double> &x, ObjFunction *inObjFunc):
    tol(inTol), iter(inIter), x0(x), objFunc(inObjFunc) {}

  ~LBFGS() {}

public:
  //void setInitialValue(std::vector<double> x);
  std::vector<double> getCurrentX();
  double getFvalueAfter();
  double lineSearch(std::vector<double> &xold, std::vector<double> &g, std::vector<double> &direction, double stpmax);
  void moveArrayToLeft(std::vector<double> a[], int r);
  bool run();

public:
  double tol;
  int iter;
  std::vector<double> x0;
  ObjFunction *objFunc;
  std::vector<double> currentX;
  double fValAfter;

private:
  int r = 20;
};

}

#endif
