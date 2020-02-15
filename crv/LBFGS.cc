#include "LBFGS.h"
#include "crvObjectiveFunctions.h"
#include <PCU.h>
#include "apfMatrix.h"

namespace crv{

//void LBFGS::setInitialValue(std::vector<double> x)
//{
//  x0 = x;
//}

std::vector<double> LBFGS::getCurrentX()
{
  return currentX;
}

static double dotP(const std::vector<double> &v1, const std::vector<double> &v2)
{
  double sum = 0.0;
  for (std::size_t i = 0; i < v1.size(); i++) {
    sum += v1[i]*v2[i];
  }
  return sum;
}

double LBFGS::getFvalueAfter()
{
  return fValAfter;
}

double LBFGS::lineSearch(std::vector<double> &xold, std::vector<double> &g, std::vector<double> &direction, double stpmax)
{
  double alpha = 1.0e-4, tolOndeltaX = 1.0e-8;
  int itrs = 2000;
  double a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
  double rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
  int n = xold.size();
  
  std::vector<double> xnew(xold.size(), 0.0);
  sum = sqrt(dotP(direction,direction));
  if (sum > stpmax) {
    for (int i = 0; i < n; i++) {
      direction[i] = direction[i]*stpmax/sum;
    }
  }
//  for (int i = 0; i < n; i++) slope += g[i]*p[i];
  slope = dotP(g,direction);
//if (slope >= 0.0) std::cout<<"slope positive"<<std::endl;

  test = 0.0;
  for (int i = 0; i < n; i++) {
    temp = std::abs(direction[i])/std::max(std::abs(xold[i]),1.0);
    if (temp > test) test = temp;
  }
  alamin = tolOndeltaX/test;
  alam = 1.0;
  
  double fold = objFunc->getValue(xold);

  for (int k = 0; k < itrs; k++) {
      for (int i =0; i < n; i++) xnew[i] = xold[i] + alam*direction[i];
      
      double fnew = objFunc->getValue(xnew);
      if (alam < alamin || fnew <= fold + alpha*alam*slope) {
        return alam;
      }
      else {
        if (alam == 1.0)
          tmplam = -slope/(2.0* (fnew-fold-slope));
        else {
          rhs1 = fnew - fold - alam*slope;
          rhs2 = f2 - fold - alam2*slope;
          a = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam - alam2);
          b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2))/(alam-alam2);
          if (a == 0.0) tmplam = -slope/(2.0*b);
          else {
            disc = b*b - 3.0*a*slope;
            if (disc < 0.0) tmplam = 0.5*alam;
            else if (b <= 0.0) tmplam = (-b + sqrt(disc))/(3.0*a);
            else tmplam = -slope/(b + sqrt(disc));
          }
        if (tmplam > 0.5*alam) tmplam = 0.5*alam;
        }  
     }
     alam2 = alam;
     alam = std::max(tmplam, 0.1*alam);
  }
  return alam;
}

void LBFGS::moveArrayToLeft(std::vector<double> a[], int r)
{
  for (int i = 0; i < r-1; i++)
      a[i] = a[i+1];   
}

bool LBFGS::run()
{
  std::vector<double> p(x0.size(), 0.0);
  std::vector<double> xs[r];
  std::vector<double> gs[r];
  std::vector<double> gdiffs[r];
  double gammas[r];

  for(int i = 0; i < r; ++i)
  {
    xs[i] = std::vector<double>(x0.size(), 0.0);
    gs[i] = std::vector<double>(x0.size(), 0.0);
    gdiffs[i] = std::vector<double>(x0.size(), 0.0);
    gammas[i] = 0.0;
  }

  xs[0] = x0;
  gs[0] = objFunc->getGrad(x0);
  
  for (std::size_t i = 0; i < xs[0].size(); i++) p[i] = -gs[0][i];

  for (int k = 0; k < iter; k++) {
    int I = 0;
    int J = 0;
    if (k+1 < r) {
      I = k+1;
      J = k;
    }
    else {
      I = r-1;
      J = I-1;
      moveArrayToLeft(xs, r);
      moveArrayToLeft(gs, r);
    }
    double stpmax = (std::max(sqrt(dotP(p,p)), double(objFunc->getSpaceDim())));
    double lambda = lineSearch(xs[J], gs[J], p, stpmax);
    
    for (std::size_t j = 0; j < xs[I].size(); j++) 
      xs[I][j] = xs[J][j] + lambda * p[j];

    gs[I] = objFunc->getGrad(xs[I]);
    
    if ( I > 0) {
      for (std::size_t jj = 0; jj < gs[I].size(); jj++)
      	gdiffs[I-1][jj] = gs[I][jj] - gs[I-1][jj];
    }

    if ((dotP(gs[I],gs[I]) < tol) || (dotP(gdiffs[I-1], gdiffs[I-1]) < tol)) {
      currentX = xs[I];
      fValAfter = objFunc->getValue(xs[I]);
      //std::cout<<"number of LBFGS iterations: "<<k<<std::endl;
      return true;
    }
    
    for (std::size_t j = 0; j < xs[I].size(); j++) p[j] = -gs[I][j];

    for (int i = I-1; i >= 0; --i) {
      std::vector<double> s(xs[I].size(), 0.0);
      std::vector<double> y(xs[I].size(), 0.0);
      for( std::size_t j = 0; j < xs[i+1].size(); j++) {
        s[j] = xs[i+1][j] - xs[i][j];
        y[j] = gs[i+1][j] - gs[i][j];
      }
      double rho = 1 / dotP(s, y);
      gammas[i] = rho * dotP(s, p);
      for(std::size_t j = 0; j < p.size(); j++)
        p[j] = p[j] - gammas[i]*y[j];      
    }
    
    for (int i = 0; i <= I-1; i++) {
      std::vector<double> s(xs[I].size(), 0.0);
      std::vector<double> y(xs[I].size(), 0.0);
      for( std::size_t j = 0; j < xs[i+1].size(); j++) {
        s[j] = xs[i+1][j] - xs[i][j];
        y[j] = gs[i+1][j] - gs[i][j];
      }
      double rho = 1 / dotP(s, y);
      double phi = rho* dotP(y, p);
      for(std::size_t j = 0; j < p.size(); j++)
        p[j] = p[j] + s[j]* (gammas[i] - phi);
    }
    
    if (dotP(p, gs[I]) > 0)
      for (std::size_t j = 0; j < p.size(); j++) p[j] = -p[j];

  }
  return false;

}

} // end of namespace     

