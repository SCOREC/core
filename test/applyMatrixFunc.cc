#include <apfMatrix.cc>
#include <iostream>

bool isClose(double a, double b)
{
  return fabs(a-b)<1E-8;
}
bool isClose(apf::Matrix3x3 a, apf::Matrix3x3 b)
{
  for (int i=0; i<3; ++i)
  {
    for(int j=0; j<3; ++j)
    {
      if(!isClose(a[i][j], b[i][j]))
        return false;
    }
  }
  return true;
}

double identityFunc(double a)
{
  return a;
}

int main()
{
  int result = 0;
  apf::Matrix3x3 tmp;

  apf::Matrix3x3 eye(1,0,0,0,1,0,0,0,1);
  apf::Matrix3x3 diag = eye*2;
  // Note we want our random matrices to be symmetric positive definite
  // so that we have all positive eigenvalues
  apf::Matrix3x3 randMat1(0.71334277, 0.23916786, 0.0265657 , 0.1655712 , 0.0240585 ,
       0.74938504, 0.86250081, 0.62461084, 0.9821244);
  apf::Matrix3x3 randMat2(0.9017294 , 0.19634398, 0.25874015, 0.82974116, 0.17699926,
       0.5478204 , 0.96691199, 0.59327342, 0.00767784);
  // make symmetric
  randMat1 = (apf::transpose(randMat1)+randMat1)*0.5;
  randMat2 = (apf::transpose(randMat2)+randMat2)*0.5;
  // make positive definite
  randMat1 = randMat1*apf::transpose(randMat1)+eye;
  randMat2 = randMat2*apf::transpose(randMat2)+eye;
  // verify with identity function
  applyMatrixFunc(eye, &identityFunc, tmp);
  result += !isClose(eye, tmp);
  applyMatrixFunc(diag, &identityFunc, tmp);
  result += !isClose(diag, tmp);
  applyMatrixFunc(randMat1, &identityFunc, tmp);
  result += !isClose(randMat1, tmp);
  applyMatrixFunc(randMat2, &identityFunc, tmp);
  result += !isClose(randMat2, tmp);

  // check that the sqrt functions give correct results
  applyMatrixFunc(eye, &sqrt, tmp);
  result += !isClose(eye, apf::transpose(tmp)*tmp);
  applyMatrixFunc(diag, &sqrt, tmp);
  result += !isClose(diag, apf::transpose(tmp)*tmp);
  applyMatrixFunc(randMat1, &sqrt, tmp);
  result += !isClose(randMat1, apf::transpose(tmp)*tmp);
  applyMatrixFunc(randMat2, &sqrt, tmp);
  result += !isClose(randMat2, apf::transpose(tmp)*tmp);

  return result;
}
