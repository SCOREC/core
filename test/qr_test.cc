#include <spr.h>

int main(int argc, char** argv)
{
  apf::DynamicMatrix A(4,3);
  A(0,0) = 1;
  A(0,1) = 1;
  A(0,2) = 1;
  A(1,0) = 1;
  A(1,1) = 1;
  A(1,2) = 1;
  A(2,0) = 1;
  A(2,1) = 1;
  A(2,2) = 1;
  A(3,0) = 1;
  A(3,1) = 1;
  A(3,2) = 1;

  apf::DynamicVector x;
  
  apf::DynamicVector b(4);
  b(0) = 1;
  b(1) = 1;
  b(3) = 1;
  b(4) = 1;
  
  spr::solveQR(A,x,b);

}

