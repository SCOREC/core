/******************************************************************************

  (c) 2004-2010 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#ifndef LIIPBMOD_H
#define LIIPBMOD_H
#include <apf.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <vector>

using namespace std;

void liipbmod_commuInt(int *ns, int **nr, map<int,int> *Neigbors, int numParts);
void liipbmod_commuDouble(double *ns, double **nr, map<int,int> *Neigbors, int numParts);

class LIIPBMod
{
private:
   double tolerance1;
   double tolerance2;
   double tolerance3;
   int IterMax;
   int numVregionMax;
   int verbose;
public:
  LIIPBMod()
  {
     tolerance1 = 1.05;
     tolerance2 = 0.02;
     tolerance3 = 1.045;
     IterMax = 300;
     numVregionMax = 5;
     verbose = 0;
  }

  LIIPBMod(int niter, int adjElm, double tol1, double tol2, double tol3, int v)
  {
     tolerance1 = tol1;
     tolerance2 = tol2;
     tolerance3 = tol3;
     IterMax = niter;
     numVregionMax = adjElm;
     verbose = v;
  }
  void setMaxIter(int niter) {IterMax = niter;}
  void setMaxAdjElm(int adjElm) {numVregionMax = adjElm;}
  void setTol1(double tol1) {tolerance1 = tol1;}
  void setTol2(double tol2) {tolerance2 = tol2;}
  void setTol3(double tol3) {tolerance3 = tol3;} 
  int run(apf::Mesh* mesh);
};

#endif
