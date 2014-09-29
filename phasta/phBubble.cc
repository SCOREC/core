#include <PCU.h>
#include "phBubble.h"
#include "phInput.h"
#include <apfMesh.h>
#include <apf.h>
#include <stdio.h>

namespace ph {

struct Bubble {
  int id;
  apf::Vector3 center;
  double radius;
};

typedef std::vector<Bubble> Bubbles;

void readBubbles(Bubbles& bubbles)
{
  char bubblefname[256];
  FILE *filebubble;
  Bubble readbubble;

  sprintf(bubblefname,"bubbles.inp");
  if (!PCU_Comm_Self())
    printf("reading bubbles info from %s\n",bubblefname);

  filebubble = fopen(bubblefname, "r");
  assert(filebubble != NULL); 
  while(1)
  {
    fscanf(filebubble, "%d %lf %lf %lf %lf", &readbubble.id, &readbubble.center[0], &readbubble.center[1], &readbubble.center[2], &readbubble.radius);
    if(feof(filebubble)) break;
    bubbles.push_back(readbubble);
  }
  fclose(filebubble);

  if (!PCU_Comm_Self())
    printf("%lu bubbles found in %s\n", bubbles.size(), bubblefname);

// Debug
/*
  for(unsigned long i=0; i<bubbles.size(); i++)
  {
    printf("%d %lf %lf %lf %lf\n", bubbles[i].id, bubbles[i].center[0], bubbles[i].center[1], bubbles[i].center[2], bubbles[i].radius);
  }
*/

}

void setBubbleScalars(apf::Mesh* m, apf::MeshEntity* v,
    Bubbles& bubbles, double* sol)
{
  apf::Vector3 v_center;
  m->getPoint(v, 0, v_center);
  /* search through bubbles, etc... */
  sol[5] = 42;
  sol[6] = 42;
}

void initBubbles(apf::Mesh* m, Input& in)
{
  Bubbles bubbles;
  readBubbles(bubbles);
  assert(in.ensa_dof >= 7);
  apf::NewArray<double> s(in.ensa_dof);
  apf::Field* f = m->findField("solution");
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::getComponents(f, v, 0, &s[0]);
    setBubbleScalars(m, v, bubbles, &s[0]);
    apf::setComponents(f, v, 0, &s[0]);
  }
  m->end(it);
}

}

