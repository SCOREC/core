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
//  unsigned long bubblecount = 0; //May be used later but the bubble id has to be read for now from the input file
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
    // File format (each line represents a bubble): x_center y_center z_center radius
    fscanf(filebubble, "%d %lf %lf %lf %lf", &readbubble.id, &readbubble.center[0], &readbubble.center[1], &readbubble.center[2], &readbubble.radius);
    if(feof(filebubble)) break;
//    bubblecount++;
//    readbubble.id = bubblecount;
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

  int bubbleid = 0; // Initialization critical here
  double distx;
  double disty;
  double distz;
  double tmpdist;
  double distance = 1e99; // Initialization critical here

  /* search through bubbles, 
     find the distance to the nearet bubble membrane (sphere) */
  for(unsigned long i=0; i<bubbles.size(); i++) 
  {
    distx = (v_center[0]-bubbles[i].center[0]);
    disty = (v_center[1]-bubbles[i].center[1]);
    distz = (v_center[2]-bubbles[i].center[2]);
    tmpdist = sqrt(distx*distx + disty*disty + distz*distz) - bubbles[i].radius;
    if(tmpdist < distance)
    {
      distance = tmpdist;
      if (distance < 0) 
      {
        bubbleid = bubbles[i].id;
        break; //if v is inside a bubble, stop searching since bubbles should not intersect each other
      }
      // A negative bubble id could be useful to know the nearest bubble of a point in the liquid phase.
      // However, phasta is not ready to handle this and the feature is not critical so ignore for now.
      // Ths can be pratical for debug purpose though.
//      else
//        bubbleid = -bubbles[i].id;    
    }
  }

//  debug
//   printf("coord: %lf %lf %lf - Dist: %lf - Bubble id: %d\n", coord[0], coord[1], coord[2], distance, bubbleid);

  sol[5] = distance;
  sol[6] = static_cast<double>(bubbleid);

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

