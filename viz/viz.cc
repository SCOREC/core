#include "viz.h"
#include <milo.h>
#include <apfMesh.h>

void Visualization::new_viz() {
  mil = milo_new(200, 200, "localhost", 4242);
  double black[3] = {0,0,0};
  milo_clear(mil, black);
  milo_zoom(mil,80);
}

void Visualization::breakpoint() {
  //milo_render(mil);
  //milo_write_png(mil,"picture");
  milo_run(mil);
}

bool Visualization::getPoint(apf::Mesh* m, apf::MeshEntity* ent, double* point) {
  if (getDimension(m,ent)!=0) {
    std::cout<<"Cannot get point of this dimension\n";
    return false;
  }
  apf::Vector3 p;
  m->getPoint(ent,0,p);
  p.toArray(point);
  return true;
}

bool Visualization::drawPoint(apf::Mesh* m, apf::MeshEntity* ent) {
  double point[3];
  getPoint(m,ent,point);
  double color[3] = {1,1,0};
  milo_dot(mil,point,color);
  return true;
}

bool Visualization::drawLine(apf::Mesh* m, apf::MeshEntity* ent) {
  double point[6];
  apf::Downward dwnVtx;
  int nDwnVtx = m->getDownward(ent,0,dwnVtx);
  for (int i=0;i<nDwnVtx;i++) {
    getPoint(m,dwnVtx[i],point+i*3);
  }
  double color[3] = {0,1,0};
  milo_line(mil,point,point+3,color);
  return true;
}

bool Visualization::drawTriangle(apf::Mesh* m, apf::MeshEntity* ent) {
  double point[9];
  /*
  apf::Downward dwnEdge;
  int nDwnEdge = m->getDownward(ent,1,dwnEdge);
  for (int i=0;i<nDwnEdge;i++) {
    apf::Downward dwnVtx;
    m->getDownward(dwnEdge[i],0,dwnVtx);
    getPoint(m,dwnVtx[0],point+i*3);
  }
  */
 
  apf::Downward dwnVtx;
  int nDwnVtx = m->getDownward(ent,0,dwnVtx);
  for (int i=0;i<nDwnVtx;i++) 
    getPoint(m,dwnVtx[i],point+i*3);

  double color[3] = {0,1,1};
  milo_triangle(mil,point,point+3,point+6,color);
  return true;
}

bool Visualization::watchEntity(apf::Mesh* m, apf::MeshEntity* ent) {
  int d = getDimension(m,ent);
  if (d==0)
    return drawPoint(m,ent);
  else if (d==1)
    return drawLine(m,ent);
  else if (d==2)
    return drawTriangle(m,ent);
  else
    std::cout<<"Cannot support this entity yet\n";
  return false;
}

bool Visualization::watchDimension(apf::Mesh* m, int d) {
  apf::MeshIterator* itr = m->begin(d);
  apf::MeshEntity* ent;
  while ((ent = m->iterate(itr))!=0) {
    watchEntity(m,ent);
  }
  return true;
}

bool Visualization::watchMesh(apf::Mesh* m) {
  return true;
}

void Visualization::end_viz() {
  milo_free(mil);
}
