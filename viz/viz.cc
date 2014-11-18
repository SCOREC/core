#include "viz.h"
#include <milo.h>
#include <apfMesh.h>

void Visualization::new_viz(int num_parts,Color color) {
  mil = milo_new("localhost", 4242);
  double color_array[3];
  getGivenColor(color_array,color);
  color_mode=0;
  max_parts = num_parts;
  milo_clear(mil, color_array);
  milo_zoom(mil,80);
}

void Visualization::breakpoint() {
  milo_run(mil);
  double black[3] = {0,0,0};
  milo_clear(mil,black);
}

bool Visualization::getPoint(apf::Mesh* m, apf::MeshEntity* ent, double* point) {
  if (getDimension(m,ent)!=0) {
    std::cout<<"Cannot get point of this dimension\n";
    return false;
  }
  apf::Vector3 p;
  m->getPoint_(ent,0,p);
  p.toArray(point);
  return true;
}
void Visualization::getPartColor(double* color, int part_num) {
  float mid = max_parts/2;
  float r,g,b;
  r=g=b=0.0;
  if (part_num==0)
    r=1;
  else if (part_num==max_parts)
    b=1;
  else if (part_num<=mid) {
    r=(mid-part_num)/mid;
    g=part_num/mid;
  }
  else {
    g=(max_parts-part_num)/mid;
    b=(part_num-mid)/mid;
  }
  color[0]=r;
  color[1]=g;
  color[2]=b;
}
void Visualization::getGivenColor(double* color_array, Color color) {
  int temp = color;
  for (int i=2;i>=0;i--) {
    color_array[i]=(temp%256)/255.0;
    temp=temp/256;
    
  }
}
bool Visualization::drawPoint(apf::Mesh* m, apf::MeshEntity* ent) {
  double point[3];
  getPoint(m,ent,point);
  double color[3] = {1,1,1};
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
  double color[3]= {1,1,1};
  milo_line(mil,point,point+3,color);
  return true;
}

bool Visualization::drawTriangle(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  double point[9];

  apf::Downward dwnVtx;
  int nDwnVtx = m->getDownward(ent,0,dwnVtx);
  for (int i=0;i<nDwnVtx;i++) 
    getPoint(m,dwnVtx[i],point+i*3);
  
  double color_array[3];
  if (color==BYPART) {
    int part = m->getOwner(ent);
    getPartColor(color_array,part);
  }
  else 
    getGivenColor(color_array,color);
  /*int part = m->getOwner(ent);
  if (part==0) {
    for (int i=0;i<3;i++)
      std::cout<<color_array[i]<<" ";
    std::cout<<std::endl;
    }*/
  milo_triangle(mil,point,point+3,point+6,color_array);
  return true;
}

bool Visualization::watchEntity(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  int d = getDimension(m,ent);
  if (d==0)
    return drawPoint(m,ent);
  else if (d==1)
    return drawLine(m,ent);
  else if (d==2)
    return drawTriangle(m,ent,color);
  else
    std::cout<<"Cannot support this entity yet\n";
  return false;
}
bool Visualization::watchDownwardEntity(apf::Mesh* m,apf::MeshEntity*ent) {
  int dim = getDimension(m,ent);
  apf::Downward dwn;
  int nDwn = m->getDownward(ent,dim-1,dwn);
  for (int i=0;i<nDwn;i++)
    watchEntity(m,dwn[i]);
  return true;
}
bool Visualization::watchDimension(apf::Mesh* m, int d) {
  apf::MeshIterator* itr = m->begin(d);
  apf::MeshEntity* ent;
  while ((ent = m->iterate(itr))!=0) {
    watchEntity(m,ent);
  }
  return true;
}

bool Visualization::watchBoundary(apf::Mesh* m,int d) {
  int mesh_dim = m->getDimension();
  apf::MeshIterator* itr = m->begin(d);
  apf::MeshEntity* ent;
  while ((ent= m->iterate(itr))!=0) {
    apf::ModelEntity* model_ent = m->toModel(ent);
    if (m->getModelType(model_ent)<mesh_dim)
      watchEntity(m,ent);
  }
}
bool Visualization::watchMesh(apf::Mesh* m) {
  for (int i=0;i<3;i++) {
    watchDimension(m,i);
  }
  return true;
}

void Visualization::end_viz() {
  milo_free(mil);
}
