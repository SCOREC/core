#include "viz.h"
#include "diffMC/maximalIndependentSet/mis.h"
#include <milo.h>
#include <apfMesh.h>

using namespace misLuby;
void Visualization::new_viz(int num_parts,Color color) {
  mil = milo_new("localhost", 4242);

  getGivenColor(color,background);
  mis_color=BLUE;
  max_parts = num_parts;
  milo_clear(mil, background);
}

void Visualization::breakpoint(std::string text) {
  milo_title(mil,text.c_str());
  milo_run(mil);
  milo_clear(mil,background);
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
void Visualization::getMISColor(double* color) {
  getGivenColor(mis_color,color);
}
void Visualization::getGivenColor(Color color, double* color_array) {
  int temp = color;
  for (int i=2;i>=0;i--) {
    color_array[i]=(temp%256)/255.0;
    temp=temp/256;
    
  }
}
void Visualization::getColor(Color color, double* color_array,int partId) {
  if (color==BYPART) {
    getPartColor(color_array,partId);
  }
  else if (color==MISCOLOR) {
    getMISColor(color_array);
  }
  else 
    getGivenColor(color,color_array);
}
bool Visualization::drawPoint(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  double point[3];
  getPoint(m,ent,point);
  double color_array[3];
  if (color==NOCOLOR)
    color=WHITE;
  int partId = m->getOwner(ent);
  getColor(color,color_array,partId);
  milo_dot(mil,point,color_array);
  return true;
}

bool Visualization::drawLine(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  double point[6];
  apf::Downward dwnVtx;
  int nDwnVtx = m->getDownward(ent,0,dwnVtx);
  for (int i=0;i<nDwnVtx;i++) {
    getPoint(m,dwnVtx[i],point+i*3);
  }
  double color_array[3];
  if (color==NOCOLOR)
    color=WHITE;
  int partId = m->getOwner(ent);
  getColor(color,color_array,partId);
  milo_line(mil,point,point+3,color_array,0);
  return true;
}

bool Visualization::drawTriangle(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  double point[9];

  apf::Downward dwnVtx;
  int nDwnVtx = m->getDownward(ent,0,dwnVtx);
  for (int i=0;i<nDwnVtx;i++) 
    getPoint(m,dwnVtx[i],point+i*3);
  
  double color_array[3];
  if (color==NOCOLOR) {
    color=BYPART;
  }
  int partId = m->getOwner(ent);
  getColor(color,color_array,partId);
  milo_triangle(mil,point,point+3,point+6,color_array);
  return true;
}

bool Visualization::watchEntity(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  int d = getDimension(m,ent);
  if (d==0)
    return drawPoint(m,ent,color);
  else if (d==1)
    return drawLine(m,ent,color);
  else if (d==2)
    return drawTriangle(m,ent,color);
  else
    std::cout<<"Cannot support this entity yet\n";
  return false;
}
bool Visualization::watchDownwardEntity(apf::Mesh* m,apf::MeshEntity*ent,Color color) {
  int dim = getDimension(m,ent);
  apf::Downward dwn;
  int nDwn = m->getDownward(ent,dim-1,dwn);
  for (int i=0;i<nDwn;i++)
    watchEntity(m,dwn[i],color);
  return true;
}
bool Visualization::watchDimension(apf::Mesh* m, int d,Color color) {
  apf::MeshIterator* itr = m->begin(d);
  apf::MeshEntity* ent;
  while ((ent = m->iterate(itr))!=0) {
    watchEntity(m,ent,color);
  }
  return true;
}

bool Visualization::watchBoundary(apf::Mesh* m,int d,Color color) {
  int mesh_dim = m->getDimension();
  apf::MeshIterator* itr = m->begin(d);
  apf::MeshEntity* ent;
  while ((ent= m->iterate(itr))!=0) {
    apf::ModelEntity* model_ent = m->toModel(ent);
    if (m->getModelType(model_ent)<mesh_dim)
      watchEntity(m,ent,color);
  }
  return true;
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
bool Visualization::setupMISColoring(apf::Mesh* m,int part_num) {
  mis_color=NOCOLOR;
  int iter=0;
  Color color_choices[9] = {RED,BLUE,GREEN,PURPLE,ORANGE,YELLOW,BROWN,PINK,GREY};
  apf::Parts neighbors;
  apf::getFacePeers(m,neighbors);

  partInfo part;
  part.id = part_num;
  for (apf::Parts::iterator itr = neighbors.begin();
       itr!=neighbors.end();itr++) {
    part.adjPartIds.push_back(*itr);
    part.net.push_back(*itr);
  }

  int randNumSeed = part_num+1;
  mis_init(randNumSeed,false);
  int isInMis = mis(part, true, true);
  if (isInMis) { 
    mis_color= color_choices[iter];
  }

  return true;
}
bool Visualization::showAxis(Color x_color,Color y_color,Color z_color) {
  double origin[3] = {0,0,0};
  double x_axis[3] = {1,0,0};
  double y_axis[3] = {0,1,0};
  double z_axis[3] = {0,0,1};
  double x_array[3];
  double y_array[3];
  double z_array[3];
  getGivenColor(x_color,x_array);
  getGivenColor(y_color,y_array);
  getGivenColor(z_color,z_array);
  milo_line(mil,origin,x_axis,x_array,1);
  milo_line(mil,origin,y_axis,y_array,1);
  milo_line(mil,origin,z_axis,z_array,1);
  return true;
}
