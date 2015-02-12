#include "viz.h"
#include <milo.h>
#include <apfMesh.h>
#include <PCU.h>

Visualization::Visualization(const char* server, unsigned int port,
    Color color, bool local) {
  mil = milo_new(server, port, local);
  getColor(color,background);
  max_parts = PCU_Comm_Peers();
  milo_clear(mil, background);
}

Visualization::Visualization(unsigned int port, Color color, bool local) {
  mil = milo_new("localhost", port, local);
  getColor(color,background);
  max_parts = PCU_Comm_Peers();
  milo_clear(mil, background);
}

void Visualization::breakpoint(std::string text) {
  if (!PCU_Comm_Self()) {
    fprintf(stdout,"Breakpoint %s\n",text.c_str());
  }
  milo_title(mil,text.c_str());
  milo_run(mil);
  milo_clear(mil,background);
}

Visualization::~Visualization() {
  milo_free(mil);
}

void Visualization::getPoint(apf::Mesh* m, apf::MeshEntity* ent, double* point) {
  if (getDimension(m,ent)!=0) {
    std::cout<<"Cannot get point of this dimension\n";
    return;
  }
  apf::Vector3 p;
  m->getPoint_(ent,0,p);
  p.toArray(point);
}

void Visualization::getPartColor(double* color) {
  float mid = max_parts/2;
  float r,g,b;
  r=g=b=0.0;
  int part_num = PCU_Comm_Self();
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

void Visualization::getGivenColor(Color color, double* color_array) {
  int temp = color;
  for (int i=2;i>=0;i--) {
    color_array[i]=(temp%256)/255.0;
    temp=temp/256;
    
  }
}

void Visualization::getColor(Color color, double* color_array) {
  if (color==BYPART) {
    getPartColor(color_array);
  }
  else 
    getGivenColor(color,color_array);
}

void Visualization::drawPoint(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  double point[3];
  getPoint(m,ent,point);
  double color_array[3];
  if (color==NOCOLOR)
    color=WHITE;
  getColor(color,color_array);
  milo_dot(mil,point,color_array);

}

void Visualization::drawLine(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  double point[6];
  apf::Downward dwnVtx;
  int nDwnVtx = m->getDownward(ent,0,dwnVtx);
  for (int i=0;i<nDwnVtx;i++) {
    getPoint(m,dwnVtx[i],point+i*3);
  }
  double color_array[3];
  if (color==NOCOLOR)
    color=WHITE;
  getColor(color,color_array);
  milo_line(mil,point,point+3,color_array,0);

}

void Visualization::drawTriangle(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  double point[9];

  apf::Downward dwnVtx;
  int nDwnVtx = m->getDownward(ent,0,dwnVtx);
  for (int i=0;i<nDwnVtx;i++) 
    getPoint(m,dwnVtx[i],point+i*3);
  
  double color_array[3];
  if (color==NOCOLOR) {
    color=BYPART;
  }
  getColor(color,color_array);
  milo_triangle(mil,point,point+3,point+6,color_array);

}

void Visualization::watchEntity(apf::Mesh* m, apf::MeshEntity* ent,Color color) {
  int d = getDimension(m,ent);
  if (d==0)
    drawPoint(m,ent,color);
  else if (d==1)
    drawLine(m,ent,color);
  else if (d==2)
    drawTriangle(m,ent,color);
  else
    std::cout<<"Cannot support this entity yet\n";

}
void Visualization::watchDownwardEntity(apf::Mesh* m,apf::MeshEntity*ent,Color color) {
  int dim = getDimension(m,ent);
  apf::Downward dwn;
  int nDwn = m->getDownward(ent,dim-1,dwn);
  for (int i=0;i<nDwn;i++)
    watchEntity(m,dwn[i],color);
}
void Visualization::watchDimension(apf::Mesh* m, int d,Color color) {
  apf::MeshIterator* itr = m->begin(d);
  apf::MeshEntity* ent;
  while ((ent = m->iterate(itr))!=0) {
    watchEntity(m,ent,color);
  }
}

void Visualization::watchBoundary(apf::Mesh* m,int d,Color color) {
  int mesh_dim = m->getDimension();
  apf::MeshIterator* itr = m->begin(d);
  apf::MeshEntity* ent;
  while ((ent= m->iterate(itr))!=0) {
    apf::ModelEntity* model_ent = m->toModel(ent);
    if (m->getModelType(model_ent)<mesh_dim)
      watchEntity(m,ent,color);
  }
}

void Visualization::watchMesh(apf::Mesh* m) {
  for (int i=0;i<3;i++) {
    watchDimension(m,i);
  }
}

void Visualization::showAxis(Color x_color,Color y_color,Color z_color) {
  double origin[3] = {0,0,0};
  double x_axis[3] = {1,0,0};
  double y_axis[3] = {0,1,0};
  double z_axis[3] = {0,0,1};
  double x_array[3];
  double y_array[3];
  double z_array[3];
  getColor(x_color,x_array);
  getColor(y_color,y_array);
  getColor(z_color,z_array);
  milo_line(mil,origin,x_axis,x_array,1);
  milo_line(mil,origin,y_axis,y_array,1);
  milo_line(mil,origin,z_axis,z_array,1);
}

void Visualization::markEnt(apf::Mesh* m, apf::MeshEntity* e,
    std::string text,Color color) {
  apf::Vector3 cen = getLinearCentroid(m,e);
  double point[3] = {0,0,0};
  cen.toArray(point);
  double color_array[3];
  getColor(color,color_array);
  milo_text(mil,point,text.c_str(),color_array);
}

void Visualization::markPart(apf::Mesh* m,std::string text,Color color) {
  int numPoints=0;
  double centroid[3]={0,0,0};
  apf::MeshIterator* itr = m->begin(0);
  apf::MeshEntity* ent;
  while((ent=m->iterate(itr))!=0) {
    double point[3] = {0,0,0};
    getPoint(m,ent,point);
    for (int i=0;i<3;i++)
      centroid[i]+=point[i];
    numPoints++;
  }
  for (int i=0;i<3;i++)
    centroid[i]/=numPoints;
  double color_array[3];
  getColor(color,color_array);
  milo_text(mil,centroid,text.c_str(),color_array);
}
