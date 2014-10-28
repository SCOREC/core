#include "viz.h"

void Visualization::new_viz() {
  mil = milo_new(200, 200, "localhost", 4242);
      double black[3] = {0,0,0};
  milo_clear(mil, black);
}

void Visualization::breakpoint() {
  
  milo_render(mil);
  //milo_write_png(mil,"picture");
}

bool Visualization::watchEntity(apf::Mesh* m, apf::MeshEntity* ent) {
  if (getDimension(m,ent)>0) {
    std::cout<<"Cannot do this dimension yet\n";
    return false;
  }
  apf::Vector3 point;
  m->getPoint(ent,0,point);
  double p[3]={0,0,0};
  //point.toArray(p);

  double color[3];
  double green[3] = {0,1,0};
  milo_dot(mil,p,green);
}

bool Visualization::watchMesh(apf::Mesh* m) {
  
}

void Visualization::end_viz() {
  milo_free(mil);
}
