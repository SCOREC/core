#include <milo.h>
#include <apfMesh.h>

class Visualization {
public:

  void new_viz() {
    mil = milo_new(200, 200, "localhost", 4242);
    milo_zoom(mil, 40);
    milo_spin(mil, -3 * M_PI / 8);
    milo_tilt(mil, M_PI / 9);
  
    double black[3] = {0,0,0};
    milo_clear(mil, black);
    
  }

  void breakpoint() {
    //milo_render(mil);
    //milo_write_png(mil,"picture");
    double point[3] = {0,0,0};
    double color[3] = {1,0,0};
    milo_text(mil, point, "debug", color);
    std::cout<<"running now\n"<<std::endl;
    milo_run(mil);
  }

bool watchEntity(apf::Mesh* m, apf::MeshEntity* ent) {
  if (getDimension(m,ent)!=0) {
    std::cout<<"Cannot do this dimension yet\n";
    return false;
  }
  apf::Vector3 point;
  m->getPoint(ent,0,point);
  double p[3];
  point.toArray(p);
  double green[3] = {0,1,0};
  milo_dot(mil,p,green);
  return true;
}

bool watchMesh(apf::Mesh* m) {
  return true;
}

void end_viz() {
  milo_free(mil);
}

private:
  milo_t mil;
};
