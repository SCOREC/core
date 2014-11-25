#ifndef VIZ_H
#define VIZ_H

#include <string>

namespace apf {
  class Mesh;
  class MeshEntity;
};

enum Color {
  RED=255*256*256,
  BLUE=255,
  GREEN=255*256,
  BLACK=0,
  WHITE=255*256*256+255*256+255,
  GREY=128*256*256+128*256+128,
  BYPART=-1,
  MISCOLOR=-2,
  NOCOLOR=-3
  
};

struct milo;

class Visualization {
public:

  //Usage of Visualization
  void new_viz(int num_parts,Color color = BLACK);
  void breakpoint(std::string text ="");
  void end_viz();

  //Display entities
  bool watchEntity(apf::Mesh* m, apf::MeshEntity* ent, Color color = NOCOLOR);
  bool watchDownwardEntity(apf::Mesh* m, apf::MeshEntity* ent,Color color=NOCOLOR);
  bool watchDimension(apf::Mesh* m, int d,Color color=NOCOLOR);
  bool watchBoundary(apf::Mesh* m,int d,Color color=NOCOLOR);
  bool watchMesh(apf::Mesh* m);
  

  //Other Methods
  bool showAxis(Color x_color=RED,Color y_color=GREEN,Color z_color=BLUE);
  bool markPart(apf::Mesh* m);
  bool markPart(apf::Mesh* m, char* text);

  
private:
  milo* mil;
  int max_parts;
  int color_mode;
  double background[3];

  bool getPoint(apf::Mesh* m, apf::MeshEntity* ent, double* point);
  void getPartColor(double* color, int part_num);
  void getMISColor(double* color, int part_num);
  void getGivenColor(Color color, double* color_array);
  void getColor(Color color, double* color_array,int partId);
  bool drawPoint(apf::Mesh* m, apf::MeshEntity* ent,Color color);
  bool drawLine(apf::Mesh* m, apf::MeshEntity* ent,Color color);
  bool drawTriangle(apf::Mesh* m, apf::MeshEntity* ent, Color color);
};

#endif
