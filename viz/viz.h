#ifndef VIZ_H
#define VIZ_H

#include <string>

namespace apf {
  class Mesh;
  class MeshEntity;
};

enum Color {
  RED=228*256*256+26*255+28, 
  BLUE=55*256*256+126*256+184, 
  GREEN=77*256*256+175*256+74, 
  PURPLE=152*256*256+78*256+163, 
  ORANGE=255*256*256+127*256, 
  YELLOW=255*256*256+255*256+51, 
  BROWN=166*256*256+86*256+40,
  PINK=247*256*256+129*256*256+191,
  GREY=153*256*256+153*256+153,
  BLACK=0,
  WHITE=255*256*256+255*256+255,
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
  bool setupMISColoring(apf::Mesh* m,int part_num);
  bool showAxis(Color x_color=RED,Color y_color=GREEN,Color z_color=BLUE);
  
private:
  milo* mil;
  int max_parts;
  Color mis_color;
  double background[3];

  bool getPoint(apf::Mesh* m, apf::MeshEntity* ent, double* point);
  void getPartColor(double* color, int part_num);
  void getMISColor(double* color);
  void getGivenColor(Color color, double* color_array);
  void getColor(Color color, double* color_array,int partId);
  bool drawPoint(apf::Mesh* m, apf::MeshEntity* ent,Color color);
  bool drawLine(apf::Mesh* m, apf::MeshEntity* ent,Color color);
  bool drawTriangle(apf::Mesh* m, apf::MeshEntity* ent, Color color);
};

#endif
