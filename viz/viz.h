#ifndef VIZ_H
#define VIZ_H

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
  BYPART=-1
};

struct milo;

class Visualization {
public:
  void new_viz(int num_parts,Color color = BLACK);
  void breakpoint();
  bool drawPoint(apf::Mesh* m, apf::MeshEntity* ent);
  bool drawLine(apf::Mesh* m, apf::MeshEntity* ent);
  bool drawTriangle(apf::Mesh* m, apf::MeshEntity* ent, Color color);
  bool watchEntity(apf::Mesh* m, apf::MeshEntity* ent, Color color = BYPART);
  bool watchDownwardEntity(apf::Mesh* m, apf::MeshEntity* ent);
  bool watchDimension(apf::Mesh* m, int d);
  bool watchBoundary(apf::Mesh* m,int d);
  bool watchMesh(apf::Mesh* m);
  void end_viz();
private:
  milo* mil;
  int max_parts;
  int color_mode;

  bool getPoint(apf::Mesh* m, apf::MeshEntity* ent, double* point);
  void getPartColor(double* color, int parts);
  void getGivenColor(double* color_array, Color color);

};

#endif
