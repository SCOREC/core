#ifndef VIZ_H
#define VIZ_H

/** \file viz.h
 * \brief Visualization tool interface
 *
 */

#include <string>

namespace apf {
  class Mesh;
  class MeshEntity;
}

/**
 * \enum Color
 * \brief Predefined Colors
 * 
 * BYPART uses linear color scale from red to blue
 */
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
  NOCOLOR=-3
};

struct milo;

/** 
 * \class Visualization 
 * \brief API between the application and viewer
 */
class Visualization {
public:

  /**
   * @brief constructs the visualization and connects to viewer
   * @param port (In) port id of viewer
   * @param color (In) background color
   * @param local (In) use COMM_SELF
   */ 
  Visualization(unsigned int port = 4242, Color color = BLACK, bool local = false);

  /**
   * @brief constructs the visualization and connects to viewer
   * @param server (In) server name of viewer
   * @param port (In) port id of viewer
   * @param color (In) background color
   * @param local (In) use COMM_SELF
   */ 
  Visualization(const char* server, unsigned int port = 4242, 
      Color color = BLACK, bool local = false);

  /**
   * @brief suspends application and sends rendering to viewer
   * @param text (In) title of breakpoint
   */
  void breakpoint(std::string text ="");

  /**
   * @brief ends the visualization and cleans up memory
   */
  ~Visualization();

  /**
   * @brief marks an entity to be visualized at the next breakpoint
   * @param m (In) the mesh
   * @param ent (In) the entity to be viewed
   * @param color (In) the color that the entity will be in the viewer
   */
  void watchEntity(apf::Mesh* m, apf::MeshEntity* ent, Color color = NOCOLOR);

  /**
   * @brief marks the downward adjancencies of an entity to be visualized at the next breakpoint
   * @param m (In) the mesh
   * @param ent (In) the entity whose downward adjancent entities will be viewed
   * @param color (In) the color that the entities will be in the viewer
   */
  void watchDownwardEntity(apf::Mesh* m, apf::MeshEntity* ent,Color color=NOCOLOR);

  /**
   * @brief marks every entity of a dimension to be visualized at the next breakpoint
   * @param m (In) the mesh
   * @param d (In) the dimension to be visualized 
   * @param color (In) the color that the entities will be in the viewer
   */
  void watchDimension(apf::Mesh* m, int d,Color color=NOCOLOR);

  /**
   * @brief marks the model boundary of a given dimension to be visualized at the next breakpoint
   * @param m (In) the mesh
   * @param d (In) the dimension to be visualized
   * @param color (In) the color that the entities will be in viewer
   */
  void watchBoundary(apf::Mesh* m,int d,Color color=NOCOLOR);

  /**
   * @brief marks all vertices, edges, and faces of the mesh to be visualized
   * @param m (In) the mesh
   */
  void watchMesh(apf::Mesh* m);
  

  //Other Methods
  /**
   * @brief creates axis from the origin
   * @param x_color (In) the color of the x axis 
   * @param y_color (In) the color of the y axis
   * @param z_color (In) the color of the z axis
   */
  void showAxis(Color x_color=RED,Color y_color=GREEN,Color z_color=BLUE);

  /** 
   * @brief writes text at the centroid of the entity
   * @param m (In) the mesh
   * @param e (In) the mesh entity
   * @param text (In) the string of text
   * @param color (In) the color of the text
   */
  void markEnt(apf::Mesh* m, apf::MeshEntity* e, std::string text,
      Color color=BLACK);

  /** 
   * @brief writes text at the centroid of the part
   * @param m (In) the mesh
   * @param text (In) the string of text
   * @param color (In) the color of the text
   */
  void markPart(apf::Mesh* m,std::string text, Color color=BLACK);

private:
  milo* mil;
  int max_parts;
  double background[3];

  void getPoint(apf::Mesh* m, apf::MeshEntity* ent, double* point);
  void getPartColor(double* color);
  void getGivenColor(Color color, double* color_array);
  void getColor(Color color, double* color_array);
  void drawPoint(apf::Mesh* m, apf::MeshEntity* ent,Color color);
  void drawLine(apf::Mesh* m, apf::MeshEntity* ent,Color color);
  void drawTriangle(apf::Mesh* m, apf::MeshEntity* ent, Color color);
};

#endif
