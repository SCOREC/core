#ifndef MA_FIRST_PROBLEM_PLANE_H
#define MA_FIRST_PROBLEM_PLANE_H

#include "maAdapt.h"
#include "maMesh.h"

namespace ma {

struct Ray
{
  Vector start;
  Vector dir;
};

enum class ProblemType
{
  TWOLARGEANGLES,
  THREELARGEANGLES,
};

class FirstProblemPlane
{
  public:
    FirstProblemPlane(Adapt* a, Vector target);
    void setVertex(Entity* v);
    void setBadElements(apf::Up& badElements);
    void getCandidateEdges(std::vector<Entity*> &edges);
    Entity* vert;
    Entity* problemFace;
    std::vector<Entity*> commEdges;
    std::vector<Entity*> problemRegions;
    Entity* problemRegion;
    Vector target;
  
  private:
    Adapt* adapter;
    Vector intersection;
    double tol;
    bool find();
    void findCandidateEdges(std::vector<Entity*> &edges);
    bool intersectRayFace(const Ray& ray, const std::vector<Vector>& coords,
    	Vector& intersection, bool& isInf);
    void findCommonEdges(apf::Up& cpRegions);
};

FirstProblemPlane* getFPP(Adapt* a, Entity* vertex, Vector target, apf::Up& invalid);
ProblemType getTetStats(Adapt* a, Entity* vert, Entity* face, Entity* region, Entity* ents[4], double area[4]);
Entity* getTetFaceOppositeVert(Mesh* m, Entity* e, Entity* v);
void getFaceCoords(Mesh* m, Entity* face, std::vector<Vector>& coords);
Vector getCenter(Mesh* mesh, Entity* face);
bool isLowInHigh(Mesh* mesh, Entity* highEnt, Entity* lowEnt);

}

#endif