/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maEdgeSwap.h"
#include "maAdapt.h"
#include "maShape.h"
#include "maShapeHandler.h"
#include <cstdio>

namespace ma {

/* most of these tables are from old MeshAdapt's
   Eswp3.cc. */
/* they represent an encoding of all possible triangulations
   of an N-vertex polygon, with N ranging from 3 to 7.
   This implies that we cannot re-triangulate polygons with
   more than 7 vertices. */

#define MAX_VERTS 7

static int triangulation_size[MAX_VERTS+1] =
{0 //0 vertices
,0 //1 vertices
,0 //2 vertices
,1 //3 vertices
,2 //4 vertices
,3 //5 vertices
,4 //6 vertices
,5 //7 vertices
};

static int triangulation_count[MAX_VERTS+1] =
{0 //0 vertices
,0 //1 vertices
,0 //2 vertices
,1 //3 vertices
,2 //4 vertices
,5 //5 vertices
,14 //6 vertices
,42 //7 vertices
};

static int unique_triangle_count[MAX_VERTS+1] =
{0 //0 vertices
,0 //1 vertices
,0 //2 vertices
,1 //3 vertices
,4 //4 vertices
,10 //5 vertices
,20 //6 vertices
,35 //7 vertices
};

static int triangles_3[1][3] = {{0,1,2}};
static int triangulations_3[1][1] = {{0}};

static int triangles_4[4][3] =
{{0,1,2}
,{0,2,3}
,{0,1,3}
,{1,2,3}
};
static int triangulations_4[2][2] =
{{0,1}
,{2,3}
};

static int triangles_5[10][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,1,4}
,{1,3,4}
,{1,2,3}
,{2,3,4}
,{0,2,4}
,{0,1,3}
,{1,2,4}
};
static int triangulations_5[5][3] =
{{0,1,2}
,{3,4,5}
,{0,6,7}
,{2,5,8}
,{3,6,9}
};

static int triangles_6[20][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,4,5}
,{0,2,5}
,{2,4,5}
,{2,3,4}
,{0,3,5}
,{3,4,5}
,{0,2,4}
,{2,3,5}
,{1,2,3}
,{0,1,3}
,{0,1,5}
,{1,4,5}
,{1,3,4}
,{0,1,4}
,{1,3,5}
,{1,2,4}
,{1,2,5}
};
static int triangulations_6[14][4] =
{{0,1,2,3}
,{0,4,5,6}
,{0,1,7,8}
,{0,3,6,9}
,{0,4,8,10}
,{2,3,11,12}
,{11,13,14,15}
,{7,8,11,12}
,{3,11,15,16}
,{8,11,13,17}
,{6,13,14,18}
,{3,6,16,18}
,{5,6,13,19}
,{8,10,13,19}
};

static int triangles_7[35][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,4,5}
,{0,5,6}
,{0,3,6}
,{3,5,6}
,{3,4,5}
,{0,4,6}
,{4,5,6}
,{0,3,5}
,{3,4,6}
,{0,2,4}
,{2,3,4}
,{0,2,6}
,{2,5,6}
,{2,4,5}
,{0,2,5}
,{2,4,6}
,{2,3,5}
,{2,3,6}
,{0,1,3}
,{1,2,3}
,{0,1,4}
,{1,3,4}
,{0,1,6}
,{1,5,6}
,{1,4,5}
,{0,1,5}
,{1,4,6}
,{1,3,5}
,{1,3,6}
,{1,2,4}
,{1,2,5}
,{1,2,6}
};
static int triangulations_7[42][5] =
{{0,1,2,3,4}
,{0,1,5,6,7}
,{0,1,2,8,9}
,{0,1,4,7,10}
,{0,1,5,9,11}
,{0,3,4,12,13}
,{0,13,14,15,16}
,{0,8,9,12,13}
,{0,4,13,16,17}
,{0,9,13,14,18}
,{0,7,14,15,19}
,{0,4,7,17,19}
,{0,6,7,14,20}
,{0,9,11,14,20}
,{2,3,4,21,22}
,{5,6,7,21,22}
,{2,8,9,21,22}
,{4,7,10,21,22}
,{5,9,11,21,22}
,{3,4,22,23,24}
,{22,24,25,26,27}
,{8,9,22,23,24}
,{4,22,24,27,28}
,{9,22,24,25,29}
,{7,22,25,26,30}
,{4,7,22,28,30}
,{6,7,22,25,31}
,{9,11,22,25,31}
,{3,4,13,23,32}
,{13,25,26,27,32}
,{8,9,13,23,32}
,{4,13,27,28,32}
,{9,13,25,29,32}
,{13,16,25,26,33}
,{4,13,16,28,33}
,{13,15,16,25,34}
,{9,13,18,25,34}
,{7,19,25,26,33}
,{4,7,19,28,33}
,{7,15,19,25,34}
,{6,7,20,25,34}
,{9,11,20,25,34}
};

/* array [8] of pointer to array [3] of int */
static int (*triangles[MAX_VERTS+1])[3] =
{0
,0
,0
,triangles_3
,triangles_4
,triangles_5
,triangles_6
,triangles_7
};

static void getTriangulation_3(int i, apf::DynamicArray<int>& t)
{
  for (int j=0; j < triangulation_size[3]; ++j)
    t[j] = triangulations_3[i][j];
}
static void getTriangulation_4(int i, apf::DynamicArray<int>& t)
{
  for (int j=0; j < triangulation_size[4]; ++j)
    t[j] = triangulations_4[i][j];
}
static void getTriangulation_5(int i, apf::DynamicArray<int>& t)
{
  for (int j=0; j < triangulation_size[5]; ++j)
    t[j] = triangulations_5[i][j];
}
static void getTriangulation_6(int i, apf::DynamicArray<int>& t)
{
  for (int j=0; j < triangulation_size[6]; ++j)
    t[j] = triangulations_6[i][j];
}
static void getTriangulation_7(int i, apf::DynamicArray<int>& t)
{
  for (int j=0; j < triangulation_size[7]; ++j)
    t[j] = triangulations_7[i][j];
}

static void getTriangulation(int loopSize, int i, apf::DynamicArray<int>& t)
{
  assert(t.getSize() == static_cast<size_t>(triangulation_size[loopSize]));
  typedef void (*GetTriangulationFunction)(int,apf::DynamicArray<int>&);
  static GetTriangulationFunction table[MAX_VERTS+1] =
  {0
  ,0
  ,0
  ,getTriangulation_3
  ,getTriangulation_4
  ,getTriangulation_5
  ,getTriangulation_6
  ,getTriangulation_7
  };
  table[loopSize](i,t);
}

class EdgeSwap2D : public EdgeSwap
{
  public:
    EdgeSwap2D(Adapt* a)
    {
      adapter = a;
      mesh = a->mesh;
      if (mesh->getDimension()==2)
        cavity.init(a);
      oldFaces.setSize(2);
    }
    virtual ~EdgeSwap2D() {}
    void findOldFaces()
    {
      apf::Up faces;
      mesh->getUp(edge,faces);
      /* filter out mesh faces in model region 
         in the case of 2D embedded in 3D */
      int n = 0;
      for (int i = 0; i < faces.n; ++i)
        if (isOnModelFace(mesh,faces.e[i]))
          oldFaces[n++]=faces.e[i];
    }
    void findQuad()
    {
      Entity* ev[2];
      mesh->getDownward(edge,0,ev);
      quad[0] = ev[0];
      quad[2] = ev[1];
      quad[1] = getTriVertOppositeEdge(mesh,oldFaces[0],edge);
      quad[3] = getTriVertOppositeEdge(mesh,oldFaces[1],edge);
    }
    void orient()
    {
      /* for surface meshes with oriented triangle normals,
         this code checks whether we are looking at the mesh "upside-down"
         and flips us right-side-up if so.
         this can be done with just adjacency orders,
         no floating point math. */
      if (!isTriEdgeAligned(mesh,oldFaces[1],edge)) {
        std::swap(oldFaces[0],oldFaces[1]);
        std::swap(quad[1],quad[3]);
      }
    }
    void getOldVerts(Entity* otv[2][3])
    {
      otv[0][0] = quad[0]; otv[0][1] = quad[1]; otv[0][2] = quad[2];
      otv[1][0] = quad[0]; otv[1][1] = quad[2]; otv[1][2] = quad[3];
    }
    void getNewVerts(Entity* ntv[2][3])
    {
      ntv[0][0] = quad[1]; ntv[0][1] = quad[3]; ntv[0][2] = quad[0];
      ntv[1][0] = quad[1]; ntv[1][1] = quad[2]; ntv[1][2] = quad[3];
    }
    bool wouldInvert()
    {
      /* again, for surface meshes, we now check
         whether we have a situation where the new normals
         would oppose the old ones */
      Entity* otv[2][3];
      getOldVerts(otv);
      Entity* ntv[2][3];
      getNewVerts(ntv);
      Vector on[2];
      on[0] = getTriNormal(mesh, otv[0]); on[1] = getTriNormal(mesh, otv[1]);
      Vector nn[2];
      nn[0] = getTriNormal(mesh, ntv[0]); nn[1] = getTriNormal(mesh, ntv[1]);
      if ((on[0] * nn[0] > 0) &&
          (on[0] * nn[1] > 0) &&
          (on[1] * nn[0] > 0) &&
          (on[1] * nn[1] > 0))
        return false;
      return true;
    }
    bool checkTopo()
    {
      Entity* ev[2];
      ev[0] = quad[1]; ev[1] = quad[3];
      /* check whether the "other" edge already exists, which
         means we are dealing with a single tet with two faces
         on the model face. we can't swap the edge between those faces. */
      return ! findElement(mesh,EDGE,ev);
    }
    bool setEdge(Entity* e)
    {
      edge = e;
      findOldFaces();
      findQuad();
      return checkTopo();
    }
    bool didImproveQuality()
    {
      return getWorstQuality(adapter,newFaces,2)
           > getWorstQuality(adapter,oldFaces);
    }
    void destroyOldFaces()
    {
      destroyElement(adapter,oldFaces[0]);
      destroyElement(adapter,oldFaces[1]);
    }
    void cancel()
    {
      destroyElement(adapter,newFaces[0]);
      destroyElement(adapter,newFaces[1]);
    }
    void makeNewFaces()
    {
      Model* c = mesh->toModel(edge);
      Entity* ntv[2][3];
      getNewVerts(ntv);
      newFaces[0] = buildElement(adapter,c,TRI,ntv[0]);
      newFaces[1] = buildElement(adapter,c,TRI,ntv[1]);
    }
    /* this function is only called when swapping
       edges on a surface triangle mesh */
    virtual bool run(Entity* e)
    {
      if (getFlag(adapter,e,DONT_SWAP))
        return false;
      if (isOnModelEdge(mesh,e))
        return false;
      if ( ! setEdge(e))
        return false;
      orient();
      if (wouldInvert())
        return false;
      cavity.beforeBuilding();
      makeNewFaces();
      cavity.afterBuilding();
      cavity.fit(oldFaces);
      if ( ! didImproveQuality())
      {
        cancel();
        return false;
      }
      cavity.transfer(oldFaces);
      destroyOldFaces();
      return true;
    }
    Entity** getOldFaces() { return &(oldFaces[0]); }
  private:
    Adapt* adapter;
    Mesh* mesh;
    Entity* edge;
    Entity* quad[4];
    EntityArray oldFaces;
    Entity* newFaces[2];
    Cavity cavity;
};

/* this class represents a loop of vertices around
   an edge to be swapped. Most of its code is dedicated
   to finding and storing this loop using adjacency searches
   from the edge and one key face attached to it */
class SwapLoop
{
  public:
    void init(Adapt* a)
    {
      mesh = a->mesh;
    }
    void setEdge(Entity* e)
    {
      edge = e;
      mesh->getDownward(e,0,edge_verts);
    }
    /* returns true iff going from vert to the other
       vertex in the tet around the edge curls in the ev[0]-ev[1]
       direction. This establishes an orientation for
       the loop and its triangles with respect to the edge
       so that the generated tets are not left-handed */
    bool isCurlOk(Entity* vert, Entity* tet)
    {
      Entity* tv0[4];
      mesh->getDownward(tet,0,tv0);
      int n = findIn(tv0,4,vert);
      Entity* tv[4];
      rotateTet(tv0,n*3,tv);
      assert(tv[0]==vert);
      int a = findIn(tv+1,3,edge_verts[0]);
      int b = findIn(tv+1,3,edge_verts[1]);
      if (b == ((a+1)%3))
        return true;
      assert(a == ((b+1)%3));
      return false;
    }
    Entity* findInitialTet(Entity* vert, Entity* face)
    {
      apf::Up tets;
      mesh->getUp(face,tets);
      for (int i=0; i < tets.n; ++i)
      {
        Entity* tet = tets.e[i];
        if (isCurlOk(vert,tet))
          return tet;
      }
      return 0;
    }
    void findFromFace(Entity* startFace)
    {
      Entity* face = startFace;
      Entity* tet = 0;
      model = 0;
      size = 0;
      while (1)
      {
        Entity* v = getTriVertOppositeEdge(mesh,face,edge);
/* there may be more than MAX_VERTS: just walk over them.
   overflow will be checked by users of SwapLoop
   (i.e. SwapCavity::findGoodTriangulation) */
        if (size < MAX_VERTS)
          verts[size] = v;
        ++size;
        if ( ! tet)
          tet = findInitialTet(v,face);
        else
          tet = getOtherTet(tet,face);
        if ( ! tet)
          break;
/* all the tets here should be in the same model region,
   record it for classifying new ones */
        if ( ! model)
          model = mesh->toModel(tet);
        else if (model != mesh->toModel(tet))
          break; /* this is analogous to (!tet) when
                    we cross a non-manifold face */
        face = getOtherFace(face,tet);
        assert(face);
        if (face == startFace)
          break;
      }
    }
    Entity* getOtherTet(Entity* oldTet, Entity* face)
    {
      apf::Up tets;
      mesh->getUp(face,tets);
      for (int i=0; i < tets.n; ++i)
      {
        Entity* newTet = tets.e[i];
        if (newTet != oldTet)
          return newTet;
      }
      return 0;
    }
    Entity* getOtherFace(Entity* oldFace, Entity* tet)
    {
      Entity* tf[4];
      mesh->getDownward(tet,2,tf);
      for (int i=0; i < 4; ++i)
      {
        if (tf[i] == oldFace) continue;
        Entity* fe[3];
        mesh->getDownward(tf[i],1,fe);
        if (findIn(fe,3,edge) != -1)
          return tf[i];
      }
      return 0;
    }
    int getSize() {return size;}
    Entity* getVert(int i) {return verts[i];}
    Entity* getEdgeVert(int i) {return edge_verts[i];}
    Model* getModel() {return model;}
  private:
    Mesh* mesh;
    Entity* edge;
    Entity* edge_verts[2];
    int size;
    Entity* verts[7];
    Model* model;
};

/* this class represents the full cavity around the
   loop of vertices. It is responsible for trying
   new triangulations as efficiently as possible
   and creating the first one that works */
class SwapCavity
{
  public:
    void init(Adapt* a)
    {
      adapter = a;
      shape = a->shape;
      mesh = a->mesh;
      loop.init(a);
      tempTet.init(a);
    }
    bool setFromEdge(Entity* edge)
    {
      return setFromEdgeAndFace(edge,mesh->getUpward(edge,0));
    }
/* returns true iff there are tets in the cavity */
    bool setFromEdgeAndFace(Entity* edge, Entity* face)
    {
      loop.setEdge(edge);
      loop.findFromFace(face);
      return loop.getSize() > 1;
    }
    bool isTetOk(Entity* tet)
    {
      double quality = shape->getQuality(tet);
      return (quality > qualityToBeat);
    }
    void getTriVerts(int tri, Entity** v)
    {
      int* ti = triangles[loop.getSize()][tri];
      for (int j=0; j < 3; ++j)
        v[j] = loop.getVert(ti[j]);
    }
    Entity* buildTopTet(Entity* triv[3])
    {
      Entity* tv[4] = {triv[0],triv[1],triv[2],loop.getEdgeVert(1)};
      return buildElement(adapter,loop.getModel(),TET,tv);
    }
    Entity* buildBottomTet(Entity* triv[3])
    {
      Entity* tv[4] = {triv[0],triv[2],triv[1],loop.getEdgeVert(0)};
      return buildElement(adapter,loop.getModel(),TET,tv);
    }
    bool checkTet(bool isTop, Entity* tv[3])
    {
      Entity* tet;
      tempTet.beforeTrying();
      if (isTop)
        tet = buildTopTet(tv);
      else
        tet = buildBottomTet(tv);
      tempTet.afterTrying();
      tempTet.fit(*oldTets);
      bool ok = isTetOk(tet);
      destroyElement(adapter,tet);
      return ok;
    }
    bool checkTriangle(int i)
    {
      Entity* tv[3];
      getTriVerts(i,tv);
      if (findElement(mesh,TRI,tv))
        return false;
      return checkTet(true,tv) && checkTet(false,tv);
    }
    void acceptTriangle(int unique_i, int local_i)
    {
      Entity* tv[3];
      getTriVerts(unique_i,tv);
      Entity* tet = buildTopTet(tv);
      this->tets[2*local_i] = tet;
      tet = buildBottomTet(tv);
      this->tets[2*local_i+1] = tet;
    }
    void cancelTriangle(int i)
    {
      Entity* tv[3];
      getTriVerts(i,tv);
/* the buildElement function just retrieves existing tets */
      destroyElement(adapter,buildTopTet(tv));
      destroyElement(adapter,buildBottomTet(tv));
    }
    bool isTriangleOk(int i)
    {
      if ( ! triangleChecked[i])
      { /* cache the expensive check */
        triangleOk[i] = checkTriangle(i);
        triangleChecked[i] = true;
      }
      return triangleOk[i];
    }
    bool tryTriangulation(int i)
    {
      getTriangulation(loop.getSize(),i,triangulation);
      for (size_t i=0; i < triangulation.getSize(); ++i)
        if ( ! isTriangleOk(triangulation[i]))
          return false;
      return true;
    }
    bool findGoodTriangulation(double q, Upward& ot)
    {
      if (loop.getSize() < 3)
        return false;
      if (loop.getSize() > MAX_VERTS)
        return false;
      qualityToBeat = q;
      oldTets = &ot;
      int unique_count = unique_triangle_count[loop.getSize()];
      triangleOk.setSize(unique_count);
      triangleChecked.setSize(unique_count);
      for (int i=0; i < unique_count; ++i)
        triangleChecked[i] = false;
      triangulation.setSize(triangulation_size[loop.getSize()]);
      for (int i=0; i < triangulation_count[loop.getSize()]; ++i)
        if (tryTriangulation(i))
          return true;
      return false;
    }
    void acceptTriangulation()
    {
      tets.setSize(2*(triangulation.getSize()));
      for (size_t i=0; i < triangulation.getSize(); ++i)
        acceptTriangle(triangulation[i],i);
    }
    EntityArray& getNewTets() {return tets;}
  private:
    Adapt* adapter;
    ShapeHandler* shape;
    Mesh* mesh;
    SwapLoop loop;
    apf::DynamicArray<int> triangleOk;
    apf::DynamicArray<int> triangleChecked;
    apf::DynamicArray<int> triangulation;
    EntityArray tets;
    double qualityToBeat;
    Cavity tempTet;
    Upward* oldTets;
};

class EdgeSwap3D : public EdgeSwap
{
  public:
    EdgeSwap3D(Adapt* a):
      adapter(a),
      mesh(a->mesh),
      swap2d(a)
    {
      halves[0].init(a);
      halves[1].init(a);
      cavity.init(a);
    }
    virtual ~EdgeSwap3D() {}
    void destroyOldTets()
    {
      for (size_t i=0; i < oldTets.getSize(); ++i)
        destroyElement(adapter,oldTets[i]);
    }
    virtual bool run(Entity* e)
    {
      if (getFlag(adapter,e,DONT_SWAP))
        return false;
      if (isOnModelEdge(mesh,e))
        return false;
      edge = e;
      mesh->getAdjacent(edge,3,oldTets);
      double oldQuality = getWorstQuality(adapter,oldTets);
      if (isOnModelFace(mesh,edge))
      {
/* note that only one model face may cross the cavity; otherwise
   the edge is on a model edge and swapping is not attempted */
        if ( ! swap2d.setEdge(edge))
          return false;
        Entity** oldFaces = swap2d.getOldFaces();
        for (int i=0; i < 2; ++i)
          cavityExists[i] = halves[i].setFromEdgeAndFace(
              edge,oldFaces[i]);
/* there must be at least one cavity */
        if (( ! cavityExists[0])&&( ! cavityExists[1]))
          return false;
        for (int i=0; i < 2; ++i)
          if (cavityExists[i])
            if ( ! halves[i].findGoodTriangulation(oldQuality,oldTets))
              return false;
        cavity.beforeBuilding();
/* if we make the mesh faces here with correct classification, they
   are not accidentally destroyed during cavity evaluation */ 
        swap2d.makeNewFaces();
/* now fill in the new tets */
        for (int i=0; i < 2; ++i)
          if (cavityExists[i])
            halves[i].acceptTriangulation();
      }
      else
      {
        cavityExists[0] = halves[0].setFromEdge(edge);
        assert(cavityExists[0]);
        if ( ! halves[0].findGoodTriangulation(oldQuality,oldTets))
          return false;
        cavity.beforeBuilding();
        halves[0].acceptTriangulation();
      }
      cavity.afterBuilding();
      cavity.fit(oldTets);
      cavity.transfer(oldTets);
      destroyOldTets();
      return true;
    }
  private:
    Adapt* adapter;
    Mesh* mesh;
    Entity* edge;
    EdgeSwap2D swap2d;
    Upward oldTets;
    SwapCavity halves[2];
    bool cavityExists[2];
    Cavity cavity;
};

EdgeSwap* makeEdgeSwap(Adapt* a)
{
  int dimension = a->mesh->getDimension();
  if (dimension==2)
    return new EdgeSwap2D(a);
  if (dimension==3)
    return new EdgeSwap3D(a);
  return 0;
}

}
