/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_ADAPT_H
#define MA_ADAPT_H

#include "maInput.h"

namespace ma {

enum {
  SPLIT         = (1<< 0),
  DONT_SPLIT    = (1<< 1),
  COLLAPSE      = (1<< 2),
  DONT_COLLAPSE = (1<< 3),
  CHECKED       = (1<< 4),
  BAD_QUALITY   = (1<< 5),
  OK_QUALITY    = (1<< 6),
  SNAP          = (1<< 7),
  DONT_SNAP     = (1<< 8),
  DONT_SWAP     = (1<< 9),
  LAYER         = (1<<10),
  LAYER_BASE    = (1<<11),
  LAYER_TOP     = (1<<12),
  DIAGONAL_1    = (1<<13),
  DIAGONAL_2    = (1<<14),
  LAYER_UNSNAP  = (1<<15)
};

class DeleteCallback;
class SolutionTransfer;
class Refine;
class ShapeHandler;

class Adapt
{
  public:
    Adapt(Input* in);
    ~Adapt();
    Input* input;
    Mesh* mesh;
    Tag* flagsTag;
    DeleteCallback* deleteCallback;
    apf::BuildCallback* buildCallback;
    SizeField* sizeField;
    SolutionTransfer* solutionTransfer;
    Refine* refine;
    ShapeHandler* shape;
    int coarsensLeft;
    int refinesLeft;
    bool hasLayer;
};

void setTolerance(Adapt* a, double t);
double getTolerance(Adapt* a, int dimension);

void setupFlags(Adapt* a);
void clearFlags(Adapt* a);
int getFlags(Adapt* a, Entity* e);
void setFlags(Adapt* a, Entity* e, int flags);
bool getFlag(Adapt* a, Entity* e, int flag);
void setFlag(Adapt* a, Entity* e, int flag);
void clearFlag(Adapt* a, Entity* e, int flag);

void clearFlagFromDimension(Adapt* a, int flag, int dimension);

void destroyElement(Adapt* a, Entity* e);

class DeleteCallback
{
  public:
    DeleteCallback(Adapt* a);
    ~DeleteCallback();
    virtual void call(Entity* e) = 0;
    Adapt* adapt;
};

bool checkFlagConsistency(Adapt* a, int dimension, int flag);

/* get the distance between two vertices in metric space */
double getDistance(Adapt* a, Entity* v0, Entity* v1);
/* given an array of entity pairs,
   return the index of the closest pair in metric space */
int getClosestPair(Adapt* a, Entity* (*pairs)[2], int n);

struct Predicate
{
  virtual bool operator()(Entity* e) = 0;
};

long markEntities(
    Adapt* a,
    int dimension,
    Predicate& predicate,
    int trueFlag,
    int falseFlag);

class NewEntities : public apf::BuildCallback
{
  public:
    void reset();
    void addEntity(Entity* e);
    virtual void call(Entity* e);
    void retrieve(EntityArray& a);
  private:
    std::vector<Entity*> entities;
};

class Cavity
{
  public:
    Cavity();
    void init(Adapt* a);
    void beforeBuilding();
    void afterBuilding();
    void beforeTrying();
    void afterTrying();
    void transfer(EntityArray& oldElements);
    void fit(EntityArray& oldElements);
    bool shouldTransfer;
    bool shouldFit;
  private:
    Adapt* adapter;
    SolutionTransfer* solutionTransfer;
    ShapeHandler* shape;
    NewEntities newEntities;
};

Entity* buildVertex(
    Adapt* a,
    Model* c,
    Vector const& point,
    Vector const& param);
Entity* buildElement(
    Adapt* a,
    Model* c,
    int type,
    Entity** verts);
Entity* rebuildElement(
    Adapt* a,
    Entity* original,
    Entity* oldVert,
    Entity* newVert);

void setBuildCallback(Adapt* a, apf::BuildCallback* cb);
void clearBuildCallback(Adapt* a);

void print(const char* format, ...) __attribute__((format(printf,1,2)));

void setFlagOnClosure(Adapt* a, Entity* e, int flag);
void syncFlag(Adapt* a, int dimension, int flag);

struct HasTag : public Predicate
{
  HasTag(Mesh* m, Tag* t);
  bool operator()(Entity* e);
  Mesh* mesh;
  Tag* tag;
};

struct HasFlag : public Predicate
{
  HasFlag(Adapt* a, int f);
  bool operator()(Entity* e);
  Adapt* adapter;
  int flag;
};

}

#endif
