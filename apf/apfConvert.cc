#include <PCU.h>
#include "apf.h"
#include "apfConvert.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include <map>

namespace apf {

class Converter
{
  public:
    Converter(Mesh *a, Mesh2 *b)
    {
      inMesh = a;
      outMesh = b;
    }
    void run()
    {
      createVertices();
      createEntities();
      for (int i = 0; i <= inMesh->getDimension(); ++i)
        createRemotes(i);
      if (inMesh->hasMatching())
        for (int i = 0; i <= inMesh->getDimension(); ++i)
          createMatches(i);
      convertQuadratic();
      convertFields();
      outMesh->acceptChanges();
    }
    ModelEntity* getNewModelFromOld(ModelEntity* oldC)
    {
      int type = inMesh->getModelType(oldC);
      int tag = inMesh->getModelTag(oldC);
      return outMesh->findModelEntity(type,tag);
    }
    void createVertices()
    {
      MeshIterator* it = inMesh->begin(0);
      MeshEntity *oldV;
      while ((oldV = inMesh->iterate(it)))
      {
        ModelEntity *oldC = inMesh->toModel(oldV);
        ModelEntity *newC = getNewModelFromOld(oldC);
        Vector3 xyz;
        inMesh->getPoint(oldV, 0, xyz);
        Vector3 param(0,0,0);
        inMesh->getParam(oldV,param);
        MeshEntity* newV = outMesh->createVertex(newC, xyz, param);
        newFromOld[oldV] = newV;
      }
      inMesh->end(it);
      assert(outMesh->count(0) == inMesh->count(0));
    }
    void createEntities()
    { 
      for (int i = 1; i < (inMesh->getDimension())+1; ++i)
        createDimension(i);
    }
    void createDimension(int dim)
    { 
      MeshIterator* it = inMesh->begin(dim);
      MeshEntity *oldE;
      while ((oldE = inMesh->iterate(it)))
      {
        int type = inMesh->getType(oldE);
        ModelEntity *oldC = inMesh->toModel(oldE);
        ModelEntity *newC = getNewModelFromOld(oldC);
        Downward down;
        int ne = inMesh->getDownward(oldE, dim-1, down);
        Downward new_down;
        for(int i=0; i<ne; ++i)
        {
          new_down[i]=newFromOld[down[i]];
        }
        MeshEntity *newE = outMesh->createEntity(type, newC, new_down); 
        newFromOld[oldE] = newE;
      }
      inMesh->end(it);
      assert(outMesh->count(dim) == inMesh->count(dim));
    }
    void createRemotes(int dim)
    {
      /*    O-------------|---|----->O
            ^old left     |   |      ^old right
                          |MPI|
            O<------------|---|------O
            ^new left     |   |      ^new right */
      /* creating links backwards is OK because
         they go both ways; as long as every link
         gets a backward copy they will all get copies */
      PCU_Comm_Begin();
      MeshIterator* it = inMesh->begin(dim);
      MeshEntity *oldLeft;
      while ((oldLeft = inMesh->iterate(it)))
      {
        MeshEntity *newLeft = newFromOld[oldLeft];
        Copies remotes;
        inMesh->getRemotes(oldLeft, remotes);
        APF_ITERATE(Copies,remotes,it) 
        {
          int rightPart = it->first;
          MeshEntity* oldRight = it->second;
          PCU_COMM_PACK(rightPart,oldRight);
          PCU_COMM_PACK(rightPart,newLeft);
        }
      }
      inMesh->end(it);
      PCU_Comm_Send();
      /* accumulation of residence sets... we could also constantly reclassify... */
      std::map<MeshEntity*, std::vector<int> > map_ent;
      while (PCU_Comm_Listen())
      {
        int leftPart = PCU_Comm_Sender();
        while ( ! PCU_Comm_Unpacked())
        {  
          MeshEntity* oldRight;
          PCU_COMM_UNPACK(oldRight);
          MeshEntity* newLeft;
          PCU_COMM_UNPACK(newLeft);
          MeshEntity *newRight = newFromOld[oldRight];
          outMesh->addRemote(newRight, leftPart, newLeft);
          map_ent[newRight].push_back(leftPart);
        }
      }
      /* assign accumulated residence sets */
      std::map <MeshEntity*, std::vector<int> >::iterator mit;
      for(mit=map_ent.begin(); mit!=map_ent.end(); ++mit)
      {
        MeshEntity* newE;
        newE = mit->first;
        std::vector<int>& vecEnt = mit->second;
        Copies remotes;
        Parts parts;
        parts.insert(vecEnt.begin(), vecEnt.end());
        parts.insert(PCU_Comm_Self());
        outMesh->setResidence(newE, parts);
      }
    }
    void convertField(Field* in, Field* out)
    {
      FieldShape* s = getShape(in);
      NewArray<double> data(countComponents(in));
      for (int d = 0; d <= 3; ++d)
      {
        if (s->hasNodesIn(d))
        {
          MeshIterator* it = inMesh->begin(d);
          MeshEntity* e;
          while ((e = inMesh->iterate(it)))
          {
            int n = s->countNodesOn(inMesh->getType(e));
            for (int i = 0; i < n; ++i)
            {
              getComponents(in,e,i,&(data[0]));
              setComponents(out,newFromOld[e],i,&(data[0]));
            }
          }
          inMesh->end(it);
        }
      }
    }
    void convertFields()
    {
      for (int i = 0; i < inMesh->countFields(); ++i) {
        Field* in = inMesh->getField(i);
        Field* out = cloneField(in, outMesh);
        convertField(in, out);
      }
    }
    void convertQuadratic()
    {
      if (inMesh->getShape() != getLagrange(2))
        return;
      if ( ! PCU_Comm_Self())
        fprintf(stderr,"transferring quadratic mesh\n");
      changeMeshShape(outMesh,getLagrange(2),/*project=*/false);
      convertField(inMesh->getCoordinateField(),outMesh->getCoordinateField());
    }
    void createMatches(int dim)
    {
      /* see createRemotes for the algorithm comments */
      PCU_Comm_Begin();
      MeshIterator* it = inMesh->begin(dim);
      MeshEntity *oldLeft;
      while ((oldLeft = inMesh->iterate(it)))
      {
        MeshEntity *newLeft = newFromOld[oldLeft];
        Matches matches;
        inMesh->getMatches(oldLeft, matches);
        for (size_t i = 0; i < matches.getSize(); ++i)
        {
          int rightPart = matches[i].peer;
          MeshEntity* oldRight = matches[i].entity;
          PCU_COMM_PACK(rightPart,oldRight);
          PCU_COMM_PACK(rightPart,newLeft);
        }
      }
      inMesh->end(it);
      PCU_Comm_Send();
      while (PCU_Comm_Listen())
      {
        int leftPart = PCU_Comm_Sender();
        while ( ! PCU_Comm_Unpacked())
        {  
          MeshEntity* oldRight;
          PCU_COMM_UNPACK(oldRight);
          MeshEntity* newLeft;
          PCU_COMM_UNPACK(newLeft);
          MeshEntity *newRight = newFromOld[oldRight];
          outMesh->addMatch(newRight, leftPart, newLeft);
        }
      }
    }
  private:
    Mesh *inMesh;
    Mesh2 *outMesh;
    std::map<MeshEntity*,MeshEntity*> newFromOld;
};

void convert(Mesh *in, Mesh2 *out)
{
  Converter c(in,out);
  c.run();
}

}
