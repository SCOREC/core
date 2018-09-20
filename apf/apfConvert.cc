#include <PCU.h>
#include "apf.h"
#include "apfConvert.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include "apfNumbering.h"
#include <map>
#include <pcu_util.h>
#include <lionPrint.h>
#include <iostream>
#include <cstdlib>

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
      convertNumberings();
      convertGlobalNumberings();
      // this must be called after anything that might create tags e.g. fields
      // or numberings to avoid problems with tag duplication
      convertTags();
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
      PCU_ALWAYS_ASSERT(outMesh->count(0) == inMesh->count(0));
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
      PCU_ALWAYS_ASSERT(outMesh->count(dim) == inMesh->count(dim));
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
    void convertNumbering(Numbering* in, Numbering* out)
    {
      FieldShape* s = getShape(in);
      int nc = countComponents(in);
      for (int d = 0; d <= 3; ++d)
      {
        if (s->hasNodesIn(d))
        {
          MeshIterator* it = inMesh->begin(d);
          MeshEntity* e;
          while ((e = inMesh->iterate(it)))
          {
            int nn = s->countNodesOn(inMesh->getType(e));
            for (int i = 0; i < nn; ++i)
            {
              for (int j = 0; j < nc; ++j) {
                number(out, newFromOld[e], i, j,
                    getNumber(in, e, i, j));
              }
            }
          }
          inMesh->end(it);
        }
      }
    }
    void convertGlobalNumbering(
        GlobalNumbering* in, GlobalNumbering* out)
    {
      FieldShape* s = getShape(in);
      PCU_DEBUG_ASSERT(countComponents(in) == 1);
      for (int d = 0; d <= 3; ++d) {
        if (s->hasNodesIn(d)) {
          MeshIterator* it = inMesh->begin(d);
          MeshEntity* e;
          while ((e = inMesh->iterate(it))) {
            int nn = s->countNodesOn(inMesh->getType(e));
            for (int i = 0; i < nn; ++i)
              number(out, newFromOld[e], i,
                  getNumber(in, e, i, 0));
          }
          inMesh->end(it);
        }
      }
    }
    void convertTag(Mesh* inMesh, MeshTag* in, Mesh* outMesh, MeshTag* out)
    {
      for (int d = 0; d <= 3; ++d) {
        int tagType = inMesh->getTagType(in);
        int tagSize = inMesh->getTagSize(in);
        PCU_DEBUG_ASSERT(tagType == outMesh->getTagType(out));
        PCU_DEBUG_ASSERT(tagSize == outMesh->getTagSize(out));
        MeshIterator* it = inMesh->begin(d);
        MeshEntity* e;
        while ((e = inMesh->iterate(it))) {
          if(inMesh->hasTag(e, in)) {
            // these initializations cannot go into the cases due to compiler
            // warnings on gcc 7.3.0
            double* dblData;
            int* intData;
            long* lngData; 
            switch (tagType) {
              case apf::Mesh::DOUBLE:
                dblData = new double[tagSize];
                inMesh->getDoubleTag(e, in, dblData);
                outMesh->setDoubleTag(newFromOld[e], out, dblData);
                break;
              case apf::Mesh::INT:
                intData = new int[tagSize];
                inMesh->getIntTag(e, in, intData);
                outMesh->setIntTag(newFromOld[e], out, intData);
                break;
              case apf::Mesh::LONG:
                lngData = new long[tagSize];
                inMesh->getLongTag(e, in, lngData);
                outMesh->setLongTag(newFromOld[e], out, lngData);
                break;
              default:
                lion_eprint(1,"Tried to convert unknown tag type\n");
                abort();
                break;
            }
        }
        }
        inMesh->end(it);
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
    void convertNumberings()
    {
      for (int i = 0; i < inMesh->countNumberings(); ++i) {
        Numbering* in = inMesh->getNumbering(i);
        Numbering* out;
        if (getField(in)) {
          // here we assume that the fields have already been copied into the
          // mesh
          Field* outField = outMesh->findField(getName(getField(in)));
          PCU_DEBUG_ASSERT(outField);
          out = createNumbering(outField);
        }
        else {
          out = createNumbering(outMesh, getName(in), getShape(in),
                                countComponents(in));
        }
        convertNumbering(in, out);
      }
    }
    void convertGlobalNumberings()
    {
      for (int i = 0; i < inMesh->countGlobalNumberings(); ++i) {
        GlobalNumbering* in = inMesh->getGlobalNumbering(i);
        GlobalNumbering* out;
        if (getField(in)) {
          // here we assume that the fields have already been copied into the
          // mesh
          Field* outField = outMesh->findField(getName(getField(in)));
          PCU_DEBUG_ASSERT(outField);
          out = createGlobalNumbering(outField);
        }
        else {
          out = createGlobalNumbering(outMesh, getName(in), getShape(in),
                                      countComponents(in));
        }
        convertGlobalNumbering(in, out);
      }
    }
    void convertTags()
    {
      DynamicArray<MeshTag*> tags;
      inMesh->getTags(tags);
      for (std::size_t i = 0; i < tags.getSize(); ++i) {
        apf::MeshTag* in = tags[i];
        PCU_DEBUG_ASSERT(in);
        // create a new tag on the outMesh
        int tagType = inMesh->getTagType(in);
        int tagSize = inMesh->getTagSize(in);
        const char* tagName = inMesh->getTagName(in);
        PCU_DEBUG_ASSERT(tagName);
        // need to make sure that the tag wasn't already created by a field or
        // numbering
        if (!outMesh->findTag(tagName)) {
          apf::MeshTag* out = NULL;
          switch (tagType) {
            case apf::Mesh::DOUBLE:
              out = outMesh->createDoubleTag(tagName, tagSize);
              break;
            case apf::Mesh::INT:
              out = outMesh->createIntTag(tagName, tagSize);
              break;
            case apf::Mesh::LONG:
              out = outMesh->createLongTag(tagName, tagSize);
              break;
            default:
              lion_eprint(1,"Tried to convert unknown tag type\n");
              abort();
          }
          PCU_DEBUG_ASSERT(out);
          // copy the tag on the inMesh to the outMesh
          convertTag(inMesh, in, outMesh, out);
        }
      }
    }
    void convertQuadratic()
    {
      if (inMesh->getShape() != getLagrange(2) && inMesh->getShape() != getSerendipity())
        return;
      if ( ! PCU_Comm_Self())
        lion_eprint(1,"transferring quadratic mesh\n");
      changeMeshShape(outMesh,inMesh->getShape(),/*project=*/false);
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
