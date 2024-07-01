#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "pcu_io.h"
#include "pcu_byteorder.h"

#include <cstdio>
#include <cstring>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>

#include <gmi.h>
#include <algorithm>
#include <iostream>

/*
read files in the AFLR3 format from Dave Marcum at Mississippi State
-little-endian if file has suffix '.lb8.ugrid'
-big-endian if file has suffix '.b8.ugrid'
-integers are 4B
-floats are 8B
*/

namespace {
  typedef std::set<int> SetInt;
  int faceTypeIdx(int type) {
    switch(type) {
      case apf::Mesh::TRIANGLE: return 0;
      case apf::Mesh::QUAD: return 1;
      default: {
                 PCU_ALWAYS_ASSERT(false);
                 return -1;
               }
    }
  }

  struct Reader {
    apf::Mesh2* mesh;
    FILE* file;
    std::map<long, apf::MeshEntity*> nodeMap;
    unsigned* faceVerts[2];
    unsigned* faceTags[2];
    unsigned* edgeVerts;
    unsigned* edgeTags;
    bool swapBytes;
  };

  struct header {
    unsigned nvtx, ntri, nquad, ntet, npyr, nprz, nhex, nbdry;
    void print() {
      lion_eprint(1,
          "nvtx %u ntri %u nquad %u ntet %u npyr %u nprz %u nhex %u\n",
          nvtx,ntri,nquad,ntet,npyr,nprz,nhex);
    }
  };

  void initReader(Reader* r, apf::Mesh2* m, const char* filename) {
    r->mesh = m;
    r->file = fopen(filename, "rb");
    if (!r->file) {
      lion_eprint(1,"ERROR couldn't open ugrid file \"%s\"\n",filename);
      abort();
    }
    unsigned endian = -1;
    if ( strstr(filename, ".b8.ugrid") ) {
      endian = PCU_BIG_ENDIAN;
    } else if ( strstr(filename, ".lb8.ugrid") ) {
      endian = PCU_LITTLE_ENDIAN;
    } else {
      lion_eprint(1,
          "ERROR file extension of \"%s\" is not supported\n", filename);
      exit(EXIT_FAILURE);
    }
    r->swapBytes = ( endian != PCU_HOST_ORDER );
  }

  void readUnsigneds(FILE* f, unsigned* v, size_t cnt, bool swap) {
    PCU_ALWAYS_ASSERT(sizeof(unsigned)==4);
    size_t read = fread(v, sizeof(unsigned), cnt, f);
    PCU_ALWAYS_ASSERT(read == cnt);
    if ( swap )
      pcu_swap_unsigneds(v,cnt);
  }

  void readDoubles(FILE* f, double* vals, size_t cnt, bool swap) {
    PCU_ALWAYS_ASSERT(sizeof(double)==8);
    size_t read = fread(vals, sizeof(double), cnt, f);
    PCU_ALWAYS_ASSERT(read == cnt);
    if ( swap )
      pcu_swap_doubles(vals,cnt);
  }

  void readHeader(Reader* r, header* h) {
    const unsigned biggest = 100*1000*1000;
    unsigned headerVals[7];
    readUnsigneds(r->file, headerVals, 7, r->swapBytes);
    for(unsigned i=0; i<7; i++) {
      PCU_ALWAYS_ASSERT(headerVals[i] < biggest);
    }
    h->nvtx = headerVals[0];
    h->ntri = headerVals[1];
    h->nquad = headerVals[2];
    h->ntet = headerVals[3];
    h->npyr = headerVals[4];
    h->nprz = headerVals[5];
    h->nhex = headerVals[6];
    h->nbdry = 0;
  }

  void readNumBdryElms(Reader* r, header* h) {
    const unsigned biggest = 100*1000*1000;
    unsigned nbdry_elem;
    readUnsigneds(r->file, &nbdry_elem, 1, r->swapBytes);
    PCU_ALWAYS_ASSERT(nbdry_elem < biggest);
    h->nbdry = nbdry_elem;
  }

  apf::MeshEntity* makeVtx(Reader* r,
      apf::Vector3& pt, apf::ModelEntity* g) {
    apf::MeshEntity* v = r->mesh->createVert(g);
    r->mesh->setPoint(v, 0, pt);
    return v;
  }

  int ftnToC(long id) {
    return --id;
  }

  apf::MeshEntity* lookupVert(Reader* r, long ftnNodeId) {
    const long cNodeId = ftnToC(ftnNodeId);
    PCU_ALWAYS_ASSERT(r->nodeMap.count(cNodeId));
    return r->nodeMap[cNodeId];
  }

  void readNodes(Reader* r, header* h) {
    const unsigned dim = 3;
    size_t cnt = h->nvtx*dim;
    double* xyz = (double*) calloc(cnt,sizeof(double));
    readDoubles(r->file, xyz, cnt, r->swapBytes);
    for(long id=0; id<h->nvtx; id++) {
      apf::Vector3 p;
      for(unsigned j=0; j<dim; j++)
        p[j] = xyz[id*dim+j];
      r->nodeMap[id] = makeVtx(r,p,0);
    }
    free(xyz);
  }

  void setNodeIds(Reader* r, header* h) {
    apf::Mesh* m = r->mesh;
    apf::MeshTag* t = m->createIntTag("ugrid-vtx-ids",1);
    for(long lid=0; lid<h->nvtx; lid++) {
      int iid = lid;
      m->setIntTag(r->nodeMap[lid],t,&iid);
    }
  }

  void readFaces(Reader* r, unsigned nfaces, int apfType) {
    const unsigned nverts = apf::Mesh::adjacentCount[apfType][0];
    size_t cnt = nfaces * nverts;
    unsigned* vtx = (unsigned*) calloc(cnt,sizeof(unsigned));
    r->faceVerts[faceTypeIdx(apfType)] = vtx;
    readUnsigneds(r->file, vtx, cnt, r->swapBytes);
  }

  void readTags(Reader* r, unsigned nfaces, int apfType) {
    unsigned* tags = (unsigned*) calloc(nfaces,sizeof(unsigned));
    r->faceTags[faceTypeIdx(apfType)] = tags;
    readUnsigneds(r->file, tags, nfaces, r->swapBytes);
  }

  void readFacesAndTags(Reader* r, header* h) {
    //read the face vertices and store them
    readFaces(r,h->ntri,apf::Mesh::TRIANGLE);
    readFaces(r,h->nquad,apf::Mesh::QUAD);
    //read the tags and store them
    readTags(r,h->ntri,apf::Mesh::TRIANGLE);
    readTags(r,h->nquad,apf::Mesh::QUAD);
  }

  void readBdryElmsAndTags(Reader* r, unsigned nbdry, int apfType) {
    const unsigned nverts = apf::Mesh::adjacentCount[apfType][0];
    size_t cnt = nbdry * (nverts + 1);
    unsigned* data = (unsigned*) calloc(cnt,sizeof(unsigned));
    readUnsigneds(r->file, data, cnt, r->swapBytes);

    unsigned* vtx = (unsigned*) calloc(nbdry*nverts,sizeof(unsigned));    
    r->edgeVerts = vtx;

    unsigned* tags = (unsigned*) calloc(nbdry,sizeof(unsigned));
    r->edgeTags = tags;

    for (size_t i = 0; i < cnt; ++i) {
      if ((i+1) % (nverts+1) != 0) {
        vtx[i - i / (nverts+1)] = data[i];
      }
      else {
        tags[i / 3] = data[i];
      }
    }
    free(data);
  }

  void readBdryElmsAndTags(Reader* r, header* h) {
    readBdryElmsAndTags(r, h->nbdry, apf::Mesh::EDGE);    
  }

  void checkFilePos(Reader* r, header* h) {
    // seven headers, vtx coords, face vertex ids, and face tags
    long expected = h->nvtx*3*sizeof(double) +
      sizeof(unsigned)*(7 + h->ntri*3 + h->nquad*4 + (h->ntri+h->nquad));
    long pos = ftell(r->file);
    PCU_ALWAYS_ASSERT( pos == expected );
  }

  void setFaceTags(Reader* r, apf::MeshTag* t, unsigned nfaces, int apfType) {
    const unsigned nverts = apf::Mesh::adjacentCount[apfType][0];
    unsigned* vtx = r->faceVerts[faceTypeIdx(apfType)];
    unsigned* tags = r->faceTags[faceTypeIdx(apfType)];
    for(unsigned id=0; id<nfaces; id++) {
      apf::Downward verts;
      for(unsigned j=0; j<nverts; j++) {
        verts[j] = lookupVert(r, vtx[id*nverts+j]);
      }
      apf::MeshEntity* f =
        apf::findElement(r->mesh, apfType, verts);
      PCU_ALWAYS_ASSERT(f);
      int val = tags[id];
      r->mesh->setIntTag(f, t, &val);
    }
    free(vtx);
    free(tags);
  }

  void setFaceTags(Reader* r, header* h) {
    apf::MeshTag* t = r->mesh->createIntTag("ugrid-face-tag", 1);
    setFaceTags(r,t,h->ntri,apf::Mesh::TRIANGLE);
    setFaceTags(r,t,h->nquad,apf::Mesh::QUAD);
  }

  void setBoundaryTags(Reader* r, apf::MeshTag* t, unsigned nbdry, int apfType)
  {
    const unsigned nverts = apf::Mesh::adjacentCount[apfType][0];
    unsigned* vtx = r->edgeVerts;
    unsigned* tags = r->edgeTags;
    for(unsigned id=0; id<nbdry; id++) {
      apf::Downward verts;
      for(unsigned j=0; j<nverts; j++) {
        verts[j] = lookupVert(r, vtx[id*nverts+j]);
      }
      apf::MeshEntity* f =
        apf::findElement(r->mesh, apfType, verts);
      PCU_ALWAYS_ASSERT(f);
      int val = tags[id];
      r->mesh->setIntTag(f, t, &val);
    }
  }

  void setBoundaryTags(Reader* r, header* h) {
    apf::MeshTag* t = r->mesh->createIntTag("ugrid-boundary-tag", 1);
    setBoundaryTags(r,t,h->nbdry,apf::Mesh::EDGE);
  }

  void classifyBoundaryElms(Reader* r, unsigned nbdry, int apfType)
  {
    const unsigned nverts = apf::Mesh::adjacentCount[apfType][0];
    unsigned* vtx = r->edgeVerts;
    unsigned* tags = r->edgeTags;
    for(unsigned id=0; id<nbdry; id++) {
      apf::Downward verts;
      for(unsigned j=0; j<nverts; j++) {
        verts[j] = lookupVert(r, vtx[id*nverts+j]);
        apf::Vector3 vtx_coord;
        r->mesh->getPoint(verts[j], 0, vtx_coord);
      }
      apf::MeshEntity* f =
        apf::findElement(r->mesh, apfType, verts);
      PCU_ALWAYS_ASSERT(f);
      apf::ModelEntity* g = r->mesh->findModelEntity(1, tags[id]);
      r->mesh->setModelEntity(f, g);
    }
    free(vtx);
    free(tags);
  }

  void classifyBoundaryElms(Reader* r, header* h) {
    classifyBoundaryElms(r,h->nbdry,apf::Mesh::EDGE);
  }

  void classifyVtx(Reader *r, header* h) {
    (void)h;
    apf::Mesh2* m = r->mesh;
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* vtx;
    while ((vtx = m->iterate(it)))
    {
      int num_up = m->countUpward(vtx);
      std::vector<int> upward_dim(num_up);
      std::vector<int> upward_id(num_up);
      for (int id = 0; id < num_up; id++)
      {
        apf::MeshEntity* ment = m->getUpward(vtx, id);
        apf::ModelEntity* gent = m->toModel(ment);
        upward_dim[id] = m->getModelType(gent);
        upward_id[id] = m->getModelTag(gent);
      }

      apf::Vector3 vtx_coord;
      m->getPoint(vtx, 0, vtx_coord);
      bool same_dim = true;
      for(size_t i=0; i<upward_dim.size(); i++) {
         if (upward_dim[0] != upward_dim[i]) {
           same_dim = false;
           break;
         }
      }
      bool same_id = true;
      for(size_t i=0; i<upward_id.size(); i++) {
         if (upward_id[0] != upward_id[i]) {
           same_id = false;
           break;
         }
      }
      if (same_dim && same_id)
      {
        /// if all edges adjacent to a vertex have the same classification
        ///   then classify the vertex on that same geometric entity
        apf::ModelEntity* gent = m->findModelEntity(upward_dim[0], upward_id[0]);
        m->setModelEntity(vtx, gent);

        if (gmi_can_get_closest_point(m->getModel()))
        {
          apf::Vector3 from, to, param;
          m->getPoint(vtx, 0, from);
          m->getClosestPoint(gent, from, to, param);
          m->setParam(vtx, param);
        }
      }
      else
      {
        /// find all the indices in the vectors where the model entity is 1D
        std::vector<int> edge_indx;
        for(size_t i = 0; i < upward_dim.size(); i++) {
          if(upward_dim[i] == 1) {
            edge_indx.push_back(i);
          }
        }
        std::vector<int> edge_id;
        for (size_t i = 0; i < edge_indx.size(); i++)
        {
          edge_id.push_back(upward_id[edge_indx[i]]);
        }
        bool same_edge_id = true;
        for(size_t i=0; i<edge_id.size(); i++) {
          if (edge_id[0] != edge_id[i]) {
            same_id = false;
            break;
          }
        }
        if (same_edge_id)
        {
          /// if the edges adjacent to a vertex that are classified on a dim 1
          ///   model entity have the same classification then classify the 
          ///   vertex on that same geometric entity.
          /// (for a vertex on the boundary)
          apf::ModelEntity* gent = m->findModelEntity(1, edge_id[0]);
          m->setModelEntity(vtx, gent);

          /// specify parametric coordinate of vertex on edge
          if (gmi_can_get_closest_point(m->getModel()))
          {
            apf::Vector3 from, to, param;
            m->getPoint(vtx, 0, from);
            m->getClosestPoint(gent, from, to, param);
            m->setParam(vtx, param);
          }
        }
        else /// a mesh vertex whose adjacent edges are classified on different
             ///   model edges must be classified on a model vertex
        {
          apf::Vector3 vtx_coord;
          m->getPoint(vtx, 0, vtx_coord);

          /// it doesn't matter which edge
          apf::ModelEntity* gent = m->findModelEntity(1, edge_id[0]);
          gmi_set* adjacent_verts = gmi_adjacent(m->getModel(), (gmi_ent*)gent, 0);
          int n_adj_verts = adjacent_verts->n;
          double p[2];
          double x[3];
          for (int j = 0; j < n_adj_verts; j++)
          {
            gmi_eval(m->getModel(), adjacent_verts->e[j], p, x);
            apf::Vector3 vec_x(x);

            /// only look at x and y dimensions of vector, model must be in x-y plane, but need not be at z=0
            double mag = pow(pow((vtx_coord[0] - vec_x[0]), 2) + pow((vtx_coord[1] - vec_x[1]), 2), 0.5);
            if (mag < 0.001)
            {
              m->setModelEntity(vtx, (apf::ModelEntity*)adjacent_verts->e[j]);
              break;
            }
          }
          gmi_free_set(adjacent_verts);
        }
      }
    }
  }

  inline unsigned ugridToMdsElmIdx(int apfType, int ugridIdx) {
    static int ugrid_to_mds_verts[6][8] = {
      {0, 1, 2,-1, -1, -1, -1, -1}, //tri
      {0, 1, 2, 3, -1, -1, -1, -1}, //quad
      {0, 1, 2, 3, -1, -1, -1, -1}, //tet
      {0, 1, 2, 3,  4,  5,  6,  7}, //hex
      {0, 1, 2, 3,  4,  5, -1, -1}, //prism
      {3, 2, 4, 0,  1, -1, -1, -1}  //pyramid
    };
    PCU_ALWAYS_ASSERT(apfType >= 2);
    return ugrid_to_mds_verts[apfType-2][ugridIdx];
  }

  void readElms(Reader* r, unsigned nelms, int apfType) {
    const unsigned nverts = apf::Mesh::adjacentCount[apfType][0];
    apf::ModelEntity* g = r->mesh->findModelEntity(3, 0);
    size_t cnt = nelms*nverts;
    unsigned* vtx = (unsigned*) calloc(cnt,sizeof(unsigned));
    readUnsigneds(r->file, vtx, cnt, r->swapBytes);
    for(unsigned i=0; i<nelms; i++) {
      apf::Downward verts;
      for(unsigned j=0; j<nverts; j++) {
        const unsigned mdsIdx = ugridToMdsElmIdx(apfType,j);
        verts[mdsIdx] = lookupVert(r, vtx[i*nverts+j]);
      }
      apf::MeshEntity* elm =
        apf::buildElement(r->mesh, g, apfType, verts);
      PCU_ALWAYS_ASSERT(elm);
    }
    free(vtx);
  }

  void readElms(Reader* r, header* h) {
    readElms(r,h->ntet,apf::Mesh::TET);
    readElms(r,h->npyr,apf::Mesh::PYRAMID);
    readElms(r,h->nprz,apf::Mesh::PRISM);
    readElms(r,h->nhex,apf::Mesh::HEX);
  }

  void read2DElms(Reader* r, unsigned nelms, int apfType) {
    const unsigned nverts = apf::Mesh::adjacentCount[apfType][0];
    size_t cnt = nelms*nverts;
    unsigned* vtx = (unsigned*) calloc(cnt,sizeof(unsigned));
    readUnsigneds(r->file, vtx, cnt, r->swapBytes);

    unsigned* elm_model_id = (unsigned*) calloc(nelms,sizeof(unsigned));
    readUnsigneds(r->file, elm_model_id, nelms, r->swapBytes);

    for(unsigned i=0; i<nelms; i++) {
      apf::Downward verts;
      for(unsigned j=0; j<nverts; j++) {
        const unsigned mdsIdx = ugridToMdsElmIdx(apfType,j);
        verts[mdsIdx] = lookupVert(r, vtx[i*nverts+j]);
      }
      apf::ModelEntity* g = r->mesh->findModelEntity(2, elm_model_id[i]);
      apf::MeshEntity* elm = apf::buildElement(r->mesh, g, apfType, verts);
      PCU_ALWAYS_ASSERT(elm);

      r->mesh->setModelEntity(elm, g);
    }
    free(vtx);
    free(elm_model_id);
  }

  void read2DElms(Reader* r, header* h) {
    read2DElms(r,h->ntri,apf::Mesh::TRIANGLE);
    read2DElms(r,h->nquad,apf::Mesh::QUAD);
  }

  void freeReader(Reader* r) {
    fclose(r->file);
  }

  void read3DUgrid(apf::Mesh2* m, const char* filename)
  {
    header hdr;
    Reader r;
    initReader(&r, m, filename);
    readHeader(&r, &hdr);
    hdr.print();
    readNodes(&r, &hdr);
    setNodeIds(&r, &hdr);
    readFacesAndTags(&r,&hdr);
    checkFilePos(&r,&hdr);
    readElms(&r,&hdr);
    setFaceTags(&r,&hdr);
    freeReader(&r);
    m->acceptChanges();
  }

  void read2DUgrid(apf::Mesh2* m, const char* filename)
  {
    header hdr;
    Reader r;
    initReader(&r, m, filename);
    readHeader(&r, &hdr);
    hdr.print();
    readNodes(&r, &hdr);
    setNodeIds(&r, &hdr);
    read2DElms(&r,&hdr);
    checkFilePos(&r,&hdr);
    readNumBdryElms(&r, &hdr);
    readBdryElmsAndTags(&r, &hdr);
    setBoundaryTags(&r,&hdr);
    classifyBoundaryElms(&r, &hdr);
    classifyVtx(&r, &hdr);
    freeReader(&r);
    m->acceptChanges();
  }

  void getMaxAndAvg(std::set<int>*& cnt, int numparts, int& max, double& avg) {
    for(int i=0; i<numparts; i++) {
      avg += cnt[i].size();
      if( (double)cnt[i].size() > max )
        max = cnt[i].size();
    }
    avg /= numparts;
  }

  template<typename T>
  void getMaxAndAvg(T* cnt, int numparts, double& max, double& avg) {
    for(int i=0; i<numparts; i++) {
      avg += cnt[i];
      if( cnt[i] > max )
        max = cnt[i];
    }
    avg /= numparts;
  }

  class ptnstats {
    public:
      int numparts;
      int* ptn;         //vertex id to part id array
      SetInt* partvtx;  //vtx ids for each part
      int* partelm;     //elm cnts
      double* partelmW; //weighted elm counts

      ptnstats() : numparts(0), ptn(NULL), partvtx(NULL),
                   partelm(NULL), partelmW(NULL) {}
      void getVtxPtn(int numVtx, const char* ptnFile) {
        FILE* f = fopen(ptnFile, "r");
        numparts = 0;
        ptn = new int[numVtx];
        for(long id=0; id<numVtx; id++) {
          int read = fscanf(f, "%d", &ptn[id]);
          PCU_ALWAYS_ASSERT(read);
          if( ptn[id] > numparts )
            numparts = ptn[id];
        }
        numparts++; //we want count, not rank
        fclose(f);
      }
      ~ptnstats() {
        delete [] ptn;
        delete [] partelm;
        delete [] partelmW;
        delete [] partvtx;
      }
      void getOwnedVtx(int nvtx) {
        //count number of non-ghosted vtx per part
        for(long id=0; id<nvtx; id++) {
          const int partid = ptn[id];
          PCU_ALWAYS_ASSERT(partid >= 0 && partid < numparts);
          PCU_ALWAYS_ASSERT(!partvtx[partid].count(id));
          partvtx[partid].insert(id);
        }
      }
      void setup() {
        partvtx = new SetInt[numparts];
        partelm = new int[numparts];
        partelmW = new double[numparts];
        for(int i=0; i<numparts; i++)
          partelm[i] = partelmW[i] = 0;
      }
  };

  void readElms(Reader* r, unsigned nelms, int apfType,
      const double weight, ptnstats& ps) {
    const unsigned nverts = apf::Mesh::adjacentCount[apfType][0];
    size_t cnt = nelms*nverts;
    unsigned* vtx = (unsigned*) calloc(cnt,sizeof(unsigned));
    readUnsigneds(r->file, vtx, cnt, r->swapBytes);
    for(unsigned i=0; i<nelms; i++) {
      SetInt elmparts;
      //determine which parts have a copy of the element
      //  based on who owns the elements bounding vertices
      for(unsigned j=0; j<nverts; j++) {
        const unsigned idx = i*nverts+j;
        PCU_ALWAYS_ASSERT(idx < cnt);
        const int vtxid = ftnToC(vtx[idx]);
        const int partid = ps.ptn[vtxid];
        //an element can only exist once on each part
        elmparts.insert(partid);
      }
      PCU_ALWAYS_ASSERT(elmparts.size() > 0 && elmparts.size() <= nverts);
      //increment the elm per part counts
      APF_ITERATE(SetInt,elmparts,ep) {
        const int partid = *ep;
        PCU_ALWAYS_ASSERT(partid >= 0 && partid < ps.numparts);
        ps.partelmW[partid] += weight;
        ps.partelm[partid]++;
        //add ghost vertices
        for(unsigned j=0; j<nverts; j++) {
          const int vtxid = ftnToC(vtx[i*nverts+j]);
          ps.partvtx[partid].insert(vtxid);
        }
      }
    }
    free(vtx);
  }

  void printPtnStats(apf::Mesh2* m, const char* ufile, const char* ptnFile,
      const double elmWeights[]) {
    header hdr;
    Reader r;
    initReader(&r, m, ufile);
    readHeader(&r, &hdr);
    hdr.print();

    ptnstats ps;
    ps.getVtxPtn(hdr.nvtx,ptnFile);
    ps.setup();
    ps.getOwnedVtx(hdr.nvtx);

    readNodes(&r,&hdr); //dummy to advance the fileptr
    readFacesAndTags(&r,&hdr); //dummy to advance the fileptr
    free(r.faceVerts[faceTypeIdx(apf::Mesh::TRIANGLE)]);
    free(r.faceVerts[faceTypeIdx(apf::Mesh::QUAD)]);
    free(r.faceTags[faceTypeIdx(apf::Mesh::TRIANGLE)]);
    free(r.faceTags[faceTypeIdx(apf::Mesh::QUAD)]);
    checkFilePos(&r,&hdr); //sanity check

    int type[] = {apf::Mesh::TET,apf::Mesh::PYRAMID,apf::Mesh::PRISM,apf::Mesh::HEX};
    unsigned cnt[] = {hdr.ntet, hdr.npyr, hdr.nprz, hdr.nhex};
    for(int i=0; i<4; i++)
      readElms(&r,cnt[i],type[i], elmWeights[ type[i] ], ps);

    //get max and avg vtx and elm per part
    double maxelm = 0; double avgelm = 0;
    double maxelmW = 0; double avgelmW = 0;
    int maxvtx = 0; double avgvtx = 0;
    getMaxAndAvg(ps.partvtx,ps.numparts,maxvtx,avgvtx);
    getMaxAndAvg(ps.partelm,ps.numparts,maxelm,avgelm);
    getMaxAndAvg(ps.partelmW,ps.numparts,maxelmW,avgelmW);
    double imbvtx = maxvtx / avgvtx;
    double imbelm = maxelm / avgelm;
    double imbelmW = maxelmW / avgelmW;
    lion_eprint(1, "imbvtx %.3f imbelmW %.3f imbelm %.3f "
        "avgvtx %.3f avgelmW %.3f avgelm %.3f\n",
        imbvtx, imbelmW, imbelm, avgvtx, avgelmW, avgelm);
  }
}

namespace apf {
  Mesh2* loadMdsFromUgrid(gmi_model* g, const char* filename, pcu::PCU *PCUObj)
  {
    Mesh2* m = makeEmptyMdsMesh(g, 0, false, PCUObj);
    header hdr;
    Reader r;
    initReader(&r, m, filename);
    readHeader(&r, &hdr);
    int dim = 2;
    if (hdr.ntet != 0 || hdr.npyr != 0 || hdr.nprz != 0 || hdr.nhex != 0)
      dim = 3;

    apf::changeMdsDimension(m, dim);
    freeReader(&r);
    if (2 == dim)
      read2DUgrid(m, filename);
    else if (3 == dim)
      read3DUgrid(m, filename);
    lion_eprint(1,"vtx %lu edge %lu face %lu rgn %lu\n",
        m->count(0), m->count(1), m->count(2), m->count(3));
    return m;
  }

  void printUgridPtnStats(gmi_model* g, const char* ufile, const char* vtxptn,
      const double elmWeights[], pcu::PCU *PCUObj) {
    Mesh2* m = makeEmptyMdsMesh(g, 0, false, PCUObj);
    apf::changeMdsDimension(m, 3);
    printPtnStats(m, ufile, vtxptn, elmWeights);
    m->destroyNative();
    apf::destroyMesh(m);
  }
}
