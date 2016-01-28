#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "pcu_io.h"
#include "pcu_byteorder.h"

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cstdlib>

/*
read files in the AFLR3 format from Dave Marcum at Mississippi State
-little-endian if file has suffix '.lb8.ugrid'
-big-endian if file has suffix '.b8.ugrid'
-integers are 4B
-floats are 8B
*/

namespace {
  int faceTypeIdx(int type) {
    switch(type) {
      case apf::Mesh::TRIANGLE: return 0;
      case apf::Mesh::QUAD: return 1;
      default: {
                 assert(false);
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
    bool swapBytes;
  };

  struct header {
    unsigned nvtx, ntri, nquad, ntet, npyr, nprz, nhex;
    void print() {
      fprintf(stderr,
          "nvtx %u ntri %u nquad %u ntet %u npyr %u nprz %u nhex %u\n",
          nvtx,ntri,nquad,ntet,npyr,nprz,nhex);
    }
  };

  void initReader(Reader* r, apf::Mesh2* m, const char* filename) {
    r->mesh = m;
    r->file = fopen(filename, "rb");
    if (!r->file) {
      fprintf(stderr,"ERROR couldn't open ugrid file \"%s\"\n",filename);
      abort();
    }
    unsigned endian = -1;
    if ( strstr(filename, ".b8.ugrid") ) {
      endian = PCU_BIG_ENDIAN;
    } else if ( strstr(filename, ".lb8.ugrid") ) {
      endian = PCU_LITTLE_ENDIAN;
    } else {
      fprintf(stderr,
          "ERROR file extension of \"%s\" is not supported\n", filename);
      exit(EXIT_FAILURE);
    }
    r->swapBytes = ( endian != PCU_HOST_ORDER );
  }

  void readUnsigneds(FILE* f, unsigned* v, size_t cnt, bool swap) {
    assert(sizeof(unsigned)==4);
    size_t read = fread(v, sizeof(unsigned), cnt, f);
    assert(read == cnt);
    if ( swap )
      pcu_swap_unsigneds(v,cnt);
  }

  void readDoubles(FILE* f, double* vals, size_t cnt, bool swap) {
    assert(sizeof(double)==8);
    size_t read = fread(vals, sizeof(double), cnt, f);
    assert(read == cnt);
    if ( swap )
      pcu_swap_doubles(vals,cnt);
  }

  void readHeader(Reader* r, header* h) {
    readUnsigneds(r->file, &h->nvtx, 1, r->swapBytes);
    readUnsigneds(r->file, &h->ntri, 1, r->swapBytes);
    readUnsigneds(r->file, &h->nquad, 1, r->swapBytes);
    readUnsigneds(r->file, &h->ntet, 1, r->swapBytes);
    readUnsigneds(r->file, &h->npyr, 1, r->swapBytes);
    readUnsigneds(r->file, &h->nprz, 1, r->swapBytes);
    readUnsigneds(r->file, &h->nhex, 1, r->swapBytes);
  }

  apf::MeshEntity* makeVtx(Reader* r,
      apf::Vector3& pt, apf::ModelEntity* g) {
    apf::MeshEntity* v = r->mesh->createVert(g);
    r->mesh->setPoint(v, 0, pt);
    return v;
  }

  apf::MeshEntity* lookupVert(Reader* r, long ftnNodeId) {
    const long cNodeId = ftnNodeId - 1;
    assert(r->nodeMap.count(cNodeId));
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
    fprintf(stderr, "read %d vtx\n", h->nvtx);
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

  void checkFilePos(Reader* r, header* h) {
    // seven headers, vtx coords, face vertex ids, and face tags
    long expected = h->nvtx*3*sizeof(double) +
      sizeof(unsigned)*(7 + h->ntri*3 + h->nquad*4 + (h->ntri+h->nquad));
    long pos = ftell(r->file);
    assert( pos == expected );
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
      assert(f);
      int val = tags[id];
      r->mesh->setIntTag(f, t, &val);
    }
    free(vtx);
    free(tags);
    fprintf(stderr, "set %d %s face tags\n",
        nfaces, apf::Mesh::typeName[apfType]);
  }

  void setFaceTags(Reader* r, header* h) {
    apf::MeshTag* t = r->mesh->createIntTag("ugrid-face-tag", 1);
    setFaceTags(r,t,h->ntri,apf::Mesh::TRIANGLE);
    setFaceTags(r,t,h->nquad,apf::Mesh::QUAD);
  }

  inline unsigned ugridToMdsElmIdx(int apfType, int ugridIdx) {
    static int ugrid_to_mds_verts[4][8] = {
      {0, 1, 2, 3, -1, -1, -1, -1}, //tet
      {0, 1, 2, 3,  4,  5,  6,  7}, //hex
      {0, 1, 2, 3,  4,  5, -1, -1}, //prism
      {3, 2, 4, 0,  1, -1, -1, -1}  //pyramid
    };
    assert(apfType >= 4);
    return ugrid_to_mds_verts[apfType-4][ugridIdx];
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
      assert(elm);
    }
    free(vtx);
    fprintf(stderr, "read %d %s\n", nelms, apf::Mesh::typeName[apfType]);
  }

  void readElms(Reader* r, header* h) {
    readElms(r,h->ntet,apf::Mesh::TET);
    readElms(r,h->npyr,apf::Mesh::PYRAMID);
    readElms(r,h->nprz,apf::Mesh::PRISM);
    readElms(r,h->nhex,apf::Mesh::HEX);
  }

  void freeReader(Reader* r) {
    fclose(r->file);
  }

  void readUgrid(apf::Mesh2* m, const char* filename)
  {
    header hdr;
    Reader r;
    initReader(&r, m, filename);
    readHeader(&r, &hdr);
    hdr.print();
    readNodes(&r, &hdr);
    readFacesAndTags(&r,&hdr);
    checkFilePos(&r,&hdr);
    readElms(&r,&hdr);
    setFaceTags(&r,&hdr);
    freeReader(&r);
    m->acceptChanges();
  }
}

namespace apf {
  Mesh2* loadMdsFromUgrid(gmi_model* g, const char* filename)
  {
    Mesh2* m = makeEmptyMdsMesh(g, 0, false);
    apf::changeMdsDimension(m, 3);
    readUgrid(m, filename);
    fprintf(stderr,"vtx %lu edge %lu face %lu rgn %lu\n",
        m->count(0), m->count(1), m->count(2), m->count(3));
    return m;
  }
}
