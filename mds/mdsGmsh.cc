#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include "gmi.h" /* this is for gmi_getline... */
#include <lionPrint.h>

#include <cstdio>
#include <cstring>
#include <pcu_util.h>
#include <cstdlib>

namespace {

bool isQuadratic(int gmshType)
{
  switch (gmshType) {
    case 8: return true;
    case 9: return true;
    case 11: return true;
    default: return false;
  }
}

int apfFromGmsh(int gmshType)
{
  switch (gmshType) {
    case 1: return apf::Mesh::EDGE;
    case 2: return apf::Mesh::TRIANGLE;
    case 3: return apf::Mesh::QUAD;
    case 4: return apf::Mesh::TET;
    case 5: return apf::Mesh::HEX;
    case 6: return apf::Mesh::PRISM;
    case 7: return apf::Mesh::PYRAMID;
    case 8: return apf::Mesh::EDGE;
    case 9: return apf::Mesh::TRIANGLE;
    case 11: return apf::Mesh::TET;
    case 15: return apf::Mesh::VERTEX;
    default: return -1;
  }
}

struct Node {
  Node() { entity = 0; }
  apf::MeshEntity* entity;
  apf::Vector3 point;
};

struct Reader {
  apf::Mesh2* mesh;
  FILE* file;
  char* line;
  char* word;
  size_t linecap;
  int major_version;
  int minor_version;
  bool isQuadratic;
  std::map<long, Node> nodeMap;
  std::map<long, apf::MeshEntity*> entMap[4];
  //the 0th vector is not used as mesh vertices don't have a 'physical entity'
  //association in the legacy 2.* gmsh format
  std::vector<int> physicalType[4];
};


void freeReader(Reader* r)
{
  free(r->line);
  fclose(r->file);
}

void getLine(Reader* r)
{
  ssize_t ret = gmi_getline(&r->line, &r->linecap, r->file);
  PCU_ALWAYS_ASSERT(ret != -1);
  r->word = r->line;
}

long getLong(Reader* r)
{
  long x;
  int nchars;
  int ret = sscanf(r->word, "%ld%n", &x, &nchars);
  PCU_ALWAYS_ASSERT(ret == 1);
  r->word += nchars;
  return x;
}

double getDouble(Reader* r)
{
  double x;
  int nchars;
  int ret = sscanf(r->word, "%lf%n", &x, &nchars);
  PCU_ALWAYS_ASSERT(ret == 1);
  r->word += nchars;
  return x;
}

bool startsWith(char const* prefix, char const* s)
{
  int ls = strlen(s);
  int lp = strlen(prefix);
  if (ls < lp)
    return false;
  return strncmp(s, prefix, lp) == 0;
}

void seekMarker(Reader* r, char const* marker)
{
  while (!startsWith(marker, r->line))
    getLine(r);
  getLine(r);
}

void checkMarker(Reader* r, char const* marker)
{
  PCU_ALWAYS_ASSERT(startsWith(marker, r->line));
}

void initReader(Reader* r, apf::Mesh2* m, const char* filename)
{
  r->mesh = m;
  r->file = fopen(filename, "r");
  if (!r->file) {
    lion_eprint(1,"couldn't open Gmsh file \"%s\"\n",filename);
    abort();
  }
  r->line = static_cast<char*>(malloc(1));
  r->line[0] = '\0';
  r->linecap = 1;
  r->isQuadratic = false;
  seekMarker(r, "$MeshFormat");
  int fileType, dataSize;
  int ret = sscanf(r->line, "%d.%d %d %d\n",
      &r->major_version, &r->minor_version, &fileType, &dataSize);
  PCU_ALWAYS_ASSERT(ret==4);
}

void readNode(Reader* r, int bm=-1)
{
  Node n;
  apf::Vector3& p = n.point;
  if(r->major_version == 4) {
    sscanf(r->line, "%lf %lf %lf", &p[0], &p[1], &p[2]);
    r->nodeMap[bm] = n;
  } else if(r->major_version == 2) {
    long id;
    sscanf(r->line, "%ld %lf %lf %lf", &id, &p[0], &p[1], &p[2]);
    r->nodeMap[id] = n;
  }
  getLine(r);
}

void readEntities(Reader* r,const char* fnameDmg) 
{
  seekMarker(r, "$Entities");
  long nlde,ilde,iud,tag,isign,nMV,nME,nMF,nMR;
  double x,y,z;
  FILE* f = fopen(fnameDmg, "w");
  sscanf(r->line, "%ld %ld %ld %ld", &nMV, &nME, &nMF, &nMR);
  fprintf(f, "%ld %ld %ld %ld \n", nMR, nMF, nME, nMV); // just reverse order 
  fprintf(f, "%f %f %f \n ", 0.0, 0.0, 0.0); // Probaby model bounding box?
  fprintf(f, "%f %f %f \n", 0.0, 0.0, 0.0); // 
   
  getLine(r); // because readNode gets the next line we need this outside  for Nodes_Block
  for (long i = 0; i < nMV; ++i){  
    sscanf(r->line, "%ld %lf %lf %lf %ld ", &tag, &x, &y, &z, &iud);
    fprintf(f, "%ld %lf %lf %lf \n",tag,x,y,z);
    getLine(r); 
  }
  for (long i = 0; i < nME; ++i){
    tag = getLong(r); 
    fprintf(f, "%ld", tag);
    for (int i=0; i< 6; ++i) x=getDouble(r);  // read past min maxes
    iud = getLong(r); 
    for(long j =0; j < iud; ++j) isign=getLong(r); // read past iud user tags
    nlde=getLong(r);  // 2 in straight edged models but...
    for(long j =0; j < nlde; ++j) {
      ilde=getLong(r);
      fprintf(f, " %ld", std::abs(ilde)); // modVerts started from 1
    }
    fprintf(f, "\n");
    getLine(r); 
  }
  for (long i = 0; i < nMF; ++i){
    tag = getLong(r); 
    fprintf(f, "%ld %d\n", tag, 1);
    for (int i=0; i< 6; ++i) x=getDouble(r);  // read past min maxes
    iud = getLong(r); 
    for(long j =0; j < iud; ++j) isign=getLong(r); // read past iud user tags
    nlde=getLong(r); 
    fprintf(f, "  %ld \n", nlde);
    for(long j =0; j < nlde; ++j) {
      ilde=getLong(r); 
      if(ilde > 0 ) 
        isign=1;
      else
        isign=0;
      fprintf(f, "    %ld %ld \n", std::abs(ilde),isign); 
    }
    getLine(r); 
  }   
  for (long i = 0; i < nMR; ++i){ 
    tag = getLong(r); 
    fprintf(f, "%ld %d \n", tag, 1);
    for (int i=0; i< 6; ++i) x=getDouble(r);  // read past min maxes
    iud = getLong(r); 
    for(long j =0; j < iud; ++j) getLong(r); // read past iud user tags
    nlde=getLong(r); 
    fprintf(f, "%ld \n", nlde);
    for(long j =0; j < nlde; ++j) {
      ilde=getLong(r); 
      if(ilde > 0 ) 
        isign=1;
      else
        isign=0;
      fprintf(f, "%ld %ld \n", std::abs(ilde),isign); 
    }
    getLine(r); 
  }   
  checkMarker(r, "$EndEntities");
  fclose(f);
}

void readNodesV2(Reader* r)
{
  PCU_ALWAYS_ASSERT(r->major_version == 2);
  seekMarker(r, "$Nodes");
  long n = getLong(r);
  getLine(r);
  for (long i = 0; i < n; ++i)
    readNode(r);
  checkMarker(r, "$EndNodes");
}

void readNodesV4(Reader* r)
{
  PCU_ALWAYS_ASSERT(r->major_version == 4);
  seekMarker(r, "$Nodes");
  long Num_EntityBlocks,Num_Nodes,Nodes_Block,edim,etag,junk1,junk2,junk3;
  sscanf(r->line, "%ld %ld %ld %ld", &Num_EntityBlocks, &Num_Nodes, &junk1, &junk2);
  getLine(r); // because readNode gets the next line we need this outside  for Nodes_Block
  for (long i = 0; i < Num_EntityBlocks; ++i){
    sscanf(r->line, "%ld %ld %ld %ld", &edim, &etag, &junk3, &Nodes_Block);
    long* blockMap = new long[Nodes_Block];
    for (long j = 0; j < Nodes_Block; ++j){
      getLine(r);
      sscanf(r->line, "%ld", &blockMap[j]);
    }
    getLine(r);
    for (long j = 0; j < Nodes_Block; ++j)
      readNode(r,blockMap[j]);  // has a genLine at end
    delete [] blockMap;
  }    
  checkMarker(r, "$EndNodes");
}

apf::MeshEntity* lookupVert(Reader* r, long nodeId, apf::ModelEntity* g)
{
  PCU_ALWAYS_ASSERT(r->nodeMap.count(nodeId));
  Node& n = r->nodeMap[nodeId];
  if (n.entity)
    return n.entity;
  n.entity = r->mesh->createVert(g);
  r->mesh->setPoint(n.entity, 0, n.point);
  return n.entity;
}

void readElement(Reader* r, long gmshType=-1, long gtag=-1)
{
  long id = getLong(r);
  if(r->major_version == 2) {
    gmshType = getLong(r);
  }
  if (isQuadratic(gmshType))
    r->isQuadratic = true;
  int apfType = apfFromGmsh(gmshType);
  PCU_ALWAYS_ASSERT(0 <= apfType);
  int nverts = apf::Mesh::adjacentCount[apfType][0];
  int dim = apf::Mesh::typeDimension[apfType];
  if(r->major_version == 2) {
    long ntags = getLong(r);
    /* The Gmsh 4.9 documentation on the legacy 2.* format states:
     * "By default, the first tag is the tag of the physical entity to which the
     * element belongs; the second is the tag of the elementary model entity to
     * which the element belongs; the third is the number of mesh partitions to
     * which the element belongs, followed by the partition ids (negative
     * partition ids indicate ghost cells). A zero tag is equivalent to no tag.
     * Gmsh and most codes using the MSH 2 format require at least the first two
     * tags (physical and elementary tags)."
     * A physical entity is a user defined grouping of elementary model entities.
     * An elementary model entity is a geometric model entity. */
    PCU_ALWAYS_ASSERT(ntags >= 2);
    const int physType = static_cast<int>(getLong(r));
    PCU_ALWAYS_ASSERT(dim>=0 && dim<4);
    r->physicalType[dim].push_back(physType);
    gtag = getLong(r);
    for (long i = 2; i < ntags; ++i)
      getLong(r); /* discard all other element tags */
  }
  apf::ModelEntity* g = r->mesh->findModelEntity(dim, gtag);
  apf::Downward verts;
  for (int i = 0; i < nverts; ++i) {
    long nid = getLong(r);
    verts[i] = lookupVert(r, nid, g);
  }
  if (dim != 0) {
    if (dim > r->mesh->getDimension())
      apf::changeMdsDimension(r->mesh, dim);
    apf::MeshEntity* ent = apf::buildElement(r->mesh, g, apfType, verts);
    if (r->isQuadratic)
      r->entMap[dim][id] = ent;
  }
  getLine(r);
}

void readElementsV2(Reader* r)
{
  PCU_ALWAYS_ASSERT(r->major_version == 2);
  seekMarker(r, "$Elements");
  long n = getLong(r);
  getLine(r);
  for (long i = 0; i < n; ++i)
    readElement(r);
  checkMarker(r, "$EndElements");
}

void readElementsV4(Reader* r)
{
  PCU_ALWAYS_ASSERT(r->major_version == 4);
  seekMarker(r, "$Elements");
  long Num_EntityBlocks,Num_Elements,Elements_Block,Edim,gtag,gmshType,junk1,junk2;
  sscanf(r->line, "%ld %ld %ld %ld", &Num_EntityBlocks, &Num_Elements, &junk1, &junk2);
  getLine(r);
  for (long i = 0; i < Num_EntityBlocks; ++i){
    sscanf(r->line, "%ld %ld %ld %ld", &Edim, &gtag, &gmshType, &Elements_Block);
    getLine(r);
    for (long j = 0; j < Elements_Block; ++j) {
       readElement(r,gmshType,gtag);
    }
  }
  checkMarker(r, "$EndElements");
}

void setElmPhysicalType(Reader* r, apf::Mesh2* m) {
  apf::MeshEntity* e;
  apf::MeshTag* tag = m->createIntTag("gmsh_physical_entity", 1);
  for(int dim=1; dim<=m->getDimension(); dim++) { //vertices don't have a physical entity ?
    if( ! r->physicalType[dim].size() ) continue;
    int* tagPtr = r->physicalType[dim].data();
    int i = 0;
    apf::MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it)))
      m->setIntTag(e, tag, &tagPtr[i++]);
    m->end(it);
  }
}

static const double gmshTet10EdgeIndices[6] = {0, 1, 2, 3, 5, 4};

static int getQuadGmshIdx(const int apfIdx, const int apfType) {
  switch (apfType) {
    case apf::Mesh::TET :
      return gmshTet10EdgeIndices[apfIdx];
      break;
    default:
      return apfIdx;
  }
}

void readQuadraticElement(Reader* r)
{
  long id = getLong(r);
  long gmshType = getLong(r);
  int apfType = apfFromGmsh(gmshType);
  if (apfType == apf::Mesh::VERTEX) return;
  PCU_ALWAYS_ASSERT(0 <= apfType);
  PCU_ALWAYS_ASSERT_VERBOSE(isQuadratic(gmshType),
      "no support for variable p-order meshes");
  int nverts = apf::Mesh::adjacentCount[apfType][0];
  int nedges = apf::Mesh::adjacentCount[apfType][1];
  int dim = apf::Mesh::typeDimension[apfType];
  long ntags = getLong(r);
  getLong(r); /* discard physical type */
  getLong(r); /* discard geometric tag */
  for (long i = 2; i < ntags; ++i)
    getLong(r); /* discard all other tags */
  for (long i = 0; i < nverts; ++i)
    getLong(r); /* discard all vertex nodes */
  apf::Downward edges;
  apf::MeshEntity* ent = r->entMap[dim][id];
  r->mesh->getDownward(ent, 1, edges);
  apf::Field* coord = r->mesh->getCoordinateField();
  std::vector<long> nids(nedges);
  for (long i = 0; i < nedges; ++i)
    nids[i] = getLong(r);
  for (long i = 0; i < nedges; ++i) {
    long nid = nids[ getQuadGmshIdx(i, apfType) ];
    apf::Vector3 point = r->nodeMap[nid].point;
    apf::setVector(coord, edges[i], 0, point);
  }
}

void readQuadratic(Reader* r, apf::Mesh2* m, const char* filename)
{
  m->changeShape(apf::getSerendipity());
  initReader(r, m, filename);
  seekMarker(r, "$Elements");
  long n = getLong(r);
  getLine(r);
  for (long i = 0; i < n; ++i) {
    readQuadraticElement(r);
    getLine(r);
  }
  checkMarker(r, "$EndElements");
  freeReader(r);
}

void readGmsh(apf::Mesh2* m, const char* filename)
{
  Reader r;
  initReader(&r, m, filename);
  if(r.major_version == 4) {
    readNodesV4(&r);
    readElementsV4(&r);
    m->acceptChanges();
  } else if(r.major_version == 2) {
    readNodesV2(&r);
    readElementsV2(&r);
    m->acceptChanges();
    setElmPhysicalType(&r,m);
  }
  if (r.isQuadratic)
    readQuadratic(&r, m, filename);
  freeReader(&r);
}
}  // closes original namespace 

namespace apf {

int gmshMajorVersion(const char* filename) {
  Reader r;
  Mesh2* m=NULL;
  initReader(&r, m,  filename);
  int version = r.major_version;
  freeReader(&r);
  return version;
}

void gmshFindDmg(const char* fnameDmg, const char* filename)
{
  Reader r;
  
  Mesh2* m=NULL;
  initReader(&r, m,  filename);
  PCU_ALWAYS_ASSERT(r.major_version==4);
  readEntities(&r, fnameDmg);
  freeReader(&r);
}


Mesh2* loadMdsFromGmsh(gmi_model* g, const char* filename, pcu::PCU *PCUObj)
{
  Mesh2* m = makeEmptyMdsMesh(g, 0, false, PCUObj);
  readGmsh(m, filename);
  return m;
}

Mesh2* loadMdsDmgFromGmsh(const char*fnameDmg, const char* filename, pcu::PCU *PCUObj)
{
  gmshFindDmg(fnameDmg, filename);  // new function that scans $Entities and writes a dmg 
  Mesh2* m = makeEmptyMdsMesh(gmi_load(fnameDmg), 0, false, PCUObj);
  readGmsh(m, filename);
  return m;
}

}
