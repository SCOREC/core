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
  bool isQuadratic;
  std::map<long, Node> nodeMap;
  std::map<long, apf::MeshEntity*> entMap[4];
  //the 0th vector is not used as mesh vertices don't have a 'physical entity'
  //association in the legacy 2.* gmsh format
  std::vector<int> physicalType[4];
};

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
}

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

void readNode(Reader* r, int bm)
{
  Node n;
  apf::Vector3& p = n.point;
  sscanf(r->line, "%lf %lf %lf", &p[0], &p[1], &p[2]);
  r->nodeMap[bm] = n;
  getLine(r);
}

void readEntities(Reader* r,const char* fnameDmg, int emap[], int* nMVskip, int pass)
{
  seekMarker(r, "$Entities");
  long nlde,ilde,iud,tag,isign;
  double x,y,z;
  FILE* f = fopen(fnameDmg, "w");
  sscanf(r->line, "%d %d %d %d", &emap[100], &emap[101], &emap[102], &emap[103]);
  if (pass==1) {
    int tmp=0;
    getLine(r); 
    for (long i = 0; i < emap[100];  ++i){ 
      sscanf(r->line, "%ld %lf %lf %lf %ld ", &tag, &x, &y, &z, &iud);
      if(iud==0) tmp++; // this is a gmsh construction model vertex that has no adjacency so has to be skipped in real read/build of dmg.
      getLine(r); 
   }
   *nMVskip=tmp;
   return;
  }
  else { // we have computed nMVskip in first pass
    emap[100]-=*nMVskip; // I can't seem to stop gmsh from writing construction vertices but we can require users to put them first 
    fprintf(f, "%d %d %d %d \n", emap[103], emap[102], emap[101], emap[100]); // just reverse order 
    fprintf(f, "%f %f %f \n ", 0.0, 0.0, 0.0); // Probaby model bounding box?
    fprintf(f, "%f %f %f \n", 0.0, 0.0, 0.0); // 
   
    getLine(r); // because readNode gets the next line we need this outside  for Nodes_Block
    int entCnt=0;
    for (long i = 0; i < emap[100] + *nMVskip; ++i){  // weird end reflects known skips but we need emap[100] correct in DMG
      sscanf(r->line, "%ld %lf %lf %lf %ld ", &tag, &x, &y, &z, &iud);
      if(iud !=0 ) {
        emap[entCnt]=tag;   // map(dmgTag)=gmshTag is backward how we need but gmshTags dup
        entCnt++; // vertex entities start from 1
        fprintf(f, "%d %lf %lf %lf \n",entCnt,x,y,z);
      }
      getLine(r); 
    }
    for (long i = 0; i < emap[101]; ++i){
 // FAIL does not advance    sscanf(r->line, "%ld %lf %lf %lf", &tag, &x, &y, &z);
      tag = getLong(r); 
      emap[entCnt]=tag;
      entCnt++;
      fprintf(f, "%d", entCnt);
      for (int i=0; i< 6; ++i) x=getDouble(r);  // read past min maxes
      iud = getLong(r); 
      for(long j =0; j < iud; ++j) isign=getLong(r); // read past iud user tags
      nlde=getLong(r);  // 2 in straight edged models but...
      for(long j =0; j < nlde; ++j) {
        ilde=getLong(r);
        int found=0;
        int k=0;
        while(!found){ // have to search since map is backwards
           if(emap[k] == abs(ilde)) found=1; 
           else k++;
        }
        fprintf(f, " %d", k+1); // modVerts started from 1
      }
      fprintf(f, "\n");
      getLine(r); 
    }   
    for (long i = 0; i < emap[102]; ++i){
      tag = getLong(r); 
      emap[entCnt]=tag; // new face tag map
      entCnt++;
      fprintf(f, "%d %d\n", entCnt, 1);
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
        int found=0;
        int k=emap[100]; // edges start at this count
        while(!found){ // have to search since map is backwards
           if(emap[k] == abs(ilde)) found=1; 
           else k++;
        }
        fprintf(f, "    %d %ld \n", k+1,isign); 
      }
      getLine(r); 
    }   
    for (long i = 0; i < emap[103]; ++i){ //not even sure that this all hangs with emap[103] > 1 but..
      tag = getLong(r); 
      emap[entCnt]=tag; // new region tag map
      entCnt++;
      fprintf(f, "%d %d \n", entCnt, 1);
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
        int found=0;
        int k=emap[100]+emap[101]; // faces start here
        while(!found){ // have to search since map is backwards
           if(emap[k] == abs(ilde)) found=1; 
           else k++;
        }
        fprintf(f, "%d %ld \n", k+1,isign); 
      }
      getLine(r); 
    }   
    checkMarker(r, "$EndEntities");
    fclose(f);
  }
}
void readNodes(Reader* r, int* emap)
{
  seekMarker(r, "$Nodes");
  long Num_EntityBlocks,Num_Nodes,Nodes_Block,edim,etag,junk1,junk2,junk3;
  sscanf(r->line, "%ld %ld %ld %ld", &Num_EntityBlocks, &Num_Nodes, &junk1, &junk2);
  getLine(r); // because readNode gets the next line we need this outside  for Nodes_Block
  for (long i = 0; i < Num_EntityBlocks; ++i){
    sscanf(r->line, "%ld %ld %ld %ld", &edim, &etag, &junk3, &Nodes_Block);
    int found=0;
    if(edim==0) {
      for (int k=0; k<emap[100]; ++k) {
         if(emap[k]==etag) found=1;
      }
    }
    if(found==1 || edim > 0){ 
      long blockMap[Nodes_Block];
      for (long j = 0; j < Nodes_Block; ++j){
        getLine(r);
        sscanf(r->line, "%ld", &blockMap[j]);
      }
      getLine(r);
      for (long j = 0; j < Nodes_Block; ++j)
        readNode(r,blockMap[j]);  // has a genLine at end
    } 
    else { // skip the construction nodes
      getLine(r); // block line only scanned
      getLine(r); // node number
      getLine(r); // coordinates
    } 
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

void readElement(Reader* r, long gmshType,long gtag)
{
  long id = getLong(r); // tag in 2 and 4
  if (isQuadratic(gmshType))
    r->isQuadratic = true;
  int apfType = apfFromGmsh(gmshType);
  PCU_ALWAYS_ASSERT(0 <= apfType);
  int nverts = apf::Mesh::adjacentCount[apfType][0];
  int dim = apf::Mesh::typeDimension[apfType];
  if(false) { // FIXME
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
//FIXME blocks compilation    long gtag = getLong(r);
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

/**
  emap[i] = gmshModelEntId where i=[0..99] is the dmg unique model entity id
  emap[100+d] = number of model verts of dimension d where d=[0..3]
  FIXME - extend the upper limit on model entity ids
*/

void readElements(Reader* r, int* emap)
{
  seekMarker(r, "$Elements");
  long Num_EntityBlocks,Num_Elements,Elements_Block,Edim,gtag,gmshType,junk1,junk2;
  sscanf(r->line, "%ld %ld %ld %ld", &Num_EntityBlocks, &Num_Elements, &junk1, &junk2);
  getLine(r); 
  int tagMapped;
  for (long i = 0; i < Num_EntityBlocks; ++i){
    sscanf(r->line, "%ld %ld %ld %ld", &Edim, &gtag, &gmshType, &Elements_Block);
    int found=0;
    if(Edim==0) { // This only determines if model vertex is Physical if yes found=1
      for (int k=0; k<emap[100]; ++k) {
         if(emap[k]==gtag) found=1;
      }
    } 
    if (found==1 || Edim > 0) {
      int kstart=0;
      for (int k=0; k<Edim; ++k) kstart+=emap[100+k]; //higher dim element start higher
      int found=0;
      int k=kstart;
      while(!found){ // have to search since map is backwards
         if(emap[k] == gtag) found=1; 
         else k++;
      }
      tagMapped= k+1; // modVerts started from 1
    }
    getLine(r);
    for (long j = 0; j < Elements_Block; ++j) {
     if(found==1 || Edim>0) 
       readElement(r,gmshType,tagMapped);
     else
       getLine(r);  // don't add non Physical nodes which will have found=0 (not in phsyical tag list)
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

void readGmsh(apf::Mesh2* m, const char* filename, int emap[])
{
  Reader r;
  initReader(&r, m, filename);
  readNodes(&r, emap); // as near as I can tell neither v2 nor my v4 mods actually USE the classification information...I suppose because V2 had none my format changes did not exploit it. 
  readElements(&r,emap);  // different story for Elements which do use the m that has the dmg info smuggle^2 though the reader r. Thuse we pass our emap down (I suppose we could put it in r too?
  freeReader(&r);
  m->acceptChanges();
  if(false) // FIXME
    setElmPhysicalType(&r,m);
  freeReader(&r);
  if (r.isQuadratic)
    readQuadratic(&r, m, filename);
}
}  // closes original namespace 

namespace apf {
void gmshFindDmg(const char* fnameDmg, const char* filename, int emap[])
{
  Reader r;
  
  Mesh2* m=NULL;
  initReader(&r, m,  filename);
  int nMVskip=0;  // first pass is a scan of model vertices to find those that don't have Physical tags  and are thus construction if we reuire users to put true model vertices into Physical Groups. 
  readEntities(&r, fnameDmg,emap,&nMVskip,1);
  freeReader(&r);
  initReader(&r, m,  filename);
  readEntities(&r, fnameDmg,emap,&nMVskip,2);
  freeReader(&r);
}


Mesh2* loadMdsFromGmsh(gmi_model* g, const char* filename)
{
  int emap[104]; //letting 0:99 (100) be the max entities until a stronger C++ coder fixes this
  Mesh2* m = makeEmptyMdsMesh(g, 0, false);
  readGmsh(m, filename,emap);
  return m;
}

Mesh2* loadMdsDmgFromGmsh(const char*fnameDmg, const char* filename)
{
  int emap[104]; //letting 99 be the max entities until a stronger C++ coder fixes this
                 // and 100+dim holds the number of model entities of a given dim
  gmshFindDmg(fnameDmg, filename, emap);  // new function that scans $Entities and writes a dmg 
  Mesh2* m = makeEmptyMdsMesh(gmi_load(fnameDmg), 0, false);
  readGmsh(m, filename,emap);
  return m;
}

}
