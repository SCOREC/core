#ifndef APF_CONVERT_TAGS_H
#define APF_CONVERT_TAGS_H

#include <PCU.h>
#include "apf.h"
#include "apfConvert.h"

namespace {
  typedef int Gid;
  static Gid getMax(const apf::GlobalToVert& globalToVert)
  {
    Gid max = -1;
    APF_CONST_ITERATE(apf::GlobalToVert, globalToVert, it)
      max = std::max(max, it->first);
    return PCU_Max_Int(max); // this is type-dependent
  }

  template <class T> inline
  apf::MeshTag* createTag(apf::Mesh*,
      const char*, const int) {
    exit(EXIT_FAILURE);
    return 0;
  }
  
  template <> inline
  apf::MeshTag* createTag<int>(apf::Mesh* m,
      const char* name, const int entries) {
    return m->createIntTag(name,entries);
  }
  
  template <> inline
  apf::MeshTag* createTag<double>(apf::Mesh* m,
      const char* name, const int entries) {
    return m->createDoubleTag(name,entries);
  }
  
  inline void setEntTag(apf::Mesh* m, apf::MeshTag* t,
      apf::MeshEntity* e, int* vals) {
    return m->setIntTag(e,t,vals);
  }
  
  inline void setEntTag(apf::Mesh* m, apf::MeshTag* t,
      apf::MeshEntity* e, double* vals) {
    return m->setDoubleTag(e,t,vals);
  }
}

namespace apf {
/** \brief Assign tag values to the mesh
  * \param m (in) mesh
  * \param tagName (in) name of returned tag
  * \param vals (in) T array of length nverts 
  * \param entries (in) number of values per vertex
  * \param nverts (in) number of vertices for this process
  * \param globalToVert (in) map from global mesh vertex ids
  *                          to local vertex * pointers
  * \details
  * See 'setCoords' for distribution details
  */
template<class T> 
apf::MeshTag* setMappedTag(Mesh2* m, const char* tagName,
    const T* vals, const int entries,
    int nverts, GlobalToVert& globalToVert)
{
  Gid max = getMax(globalToVert);
  Gid total = max + 1;
  int peers = PCU_Comm_Peers();
  int quotient = total / peers;
  int remainder = total % peers;
  int mySize = quotient;
  int self = PCU_Comm_Self();
  if (self == (peers - 1))
    mySize += remainder;
  int myOffset = self * quotient;

  /* Force each peer to have exactly mySize verts.
     This means we might need to send and recv some coords */
  T* c = new T[mySize*entries];

  int start = PCU_Exscan_Int(nverts);

  PCU_Comm_Begin();
  int to = std::min(peers - 1, start / quotient);
  int n = std::min((to+1)*quotient-start, nverts);
  while (nverts > 0) {
    PCU_COMM_PACK(to, start);
    PCU_COMM_PACK(to, n);
    PCU_Comm_Pack(to, vals, n*entries*sizeof(T));

    nverts -= n;
    start += n;
    vals += n*entries;
    to = std::min(peers - 1, to + 1);
    n = std::min(quotient, nverts);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    PCU_COMM_UNPACK(start);
    PCU_COMM_UNPACK(n);
    PCU_Comm_Unpack(&c[(start - myOffset) * entries], n*entries*sizeof(T));
  }

  /* Tell all the owners of the data what we need */
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  PCU_Comm_Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) {
    int gid = it->first;
    int to = std::min(peers - 1, gid / quotient);
    PCU_COMM_PACK(to, gid);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    int gid;
    PCU_COMM_UNPACK(gid);
    int from = PCU_Comm_Sender();
    tmpParts.at(gid - myOffset).push_back(from);
  }

  /* Send the data to everybody who want them */
  PCU_Comm_Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      int gid = i + myOffset;
      PCU_COMM_PACK(to, gid);
      PCU_Comm_Pack(to, &c[i*entries], entries*sizeof(T));
    }
  }
  PCU_Comm_Send();
  apf::MeshTag* t = createTag<T>(m,tagName,entries);
  T* v = new T[entries];
  while (PCU_Comm_Receive()) {
    int gid;
    PCU_COMM_UNPACK(gid);
    PCU_Comm_Unpack(v, entries*sizeof(T));
    setEntTag(m,t,globalToVert[gid],v);
  }
  delete [] v;
  delete [] c;
  return t;
}

}
#endif
