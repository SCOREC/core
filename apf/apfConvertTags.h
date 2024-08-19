#ifndef APF_CONVERT_TAGS_H
#define APF_CONVERT_TAGS_H

#include "apf.h"
#include "apfConvert.h"

namespace {
  static apf::Gid getMax(const apf::GlobalToVert& globalToVert, pcu::PCU *pcu_obj)
  {
    apf::Gid max = -1;
    APF_CONST_ITERATE(apf::GlobalToVert, globalToVert, it)
      max = std::max(max, it->first);
    return pcu_obj->Max<long>(max); // this is type-dependent
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
  apf::MeshTag* createTag<long>(apf::Mesh* m,
      const char* name, const int entries) {
    return m->createLongTag(name,entries);
  }
  
  template <> inline
  apf::MeshTag* createTag<double>(apf::Mesh* m,
      const char* name, const int entries) {
    return m->createDoubleTag(name,entries);
  }
  
  inline void setEntTag(apf::Mesh* m, apf::MeshTag* t,
      apf::MeshEntity* e, long* vals) {
    return m->setLongTag(e,t,vals);
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
  apf::Gid max = getMax(globalToVert, m->getPCU());
  apf::Gid total = max + 1;
  int peers = m->getPCU()->Peers();
  int quotient = total / peers;
  int remainder = total % peers;
  int mySize = quotient;
  int self = m->getPCU()->Self();
  if (self == (peers - 1))
    mySize += remainder;
  apf::Gid myOffset = (long)self * quotient;

  /* Force each peer to have exactly mySize verts.
     This means we might need to send and recv some coords */
  T* c = new T[mySize*entries];

  apf::Gid start = m->getPCU()->Exscan(nverts);

  m->getPCU()->Begin();
  apf::Gid tmpL=start / quotient;
  int tmpI=tmpL;
  int to = std::min(peers - 1, tmpI);
  tmpL=(to+1)*(long)quotient-start;
  tmpI=tmpL;
  int n = std::min(tmpI, nverts);
  while (nverts > 0) {
    m->getPCU()->Pack(to, start);
    m->getPCU()->Pack(to, n);
    m->getPCU()->Pack(to, vals, n*entries*sizeof(T));

    nverts -= n;
    start += n;
    vals += n*entries;
    to = std::min(peers - 1, to + 1);
    n = std::min(quotient, nverts);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    m->getPCU()->Unpack(start);
    m->getPCU()->Unpack(n);
    m->getPCU()->Unpack(&c[(start - myOffset) * entries], n*entries*sizeof(T));
  }

  /* Tell all the owners of the data what we need */
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  m->getPCU()->Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) {
    apf::Gid gid = it->first;
    tmpL=gid / quotient;
    tmpI=tmpL;
    int to = std::min(peers - 1, tmpI);
// replaced     int to = std::min(peers - 1, gid / quotient);
    m->getPCU()->Pack(to, gid);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    apf::Gid gid;
    m->getPCU()->Unpack(gid);
    int from = m->getPCU()->Sender();
    tmpParts.at(gid - myOffset).push_back(from);
  }

  /* Send the data to everybody who want them */
  m->getPCU()->Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      apf::Gid gid = i + myOffset;
      m->getPCU()->Pack(to, gid);
      m->getPCU()->Pack(to, &c[i*entries], entries*sizeof(T));
    }
  }
  m->getPCU()->Send();
  apf::MeshTag* t = createTag<T>(m,tagName,entries);
  T* v = new T[entries];
  while (m->getPCU()->Receive()) {
    apf::Gid gid;
    m->getPCU()->Unpack(gid);
    m->getPCU()->Unpack(v, entries*sizeof(T));
    setEntTag(m,t,globalToVert[gid],v);
  }
  delete [] v;
  delete [] c;
  return t;
}

}
#endif
