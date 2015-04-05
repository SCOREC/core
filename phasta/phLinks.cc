#include <PCU.h>
#include "phLinks.h"
#include "phAdjacent.h"
#include <apf.h>

namespace ph {

bool LinkKey::operator<(LinkKey const& other) const
{
  if (send != other.send)
    return send;
  return peer < other.peer;
}

/* the PhastaSharing class is responsible for ensuring that
   ILWORK links matched entities correctly.
   Only on-part owners are treated as remotes of one another,
   where on-part ownership is defined by smallest entity pointer */

struct PhastaSharing : public apf::Sharing {
  PhastaSharing(apf::Mesh* m)
  {
    mesh = m;
    helper = apf::getSharing(m);
  }
  ~PhastaSharing()
  {
    delete helper;
  }
  bool isOwned(apf::MeshEntity* e)
  {
    return helper->isOwned(e);
  }
  /* this will only be called for global masters */
  void getCopies(apf::MeshEntity* e,
      apf::CopyArray& copies)
  {
    helper->getCopies(e, copies);
    if ( ! mesh->hasMatching())
      return;
    /* filter out matches which are on the same part as the global master */
    int self = PCU_Comm_Self();
    size_t i = 0;
    for (size_t j = 0; j < copies.getSize(); ++j)
      if (copies[j].peer != self)
        copies[i++] = copies[j];
    copies.setSize(i);
  }
  apf::Mesh* mesh;
  apf::Sharing* helper;
};

/* this algorithm is essential to parallel
   scalability: generate local inter-part 
   communication arrays describing shared entities.

   for every partition model face, each part
   will store an array of all entities it has
   classified on that partition model face.
   the two arrays are aligned such that
   the same index in both arrays denotes the same vertex.

   phParAdapt had an implementation which 
   used arrays of size equal to the number of processors.

   This version is copied from code in MDS
   (and/or PUMI) that does the exact same thing
   to generate SMB files, but uses space and
   time proportional only to the neighborhood size,
   thanks in part to PCU algorithms

   In a perfect world, we should only have one copy
   of this code.
*/

void getLinks(apf::Mesh* m, int dim, Links& links)
{
  PhastaSharing shr(m);
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
/* the alignment is such that the owner part's
   array follows the order of its vertex iterator
   traversal. The owner dictates the order to the
   other part by sending remote copies */
    if ( ! shr.isOwned(v))
      continue;
    apf::CopyArray remotes;
    shr.getCopies(v, remotes);
    for (size_t i = 0; i < remotes.getSize(); ++i) {
      /* in matching we may accumulate multiple occurrences
         of the same master in the outgoing links array
         to a part that contains multiple copies of it. */
      links[LinkKey(1, remotes[i].peer)].push_back(v);
      PCU_COMM_PACK(remotes[i].peer, remotes[i].entity);
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen()) {
    int peer = PCU_Comm_Sender();
    while (!PCU_Comm_Unpacked()) {
      apf::MeshEntity* v;
      PCU_COMM_UNPACK(v);
      links[LinkKey(0, peer)].push_back(v);
    }
  }
}

/* encode the local links into a big array of integers
   per phParAdapt's ILWORK array format.

   ilwork[0] = total number of links
   following that, each link is serialized as:

   tag (used to be something like (peer*peers + self), now zero)
   type (0 for send, 1 for receive)
   peer
   #entities

   then for each entity there are two integers

   local id
   #dofs on entity
*/

void encodeILWORK(apf::Numbering* n, Links& links, int& size, int*& a)
{
  size = 1; // total links
  APF_ITERATE(Links, links, it) {
    size += 4; //link header
    size += it->second.size() * 2; //entity entries
  }
  a = new int[size];
  a[0] = links.size();
  int i = 1;
  APF_ITERATE(Links, links, it) {
    LinkKey k = it->first;
    Link& l = it->second;
    a[i++] = 0;
    a[i++] = k.send;
    a[i++] = k.peer + 1; /* peers numbered from 1 */
    a[i++] = l.size();
    APF_ITERATE(Link, l, lit) {
      /* entities also numbered from 1 */
      a[i++] = apf::getNumber(n, *lit, 0, 0) + 1;
      a[i++] = 1;
    }
  }
  assert(i == size);
}

static apf::MeshEntity* getSideElement(apf::Mesh* m, apf::MeshEntity* s)
{
  return m->getUpward(s, 0);
}

void encodeILWORKF(apf::Numbering* n, Links& links, int& size, int*& a)
{
  apf::Mesh* m = apf::getMesh(n);
  size = 1; // total links
  APF_ITERATE(Links, links, it) {
    size += 2; //link header
    size += it->second.size(); //entity entries
  }
  a = new int[size];
  a[0] = links.size();
  int i = 1;
  APF_ITERATE(Links, links, it) {
    LinkKey k = it->first;
    Link& l = it->second;
    a[i++] = k.peer + 1; /* peers numbered from 1 */
    a[i++] = l.size();
    APF_ITERATE(Link, l, lit) {
      /* entities also numbered from 1 */
      apf::MeshEntity* e = getSideElement(m, *lit);
      a[i++] = apf::getNumber(n, e, 0, 0) + 1;
    }
  }
  assert(i == size);
}

static apf::MeshEntity* getOtherElem(apf::Mesh* m, apf::MeshEntity* elem,
    apf::MeshEntity* face)
{
  apf::Up up;
  m->getUp(face, up);
  if (up.n == 2)
    return up.e[1 - apf::findIn(up.e, 2, elem)];
  if (!m->hasMatching())
    return 0;
  apf::Matches matches;
  m->getMatches(face, matches);
  int self = PCU_Comm_Self();
  for (size_t i = 0; i < matches.getSize(); ++i)
    if (matches[i].peer == self)
      return m->getUpward(matches[i].entity, 0);
  return 0;
}

int* formIENNEIGH(apf::Numbering* ln)
{
  apf::Mesh* m = getMesh(ln);
  int dim = m->getDimension();
  int sideDim = dim - 1;
  int type = getFirstType(m, dim);
  int nsides = apf::Mesh::adjacentCount[type][sideDim];
  size_t nelem = m->count(dim);
  int* ienneigh = new int[nelem * nsides];
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  int i = 0;
  while ((e = m->iterate(it))) {
    assert(m->getType(e) == type);
    assert(face_apf2ph[type]);
    apf::Downward sides;
    m->getDownward(e, sideDim, sides);
    for (int j = 0; j < nsides; ++j) {
      apf::MeshEntity* oe = getOtherElem(m, e, sides[j]);
      int oj = face_apf2ph[type][j];
      int oi = oe ? (getNumber(ln, oe, 0, 0) + 1) : 0;
      ienneigh[oj * nelem + i] = oi;
    }
    ++i;
  }
  m->end(it);
  return ienneigh;
}

}
