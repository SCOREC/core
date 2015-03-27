#include <PCU.h>
#include "phLinks.h"
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
   will store an array of all vertices it has
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

void getVertexLinks(apf::Numbering* n, Links& links)
{
  apf::Mesh* m = apf::getMesh(n);
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  PhastaSharing shr(m);
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
      LinkPair p;
      p.local = apf::getNumber(n, v, 0, 0);
      /* in matching we may accumulate multiple occurrences
         of the same master in the outgoing links array
         to a part that contains multiple copies of it. */
      links[LinkKey(1, remotes[i].peer)].push_back(p);
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
      LinkPair p;
      p.local = apf::getNumber(n, v, 0, 0);
      links[LinkKey(0, peer)].push_back(p);
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

void encodeILWORK(Links& links, int& size, int*& a)
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
      a[i++] = lit->local + 1;
      a[i++] = 1;
    }
  }
  assert(i == size);
}

void encodeILWORKF(Links& links, int& size, int*& a)
{
  size = 1; // total links
  APF_ITERATE(Links, links, it) {
    size += 2; //link header
    size += it->second.size() * 2; //entity entries
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
      a[i++] = lit->local + 1;
      a[i++] = lit->remote + 1;
    }
  }
  assert(i == size);
}

void separateElementGraph(apf::Mesh* m, apf::LocalCopy* e2e,
    Links& links, int*& ienneigh)
{
  int dim = m->getDimension();
  int type = getFirstType(m, dim);
  int nsides = apf::Mesh::adjacentCount[type][dim - 1];
  size_t nelem = m->count(dim);
  int self = PCU_Comm_Self();
  ienneigh = new int[nelem * nsides];
  for (size_t i = 0; i < nelem; ++i)
    for (int j = 0; j < nsides; ++j) {
      apf::LocalCopy& lc = e2e[i * nsides + j];
      if (lc.isNull())
        ienneigh[j * nelem + i] = 0;
      else if (lc.peer == self)
        ienneigh[j * nelem + i] = lc.localNumber;
      else {
        ienneigh[j * nelem + i] = 0;
        LinkKey k(0, lc.peer);
        LinkPair p;
        p.local = i;
        p.remote = lc.localNumber;
        /* note: this does not match the order of the
           tasks on the two ranks the way we did for
           the nodes... */
        links[k].push_back(p);
      }
    }
}

}
