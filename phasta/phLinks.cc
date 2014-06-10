#include "phLinks.h"
#include <PCU.h>
#include <apf.h>

namespace ph {

bool LinkKey::operator<(LinkKey const& other) const
{
  if (send != other.send)
    return send;
  return peer < other.peer;
}

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

void getVertexLinks(apf::Mesh* m, Links& links)
{
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    if (!m->isShared(v))
      continue;
/* the alignement is such that the owner part's
   array follows the order of its vertex iterator
   traversal. The owner dictates the order to the
   other part by sending remote copies */
    if (!m->isOwned(v))
      continue;
    apf::Copies remotes;
    m->getRemotes(v, remotes);
    APF_ITERATE(apf::Copies, remotes, rit) {
      int peer = rit->first;
      apf::MeshEntity* c = rit->second;
      links[LinkKey(true, peer)].push_back(v);
      PCU_COMM_PACK(peer, c);
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen()) {
    int peer = PCU_Comm_Sender();
    while (!PCU_Comm_Unpacked()) {
      apf::MeshEntity* v;
      PCU_COMM_UNPACK(v);
      links[LinkKey(false, peer)].push_back(v);
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

void encodeLinks(apf::Numbering* n, Links& links, size_t& size, int*& a)
{
  size = 1; // total links
  APF_ITERATE(Links, links, it) {
    size += 4; //link header
    size += it->second.size() * 2; //entity entries
  }
  a = new int[size];
  a[0] = links.size();
  size_t i = 1;
  APF_ITERATE(Links, links, it) {
    LinkKey k = it->first;
    Link& l = it->second;
    a[i++] = 0;
    a[i++] = k.send ? 0 : 1;
    a[i++] = k.peer;
    a[i++] = l.size();
    APF_ITERATE(Link, l, lit) {
      int id = apf::getNumber(n, *lit, 0, 0);
      a[i++] = id;
      a[i++] = 1;
    }
  }
  assert(i == size);
}

}
