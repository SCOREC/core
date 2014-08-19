#ifndef PH_LINKS_H
#define PH_LINKS_H

#include <vector>
#include <map>
#include <apfNumbering.h>

namespace ph {

struct LinkKey
{
  LinkKey(int s, int p):send(s),peer(p) {}
  int send;
  int peer;
  bool operator<(LinkKey const& other) const;
};

typedef std::vector<apf::MeshEntity*> Link;

typedef std::map<LinkKey, Link> Links;

void getVertexLinks(apf::Mesh* m, Links& links);

void encodeLinks(apf::Numbering* n, Links& links, size_t& size, int*& a);

}

#endif
