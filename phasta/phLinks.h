#ifndef PH_LINKS_H
#define PH_LINKS_H

#include <vector>
#include <map>
#include <apfNumbering.h>
#include "apfZoltan.h"

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

void getLinks(apf::Mesh* m, int dim, Links& links);

void encodeILWORK(apf::Numbering* n, Links& links, int& size, int*& a);

void encodeILWORKF(apf::Numbering* n, Links& links, int& size, int*& a);

int* formIENNEIGH(apf::Numbering* ln);

}

#endif
