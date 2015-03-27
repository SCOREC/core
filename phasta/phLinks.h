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

struct LinkPair {
  int local;
  int remote;
};

typedef std::vector<LinkPair> Link;

typedef std::map<LinkKey, Link> Links;

void getVertexLinks(apf::Numbering* n, Links& links);

void encodeILWORK(Links& links, int& size, int*& a);

void encodeILWORKF(Links& links, int& size, int*& a);

void separateElementGraph(apf::Mesh* m, apf::LocalCopy* e2e,
    Links& links, int*& ienneigh);

}

#endif
