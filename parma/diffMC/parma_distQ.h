#ifndef PARMA_DISTQ_H_
#define PARMA_DISTQ_H_

#include <apfMesh.h>
#include <map>

namespace parma {
  struct Greater {
    bool operator() (const int& l, const int& r) const {
      return (l > r);
    }
  };
  struct Less {
    bool operator() (const int& l, const int& r) const {
      return (l < r);
    }
  };

  template <class Compare> class DistanceQueue {
    typedef typename std::multimap<int, apf::MeshEntity*, Compare> DistanceQ;
    typedef typename DistanceQ::iterator DistanceQIter;

    public:
    DistanceQueue(apf::Mesh* mesh) : m(mesh) 
    {
      t = m->createIntTag("parmaDistanceQueue",1);
    }

    ~DistanceQueue()
    {
      apf::removeTagFromDimension(m,t,0);
      m->destroyTag(t);
    }

    void push(apf::MeshEntity* e, int dist)
    {
      DistanceQIter it = q.begin();
      if ( m->hasTag(e, t) ) {
        it = erase(dist, e);
      }
      int one = 1;
      m->setIntTag(e, t, &one);
      q.insert(it, std::make_pair(dist, e));
    }

    apf::MeshEntity* pop()
    {
      DistanceQIter it;
      it = q.begin();
      apf::MeshEntity* e = it->second;
      q.erase(it);
      return e;
    }

    bool empty() 
    {
      return q.empty();
    }

    private:
    apf::Mesh* m;
    apf::MeshTag* t;
    DistanceQ q;

    DistanceQIter erase(int dist, apf::MeshEntity* e) 
    {
      assert( m->hasTag(e, t) );
      DistanceQIter it = q.find(dist);
      DistanceQIter rit = ( it != q.begin() ) ? it-- : q.end();
      while( it != q.end() ) {
        if( it->second == e ) {
          q.erase(it);
          break;
        }
        rit = it;
        it++;
      }
      return rit;
    }
  };
}

#endif

/*  Angry Programmer:
    ~~  ~~ 
    {. {.      ____________________________
      /       /                            |
    ======   /  Arghhh... why templates!?  |
    |WWWW|  <______________________________|
    ======
*/  

