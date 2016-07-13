#ifndef MA_CRAWLER_H
#define MA_CRAWLER_H

#include <vector>
#include "maAdapt.h"

namespace ma {

struct Crawler
{
  Crawler(Mesh* m):mesh(m) {}
  virtual ~Crawler() {}
  typedef std::vector<Entity*> Layer;
  virtual void begin(Layer& first) = 0;
  virtual void end() {}
  virtual Entity* crawl(Entity* e) = 0;
  virtual void send(Entity* e, int to) = 0;
  virtual bool recv(Entity* e, int from) = 0;
  Mesh* mesh;
};

void crawlLayers(Crawler* c);
void crawlLayer(Crawler* c, Crawler::Layer& layer);
void syncLayer(Crawler* c, Crawler::Layer& layer);
void getDimensionBase(Adapt* a, int d, Crawler::Layer& base);
Entity* getOtherVert(Mesh* m, Entity* v, Predicate& visited);
Entity* getOtherEdge(Mesh* m, Entity* e, Predicate& visited);

void flagLayerTop(Adapt* a);

}

#endif
