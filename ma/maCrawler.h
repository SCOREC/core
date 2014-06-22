#ifndef MA_CRAWLER_H
#define MA_CRAWLER_H

#include <vector>
#include "maMesh.h"

namespace ma {

class Adapt;

struct Crawler
{
  Crawler(Adapt* a):adapter(a) {}
  typedef std::vector<Entity*> Layer;
  virtual void begin(Layer& first) = 0;
  virtual void end() = 0;
  virtual Entity* crawl(Entity* e) = 0;
  virtual void send(Entity* e, int to) = 0;
  virtual bool recv(Entity* e, int from) = 0;
  Adapt* adapter;
};

class QuadFlagger : public Crawler
{
  public:
    QuadFlagger(Adapt* a):Crawler(a) {}
    virtual void begin(Layer& first);
    virtual void end();
    virtual Entity* crawl(Entity* e);
    virtual void send(Entity* e, int to);
    virtual bool recv(Entity* e, int from);
};

void crawlLayers(Crawler* c);

}

#endif
