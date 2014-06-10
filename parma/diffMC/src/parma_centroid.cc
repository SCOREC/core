#include <apfPartition.h>
#include <apf.h>
#include <PCU.h>

namespace parma {

class CentroidStepper
{
  public:
    typedef std::map<int,int> CountMap;
    typedef std::map<int,double> WeightMap;
    typedef std::multimap<double,apf::MeshEntity*> DistanceQueue;
    typedef std::map<int,apf::Vector3> CentroidMap;

    apf::Mesh* m;
    int dim;
    apf::MeshTag* weightTag;
    apf::MeshTag* sideTag;
    double selfWeight;
    CountMap sides;
    int totalSides;
    WeightMap peers;
    WeightMap targets;
    WeightMap sending;
    apf::Migration* migration;
    apf::Vector3 centroid;
    apf::MeshTag* sendTag;
    DistanceQueue distanceQueue;
    CentroidMap centroids;
    double magicFactor;

    CentroidStepper(apf::Mesh* m_in, double f)
    {
      m = m_in;
      dim = m->getDimension();
      magicFactor = f;
    }

    void tagSide(apf::MeshEntity* s, int peer)
    {
      m->setIntTag(s,sideTag,&peer);
    }

    void findBoundaries()
    {
      sideTag = m->createIntTag("ma_side",1);
      apf::MeshEntity* s;
      apf::MeshIterator* it = m->begin(dim-1);
      totalSides = 0;
      while ((s = m->iterate(it)))
        if (m->countUpward(s)==1)
        {
          int peer;
          if (m->isShared(s))
          {
            peer = apf::getOtherCopy(m,s).first;
            ++(sides[peer]);
            ++totalSides;
          }
          else
            peer = -1; //geometric boundary
          tagSide(s,peer);
        }
      m->end(it);
    }

    void clearSideTag()
    {
      apf::removeTagFromDimension(m,sideTag,dim-1);
      m->destroyTag(sideTag);
    }

    double getElementWeight(apf::MeshEntity* e)
    {
      double weight;
      m->getDoubleTag(e,weightTag,&weight);
      return weight;
    }

    void getSelfWeight()
    {
      apf::MeshIterator* it = m->begin(dim);
      apf::MeshEntity* e;
      double sum = 0;
      while ((e = m->iterate(it)))
        sum += getElementWeight(e);
      m->end(it);
      selfWeight = sum;
    }

    void exchangeWeights()
    {
      PCU_Comm_Begin();
      APF_ITERATE(CountMap,sides,it)
        PCU_COMM_PACK(it->first,selfWeight);
      PCU_Comm_Send();
      while (PCU_Comm_Listen())
      {
        double otherWeight;
        PCU_COMM_UNPACK(otherWeight);
        peers[PCU_Comm_Sender()] = otherWeight;
      }
    }

    void getTargets()
    {
      APF_ITERATE(WeightMap,peers,it)
        if (it->second < selfWeight) {
          double difference = selfWeight - it->second;
          double sideFraction = sides[it->first];
          sideFraction /= totalSides;
          targets[it->first] =
            difference * sideFraction * magicFactor;
        }
    }

    int getSidePeer(apf::MeshEntity* side)
    {
      int peer;
      m->getIntTag(side,sideTag,&peer);
      return peer;
    }

    bool isSide(apf::MeshEntity* e)
    {
      return m->hasTag(e,sideTag);
    }

    apf::MeshEntity* getOtherSide(
        apf::MeshEntity* side,
        apf::MeshEntity* hinge)
    {
      apf::Up up;
      m->getUp(hinge,up);
      for (int i=0; i < up.n; ++i)
      {
        apf::MeshEntity* s = up.e[i];
        if ((s != side)&&(isSide(s)))
          return s;
      }
      return 0;
    }

    void getCentroid()
    {
      apf::Vector3 x(0,0,0);
      apf::MeshEntity* e;
      apf::MeshIterator* it = m->begin(dim);
      while ((e = m->iterate(it)))
        x = x + (apf::getLinearCentroid(m, e) * getElementWeight(e));
      m->end(it);
      centroid = x / selfWeight;
    }

    void tagSend(apf::MeshEntity* e, int to)
    {
      m->setIntTag(e,sendTag,&to);
    }

    int getTaggedSend(apf::MeshEntity* e)
    {
      int to;
      m->getIntTag(e,sendTag,&to);
      return to;
    }

    apf::MeshEntity* popDistanceQueue()
    {
      DistanceQueue::iterator it = distanceQueue.begin();
      apf::MeshEntity* e = it->second;
      distanceQueue.erase(it);
      return e;
    }

    void clearSendTag()
    {
      apf::removeTagFromDimension(m,sendTag,dim);
      m->destroyTag(sendTag);
    }

    apf::MeshEntity* getOtherElement(apf::MeshEntity* e, apf::MeshEntity* s)
    {
      apf::Up up;
      m->getUp(s,up);
      for (int i=0; i < up.n; ++i)
      {
        apf::MeshEntity* o = up.e[i];
        if (o != e)
          return o;
      }
      return 0;
    }

    int getAdjacentElements(apf::MeshEntity* e, apf::Downward oe)
    {
      int n = 0;
      apf::Downward s;
      int ns = m->getDownward(e,dim-1,s);
      for (int i=0; i < ns; ++i)
      {
        oe[n] = getOtherElement(e,s[i]);
        if (oe[n]) ++n;
      }
      return n;
    }

    void exchangeCentroids()
    {
      PCU_Comm_Begin();
      APF_ITERATE(WeightMap,peers,it)
        PCU_COMM_PACK(it->first,centroid);
      PCU_Comm_Send();
      while (PCU_Comm_Listen())
      {
        apf::Vector3 otherCentroid;
        PCU_COMM_UNPACK(otherCentroid);
        centroids[PCU_Comm_Sender()] = otherCentroid;
      }
    }

    double getDistanceTo(apf::MeshEntity* e, int to)
    {
      return (apf::getLinearCentroid(m, e) - centroids[to]).getLength();
    }

    void pushDistanceQueue(apf::MeshEntity* e, int to)
    {
      if (m->hasTag(e,sendTag))
        return;
      double d = getDistanceTo(e,to);
      tagSend(e,to);
      distanceQueue.insert(std::make_pair(d,e));
    }

    void initDistanceQueue()
    {
      apf::MeshEntity* s;
      sendTag = m->createIntTag("ma_send",1);
      apf::MeshIterator* it = m->begin(dim-1);
      typedef std::map<int,std::pair<double,apf::MeshEntity*> > ClosestMap;
      ClosestMap closest;
      while ((s = m->iterate(it)))
        if (isSide(s))
        {
          int peer = getSidePeer(s);
          if ( ! targets.count(peer))
            continue;
          apf::MeshEntity* e = m->getUpward(s,0);
          double d = getDistanceTo(e,peer);
          if (( ! closest.count(peer)) ||
              (d < closest[peer].first))
            closest[peer] = std::make_pair(d,e);
        }
      m->end(it);
      assert(closest.size() == targets.size());
      APF_ITERATE(ClosestMap,closest,it)
        pushDistanceQueue(it->second.second,it->first);
    }

    void trySending(apf::MeshEntity* e)
    {
      assert( ! migration->has(e));
      int peer = getTaggedSend(e);
      if (sending[peer] >= targets[peer])
        return;
      migration->send(e,peer);
      sending[peer] += getElementWeight(e);
      apf::Downward oe;
      int noe = getAdjacentElements(e,oe);
      for (int i=0; i < noe; ++i)
        pushDistanceQueue(oe[i],peer);
    }

    void migrate()
    {
      migration = new apf::Migration(m);
      while ( ! distanceQueue.empty())
        trySending(popDistanceQueue());
      clearSideTag();
      clearSendTag();
      m->migrate(migration);
    }

    double getImbalance()
    {
      double totalWeight = selfWeight;
      PCU_Add_Doubles(&totalWeight,1);
      double averageWeight = totalWeight / PCU_Comm_Peers();
      double maxWeight = selfWeight;
      PCU_Max_Doubles(&maxWeight, 1);
      return maxWeight / averageWeight;
    }

/* send to other centroid strategy.
   great for the most part, only issue is draining
   part with many boundaries with different targets.
   mostly fixed that by only initializing the queue with
   the one closest element */
    bool run(apf::MeshTag* w, double imbalance)
    {
      weightTag = w;
      getSelfWeight();
      if (getImbalance() < imbalance)
        return false;
      findBoundaries();
      exchangeWeights();
      getTargets();
      getCentroid();
      exchangeCentroids();
      initDistanceQueue();
      migrate();
      return true;
    }
    
};

/* Stepper above and Diffuser could be one
   class. You'd just have to be careful to
   destroy certain structures between steps.
   This is the lazy solution */
class CentroidDiffuser : public apf::Balancer
{
  public:
    CentroidDiffuser(apf::Mesh* m, double f)
    {
      mesh = m;
      factor = f;
    }
    bool runStep(apf::MeshTag* weights, double tolerance)
    {
      CentroidStepper d(mesh,factor);
      return d.run(weights,tolerance);
    }
    virtual void balance(apf::MeshTag* weights, double tolerance)
    {
      double t0 = MPI_Wtime();
      while (runStep(weights,tolerance));
      double t1 = MPI_Wtime();
      if (!PCU_Comm_Self())
        printf("centroid balanced to %f in %f seconds\n", tolerance, t1-t0);
    }
  private:
    apf::Mesh* mesh;
    double factor;
};

}

apf::Balancer* Parma_MakeCentroidDiffuser(apf::Mesh* m, double stepFactor)
{
  return new parma::CentroidDiffuser(m,stepFactor);
}
