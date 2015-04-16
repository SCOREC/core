#include "dspSmoothers.h"
#include <math.h>
#include <apf.h>

namespace dsp {
  
  Smoother::~Smoother()
  {
  }
  
  class LaplacianSmoother : public Smoother {
  public:
    void smooth(apf::Field* df, Boundary& fixed, Boundary& moving)
    {
      apf::Mesh* m = apf::getMesh(df);
      /* start Fan's code */
      apf::MeshIterator* it;
      apf::MeshEntity* v;
      apf::ModelEntity* me;
      apf::Vector3 p;
      
      //---------------------------------------------------------
      int n_mb = 0; int n_in = 0; int n_fb = 0;
      
      //iterate vertex to count the number of each type of vertex
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        me = m->toModel(v);
        if (moving.count(me)) n_mb++;
        else if (fixed.count(me)) n_fb++;
        else n_in++;
      }
      m->end(it);
      
      //----------------------------------------------------------
      int MB_begin = 0;            int MB_end = n_mb - 1;
      int IN_begin = n_mb;         int IN_end = n_mb + n_in - 1;
      int FB_begin = n_mb + n_in;  int FB_end = n_mb + n_in + n_fb - 1;
      int ithMB = 0; int ithIN =0; int ithFB = 0;
      
      vector < apf::MeshEntity* > V_total(n_mb + n_in + n_fb);
      vector < apf::Vector3 > P_total(n_mb + n_in + n_fb);
      apf::Numbering* numbers = numberOwnedDimension(m, "my_numbering", 0);
      
      //iterate to store vertices and points
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        m->getPoint(v, 0, p);
        me = m->toModel(v);
        if (moving.count(me)) {
          V_total[MB_begin + ithMB] = v;
          P_total[MB_begin + ithMB] = p;
          apf::number(numbers, v, 0, 0, MB_begin + ithMB);
          ithMB++;
        }
        else if (fixed.count(me)) {
          V_total[FB_begin + ithFB] = v;
          P_total[FB_begin + ithFB] = p;
          apf::number(numbers, v, 0, 0, FB_begin + ithFB);
          ithFB++;
        }
        else {
          V_total[IN_begin + ithIN] = v;
          P_total[IN_begin + ithIN] = p;
          // apf::number(numbers, v, 0, 0, IN_begin + ithIN);
          ithIN++;
        }
      }
      m->end(it);
      
      //----------------------------------------------------------
      ithIN = 0;
      apf::Adjacent adj;
      
      //make a queue and put all MB vertex in it
      std::queue < apf::MeshEntity* > q;
      for (int i = MB_begin ; i < MB_end + 1 ; i++) {
        q.push(V_total[i]);
      }
      
      //tag = 1, indicates this is in queue before AND this is a interior vertex
      apf::MeshTag* in_queue_tag = m->createIntTag("In_Queue_Tag", 1);
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        m->setIntTag(v, in_queue_tag, 0);
      }
      m->end(it);
      
      //find its adj, if it is never in queue, put it in queue
      while (!q.empty()) {
        v = q.front();
        me = m->toModel(v);
        m->getPoint(v, 0, p);
        m->getAdjacent(v, 1, adj);
        int num_adj = adj.getSize();
        for (int i = 0 ; i < num_adj ; i++) {
          apf::MeshEntity* adj_v = apf::getEdgeVertOppositeVert(m, adj[i], v);
          apf::ModelEntity* adj_v_me = m->toModel(adj_v);
          if ((!moving.count(adj_v_me)) & (!fixed.count(adj_v_me))) {
            int tag;
            m->getIntTag(adj_v, in_queue_tag, tag);
            if (tag == 0) {
              q.push(adj_v);
              m->setIntTag(adj_v, in_queue_tag, 1);
            }
          }
        }
        if ((!moving.count(me)) & (!fixed.count(me))) {
          V_total[IN_begin + ithIN] = v;
          P_total[IN_begin + ithIN] = p;
          apf::number(numbers, v, 0, 0, IN_begin + ithIN);
          ithIN++;
        }
        q.pop();
      }
      
      //----------------------------------------------------------
      double tol = 1.0E-5; //tolerance
      vector < apf::Vector3 > delta_P(n_mb + n_in + n_fb);
      apf::Vector3 P_temp = apf::Vector3(0, 0, 0);
      
      for (int i = 0 ; i < n_mb + n_in + n_fb ; i++) {
        delta_P[i] = apf::Vector3(0, 0, 0);
      }
      
      // average nodal position = sum(all adj_V's position)/num of adj_V
      // check max, stop until it is less the tolerance
      double max = 1.0;
      while (max > tol) {
        max = 0.0;
        for (int i = IN_begin ; i < IN_end + 1 ; i++) {
          m->getAdjacent(V_total[i], 1, adj);
          int num_adj = adj.getSize();
          for (int j = 0 ; j < num_adj ; j++) {
            apf::MeshEntity* adj_v = apf::getEdgeVertOppositeVert(m, adj[j], V_total[i]);
            int adj_v_id = apf::getNumber(numbers, adj_v, 0, 0);
            P_temp = P_temp + P_total[adj_v_id];
          }
          P_temp[0] = P_temp[0]/num_adj;
          P_temp[1] = P_temp[1]/num_adj;
          P_temp[2] = P_temp[2]/num_adj;
          delta_P[i] = P_temp - P_total[i];
          
          double temp_max = sqrt(pow(delta_P[i][0], 2) + pow(delta_P[i][1], 2) + pow(delta_P[i][2], 2));
          if (max < temp_max) {
            max = temp_max;
          }
          //update IN
          P_total[i] = P_temp;
          mesh->setPoint(V_total[i], 0, P_total[i]);
        }
      }
      
      /* end Fan's code */
      (void)m;
      (void)df;
      (void)fixed;
      (void)moving;
    }
  };
  
  class EmptySmoother : public Smoother {
  public:
    void smooth(apf::Field* df, Boundary& fixed, Boundary& moving)
    {
      (void)df;
      (void)fixed;
      (void)moving;
    }
  };
  
  Smoother* Smoother::makeLaplacian()
  {
    return new LaplacianSmoother();
  }
  
  Smoother* Smoother::makeEmpty()
  {
    return new EmptySmoother();
  }
  
}
