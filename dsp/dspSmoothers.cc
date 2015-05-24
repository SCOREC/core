#include "dspSmoothers.h"
#include <math.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <cmath>
#include <vector>
#include <queue>
#include <time.h>
#include <apfShape.h>

using namespace std;

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
      apf::Vector3 d;
      //erase numbering and tag
      apf::Numbering* numbers;
      apf::MeshTag* in_queue_tag;
      apf::Field* qfield;
      apf::destroyNumbering(numbers);
      m->destroyTag(in_queue_tag);
      apf::destroyField(qfield);
      
      //---------------------------------------------------------
      //data structure
      int mb_0 = 0; int in_0 = 0; int fb_0 = 0;
      vector < apf::MeshEntity* > V_total;
      vector < apf::Vector3 > D_total;
      numbers = apf::createNumbering(m, "my_numbers", m->getShape(), 1);
      
      //iterate vertex to count the number of each type of vertex
      int id = 0;
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        me = m->toModel(v);
        if (id == 0)
          if (apf::isNumbered(numbers, v, 0, 0)) {
            apf::destroyNumbering(numbers);
            m->destroyTag(in_queue_tag);
            apf::destroyField(qfield);
          }
        if (moving.count(me)) {
          V_total.push_back(v);
          apf::getVector(df, v, 0, d);
          D_total.push_back(d);
          apf::number(numbers, v, 0, 0, id);
          id++;
        }
      }
      m->end(it);
      
      in_0 = id;
      
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        me = m->toModel(v);
        if ((!fixed.count(me))&(!moving.count(me))) {
          V_total.push_back(v);
          apf::getVector(df, v, 0, d);
          D_total.push_back(d);
          apf::number(numbers, v, 0, 0, id);
          id++;
        }
      }
      m->end(it);
      
      fb_0 = id;
      
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        me = m->toModel(v);
        if (fixed.count(me)) {
          V_total.push_back(v);
          apf::getVector(df, v, 0, d);
          D_total.push_back(d);
          apf::number(numbers, v, 0, 0, id);
          id++;
        }
      }
      m->end(it);
      
      //----------------------------------------------------------
      clock_t t;
      t = clock();
      //reordering
      //data structure
      //tag = 1, indicates this is in queue before AND this is a interior vertex
      in_queue_tag = m->createIntTag("In_Queue_Tag", 1);
      int zero = 0; int one = 1; int tag;
      apf::Adjacent adj;
      
      id = in_0;
      
      //make a queue and put all MB vertex in it
      queue < apf::MeshEntity* > q;
      for (int i = mb_0 ; i < in_0 ; i++) {
        q.push(V_total[i]);
      }
      
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        m->setIntTag(v, in_queue_tag, &zero);
      }
      m->end(it);
      
      //find its adj, if it is never in queue, put it in queue
      while (!q.empty()) {
        v = q.front();
        me = m->toModel(v);
        apf::getVector(df, v, 0, d);
        m->getAdjacent(v, 1, adj);
        int num_adj = adj.getSize();
        for (int i = 0 ; i < num_adj ; i++) {
          apf::MeshEntity* adj_v = apf::getEdgeVertOppositeVert(m, adj[i], v);
          apf::ModelEntity* adj_v_me = m->toModel(adj_v);
          if ((!moving.count(adj_v_me)) & (!fixed.count(adj_v_me))) {
            m->getIntTag(adj_v, in_queue_tag, &tag);
            if (tag == 0) {
              q.push(adj_v);
              m->setIntTag(adj_v, in_queue_tag, &one);
            }
          }
        }
        if ((!moving.count(me)) & (!fixed.count(me))) {
          V_total[id] = v;
          D_total[id] = d;
          apf::number(numbers, v, 0, 0, id);
          id++;
        }
        q.pop();
      }
      t = clock() - t;
      cout << "Reordering time = " << ((float)t)/CLOCKS_PER_SEC << endl;
      
      //----------------------------------------------------------
      t = clock();
      double tol = 1.0E-5; //tolerance
      apf::Vector3 D_temp = apf::Vector3(0.0, 0.0, 0.0);
      
      // average nodal displacement = sum(all adj_V's displacement)/num of adj_V
      // check max, stop until it is less the tolerance
      double max = 1.0; int loop_times = 0;
      while (max > tol) {
        max = 0.0;
        for (int i = in_0 ; i < fb_0 ; i++) {
          m->getAdjacent(V_total[i], 1, adj);
          int num_adj = adj.getSize();
          for (int j = 0 ; j < num_adj ; j++) {
            apf::MeshEntity* adj_v = apf::getEdgeVertOppositeVert(m, adj[j], V_total[i]);
            int adj_v_id = apf::getNumber(numbers, adj_v, 0, 0);
            D_temp = D_temp + D_total[adj_v_id];
          }
          D_temp[0] = D_temp[0]/num_adj;
          D_temp[1] = D_temp[1]/num_adj;
          D_temp[2] = D_temp[2]/num_adj;
          D_total[i] = D_temp;
          apf::getVector(df, V_total[i], 0, D_temp);
          double d1 = D_total[i][0] - D_temp[0];
          double d2 = D_total[i][1] - D_temp[1];
          double d3 = D_total[i][2] - D_temp[2];
          apf::setVector(df, V_total[i], 0, D_total[i]);
          
          double temp_max = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
          if (max < temp_max) {
            max = temp_max;
          }
        }
        loop_times++;
      }
      
      t = clock() - t;
      cout << "Loop times = " << loop_times << endl;
      cout << "CPU time = " << ((float)t)/CLOCKS_PER_SEC << endl;
      
      apf::Downward down;
      double quality;
      int badTetNum = 0;
      qfield = apf::createField(m, "quality", apf::SCALAR, apf::getConstant(3));
      
      it = m->begin(3);
      while ((v = m->iterate(it))) {
        int num_down = m->getDownward(v,1,down);
        if (num_down != 6)
          cout << "WARNING! NOT A TET!" << endl;
        double l[6];
        for (int i=0; i < 6; ++i) {
          apf::MeshElement* melm = apf::createMeshElement(m,down[i]);
          l[i] = apf::measure(melm);
          apf::destroyMeshElement(melm);
        }
        apf::MeshElement* melm = apf::createMeshElement(m,v);
        double V = apf::measure(melm);
        apf::destroyMeshElement(melm);
        double s=0;
        for (int i=0; i < 6; ++i)
          s += l[i]*l[i];
        if (V < 0)
          quality = -15552.0*(V*V)/(s*s*s);
        quality = 15552.0*(V*V)/(s*s*s);
        if (quality <= 0.027)
          badTetNum++;
        apf::setScalar(qfield, v, 0, quality);
      }
      m->end(it);
      
      cout << "Number of bad tets = " << badTetNum << endl;
      
      /* end Fan's code */
      (void)m;
      (void)df;
      (void)fixed;
      (void)moving;
    }
  };
  
  class SemiSpringSmoother : public Smoother {
  public:
    void smooth(apf::Field* df, Boundary& fixed, Boundary& moving)
    {
      apf::Mesh* m = apf::getMesh(df);
      /* start Fan's code */
      apf::MeshIterator* it;
      apf::MeshEntity* v;
      apf::ModelEntity* me;
      apf::Vector3 d;
      //erase numbering and tag
      apf::Numbering* numbers;
      apf::MeshTag* in_queue_tag;
      apf::Field* qfield;
      

      //---------------------------------------------------------
      //data structure
      int mb_0 = 0; int in_0 = 0; int fb_0 = 0;
      vector < apf::MeshEntity* > V_total;
      numbers = apf::createNumbering(m, "my_numbers", m->getShape(), 1);
      
      //iterate vertex to count the number of each type of vertex
      int id = 0;
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        me = m->toModel(v);
        if (id == 0)
          if (apf::isNumbered(numbers, v, 0, 0)) {
            apf::destroyNumbering(numbers);
            m->destroyTag(in_queue_tag);
            apf::destroyField(qfield);
          }
        if (moving.count(me)) {
          V_total.push_back(v);
          apf::number(numbers, v, 0, 0, id);
          id++;
        }
      }
      m->end(it);
      
      in_0 = id;
      
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        me = m->toModel(v);
        if ((!fixed.count(me))&(!moving.count(me))) {
          V_total.push_back(v);
          apf::number(numbers, v, 0, 0, id);
          id++;
        }
      }
      m->end(it);
      
      fb_0 = id;
      
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        me = m->toModel(v);
        if (fixed.count(me)) {
          V_total.push_back(v);
          apf::number(numbers, v, 0, 0, id);
          id++;
        }
      }
      m->end(it);
      
      //----------------------------------------------------------
      clock_t t;
      t = clock();
      //reordering
      //data structure
      //tag = 1, indicates this is in queue before AND this is a interior vertex
      in_queue_tag = m->createIntTag("In_Queue_Tag", 1);
      int zero = 0; int one = 1; int tag;
      apf::Adjacent adj;
      
      id = in_0;
      
      //make a queue and put all MB vertex in it
      queue < apf::MeshEntity* > q;
      for (int i = mb_0 ; i < in_0 ; i++) {
        q.push(V_total[i]);
      }
      
      it = m->begin(0);
      while ((v = m->iterate(it))) {
        m->setIntTag(v, in_queue_tag, &zero);
      }
      m->end(it);
      
      //find its adj, if it is never in queue, put it in queue
      while (!q.empty()) {
        v = q.front();
        me = m->toModel(v);
        m->getAdjacent(v, 1, adj);
        int num_adj = adj.getSize();
        for (int i = 0 ; i < num_adj ; i++) {
          apf::MeshEntity* adj_v = apf::getEdgeVertOppositeVert(m, adj[i], v);
          apf::ModelEntity* adj_v_me = m->toModel(adj_v);
          if ((!moving.count(adj_v_me)) & (!fixed.count(adj_v_me))) {
            m->getIntTag(adj_v, in_queue_tag, &tag);
            if (tag == 0) {
              q.push(adj_v);
              m->setIntTag(adj_v, in_queue_tag, &one);
            }
          }
        }
        if ((!moving.count(me)) & (!fixed.count(me))) {
          V_total[id] = v;
          apf::number(numbers, v, 0, 0, id);
          id++;
        }
        q.pop();
      }
      t = clock() - t;
      cout << "Reordering time = " << ((float)t)/CLOCKS_PER_SEC << endl;
      
      //----------------------------------------------------------
      t = clock();
      double tol = 1.0E-5; //tolerance
      apf::Downward down;
      double stiffness_temp; apf::Vector3 force_temp;
      vector < apf::Vector3 > tet_OP(3);
      vector < apf::Vector3 > tet_DP(3);
      apf::Vector3 P_temp; //coordinate
      apf::Vector3 D_temp; //displcament
      apf::Vector3 D_new; //displcament
      double d1, d2, d3;
      
      // semi-torsional spring: stiffness = 1/length^2 + sum(1/sin(theta)^2)
      // check max, stop until it is less the tolerance
      double max = 1.0; int loop_times = 0;
      while (max > tol) {
        max = 0.0;
        for (int i = in_0 ; i < fb_0 ; i++) {
          double stiffness_sum = 0.0;
          apf::Vector3 force_sum = apf::Vector3(0.0, 0.0, 0.0);
          m->getPoint(V_total[i], 0, P_temp);
          apf::getVector(df, V_total[i], 0, D_temp);
          m->getAdjacent(V_total[i], 3, adj);
          int num_adj = adj.getSize();
          for (int j = 0 ; j < num_adj ; j++) {
            int num_down; //num_down is supposed to be 4
            num_down = m->getDownward(adj[j], 0, down);
            int tet_id = 0;
            for (int k = 0 ; k < num_down ; k++) {
              if (down[k] != V_total[i]) {
                m->getPoint(down[k], 0, tet_OP[tet_id]);
                apf::getVector(df, down[k], 0, tet_DP[tet_id]);
                tet_id++;
              }
            }
            if (tet_id != 3)
              cout << "Find more than 3 adjacnet vertices in one tet!" << endl;
            //----------------------------------------------------
            apf::Vector3 n_1; apf::Vector3 n_2;
            for (int K = 0 ; K < 3 ; K++) {
              int A; int B;
              if (K == 0)      { A = 1; B = 2; }
              else if (K == 1) { A = 2; B = 0; }
              else if (K == 2) { A = 0; B = 1; }
              
              n_1[0] = (tet_OP[B][1] - tet_OP[K][1]) * (tet_OP[A][2] - tet_OP[K][2])
              - (tet_OP[B][2] - tet_OP[K][2]) * (tet_OP[A][1] - tet_OP[K][1]);
              
              n_1[1] = (tet_OP[B][2] - tet_OP[K][2]) * (tet_OP[A][0] - tet_OP[K][0])
              - (tet_OP[B][0] - tet_OP[K][0]) * (tet_OP[A][2] - tet_OP[K][2]);
              
              n_1[2] = (tet_OP[B][0] - tet_OP[K][0]) * (tet_OP[A][1] - tet_OP[K][1])
              - (tet_OP[B][1] - tet_OP[K][1]) * (tet_OP[A][0] - tet_OP[K][0]);
              
              n_2[0] = (tet_OP[B][1] - P_temp[1]) * (tet_OP[A][2] - P_temp[2])
              - (tet_OP[B][2] - P_temp[2]) * (tet_OP[A][1] - P_temp[1]);
              
              n_2[1] = (tet_OP[B][2] - P_temp[2]) * (tet_OP[A][0] - P_temp[0])
              - (tet_OP[B][0] - P_temp[0]) * (tet_OP[A][2] - P_temp[2]);
              
              n_2[2] = (tet_OP[B][0] - P_temp[0]) * (tet_OP[A][1] - P_temp[1])
              - (tet_OP[B][1] - P_temp[1]) * (tet_OP[A][0] - P_temp[0]);
              
              double cos_squ_up = (n_1[0] * n_2[0] + n_1[1] * n_2[1] + n_1[2] * n_2[2])
              * (n_1[0] * n_2[0] + n_1[1] * n_2[1] + n_1[2] * n_2[2]);
              
              double cos_squ_dw = (n_1[0] * n_1[0] + n_1[1] * n_1[1] + n_1[2] * n_1[2])
              * (n_2[0] * n_2[0] + n_2[1] * n_2[1] + n_2[2] * n_2[2]);
              
              double cos_squ = cos_squ_up / cos_squ_dw;
              
              d1 = tet_OP[K][0] - P_temp[0];
              d2 = tet_OP[K][1] - P_temp[1];
              d3 = tet_OP[K][2] - P_temp[2];
              
              double length_squ = d1 * d1 + d2 * d2 + d3 * d3;
              
              stiffness_temp = 1/(1 - cos_squ) + 1/8/length_squ;
              force_temp[0] = stiffness_temp * tet_DP[K][0];
              force_temp[1] = stiffness_temp * tet_DP[K][1];
              force_temp[2] = stiffness_temp * tet_DP[K][2];
              
              stiffness_sum = stiffness_sum + stiffness_temp;
              force_sum = force_sum + force_temp;
            }
            //----------------------------------------------------
          }
          D_new[0] = force_sum[0]/stiffness_sum;
          D_new[1] = force_sum[1]/stiffness_sum;
          D_new[2] = force_sum[2]/stiffness_sum;
          
          d1 = D_new[0] - D_temp[0];
          d2 = D_new[1] - D_temp[1];
          d3 = D_new[2] - D_temp[2];
          apf::setVector(df, V_total[i], 0, D_new);
          
          double temp_max = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
          if (max < temp_max) {
            max = temp_max;
          }
        }
        loop_times++;
      }
      
      t = clock() - t;
      cout << "Loop times = " << loop_times << endl;
      cout << "CPU time = " << ((float)t)/CLOCKS_PER_SEC << endl;
      
      //apf::Downward down;
      double quality;
      int badTetNum = 0;
      qfield = apf::createField(m, "quality", apf::SCALAR, apf::getConstant(3));

      it = m->begin(3);
      while ((v = m->iterate(it))) {
        int num_down = m->getDownward(v,1,down);
        if (num_down != 6)
          cout << "WARNING! NOT A TET!" << endl;
        double l[6];
        for (int i=0; i < 6; ++i) {
          apf::MeshElement* melm = apf::createMeshElement(m,down[i]);
          l[i] = apf::measure(melm);
          apf::destroyMeshElement(melm);
        }
        apf::MeshElement* melm = apf::createMeshElement(m,v);
        double V = apf::measure(melm);
        apf::destroyMeshElement(melm);
        double s=0;
        for (int i=0; i < 6; ++i)
          s += l[i]*l[i];
        if (V < 0)
          quality = -15552.0*(V*V)/(s*s*s);
        quality = 15552.0*(V*V)/(s*s*s);
        if (quality <= 0.027)
          badTetNum++;
        apf::setScalar(qfield, v, 0, quality);
      }
      m->end(it);
      
      cout << "Number of bad tets = " << badTetNum << endl;
      
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
  
  Smoother* Smoother::makeSemiSpring()
  {
    return new SemiSpringSmoother();
  }
  
  Smoother* Smoother::makeEmpty()
  {
    return new EmptySmoother();
  }
  
}
