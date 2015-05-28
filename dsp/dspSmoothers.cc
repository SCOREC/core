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
  
  void Smoother::preprocess(apf::Mesh* m, Boundary& fixed, Boundary& moving, vector < apf::MeshEntity* >& V_total, int& in_0, int& fb_0)
  {
    /* start Fan's code */
    //data structure
    apf::MeshIterator* it;
    apf::MeshEntity* v;
    apf::ModelEntity* me;
    apf::Vector3 d;
    apf::Adjacent adj;
    apf::Numbering* numbers;
    apf::MeshTag* in_queue_tag;
    int mb_0 = 0; in_0 = 0; fb_0 = 0;
    numbers = apf::createNumbering(m, "my_numbers", m->getShape(), 1);
    
    //iterate vertex to count the number of each type of vertex
    int id = 0;
    it = m->begin(0);
    while ((v = m->iterate(it))) {
      me = m->toModel(v);
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
    
    //reordering
    //data structure
    //tag = 1, indicates this is in queue before AND this is a interior vertex
    in_queue_tag = m->createIntTag("In_Queue_Tag", 1);
    int zero = 0; int one = 1; int tag;
    
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
    m->destroyTag(in_queue_tag);
  }
  
  void Smoother::cleanup(apf::Mesh* m)
  {
    //apf::Numbering* numbers = m->findNumbering("my_numbers");
    //apf::destroyNumbering(numbers);

    apf::Field* qfield = m->findField("quality");
    apf::destroyField(qfield);
  }
  
  class LaplacianSmoother : public Smoother {
  public:
    void smooth(apf::Field* df, vector < apf::MeshEntity* >& V_total, int in_0, int fb_0)
    {
      apf::Mesh* m = apf::getMesh(df);
      /* start Fan's code */
      //----------------------------------------------------------
      //data structure
      apf::MeshIterator* it;
      apf::MeshEntity* v;
      apf::Vector3 d;
      apf::Field* qfield;
      apf::Adjacent adj;
      clock_t t;
      //apf::Numbering* numbers = m->findNumbering("my_numbers");
      
      //----------------------------------------------------------
      //update mesh
      t = clock();
      double tol = 1.0E-5; //tolerance
      apf::Vector3 D_temp = apf::Vector3(0.0, 0.0, 0.0);
      apf::Vector3 D_adj;
      apf::Vector3 D_new;
      
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
            apf::getVector(df, adj_v, 0, D_adj);
            D_temp = D_temp + D_adj;
          }
          D_temp[0] = D_temp[0]/num_adj;
          D_temp[1] = D_temp[1]/num_adj;
          D_temp[2] = D_temp[2]/num_adj;
          D_new = D_temp;
          apf::getVector(df, V_total[i], 0, D_temp);
          double d1 = D_new[0] - D_temp[0];
          double d2 = D_new[1] - D_temp[1];
          double d3 = D_new[2] - D_temp[2];
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
      
      //----------------------------------------------------------
      //print out quality of elements
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
      //----------------------------------------------------------
      /* end Fan's code */
      (void)m;
      (void)df;
    }
  };
  
  class SemiSpringSmoother : public Smoother {
  public:
    void smooth(apf::Field* df, vector < apf::MeshEntity* >& V_total, int in_0, int fb_0)
    {
      apf::Mesh* m = apf::getMesh(df);
      /* start Fan's code */
      //---------------------------------------------------------
      //data structure
      apf::MeshIterator* it;
      apf::MeshEntity* v;
      apf::Vector3 d;
      apf::Field* qfield;
      apf::Adjacent adj;
      clock_t t;
      //apf::Numbering* numbers = m->findNumbering("my_numbers");
      
      //----------------------------------------------------------
      //update mesh
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
      
      //----------------------------------------------------------
      //print out quality
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
      //----------------------------------------------------------
      /* end Fan's code */
      (void)m;
      (void)df;
    }
  };
  
  class ElasticSmoother : public Smoother {
  public:
    void smooth(apf::Field* df, vector < apf::MeshEntity* >& V_total, int in_0, int fb_0)
    {
      cout << "Elastic" << endl;
      apf::Mesh* m = apf::getMesh(df);
      /* start Fan's code */
      //---------------------------------------------------------
      //data structure
      apf::MeshIterator* it;
      apf::MeshEntity* v;
      apf::Vector3 d;
      apf::Field* qfield;
      apf::Adjacent adj;
      apf::Downward down;
      clock_t t;
      apf::Numbering* numbers = m->findNumbering("my_numbers");
      //----------------------------------------------------------
      //update mesh
      t = clock();
      double tol = 1.0E-5; //tolerance
      //----------------------------------------------------------
      //start elastic
      int num_elm = 0;
      it = m->begin(3);
      while ((v = m->iterate(it))) {
        num_elm++;
      }
      m->end(it);
      //----------------------------------------------------------
      //ien and id
      int** ien = new int*[num_elm];
      for (int i = 0 ; i < num_elm ; i++)
      {
        ien[i] = new int[4]; //only works for tet
      }
      
      int elm_id = 0;
      it = m->begin(3);
      while ((v = m->iterate(it))) {
        int numVertices = m->getDownward(v, 0, down);
        if (numVertices != 4)
          cout << "ERROR! It is not a tet. " << endl;
        for (int i = 0 ; i < numVertices ; i++) {
          int vertex_id = apf::getNumber(numbers, down[i], 0, 0);
          ien[elm_id][i] = vertex_id;
        }
        elm_id++;
      }
      m->end(it);
      
      int num_v = V_total.size();
      
      int** id = new int*[num_v];
      for (int i = 0 ; i < num_v ; i++)
      {
        id[i] = new int[3];
      }
      
      int id_temp = 1;
      for (int i = 0 ; i < in_0 ; i++) {
        id[i][0] = 0;
        id[i][1] = 0;
        id[i][2] = 0;
      }
      for (int i = in_0 ; i < fb_0 ; i++) {
        id[i][0] = id_temp++;
        id[i][1] = id_temp++;
        id[i][2] = id_temp++;
      }
      for (int i = fb_0 ; i < num_v ; i++) {
        id[i][0] = 0;
        id[i][1] = 0;
        id[i][2] = 0;
      }
      int id_max = id_temp - 1;
      //----------------------------------------------------------
      //global and local setup
      //global setup
      double** K_global = new double*[id_max];
      for (int i = 0 ; i < id_max ; i++) {
        K_global[i] = new double[id_max];
        for (int j = 0 ; j < id_max ; j++) {
          K_global[i][j] = 0.0;
        }
      }
      double* F_global = new double[id_max];
      for (int i = 0 ; i < id_max ; i++) {
        F_global[i] = 0.0;
      }
      
      //local setup
      double** K = new double*[3];
      for (int i = 0 ; i < 3 ; i++) {
        K[i] = new double[3];
        for (int j = 0 ; j < 3 ; j++)
          K[i][j] = 0.0;
      }
      double Na1,Na2,Na3,Nb1,Nb1,Nb1;
      double preDisp;
      
      double shapeFunction[4][3] = {0.0};
      shapeFunction[0][0] = -1; shapeFunction[0][1] = -1; shapeFunction[0][2] = -1;
      shapeFunction[1][0] =  1; shapeFunction[2][1] =  1; shapeFunction[3][2] =  1;

      //define physics parameter
      double E = 4.0E8;
      double nu = 0.3;
      double lamda = nu * E/((1.0 + nu) * (1.0 - 2.0 * nu));
      double mu = E/(2*(1.0 + nu));
      //----------------------------------------------------------
      //assembly
      for (int elm_id = 0 ; elm_id < num_elm; elm_id++) {
        apf::Vector3 p0, p1, p2, p3;
        m->getPoint(V_total[ien[elm_id][0]], 0, p0);
        m->getPoint(V_total[ien[elm_id][1]], 0, p1);
        m->getPoint(V_total[ien[elm_id][2]], 0, p2);
        m->getPoint(V_total[ien[elm_id][3]], 0, p3);
        
        double J11, J12, J13, J21, J22, J23, J31, J32, J33,
        double JacoDeter;
        J11 = p1[0] - p0[0]; J12 = p1[1] - p0[1]; J13 = p1[2] - p0[2];
        J21 = p2[0] - p0[0]; J22 = p2[1] - p0[1]; J23 = p2[2] - p0[2];
        J31 = p3[0] - p0[0]; J32 = p3[1] - p0[1]; J33 = p3[2] - p0[2];
        JacoDeter = 1.0/24.0 * J11 * (J22 * J33 - J23 * J32) - J12 * (J21 * J33 - J23 * J31) + J13 * (J21 * J32 - J22 * J31);

        //loop over Gaussian Quadrature points
        for (int GqPtId = 0 ; GqPtId < 4 ; GqPtId++) {
          //loop over nodes
          for (int i = 0 ; i < 3 ; i++) {
            for (int j = 0 ; j < 3 ; j++) {
              
              Nb1 = shapeFunction[i][0]; Nb2 = shapeFunction[i][1]; Nb3 = shapeFunction[i][2];
              Na1 = shapeFunction[j][0]; Na2 = shapeFunction[j][1]; Na3 = shapeFunction[j][2];
              
              K[0][0] = JacoDeter * (Na2*Nb2*mu + Na3*Nb3*mu + Na1*Nb1*(lamda + 2*mu));
              K[0][1] = JacoDeter * (Na2*Nb1*lamda + Na1*Nb2*mu);
              K[0][2] = JacoDeter * (Na3*Nb1*lamda + Na1*Nb3*mu);
              K[1][0] = JacoDeter * (Na1*Nb2*lamda + Na2*Nb1*mu);
              K[1][1] = JacoDeter * (Na1*Nb1*mu + Na3*Nb3*mu + Na2*Nb2*(lamda + 2*mu));
              K[1][2] = JacoDeter * (Na3*Nb2*lamda + Na2*Nb3*mu);
              K[2][0] = JacoDeter * (Na1*Nb3*lamda + Na3*Nb1*mu);
              K[2][1] = JacoDeter * (Na2*Nb3*lamda + Na3*Nb2*mu);
              K[2][2] = JacoDeter * (Na1*Nb1*mu + Na2*Nb2*mu + Na3*Nb3*(lamda + 2*mu));
              
              int LMi0 = id[ien[elm_id][i]][0]; int LMj0 = id[ien[elm_id][j]][0];
              int LMi1 = id[ien[elm_id][i]][1]; int LMj1 = id[ien[elm_id][j]][1];
              int LMi2 = id[ien[elm_id][i]][2]; int LMj2 = id[ien[elm_id][j]][2];

              if (LMi0 != 0) {
                if (LMj0 != 0) {
                  K_global[LMi0 - 1][LMj0 - 1] = K_global[LMi0 - 1][LMj0 - 1] + K[0][0];
                }
                else if (LMj0 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi0 - 1] = F_global[LMi0 - 1] - K[0][0] * preDisp[0];
                }
                if (LMj1 != 0) {
                  K_global[LMi0 - 1][LMj1 - 1] = K_global[LMi0 - 1][LMj1 - 1] + K[0][1];
                }
                else if (LMj1 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi0 - 1] = F_global[LMi0 - 1] - K[0][1] * preDisp[1];
                }
                if (LMj2 != 0) {
                  K_global[LMi0 - 1][LMj2 - 1] = K_global[LMi0 - 1][LMj2 - 1] + K[0][2];
                }
                else if (LMj2 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi0 - 1] = F_global[LMi0 - 1] - K[0][2] * preDisp[2];
                }
              }
              if (LMi1 != 0) {
                if (LMj0 != 0) {
                  K_global[LMi1 - 1][LMj0 - 1] = K_global[LMi1 - 1][LMj0 - 1] + K[1][0];
                }
                else if (LMj0 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi1 - 1] = F_global[LMi1 - 1] - K[1][0] * preDisp[0];
                }
                if (LMj1 != 0) {
                  K_global[LMi1 - 1][LMj1 - 1] = K_global[LMi1 - 1][LMj1 - 1] + K[1][1];
                }
                else if (LMj1 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi1 - 1] = F_global[LMi1 - 1] - K[1][1] * preDisp[1];
                }
                if (LMj2 != 0) {
                  K_global[LMi1 - 1][LMj2 - 1] = K_global[LMi1 - 1][LMj2 - 1] + K[1][2];
                }
                else if (LMj2 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi1 - 1] = F_global[LMi1 - 1] - K[1][2] * preDisp[2];
                }
              }
              if (LMi2 != 0) {
                if (LMj0 != 0) {
                  K_global[LMi2 - 1][LMj0 - 1] = K_global[LMi2 - 1][LMj0 - 1] + K[2][0];
                }
                else if (LMj0 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi2 - 1] = F_global[LMi2 - 1] - K[2][0] * preDisp[0];
                }
                if (LMj1 != 0) {
                  K_global[LMi2 - 1][LMj1 - 1] = K_global[LMi2 - 1][LMj1 - 1] + K[2][1];
                }
                else if (LMj1 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi2 - 1] = F_global[LMi2 - 1] - K[2][1] * preDisp[1];
                }
                if (LMj2 != 0) {
                  K_global[LMi2 - 1][LMj2 - 1] = K_global[LMi2 - 1][LMj2 - 1] + K[2][2];
                }
                else if (LMj2 == 0) {
                  apf::getVector(df, V_total[ien[elm_id][j]], 0, preDisp);
                  F_global[LMi2 - 1] = F_global[LMi2 - 1] - K[2][2] * preDisp[2];
                }
              }
            } //end for j
          } //end for i
        } //end for GqPtId
      } //end for elm_id
      //----------------------------------------------------------
      //solver
      //conjugate Gradient Method
      double* d_global = new double[id_max];
      for (int i = 0 ; i < id_max ; i++) {
        d_global[i] = 0.0;
      }
      
      double r[id_max];
      double d[id_max];
      double Ad[id_max];
      double alpha, beta, rTr, dTAd, rTr_new;
      for (int i = 0 ; i < id_max ; i++) {
        d[i] = F_global[i];
        r[i] = d[i];
      }
      
      for (int i = 0 ; i < id_max ; i++) {
        double r_max = *max_element(r,r+id_max)
        if (r_max <= 1e-16)
          break;
        rTr = 0.0; dTAd = 0.0; rTr_new = 0.0;
        for (int j = 0 ; j < id_max ; j++)
          rTr += r[j] * r[j];
        for (int j = 0 ; j < id_max ; j++)
          Ad[j] = 0.0;
          for (int k = 0 ; k < id_max ; k++)
            Ad[j] += K_global[j][k] * d[k];
        for (int j = 0 ; j < id_max ; j++)
          dTAd += d[j] * Ad[j];
        alpha = rTr / dTAd;
        for (int j = 0 ; j < id_max ; j++)
          d_global[j] += alpha * d[j];
        for (int j = 0 ; j < id_max ; j++)
          r[j] -= alpha * Ad[j];
        for (int j = 0 ; j < id_max ; j++)
          rTr_new += r[j] * r[j];
        beta = rTr_new / rTr;
        for (int j = 0 ; j < id_max ; j++)
          d[j] = r[j] + beta * d[j];
      }
      //----------------------------------------------------------
      //post-process
      apf::Vector3 D_new;
      for (int i = 0 ; i < num_v ; i++) {
        if (id[i][0] != 0) {
          apf::getVector(df, V_total[i], 0, D_new);
          D_new[0] = d_global[id[i][0] - 1];
          apf::setVector(df, V_total[i], 0, D_new);
        }
        if (id[i][1] != 0) {
          apf::getVector(df, V_total[i], 0, D_new);
          D_new[1] = d_global[id[i][1] - 1];
          apf::setVector(df, V_total[i], 0, D_new);
        }
        if (id[i][2] != 0) {
          apf::getVector(df, V_total[i], 0, D_new);
          D_new[2] = d_global[id[i][2] - 1];
          apf::setVector(df, V_total[i], 0, D_new);
        }
      }
      //end elastic
      //----------------------------------------------------------
      t = clock() - t;
      cout << "CPU time = " << ((float)t)/CLOCKS_PER_SEC << endl;
      
      //----------------------------------------------------------
      //print out quality
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
      //----------------------------------------------------------
      /* end Fan's code */
      (void)m;
      (void)df;
    }
  };

  class EmptySmoother : public Smoother {
  public:
    void smooth(apf::Field* df, vector < apf::MeshEntity* >& V_total, int in_0, int fb_0)
    {
      (void)df;
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
  
  Smoother* Smoother::makeElastic()
  {
    return new ElasticSmoother();
  }
  
  Smoother* Smoother::makeEmpty()
  {
    return new EmptySmoother();
  }
  
}
