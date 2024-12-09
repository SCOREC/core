%module pyCore
%{
#include <mpi.h>
#include <vector>
#include <PCU_C.h>
#include <pcu_util.h>
#include <gmi.h>
#include <gmi_mesh.h>
#include <gmi_null.h>

#ifdef HAVE_SIMMETRIX
  #include <gmi_sim.h>
#endif

#include <lionPrint.h>

#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apf.h>
#include <apfNumbering.h>
#include <apfShape.h>
/* #include <apfVector.h> */
/* #include <apfMatrix.h> */
/* #include <apfDynamicVector.h> */
/* #include <canArray.h> */
#include <spr.h>
#include <maInput.h>
#include <ma.h>
#include <crv.h>

#ifdef HAVE_SIMMETRIX
  #include <sim_helper.h>
#endif
%}


%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

/* PCU RELATED WRAPPERS */
/* ==== FROM PCU.h ====*/
MPI_Comm PCU_Get_Comm(void);
int PCU_Comm_Init(void);
int PCU_Comm_Free(void);

int PCU_Comm_Self(void);
int PCU_Comm_Peers(void);
double PCU_Time(void);
bool PCU_Comm_Initialized(void);

/* ==== FROM pcu_util.h ====*/
void PCU_Assert_Fail(const char* msg);

/* This are defined as macros in the .h file. Apparaently, it is OK to Lie
to SWIG that these are functions ;)
*/
void PCU_ALWAYS_ASSERT(int cond);
void PCU_ALWAYS_ASSERT_VERBOSE(int cond, const char* msg);

/* These are helper functions used to inintialize/finalize SimModSuite */
#ifdef HAVE_SIMMETRIX
  void start_sim(const char* logfile = 0);
  void stop_sim();
  bool is_sim_started();
#endif

/* GMI RELATED WRAPPERS */
void gmi_register_mesh(void);
void gmi_register_null(void);

#ifdef HAVE_SIMMETRIX
  void gmi_register_sim(void);
  void gmi_sim_start(void);
  void gmi_sim_stop(void);
  void gmi_sim_stop(void);
  gmi_model* gmi_sim_load(const char* nativefile, const char* smdfile);
#endif


/* LIONPRINT RELATED WRAPPER */
void lion_set_verbosity(int lvl);

/* APF RELATED WRAPPERS */
%ignore can::Array::operator=;
%ignore can::Array::operator[];
%include<canArray.h>
%template(dcanarry) can::Array<double,0>;
%template(dcanarry3) can::Array<double,3>;


%include<apfDynamicArray.h>
/* %template(ddynamicarry) apf::DynamicArray<double>; */


%include<apfArray.h>
/* %template(darry3) apf::Array<double,3>; */
%include<apfNew.h>

%include<apfVector.h>
/* %template(dvector3) apf::Vector<3>; */
%include<apfMatrix.h>
/* %ignore apf::Matrix; */

%ignore apf::ElementVertOp;
%ignore apf::BuildCallback;
%include<apfMesh.h>
%include<apfMesh2.h>
%extend apf::Mesh2{
  void setVertScalarField(apf::Field* f, apf::MeshEntity* e, int downId, int downNode, double value)
  {
    apf::MeshEntity* downs[12];
    int nd = self->getDownward(e, 0, downs);
    PCU_ALWAYS_ASSERT(downId < nd);
    apf::setScalar(f, downs[downId], downNode, value);
  }
  void setVertVectorField(apf::Field* f, apf::MeshEntity* e, int downId, int downNode,
    double v1, double v2, double v3)
  {
    apf::MeshEntity* downs[12];
    int nd = self->getDownward(e, 0, downs);
    PCU_ALWAYS_ASSERT(downId < nd);
    apf::setVector(f, downs[downId], downNode, apf::Vector3(v1, v2, v3));
  }
  int getVertNumbering(apf::Numbering* n, apf::MeshEntity* e, int downId, int downNode, int downComponent)
  {
    apf::MeshEntity* downs[12];
    int nd = self->getDownward(e, 0, downs);
    PCU_ALWAYS_ASSERT(downId < nd);
    return apf::getNumber(n, downs[downId], downNode, downComponent);
  }
  double getVertScalarField(apf::Field* f, apf::MeshEntity* e, int downId, int downNode)
  {
    apf::MeshEntity* downs[12];
    int nd = self->getDownward(e, 0, downs);
    PCU_ALWAYS_ASSERT(downId < nd);
    return apf::getScalar(f, downs[downId], downNode);
  }
  apf::Vector3 getVertVectorField(apf::Field* f, apf::MeshEntity* e, int downId, int downNode)
  {
    apf::MeshEntity* downs[12];
    int nd = self->getDownward(e, 0, downs);
    PCU_ALWAYS_ASSERT(downId < nd);
    apf::Vector3 p;
    apf::getVector(f, downs[downId], downNode, p);
    return p;
  }
  void setIPxyz(apf::Field* f, apf::MeshEntity* e, int node, double x, double y, double z)
  {
    int dim = self->getDimension();
    PCU_ALWAYS_ASSERT(dim >= 2);
    double vec[3] = {x, y, z};
    apf::setComponents(f, e, node, vec);
  }
  double measureSize(apf::MeshEntity* v)
  {
    double sum = 0.0;
    int count = 1;
    int upCount = self->countUpward(v);
    for(int i = 0; i < upCount; ++i) {
      sum += measure(self, self->getUpward(v, i));
      count++;
    }
    return sum/count;
  }
  apf::Field* getCurrentIsoSize(const char* name)
  {
    apf::Field* currentSize = apf::createField(
        self, name, apf::SCALAR, apf::getLagrange(1));

    apf::Field* cnt = apf::createField(
        self, "current_size_cnt", apf::SCALAR, apf::getLagrange(1));

    apf::zeroField(currentSize);
    apf::zeroField(cnt);


    apf::MeshEntity* e;
    apf::MeshIterator* it = self->begin(0);
    while ( (e = self->iterate(it)) ) {
      double local_sum = 0.;
      double local_cnt = 0.;
      for (int i = 0; i < self->countUpward(e); i++) {
        local_sum += apf::measure(self, self->getUpward(e, i));
        local_cnt += 1.0;
      }
      apf::setScalar(currentSize, e, 0, local_sum);
      apf::setScalar(cnt, e, 0, local_cnt);

    }
    self->end(it);

    apf::accumulate(currentSize);
    apf::accumulate(cnt);

    it = self->begin(0);
    while ( (e = self->iterate(it)) ) {
      if (!self->isOwned(e)) continue;
      double sum = apf::getScalar(currentSize, e, 0);
      double count = apf::getScalar(cnt, e, 0);
      apf::setScalar(currentSize, e, 0, sum/count);
    }
    self->end(it);

    apf::synchronize(currentSize);
    self->removeField(cnt);
    apf::destroyField(cnt);
    return currentSize;
  }
  double getMinOfScalarField(apf::Field* field)
  {
    PCU_ALWAYS_ASSERT(apf::getValueType(field) == apf::SCALAR);
    double local_min = 1.0e32;
    apf::MeshEntity* e;
    apf::MeshIterator* it = self->begin(0);
    while ( (e = self->iterate(it)) ) {
      if (!self->isOwned(e))
        continue;
      double val = apf::getScalar(field, e, 0);
      if (val < local_min)
        local_min = val;
    }
    self->end(it);
    PCU_Min_Doubles(&local_min, 1);
    return local_min;
  }
  double getMaxOfScalarField(apf::Field* field)
  {
    PCU_ALWAYS_ASSERT(apf::getValueType(field) == apf::SCALAR);
    double local_max = -1.0e32;
    apf::MeshEntity* e;
    apf::MeshIterator* it = self->begin(0);
    while ( (e = self->iterate(it)) ) {
      if (!self->isOwned(e))
        continue;
      double val = apf::getScalar(field, e, 0);
      if (val > local_max)
        local_max = val;
    }
    self->end(it);
    PCU_Max_Doubles(&local_max, 1);
    return local_max;
  }
  bool isBoundingModelRegion(int rtag, int dim, int tag)
  {
    if (dim != 2) return false;
    gmi_model* gmodel = self->getModel();
    gmi_ent* gregion = gmi_find(gmodel, 3, rtag);
    gmi_set* adj = gmi_adjacent(gmodel, gregion, dim);

    for(int i = 0; i < adj->n; i++) {
      int adj_g_tag = gmi_tag(gmodel, adj->e[i]);
      if (adj_g_tag == tag) {
        gmi_free_set(adj);
        return true;
      }
    }
    gmi_free_set(adj);
    return false;
  }
}

#define __attribute__(x)
%ignore apf::fail;
%include<apf.h>
%include<apfNumbering.h>
%include<apfShape.h>



namespace apf {
  apf::Mesh2* makeEmptyMdsMesh(gmi_model* model, int dim, bool isMatched);
  apf::Mesh2* loadMdsMesh(const char* modelfile, const char* meshfile);
  apf::Mesh2* loadMdsMesh(gmi_model* model, const char* meshfile);
  void writeASCIIVtkFiles(const char* prefix, apf::Mesh2* m);
  /* void writeVtkFiles(const char* prefix, apf::Mesh* m, int cellDim = -1); */
  /* void writeVtkFiles(const char* prefix, apf::Mesh* m, */
  /*   std::vector<std::string> writeFields, int cellDim = -1); */
}
/* void getAlignment(Mesh* m, MeshEntity* elem, MeshEntity* boundary, */
/*     int& which, bool& flip, int& rotate); */

/* SPR RELATED WRAPPERS */
%include<spr.h>


/* MA RELATED WRAPPERS */
/* let swig know about the typedefs */
namespace ma {
  typedef apf::Vector3 Vector;
  typedef apf::Matrix3x3 Matrix;
  typedef apf::Mesh2 Mesh;
  typedef apf::MeshEntity Entity;
  typedef apf::MeshIterator Iterator;
  typedef apf::MeshTag Tag;
  typedef apf::DynamicArray<Entity*> EntityArray;
  typedef std::set<Entity*> EntitySet;
  typedef EntityArray Upward;
  typedef apf::Downward Downward;
  typedef apf::ModelEntity Model;
  Vector getPosition(Mesh* m, Entity* vertex);
  typedef apf::Copies Remotes;
  typedef apf::Parts Parts;
}

%include <maSolutionTransfer.h>
%include <maSize.h>
%include <maInput.h>

namespace ma {
  void adapt(Input* in);
  void adaptVerbose(Input* in, bool verbosef = false);
}

/* CRV RELATED WRAPPERS */
%rename(crvadapt) crv::adapt;
namespace crv {
  void adapt(ma::Input* in);
  class BezierCurver : public MeshCurver
  {
    public:
      BezierCurver(apf::Mesh2* m, int P, int B) : MeshCurver(m,P)
      {
        setBlendingOrder(apf::Mesh::TYPES,B);
      };
      virtual bool run();
  };
  void writeCurvedVtuFiles(apf::Mesh* m, int type, int n, const char* prefix);
  void writeCurvedWireFrame(apf::Mesh* m, int n, const char* prefix);
}
