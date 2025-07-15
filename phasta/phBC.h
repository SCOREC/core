#ifndef PH_BC_H
#define PH_BC_H

#include <map>
#include <set>
#include <string>
#include <apfVector.h>
#include <gmi.h>
#include "phInput.h"

/* full names and abbreviations for boundary conditions:

   Essential boundary conditions:

   density         D   rho
   temperature     T
   pressure        P
   comp1           C1
   comp3           C3
   scalar_1        S1  sc1
   scalar_2        S2  sc1
   scalar_3        S3  sc1
   scalar_4        S4  sc1
   comp1_elas      EC1
   comp3_elas      EC3
 
   Natural boundary conditions:

   mass flux        MF
   natural pressure NP
   traction vector  TV
   traction vector melas TVM
   heat flux        HF
 * turbulence wall  TW
   scalar_1 flux    F1
   scalar_2 flux    F2
   scalar_3 flux    F3
   scalar_4 flux    F4
 * surf ID          SID

   Initial conditions:

   initial pressure
   initial velocity
   initial temperature
   initial scalar_1
   initial scalar_2
   initial scalar_3
   initial scalar_4

   A .spj file is composed of two kinds of lines:

   # comments begin with a pound sign
   condition_name: model_id model_dim c0 c1 c2 c3

   The model_id and model_dim specify a geometric
   model entity. Typically model_dim=2 for boundary conditions
   on model faces, or model_dim=3 for initial conditions
   on model regions.

   All conditions have one component except for the following:
   comp1            4
   comp3            4
   initial velocity 3
   traction vector  3
   comp1_elas       4
   comp3_elas       4 
   traction vector melas  3

   A bit about comp1 and comp3: both of them are essential
   boundary conditions affecting velocity.
   comp3 is a strict constraint, velocity = m*(x,y,z)
   comp1 is a planar constraint, dot(velocity,(x,y,z)) = m;
   Both are listed in the spj file as: m x y z

   Inital velocity and traction vector are both listed
   as: x y z

   The scalars 1 through 4 are not always present.
   The number of scalars is computed as (ensa_dof - 5),
   and the BC and BCB arrays contain only enough entries
   for the scalars present. that is why the scalars are at
   the end of said arrays.
*/

namespace apf {
class Mesh;
class MeshEntity;
class ModelEntity;
}

namespace ph {

struct BC
{
  virtual ~BC() {}
  int tag;
  int dim;
  virtual double* eval(apf::Vector3 const& x) = 0;
  bool operator<(const BC& other) const;
};

struct ConstantBC : public BC
{
  ConstantBC();
  ~ConstantBC();
  virtual double* eval(apf::Vector3 const& x);
  double* value;
};

struct BCPointerLess
{
  bool operator()(BC const* a, BC const* b) const { return *a < *b; };
};

struct FieldBCs
{
  ~FieldBCs();
  typedef std::set<BC*, BCPointerLess> Set;
  Set bcs;
};

struct BCs
{
  typedef std::map<std::string, FieldBCs> Map;
  Map fields;
};

void readBCs(gmi_model* m, const char* attFile, bool axisymmetry, BCs& bcs);
void loadModelAndBCs(ph::Input& in, gmi_model*& m, BCs& bcs, pcu::PCU *PCUObj);

bool applyNaturalBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, apf::Vector3 const& x, double* values, int* bits);
bool applyEssentialBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, apf::Vector3 const& x, double* values, int* bits);
bool applySolutionBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, apf::Vector3 const& x, double* values);

bool applyVelocityConstaints(gmi_model* gm, BCs& bcs, gmi_ent* e,
    apf::Vector3 const& x, double* BC, int* iBC);

bool applyElasticConstaints(gmi_model* gm, BCs& bcs, gmi_ent* e,
    apf::Vector3 const& x, double* BC, int* iBC);

double* getBCValue(gmi_model* gm, FieldBCs& bcs, gmi_ent* ge,
    apf::Vector3 const& x = apf::Vector3(0,0,0));

bool haveBC(BCs& bcs, std::string const& name);

ConstantBC* makeConstantBC(BCs& bcs, std::string const& name, int dim, int tag,
    int nvals);

}

#endif
