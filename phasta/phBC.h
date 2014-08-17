
#ifndef PH_BC_H
#define PH_BC_H

#include <map>
#include <set>
#include <string>

#include <gmi.h>

/* full names and abbreviations for boundary conditions:

   Essential boundary conditions:

   density         D   rho
   temperature     T
   pressure        P
   comp 1          C1
   comp 3          C3
   scalar_1        S1  sc1
   scalar_2        S2  sc1
   scalar_3        S3  sc1
   scalar_4        S4  sc1
 * take bc from ic TBI

   Natural boundary conditions:

   mass flux        MF
   natural pressure NP
   traction vector  TV
   heat flux        HF
 * turbulence wall  TW
   scalar_1 flux    F1
   scalar_2 flux    F2
   scalar_3 flux    F3
   scalar_4 flux    F4
 * surf ID          SID

   A bit about comp1 and comp3: seems like both affect the
   first (direction, magnitude pair in the BC array).
   This is fluid velocity.
   Either comp1 or comp3 is applied, not both.
   comp3 implies that all 3 components of the direction are
   fixed, while comp1 means only one of them is set.
   The old code chooses the largest absolute value.

   Although the .spj file has four reals for the comp1/3
   condition, in practice the last is always zero and the
   code decomposes the vector into direction and magnitude
   based on the first three reals.

   TBI means to inherit the conditions placed on the model
   region on the model boundary as well.

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
  BC();
  ~BC();
  int tag;
  int dim;
  double* values;
  bool operator<(const BC& other) const;
};

struct FieldBCs
{
  int size;
  typedef std::set<BC> Set;
  Set bcs;
};

struct BCs
{
  typedef std::map<std::string, FieldBCs> Map;
  Map fields;
};

void readBCs(const char* filename, BCs& bcs);

bool applyNaturalBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs,
    double* values, int* bits);
bool applyEssentialBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs,
    double* values, int* bits);
bool applySolutionBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, double* values);

void getBCFaces(apf::Mesh* m, BCs& bcs, std::set<apf::ModelEntity*>& faces);

bool applyVelocityConstaints(gmi_model* gm, BCs& bcs, gmi_ent* e,
    double* BC, int* iBC);

double* getValuesOn(gmi_model* gm, FieldBCs& bcs, gmi_ent* ge);

bool hasBC(BCs& bcs, std::string const& name);

}

#endif
