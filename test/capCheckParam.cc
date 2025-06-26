#include <iostream>

#include <PCU.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfMesh2.h>
#include <gmi_cap.h>
#include <pcu_util.h>

#include <CreateMG_Framework_Geometry.h>
#include <CreateMG_Framework_Mesh.h>

void checkParametrization(CreateMG::MDBI* mdb, CreateMG::GDBI* gdb);

int main(int argc, char** argv)
{
  pcu::Init(&argc, &argv);
  { // pcu object scope
  pcu::PCU PCUobj;

  if (argc != 2) {
    if (PCUobj.Self() == 0)
      std::cerr << "usage: " << argv[0] << " <cre file .cre>\n";
    return EXIT_FAILURE;
  }

  const char* creFileName = argv[1];

  gmi_cap_start();
  gmi_register_cap();

  gmi_model* model = gmi_cap_load(creFileName);
  apf::Mesh2* mesh = apf::createCapMesh(model, &PCUobj);

  // check parametrization using capstone apis
  checkParametrization(apf::getCapNative(mesh), gmi_export_cap(model));
  apf::destroyMesh(mesh);
  gmi_cap_stop();
  } // pcu object scope
  pcu::Finalize();
}

void checkParametrization(CreateMG::MDBI* mdb, CreateMG::GDBI* gdb) {
  using namespace CreateMG;
  using namespace CreateMG::Mesh;
  MeshSmartIterator miter(mdb);
  mdb->get_topo_iterator(TOPO_VERTEX, miter);
  int count = 0;
  double sum = 0.0;
  for(mdb->iterator_begin(miter); !mdb->iterator_end(miter); mdb->iterator_next(miter)) {
    M_MTopo vert = mdb->iterator_value(miter);
    M_GTopo geom;
    GeometryTopoType gtype;
    mdb->get_geom_entity(vert, gtype, geom);
    if (!gdb->is_face(geom)) continue;
    double range_u[2];
    double range_v[2];
    gdb->get_parametrization_range(geom, 0, range_u[0], range_u[1]);
    gdb->get_parametrization_range(geom, 1, range_v[0], range_v[1]);
    GeometryTopoType gtype1;
     double u,v;
    mdb->get_vertex_uv_parameters(vert, u, v, gtype1);
    PCU_ALWAYS_ASSERT(gtype1 == gtype);

    // coordinate from mesh
    apf::Vector3 coord;
    mdb->get_vertex_coord(vert, &(coord[0]));

    // coordinate from surface
    vec3d x;
    gdb->get_point(geom, vec3d(u, v, 0.0), x);
    apf::Vector3 pcoord(x[0], x[1], x[2]);

    if (count < 50) {
      std::cout << count << ", "
        << u << ", " << v << ", "
        << range_u[0] << ", " << range_u[1]
        << ", " << range_v[0] << ", " << range_v[1]
        << ", " << (coord - pcoord).getLength() << std::endl;
    }
    sum += (coord-pcoord) * (coord-pcoord);
    count++;
  }
  std::cout << "norm of the difference vector is " << std::sqrt(sum)
    << std::endl;
}
