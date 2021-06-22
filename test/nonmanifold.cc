#include "apfMesh2.h"
#include "apfMDS.h"
#include "gmi.h"
#include "gmi_null.h"
#include "gmi_mesh.h"
#include "PCU.h"
#include <lionPrint.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);                                                                                              
    PCU_Comm_Init(); 
    lion_set_verbosity(1);  //write basic info

    int dim=2;

    gmi_register_null();
    gmi_model* gm = gmi_load(".null");
    apf::Mesh2 *m=apf::makeEmptyMdsMesh(gm, dim, false);

    // > create edge
    apf::Vector3 coord0(0., 0., 0.);
    apf::Vector3 coord1(1., 0., 0.);
    apf::Vector3 param(0,0,0);

    apf::MeshEntity *v0=m->createVertex(NULL, coord0, param);
    apf::MeshEntity *v1=m->createVertex(NULL, coord1, param);

    std::vector<apf::MeshEntity *> verts={v0, v1};
    apf::buildElement(m, NULL, apf::Mesh::EDGE, verts.data());
    // <

    // > create quad
    apf::Vector3 coord2(0., 0., 1.);
    apf::Vector3 coord3(1., 0., 1.);
    apf::Vector3 coord4(1., 1., 1.);
    apf::Vector3 coord5(0., 1., 1.);

    apf::MeshEntity *v2=m->createVertex(NULL, coord2, param);
    apf::MeshEntity *v3=m->createVertex(NULL, coord3, param);
    apf::MeshEntity *v4=m->createVertex(NULL, coord4, param);
    apf::MeshEntity *v5=m->createVertex(NULL, coord5, param);

    verts={v2, v3, v4, v5};
    apf::buildElement(m, NULL, apf::Mesh::QUAD, verts.data());
    // <

    // apf::alignMdsRemotes(m);
    apf::deriveMdsModel(m);
    m->acceptChanges();
    m->verify();
    
    // > test some adjaceny queries 
    apf::MeshEntity *e;
    apf::Adjacent adj;
    for (int d=0; d<=dim; ++d) {
        auto mit = m->begin(d);
        while ((e=m->iterate(mit))) {
            std::cout << "Entitiy #" << apf::getMdsIndex(m, e) << " of dim " << d << std::endl;           
            for (int dd=0; dd<=dim; ++dd) {
                m->getAdjacent(e, dd, adj);
                std::cout << "  dim " << dd << " adjacency: [";
                for (auto ee: adj) {
                    std::cout << " " << apf::getMdsIndex(m, ee);
                }
                std::cout << " ]\n";
            }
        }
        m->end(mit);
    }
    // <

    return 0;
}
