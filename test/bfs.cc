#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <gmi_mesh.h>
#include <PCU.h>

apf::MeshEntity* grabFirstVertex(apf::Mesh* m)
{
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v = m->iterate(it);
  m->end(it);
  return v;
}

void renderIntTag(apf::Mesh* m, apf::MeshTag* tag, const char* filename)
{
  apf::Numbering* n = apf::createNumbering(m,
      m->getTagName(tag), m->getShape(), 1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    if (m->hasTag(v, tag)) {
      int x;
      m->getIntTag(v, tag, &x);
      apf::number(n, v, 0, 0, x);
    } else {
      apf::number(n, v, 0, 0, 0);
    }
  }
  m->end(it);
  apf::writeVtkFiles(filename, m);
  apf::destroyNumbering(n);
}

void visit(apf::Mesh* m, apf::MeshTag* visited, apf::MeshEntity* vertex,
    int layer)
{
  int yes = layer;
  m->setIntTag(vertex, visited, &yes);
}

void runBFS(apf::Mesh* m, apf::MeshEntity* startVertex)
{
  apf::MeshTag* visited = m->createIntTag("visited",1);
  int layer=2;
  std::vector<apf::MeshEntity*> current;
  std::vector<apf::MeshEntity*> next;
  current.push_back(startVertex);
  visit(m,visited,startVertex,1);
  for (int i=0;i<layer;i++) {
    for (unsigned int j=0;j<current.size();j++) {
      apf::MeshEntity* vertex = current[j];
      apf::Up edges;
      m->getUp(vertex,edges);
      for (int k=0;k<edges.n;k++) {
        apf::MeshEntity* edge = edges.e[k];
        apf::Downward vertices;
        int nvertices = m->getDownward(edge,getDimension(m,edge)-1,vertices);
        assert(nvertices==2);
        apf::MeshEntity* v = vertices[0];
        if (v==vertex)
          v=vertices[1];
        if (m->hasTag(v,visited)) continue;
        next.push_back(v);
        visit(m,visited,v,i+2);
      }
      
    }
    current=next;
    next.clear();
  }
  renderIntTag(m, visited, "after");
}

int main(int argc, char** argv)
{
  assert(argc == 3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  runBFS(m, grabFirstVertex(m));
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

