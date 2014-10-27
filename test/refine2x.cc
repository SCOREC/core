#include <ma.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

class AnisotropicX: public ma::AnisotropicFunction {
  public:
    AnisotropicX(ma::Mesh* mesh, int splitDir) {
      m = mesh;
      dir = splitDir;
      assert(dir>=0 && dir<3);
      isf = new ma::IdentitySizeField(m);
      fLen = createFieldOn(m, "incidentEdgeLength", apf::SCALAR);
      fCnt = createFieldOn(m, "incidentEdgeCount", apf::SCALAR);
      getVtxSize();
    }
    ~AnisotropicX() {
      apf::destroyField(fLen);
      apf::destroyField(fCnt);
      delete isf;
    }
    void getValue(ma::Entity* v, ma::Matrix& r, ma::Vector& h) {
      r = ma::Matrix(1,0,0,
                     0,1,0,
                     0,0,1);
      int cnt = static_cast<int>(apf::getScalar(fCnt, v, 0));
      double l = apf::getScalar(fLen, v, 0) / cnt;
      const double f = 1.8;
      double sz[3] = {l,l,l};
      for(int i=0; i<3; i++) 
        if( i == dir ) 
          sz[i] /= f;
      h = ma::Vector(sz);
    }
  private:
    ma::Mesh* m;
    int dir;
    ma::SizeField* isf;
    apf::Field* fLen;
    apf::Field* fCnt;

    void getVtxSize() {
      double len;
      int cnt;
      ma::Entity* vtx;
      apf::MeshIterator* itr = m->begin(0);
      while( (vtx = m->iterate(itr)) ) {
        getEdgeLenAndCnt(vtx, len, cnt);
        apf::setScalar(fLen, vtx, 0, len);
        apf::setScalar(fCnt, vtx, 0, cnt);
      }
      m->end(itr);
      apf::accumulate(fLen);
      apf::accumulate(fCnt);
      apf::synchronize(fLen);
      apf::synchronize(fCnt);
    }
    void getEdgeLenAndCnt(ma::Entity* v, double& len, int& cnt) {
      len = 0;
      cnt = 0;
      apf::Up edges;
      m->getUp(v, edges);
      for(int eIdx=0; eIdx < edges.n; eIdx++) {
        if( m->isOwned(edges.e[eIdx])) {
          cnt++;
          len += isf->measure(edges.e[eIdx]);
        }
      }
    }
};

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if (argc != 5) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0] 
        << " <model file> <in mesh> <split direction=[0-2]> <out mesh> \n";
    return EXIT_FAILURE;
  }
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(argv[1],argv[2]);
  AnisotropicX* ansx = new AnisotropicX(m, atoi(argv[3]));
  ma::Input* in = ma::configure(m, ansx);
  in->shouldRunPreParma = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->maximumIterations = 10;
  if (in->shouldSnap) {
    in->shouldSnap = false;
    assert(in->shouldTransferParametric);
  }
  in->shouldFixShape = false;
  ma::adapt(in);
  m->verify();
  delete ansx;
  m->writeNative(argv[4]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
