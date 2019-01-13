#include <PCU.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <samSz.h>
#include <queue>
#include <ma.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>


#include "CapstoneModule.h"
#include "CreateMG_Framework_Core.h"
#include "CreateMG_Framework_Analysis.h"
#include "CreateMG_Framework_Application.h"
#include "CreateMG_Framework_Attributes.h"
#include "CreateMG_Framework_Core.h"
#include "CreateMG_Framework_Geometry.h"
#include "CreateMG_Framework_Mesh.h"

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;


typedef std::vector<double> row;

std::vector<row> readTable(const char* name);

std::vector<double> extractSurfaceData(const std::vector<row> &table,
    const int strandSize, const int col);

void removeUnusedVerts(apf::Mesh2* m, int offset);

apf::Mesh2* createVolumeMesh(apf::Mesh2* m, const std::vector<row> &t, int s);

apf::Field* addScalarField(apf::Mesh2* m, const std::vector<row> t, const char* name, int col, int strandSize);

apf::Field* addVector3Field(apf::Mesh2* m, const std::vector<row> t, const char* name,
    int col0, int col1, int col2, int strandSize);

struct SortingStruct
{
  apf::Vector3 v;
  double wm;
  bool operator<(const SortingStruct &other) const
  {
    return wm < other.wm;
  }
};

//gradation routines from Proteus

int gradeSizeModify(apf::Mesh* m, double gradingFactor, 
    double size[2], apf::Adjacent edgAdjVert, 
    apf::Adjacent vertAdjEdg,
    std::queue<apf::MeshEntity*> &markedEdges,
    apf::MeshTag* isMarked,
    int fieldType,
    int vecPos, //which idx of sizeVec to modify
    int idxFlag)

//General function to actually modify sizes
{
    //Determine a switching scheme depending on which vertex needs a modification
    int idx1,idx2;
    if(idxFlag == 0){
      idx1=0;
      idx2=1;
    } 
    else{
      idx1=1;
      idx2 = 0;
    } 
    
    int marker[3] = {0,1,0}; 
    double marginVal = 0.01;
    int needsParallel=0;

    if(fieldType == apf::SCALAR){
      apf::Field* size_iso = m->findField("size");

      if(size[idx1]>(gradingFactor*size[idx2])*(1+marginVal))
      {
        if(m->isOwned(edgAdjVert[idx1]))
        {
          size[idx1] = gradingFactor*size[idx2];
          apf::setScalar(size_iso,edgAdjVert[idx1],0,size[idx1]);
          m->getAdjacent(edgAdjVert[idx1], 1, vertAdjEdg);
          for (std::size_t i=0; i<vertAdjEdg.getSize();++i){
            m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
            //if edge is not already marked
            if(!marker[2]){
              m->setIntTag(vertAdjEdg[i],isMarked,&marker[1]);
              markedEdges.push(vertAdjEdg[i]);
            }
          }
        } //end isOwned
        else
        { //Pack information to owning processor
          needsParallel=1;
          apf::Copies remotes;
          m->getRemotes(edgAdjVert[idx1],remotes);
          double newSize = gradingFactor*size[idx2];
          int owningPart=m->getOwner(edgAdjVert[idx1]);
          PCU_COMM_PACK(owningPart, remotes[owningPart]);
          PCU_COMM_PACK(owningPart,newSize);
        }
      }

    }//end if apf::SCALAR
  return needsParallel;
}

void markEdgesInitial(apf::Mesh* m, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
//Function used to initially determine which edges need to be considered for gradation
{
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  double size[2];
  apf::MeshTag* isMarked = m->findTag("isMarked");
  apf::Field* size_iso = m->findField("size");
  apf::Adjacent edgAdjVert;
  apf::MeshEntity* edge;
  apf::MeshIterator* it = m->begin(1);
  while((edge=m->iterate(it))){
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      size[i]=apf::getScalar(size_iso,edgAdjVert[i],0);
    }
    if( (size[0] > gradingFactor*size[1]) || (size[1] > gradingFactor*size[0]) ){
      //add edge to a queue 
      markedEdges.push(edge);
      //tag edge to indicate that it is part of queue 
      m->setIntTag(edge,isMarked,&marker[1]); 
    }
    else{
      m->setIntTag(edge,isMarked,&marker[0]); 
    }
  }
  m->end(it); 
}

int serialGradation(apf::Mesh* m, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
//Function used loop over the mesh edge queue for gradation and modify the sizes
{
  double size[2];
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 
  apf::MeshTag* isMarked = m->findTag("isMarked");
  apf::Field* size_iso = m->findField("size");
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  apf::MeshEntity* edge;
  apf::MeshIterator* it = m->begin(1);
  int needsParallel=0;

  //perform serial gradation while packing necessary info for parallel
  while(!markedEdges.empty()){ 
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      size[i] = apf::getScalar(size_iso,edgAdjVert[i],0);
    }

    needsParallel+=gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::SCALAR,0, 0);
    needsParallel+=gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::SCALAR,0, 1);

    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }
  return needsParallel;
}

int gradeMesh(apf::Mesh* m)
//Function to grade isotropic mesh through comparison of edge vertex size ratios
//This implementation accounts for parallel meshes as well
//First do serial gradation. 
//If a shared entity has its size modified, then send new size to owning copy.
//After full loop over entities, have owning copy take minimum of all sizes received
//Flag adjacent entities to owning copy.
//Communicate to remote copies that a size was modified, and flag adjacent edges to remote copies for further gradation
{
  //unique to this code
  apf::Field* size_iso = m->findField("size");
  double gradingFactor = 1.2;
  //

  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  double size[2];
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  apf::MeshIterator* it;
  markEdgesInitial(m,markedEdges,gradingFactor);

  int needsParallel=1;
  int nCount=1;
  while(needsParallel)
  {
    PCU_Comm_Begin();
    needsParallel = serialGradation(m,markedEdges,gradingFactor);

    PCU_Add_Ints(&needsParallel,1);
    PCU_Comm_Send(); 

    apf::MeshEntity* ent;
    double receivedSize;
    double currentSize;
    double newSize;

    //Need a container to get all entitites that need to be updated on remotes
    std::queue<apf::MeshEntity*> updateRemoteVertices;

    apf::Copies remotes;
    //owning copies are receiving
    while(PCU_Comm_Receive())
    {
      PCU_COMM_UNPACK(ent);
      PCU_COMM_UNPACK(receivedSize);

      if(!m->isOwned(ent)){
        std::cout<<"THERE WAS AN ERROR"<<std::endl;
        std::exit(1);
      }

      currentSize = apf::getScalar(size_iso,ent,0);
      newSize = std::min(receivedSize,currentSize);
      apf::setScalar(size_iso,ent,0,newSize);
      
      //add adjacent edges into Q
      m->getAdjacent(ent, 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i)
      {
        edge = vertAdjEdg[i];
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        if(!marker[2])
        {
          markedEdges.push(edge);
          //tag edge to indicate that it is part of queue 
          m->setIntTag(edge,isMarked,&marker[1]);
        }
      }
      updateRemoteVertices.push(ent);
    }

    PCU_Comm_Begin();

    while(!updateRemoteVertices.empty())
    { 
      ent = updateRemoteVertices.front();
      //get remote copies and send updated mesh sizes
      m->getRemotes(ent,remotes);
      currentSize = apf::getScalar(size_iso,ent,0);
      for(apf::Copies::iterator iter=remotes.begin(); iter!=remotes.end();++iter)
      {
        PCU_COMM_PACK(iter->first, iter->second);
      }
      updateRemoteVertices.pop();
    }

    PCU_Comm_Send();
    //while remote copies are receiving
    while(PCU_Comm_Receive())
    {
      //unpack
      PCU_COMM_UNPACK(ent);
      //PCU_COMM_UNPACK(receivedSize);
      assert(!m->isOwned(ent));

      if(m->isOwned(ent)){
        std::cout<<"Problem occurred\n";
        std::exit(1);
      }

      //add adjacent edges into Q
      m->getAdjacent(ent, 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i)
      {
        edge = vertAdjEdg[i];
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        if(!marker[2])
        {
          markedEdges.push(edge);
          //tag edge to indicate that it is part of queue 
          m->setIntTag(edge,isMarked,&marker[1]);
        }
      }
    }
    apf::synchronize(size_iso);

  } //end outer while

  //Cleanup of edge marker field
  it = m->begin(1);
  while((edge=m->iterate(it))){
    m->removeTag(edge,isMarked);
  }
  m->end(it); 
  m->destroyTag(isMarked);

  //apf::synchronize(size_iso);
  return needsParallel;
}


//



int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  if (argc != 7) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0]
        << " <cre file .cre> <data file .txt> <data offset> <strand size> <desired max size> <error reduction factor>\n";
    return EXIT_FAILURE;
  }


  gmi_register_mesh();
  gmi_register_null();

  const char* creFileName = argv[1];
  const char* dataFileName = argv[2];
  int offset = atoi(argv[3]);
  int strandSize = atoi(argv[4]);
  double h_global = atof(argv[5]);
  double factor = atof(argv[6]);

  // load capstone mesh
  // create an instance of the Capstone Module activating CREATE/CREATE/CREATE
  // for the Geometry/Mesh/Attribution databases
  /* const std::string gdbName("Geometry Database : Create");// Switch Create with SMLIB for CAD */
  const std::string gdbName("Geometry Database : SMLIB");// Switch Create with SMLIB for CAD
  const std::string mdbName("Mesh Database : Create");
  const std::string adbName("Attribution Database : Create");

  CapstoneModule  cs("the_module", gdbName.c_str(), mdbName.c_str(), adbName.c_str());

  GeometryDatabaseInterface     *g = cs.get_geometry();
  MeshDatabaseInterface         *m = cs.get_mesh();
  AppContext                    *c = cs.get_context();


  PCU_ALWAYS_ASSERT(g);
  PCU_ALWAYS_ASSERT(m);
  PCU_ALWAYS_ASSERT(c);

  v_string filenames;
  filenames.push_back(creFileName);

  M_GModel gmodel = cs.load_files(filenames);

  int numbreps = 0;
  MG_CALL(g->get_num_breps(numbreps));
  std::cout << "number of b reps is " << numbreps << std::endl;
  if(numbreps == 0)
      error(HERE, ERR_INVALID_INPUT, "Model is empty");

  M_MModel mmodel;
  // Pick the volume mesh model from associated mesh models to this geom model
  std::vector<M_MModel> mmodels;
  MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));
  for(std::size_t i = 0; i < mmodels.size(); ++i)
  {
      M_MModel ammodel = mmodels[i];
      std::size_t numregs = 0;
      std::size_t numfaces = 0;
      std::size_t numedges = 0;
      std::size_t numverts = 0;
      MG_API_CALL(m, set_current_model(ammodel));
      MG_API_CALL(m, get_num_topos(TOPO_REGION, numregs));
      MG_API_CALL(m, get_num_topos(TOPO_FACE, numfaces));
      MG_API_CALL(m, get_num_topos(TOPO_EDGE, numedges));
      MG_API_CALL(m, get_num_topos(TOPO_VERTEX, numverts));
      std::cout << "num regions is " << numregs << std::endl;
      std::cout << "num faces   is " << numfaces << std::endl;
      std::cout << "num edges   is " << numedges << std::endl;
      std::cout << "num verts   is " << numverts << std::endl;
      std::cout << "-----------" << std::endl;
      if(numregs > 0)
      {
	  mmodel = ammodel;
	  break;
      }
  }

  /* SET THE ADJACENCIES */
  MG_API_CALL(m, set_adjacency_state(REGION2FACE|
                                     REGION2EDGE|
                                     REGION2VERTEX|
                                     FACE2EDGE|
                                     FACE2VERTEX));
  MG_API_CALL(m, set_reverse_states());
  MG_API_CALL(m, set_adjacency_scope(TOPO_EDGE, SCOPE_FULL));
  MG_API_CALL(m, set_adjacency_scope(TOPO_FACE, SCOPE_FULL));
  MG_API_CALL(m, compute_adjacency());


  gmi_cap_start();
  gmi_register_cap();

  printf("---- CapStone Mesh Loaded. \n");

  // read the data into a vector array for now
  printf("\n---- Reading Strand Data. \n");
  std::vector<row> table = readTable(dataFileName);
  printf("---- Reading Strand Data: Done. \n");

  // create the mesh object (this one is CapStone underneath)
  printf("\n---- Creating Mesh Wrapper Object. \n");
  apf::Mesh2* mesh = apf::createMesh(m,g);
  printf("---- Creating Mesh Wrapper Object: Done. \n");

  // remove unused verts
  printf("\n---- Removing Extra Verts. \n");
  removeUnusedVerts(mesh, offset);
  printf("---- Removing Extra Verts: Done. \n");

  // make the volume mesh (this one is MDS underneath)
  printf("\n---- Creating Volume Mesh. \n");
  apf::Mesh2* volMesh = createVolumeMesh(mesh, table, strandSize);
  printf("---- Creating Volume Mesh: Done. \n");


  // get all fields and add them to the mesh
  printf("\n---- Adding Fields to Volume Mesh. \n");
  apf::Field* rhoField	   = addScalarField(volMesh, table, "rho", 3, strandSize);
  apf::Field* rho_uvwField = addVector3Field(volMesh, table, "rho_uvw", 4, 5, 6, strandSize);
  apf::Field* eField       = addScalarField(volMesh, table, "e", 7, strandSize);
  apf::Field* nuField      = addScalarField(volMesh, table, "nu", 8, strandSize);
  printf("---- Adding Fields to Volume Mesh: Done. \n");


  printf("\n---- Writing Volume Mesh/Fields to VTK. \n");
  apf::writeVtkFiles("v_mesh", volMesh);
  printf("---- Writing Volume Mesh/Fields to VTK: Done. \n");

  printf("\n---- Printing Mesh Stats. \n");
  printf("number of mesh verts: %d\n", volMesh->count(0));
  printf("number of mesh regions: %d\n", volMesh->count(3));
  printf("---- Printing Mesh Stats: Done. \n");


  //Get Size Field for Adapt
  apf::Field* gradRhoField = apf::recoverGradientByVolume(rhoField); 
  apf::Field* hessianRhoField = apf::recoverGradientByVolume(gradRhoField); 
  apf::Field* gradEField  = apf::recoverGradientByVolume(eField); 
  apf::Field* hessianEField  = apf::recoverGradientByVolume(gradEField); 

  apf::Field* metricField = apf::createLagrangeField(volMesh,"metric",apf::MATRIX,1);
  apf::Field* lambdaMaxField = apf::createLagrangeField(volMesh,"lambdaMax",apf::SCALAR,1);
  apf::Field* sizeField = apf::createLagrangeField(volMesh,"size",apf::SCALAR,1);
  apf::Field* surfaceSizeField = apf::createLagrangeField(mesh,"surface_size",apf::SCALAR,1);
  
  //get current size field
  apf::Field* currentSize = samSz::isoSize(mesh);
  

  apf::MeshIterator* it;
  apf::MeshEntity* vert;

  //transfer rho and e field to surface mesh
  apf::Field* surfaceRhoField = apf::createLagrangeField(mesh,"surfaceDensity",apf::SCALAR,1);
  apf::Field* surfaceEField = apf::createLagrangeField(mesh,"surfaceE",apf::SCALAR,1);

  it = mesh->begin(0);
  apf::MeshIterator* itVol = volMesh->begin(0); 

  while( (vert = mesh->iterate(it)) )
  {
    apf::MeshEntity* vertVol = volMesh->iterate(itVol);
    double rhoVal = apf::getScalar(rhoField,vertVol,0); 
    double eVal = apf::getScalar(eField,vertVol,0); 
    apf::setScalar(surfaceRhoField,vert,0,rhoVal);
    apf::setScalar(surfaceEField,vert,0,eVal);
  }
  mesh->end(it);
  volMesh->end(itVol);

  //getSpeed
  apf::Field* speedField = apf::createLagrangeField(volMesh,"speed",apf::SCALAR,1);
  it = volMesh->begin(0);
  double minNu = 1.0;
  while( (vert = volMesh->iterate(it)) )
  {
    double rhoVal = apf::getScalar(rhoField,vert,0); 
    apf::Vector3 rhoVelocity;
    apf::getVector(rho_uvwField,vert,0,rhoVelocity); 
    double speed = rhoVelocity.getLength()/rhoVal;
    apf::setScalar(speedField,vert,0,speed);
    if(apf::getScalar(nuField,vert,0)<minNu)
      minNu = apf::getScalar(nuField,vert,0);
  }
  volMesh->end(it);
  apf::Field* gradSpeedField = apf::recoverGradientByVolume(speedField); 
  apf::Field* hessianSpeedField = apf::recoverGradientByVolume(gradSpeedField); 

  //get eigenvalues
  it = volMesh->begin(0);
  while( (vert = volMesh->iterate(it)) )
  {
      apf::Matrix3x3 metric;
      //apf::getMatrix(hessianRhoField, vert, 0, metric);
      apf::getMatrix(hessianEField, vert, 0, metric);
      //apf::getMatrix(hessianSpeedField, vert, 0, metric);

      apf::Vector3 eigenVectors[3];
      double eigenValues[3];
      metric = (metric+apf::transpose(metric))/2.0;
      apf::setMatrix(metricField,vert,0,metric);
      apf::eigen(metric, eigenVectors, eigenValues);

      // Sort the eigenvalues and corresponding vectors
      // Larger eigenvalues means a need for a finer mesh
      SortingStruct ssa[3];
      for (int i = 0; i < 3; ++i)
      {
        ssa[i].v = eigenVectors[i];
        ssa[i].wm = std::fabs(eigenValues[i]);
      }
      std::sort(ssa, ssa + 3);

      assert(ssa[2].wm >= ssa[1].wm);
      assert(ssa[1].wm >= ssa[0].wm);
      apf::setScalar(lambdaMaxField,vert,0,ssa[2].wm);
  }
  volMesh->end(it);

 
  //find global lambda max, h_lambdamax
  //on the surface mesh, based on the results of the volume mesh
  it = mesh->begin(0);
  itVol = volMesh->begin(0); 
  double lambda_max = 0.0;
  double h_lambdamax =0.0;
  while( (vert = mesh->iterate(it)) )
  {
    apf::MeshEntity* vertVol = volMesh->iterate(itVol);
    if(lambda_max < apf::getScalar(lambdaMaxField,vertVol,0))
    {
      lambda_max = apf::getScalar(lambdaMaxField,vertVol,0);
      h_lambdamax = apf::getScalar(currentSize,vert,0);
    }
  }
  volMesh->end(it);

  double lambda_cutoff = lambda_max*1e-10;
  
  //set size field
  it=volMesh->begin(0);
  int counter3 = 0;
  double h_special = -1;
  while( (vert = volMesh->iterate(it)) )
  {    
    double h_v = h_special; //set default size as user specified input
    double lambda_vert = apf::getScalar(lambdaMaxField,vert,0);
    if(lambda_vert > lambda_cutoff)
    {
      h_v = sqrt(lambda_max/lambda_vert/factor)*h_lambdamax;
    }
    else
      counter3++;

    //need to expose this as user input
    if(h_v > h_global) //maximum value for size field
      h_v = h_global;
    apf::setScalar(sizeField,vert,0,h_v);
  }
  volMesh->end(it);

  gradeMesh(volMesh);
  
  apf::writeVtkFiles("initialSurfaceMesh", mesh);
  apf::writeVtkFiles("volumeMeshSizeField", volMesh);

  //get the size field on the surface
  //assumes that the vertices are ordered exactly the same way
 
  it = mesh->begin(0);
  itVol = volMesh->begin(0); 
  while( (vert = mesh->iterate(it)) )
  {
    apf::MeshEntity* vertVol = volMesh->iterate(itVol);
    apf::setScalar(surfaceSizeField,vert,0,apf::getScalar(sizeField,vertVol,0));
  }
  mesh->end(it);
  volMesh->end(itVol);

  //adapt the mesh
  ma::Input* in;
  apf::Field* adaptSize  = apf::createFieldOn(mesh, "adapt_size", apf::SCALAR);
  apf::copyData(adaptSize, surfaceSizeField);
  in = ma::configure(mesh, adaptSize);
  ma::validateInput(in);
  in->shouldRunPreZoltan = true;
  in->shouldRunMidZoltan = true;
  in->shouldRunPostZoltan = true;
  in->maximumImbalance = 1.05;
  in->maximumIterations = 10;
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  //in->shouldTransferToClosestPoint = true;
  in->shouldForceAdaptation = true;
  //in->debugFolder = "./debug";
  ma::adaptVerbose(in);
  //ma::adaptVerbose(in, true);
  //
  apf::writeVtkFiles("adaptedMesh", mesh);


  apf::destroyField(gradRhoField);
  apf::destroyField(hessianRhoField);
  apf::destroyField(gradEField);
  apf::destroyField(hessianEField);
  apf::destroyField(surfaceRhoField);
  apf::destroyField(surfaceEField);
  apf::destroyField(metricField);
  apf::destroyField(lambdaMaxField);
  apf::destroyField(sizeField);
  apf::destroyField(surfaceSizeField);
  apf::destroyField(currentSize);
  apf::destroyField(speedField);
  apf::destroyField(gradSpeedField);
  apf::destroyField(hessianSpeedField);

  //End TODO

  gmi_cap_stop();

  PCU_Comm_Free();
  MPI_Finalize();
}

std::vector<row> readTable(const char* name)
{
  std::vector<row> table;
  std::ifstream file(name);
  while(!file.eof()){
    std::string str;
    std::getline(file, str);
    std::stringstream ss(str);
    double value;
    row r;
    while (ss >> value)
      r.push_back(value);
    if (r.size() > 0)
      table.push_back(r);
  }
  return table;
}

static std::vector<double> extractNthStrandData(const std::vector<row> &table,
    const int strandSize, const int n, const int col)
{
  std::vector<double> outData;
  int i = 0;
  while (i*strandSize+(strandSize-n) < int(table.size())) {
    outData.push_back(table[i*strandSize+(strandSize - n)][col]);
    i++;
  }
  return outData;
}

std::vector<double> extractSurfaceData(const std::vector<row> &table,
    const int strandSize, const int col)
{
  // surface data corresponds to strand 1;
  return extractNthStrandData(table, strandSize, 1, col);
}

apf::Field* addScalarField(apf::Mesh2* m, const std::vector<row> t, const char* name, int col, int strandSize)
{
  apf::Numbering* nLayer = m->findNumbering("layer_num");
  apf::Numbering* nVid   = m->findNumbering("vid_num");
  PCU_ALWAYS_ASSERT(nLayer);
  PCU_ALWAYS_ASSERT(nVid);

  apf::Field* f = apf::createFieldOn(m, name, apf::SCALAR);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  int count = 0;
  while ( (e = m->iterate(it)) ) {
    int layer = apf::getNumber(nLayer, e, 0, 0);
    int vid   = apf::getNumber(nVid, e, 0, 0);
    double value = t[vid*strandSize+(strandSize - layer)][col];
    apf::setScalar(f, e, 0, value);
    count++;
  }
  m->end(it);
  return f;
}

apf::Field* addVector3Field(apf::Mesh2* m, const std::vector<row> t, const char* name,
    int col0, int col1, int col2, int strandSize)
{
  apf::Numbering* nLayer = m->findNumbering("layer_num");
  apf::Numbering* nVid   = m->findNumbering("vid_num");
  PCU_ALWAYS_ASSERT(nLayer);
  PCU_ALWAYS_ASSERT(nVid);

  apf::Field* f = apf::createFieldOn(m, name, apf::VECTOR);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  int count = 0;
  while ( (e = m->iterate(it)) ) {
    int layer = apf::getNumber(nLayer, e, 0, 0);
    int vid   = apf::getNumber(nVid, e, 0, 0);
    double value0 = t[vid*strandSize+(strandSize - layer)][col0];
    double value1 = t[vid*strandSize+(strandSize - layer)][col1];
    double value2 = t[vid*strandSize+(strandSize - layer)][col2];
    apf::setVector(f, e, 0, apf::Vector3(value0, value1, value2));
    count++;
  }
  m->end(it);
  return f;
}



void removeUnusedVerts(apf::Mesh2* m, int offset)
{
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  int count = 0;
  while ( (e = m->iterate(it)) ) {
    if (count < offset) {
      PCU_ALWAYS_ASSERT(m->countUpward(e) == 0);
      m->destroy(e);
    }
    else
      PCU_ALWAYS_ASSERT(m->countUpward(e) != 0);
    count++;
  }
  m->end(it);
}

static std::vector<apf::Vector3>
readLayerCoordinates(int layer, const std::vector<row> &t, int s)
{
  std::vector<double> xs = extractNthStrandData(t, s, layer, 0);
  std::vector<double> ys = extractNthStrandData(t, s, layer, 1);
  std::vector<double> zs = extractNthStrandData(t, s, layer, 2);

  std::vector<apf::Vector3> res;

  for (std::size_t i = 0; i < xs.size(); i++)
    res.push_back(apf::Vector3(xs[i], ys[i], zs[i]));

  return res;
}


apf::Mesh2* createVolumeMesh(apf::Mesh2* m, const std::vector<row> &t, int s)
{
  // add numbering to the surface mash
  apf::Numbering* n = apf::createNumbering(m, "surf_verts_num", apf::getLagrange(1) , 1);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  int count = 0;
  while ( (e = m->iterate(it)) ) {
    apf::number(n, e, 0, 0, count);
    count++;
  }
  m->end(it);


  // get the original faces
  std::vector<apf::MeshEntity*> faces;
  it = m->begin(2);
  while ( (e = m->iterate(it)) )
    faces.push_back(e);
  m->end(it);

  // create and empty mds mesh
  apf::Mesh2* vMesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 3, false);
  apf::Numbering* nLayer = apf::createNumbering(vMesh, "layer_num", apf::getLagrange(1) , 1);
  apf::Numbering* nVID = apf::createNumbering(vMesh, "vid_num", apf::getLagrange(1) , 1);

  // process each layer
  // layer = 1 corresponds to the surface nodes
  // treat it out side of the loop
  std::vector<apf::Vector3> coords = readLayerCoordinates(1, t, s);
  std::vector<apf::MeshEntity*> lastVs;
  for (std::size_t i = 0; i < coords.size(); i++) {
    apf::MeshEntity* v = vMesh->createVertex(0, coords[i], apf::Vector3(0, 0, 0));
    apf::number(nLayer, v, 0, 0, 1);
    apf::number(nVID, v, 0, 0, i);
    lastVs.push_back(v);
  }

  PCU_ALWAYS_ASSERT(coords.size() == lastVs.size());

  // now process the remaining layers
  for (int layer = 2; layer < s+1; layer++) {
    // read the coords first
    std::vector<apf::Vector3> coords = readLayerCoordinates(layer, t, s);
    // create the new verts
    std::vector<apf::MeshEntity*> vs;
    for (std::size_t i = 0; i < coords.size(); i++) {
      apf::MeshEntity* v = vMesh->createVertex(0, coords[i], apf::Vector3(0, 0, 0));
      apf::number(nLayer, v, 0, 0, layer);
      apf::number(nVID, v, 0, 0, i);
      vs.push_back(v);
    }

    // create the prisms
    for (std::size_t i = 0; i < faces.size(); i++) {
      apf::MeshEntity* f = faces[i];
      apf::MeshEntity* dv[3];
      m->getDownward(f, 0, dv);
      int id0 = apf::getNumber(n, dv[0], 0, 0);
      int id1 = apf::getNumber(n, dv[1], 0, 0);
      int id2 = apf::getNumber(n, dv[2], 0, 0);

      apf::MeshEntity* prismVs[6];

      prismVs[0] = lastVs[id0];
      prismVs[1] = lastVs[id1];
      prismVs[2] = lastVs[id2];

      prismVs[3] = vs[id0];
      prismVs[4] = vs[id1];
      prismVs[5] = vs[id2];

      apf::buildElement(vMesh, 0, apf::Mesh::PRISM, prismVs);
    }
    // set the lastVs to Vs
    for (std::size_t i = 0; i < vs.size(); i++)
      lastVs[i] = vs[i];
  }
  vMesh->acceptChanges();
  apf::deriveMdsModel(vMesh);
  apf::verify(vMesh);
  return vMesh;
}
