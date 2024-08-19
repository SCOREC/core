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
#include <lionPrint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <memory>


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

struct FieldOfInterest{
  FieldOfInterest()
  {
    lambda_max = 0.0;
    h_lambdamax = 0.0;
  }
  double lambda_max;
  double h_lambdamax;
  apf::Field* lambdaMaxField;
  apf::Field* lambdaStrandMax;
  apf::Field* sizeField;
  double lambda_cutoff()
  {
    return lambda_max*1e-10;
  }
};

std::vector<row> readTable(const char* name);

std::vector<double> extractSurfaceData(const std::vector<row> &table,
    const int strandSize, const int col);

apf::Mesh2* createVolumeMesh(apf::Mesh2* m, const std::vector<row> &t, int s,
    std::vector<std::vector<apf::MeshEntity*> > &surfToStrandMap);

apf::Field* addScalarField(apf::Mesh2* m, const std::vector<row> t, const char* name, int col, int strandSize);

apf::Field* addVector3Field(apf::Mesh2* m, const std::vector<row> t, const char* name,
    int col0, int col1, int col2, int strandSize);

void writeCre(CapstoneModule& cs, const std::string& fileName);

void writeMdsMesh(apf::Mesh2* m, const char* name, const char* fieldName);

void computeSizeDistribution(apf::Mesh2* m, int factor,
    std::vector<int>& binCount, std::vector<double>& binArea);


void adjustRefinementLevel(apf::Mesh2* m, apf::Field* finalSize,
    apf::Field* currentSize, int maxLevel);

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

int gradeSizeModify(apf::Mesh* m, apf::Field* size_iso,double gradingFactor, 
    double size[2], apf::Adjacent edgAdjVert, 
    apf::Adjacent vertAdjEdg,
    std::queue<apf::MeshEntity*> &markedEdges,
    apf::MeshTag* isMarked,
    int fieldType,
    int vecPos, //which idx of sizeVec to modify
    int idxFlag)

//General function to actually modify sizes
{
    (void)vecPos;
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
      //apf::Field* size_iso = m->findField("size");

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
          m->getPCU()->Pack(owningPart, remotes[owningPart]);
          m->getPCU()->Pack(owningPart,newSize);
        }
      }

    }//end if apf::SCALAR
  return needsParallel;
}

void markEdgesInitial(apf::Mesh* m, apf::Field* size_iso, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
//Function used to initially determine which edges need to be considered for gradation
{
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  double size[2];
  apf::MeshTag* isMarked = m->findTag("isMarked");
//  apf::Field* size_iso = m->findField("size");
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

int serialGradation(apf::Mesh* m,apf::Field* size_iso, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
//Function used loop over the mesh edge queue for gradation and modify the sizes
{
  double size[2];
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 
  apf::MeshTag* isMarked = m->findTag("isMarked");
//  apf::Field* size_iso = m->findField("size");
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  apf::MeshEntity* edge;
  int needsParallel=0;

  //perform serial gradation while packing necessary info for parallel
  while(!markedEdges.empty()){ 
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      size[i] = apf::getScalar(size_iso,edgAdjVert[i],0);
    }

    needsParallel+=gradeSizeModify(m,size_iso, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::SCALAR,0, 0);
    needsParallel+=gradeSizeModify(m,size_iso, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::SCALAR,0, 1);

    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }
  return needsParallel;
}

int gradeMesh(apf::Mesh* m,apf::Field* size_iso)
//Function to grade isotropic mesh through comparison of edge vertex size ratios
//This implementation accounts for parallel meshes as well
//First do serial gradation. 
//If a shared entity has its size modified, then send new size to owning copy.
//After full loop over entities, have owning copy take minimum of all sizes received
//Flag adjacent entities to owning copy.
//Communicate to remote copies that a size was modified, and flag adjacent edges to remote copies for further gradation
{
  //unique to this code
  //apf::Field* size_iso = m->findField("size");
  double gradingFactor = 1.2;
  //

  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  apf::MeshIterator* it;
  markEdgesInitial(m,size_iso,markedEdges,gradingFactor);

  int needsParallel=1;
  while(needsParallel)
  {
    m->getPCU()->Begin();
    needsParallel = serialGradation(m,size_iso,markedEdges,gradingFactor);

    m->getPCU()->Add<int>(&needsParallel,1);
    m->getPCU()->Send(); 

    apf::MeshEntity* ent;
    double receivedSize;
    double currentSize;
    double newSize;

    //Need a container to get all entitites that need to be updated on remotes
    std::queue<apf::MeshEntity*> updateRemoteVertices;

    apf::Copies remotes;
    //owning copies are receiving
    while(m->getPCU()->Receive())
    {
      m->getPCU()->Unpack(ent);
      m->getPCU()->Unpack(receivedSize);

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

    m->getPCU()->Begin();

    while(!updateRemoteVertices.empty())
    { 
      ent = updateRemoteVertices.front();
      //get remote copies and send updated mesh sizes
      m->getRemotes(ent,remotes);
      currentSize = apf::getScalar(size_iso,ent,0);
      for(apf::Copies::iterator iter=remotes.begin(); iter!=remotes.end();++iter)
      {
        m->getPCU()->Pack(iter->first, iter->second);
      }
      updateRemoteVertices.pop();
    }

    m->getPCU()->Send();
    //while remote copies are receiving
    while(m->getPCU()->Receive())
    {
      //unpack
      m->getPCU()->Unpack(ent);
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

//Eigenvalue routines for volume and surface meshes
//
  //function to get volume-based eigenvalues
//apf::Field* getLambdaMax(mesh,hessianField,)
void getLambdaMax(apf::Mesh* mesh,apf::Field* hessianField,apf::Field* lambdaMaxField)
{
  apf::MeshEntity* vert;
  apf::MeshIterator* it;
  it = mesh->begin(0);
  //apf::Field* lambdaMaxField = apf::createLagrangeField(mesh,"lambdaMax",apf::SCALAR,1);
  while( (vert = mesh->iterate(it)) )
  {
      apf::Matrix3x3 metric;
      apf::getMatrix(hessianField, vert, 0, metric);
      apf::Vector3 eigenVectors[3];
      double eigenValues[3];
      metric = (metric+apf::transpose(metric))/2.0;
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
  mesh->end(it);
}   

void getVolMaxPair(apf::Mesh* mesh,std::vector<std::vector<apf::MeshEntity*> > surfToStrandMap, apf::Field* lambdaMaxField,apf::Field* lambdaStrandMax, apf::Field* currentSize,double &lambda_max,double &h_lambdamax)
{
  apf::MeshIterator* it = mesh->begin(0);
  apf::MeshEntity* vert, *vertVol;
  int counter = 0;
  while( (vert = mesh->iterate(it)) )
  {
    //the goal of this loop is to identify lambda max global, find corresponding surface h-max global, and then set the surface lambda field
    int strandSize = surfToStrandMap[counter].size();
    for(int i=0; i<strandSize;i++)
    {
      vertVol = surfToStrandMap[counter][i];
      if(lambda_max < apf::getScalar(lambdaMaxField,vertVol,0))
      {
        lambda_max = apf::getScalar(lambdaMaxField,vertVol,0);
        h_lambdamax = apf::getScalar(currentSize,vert,0);
      }
    }

    //find local strand lambda max and set the field
    double lambda_local = 0.0;
    for(int i=0; i<strandSize;i++)
    {
      vertVol = surfToStrandMap[counter][i];
      if(lambda_local < apf::getScalar(lambdaMaxField,vertVol,0))
      {
        lambda_local = apf::getScalar(lambdaMaxField,vertVol,0);
      }
    }
    apf::setScalar(lambdaStrandMax,vert,0,lambda_local);

    counter++;
  }
  mesh->end(it);
}

void getSurfMaxPair(apf::Mesh* mesh,apf::Field* lambdaMaxField,apf::Field* currentSize,double &lambda_max,double &h_lambdamax)
{
  apf::MeshIterator* it = mesh->begin(0);
  apf::MeshEntity* vert;
  while( (vert = mesh->iterate(it)) )
  {
    if(lambda_max < apf::getScalar(lambdaMaxField,vert,0))
    {
      lambda_max = apf::getScalar(lambdaMaxField,vert,0);
      h_lambdamax = apf::getScalar(currentSize,vert,0);
    }
  }
  mesh->end(it);
}

void setSizeField(apf::Mesh* mesh, apf::Field* lambdaMaxField,apf::Field* sizeField,double lambda_max,double lambda_cutoff, double h_lambdamax,double h_global,double factor)
{
  apf::MeshIterator* it=mesh->begin(0);
  apf::MeshEntity* vert;
  int counter3 = 0;
  double h_special = -1;
  while( (vert = mesh->iterate(it)) )
  {    
    double h_v = h_special; //set default size as user specified input
    double lambda_vert = apf::getScalar(lambdaMaxField,vert,0);
    if(lambda_vert > lambda_cutoff)
    {
      h_v = sqrt(lambda_max/lambda_vert/factor)*h_lambdamax;
    }
    else
      counter3++;
    //h_global is user-specified max size
    if(h_v > h_global) //maximum value for size field
      h_v = h_global;
    apf::setScalar(sizeField,vert,0,h_v);
  }
  mesh->end(it);

}

static std::vector<bool> decodeBitFields(const char* bitFields)
{
  std::vector<bool> res;
  res.resize(strlen(bitFields));
  printf("length is %d\n", strlen(bitFields));
  for (int i = 0; i < strlen(bitFields); i++) {
    if (bitFields[i] == '0')
      res[i] = false;
    else if (bitFields[i] == '1')
      res[i] = true;
    else
      PCU_ALWAYS_ASSERT(0);
  }
  return res;
}

void isotropicIntersect(apf::Mesh* m, std::queue<apf::Field*> sizeFieldList, const char* bitFields, apf::Field* finalSizeField,apf::Field* finalChoiceField)
{
  std::vector<bool> userInput = decodeBitFields(bitFields);
  PCU_ALWAYS_ASSERT(userInput.size() == sizeFieldList.size());
  printf("HERE\n");
  apf::MeshEntity *vert;
  apf::MeshIterator *it = m->begin(0);

  apf::Field *field = sizeFieldList.front();
  apf::copyData(finalSizeField,field);
  sizeFieldList.pop();
  int choiceIdx = 1; //assumes the initial field was set to choice 0
  while(!sizeFieldList.empty())
  {
    field = sizeFieldList.front();
    while( (vert = m->iterate(it)) )
    {
      double value1 = apf::getScalar(finalSizeField,vert,0);
      double value2 = apf::getScalar(field,vert,0);
      if (!userInput[0]) value1 *= 1.e16;
      if (!userInput[1]) value2 *= 1.e16;
      double minValue = std::min(value1,value2);
      apf::setScalar(finalSizeField,vert,0,minValue);
      if(value1 > value2)
      {
        apf::setScalar(finalChoiceField,vert,0,choiceIdx);  
      }
    }
    sizeFieldList.pop();
    choiceIdx++;
  }
  m->end(it);
}


//

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
  double initialTime = pcu::Time();

  if (argc != 7) {
    if(0==PCUObj.get()->Self()) {
      std::cerr << "usage: " << argv[0]
        << " <data file .txt> <target field(s)> <strand size> <desired max size> <error reduction factor> <max refinement level>\n";
      std::cerr << "*target field(s) is a bit string to select which field(s) are used for error estimation\n";
      std::cerr << "1st bit --- pressure\n";
      std::cerr << "2nd bit --- skin friction\n";
    }
    return EXIT_FAILURE;
  }


  gmi_register_mesh();
  gmi_register_null();

  const char* dataFileName = argv[1];
  const char* bitFields = argv[2];
  const int strandSize = atoi(argv[3]);
  double h_global = atof(argv[4]);
  double factor = atof(argv[5]);
  const int maxLevel = atoi(argv[6]);

  gmi_cap_start();
  gmi_register_cap();

  // read the data into a vector array for now
  printf("\n---- Reading Strand Data. \n");
  std::vector<row> table = readTable(dataFileName);
  printf("---- Reading Strand Data: Done. \n");

  // create the mesh object (this one is CapStone underneath)
  printf("\n---- Creating Mesh Wrapper Object. \n");
  apf::Mesh2* mesh = apf::createMesh(m,g,PCUObj.get());
  printf("---- Creating Mesh Wrapper Object: Done. \n");

  // make the volume mesh (this one is MDS underneath)
  printf("\n---- Creating Volume Mesh. \n");
  // create an empty array with the correct sizes to hold the surface 2 strand map
  int surfVertCount = 
  std::vector<std::vector<apf::MeshEntity*> > surfToStrandMap(surfVertCount,
      std::vector<apf::MeshEntity*>(strandSize, NULL));
  apf::Mesh2* volMesh = createVolumeMesh(mesh, table, strandSize, surfToStrandMap);
  printf("---- Creating Volume Mesh: Done. \n");

  // get all fields and add them to the mesh
  printf("\n---- Adding Fields to Volume Mesh. \n");
  apf::Field* rhoField	   = addScalarField(volMesh, table, "rho", 3, strandSize);
  apf::Field* rho_uvwField = addVector3Field(volMesh, table, "rho_uvw", 4, 5, 6, strandSize);
  apf::Field* eField       = addScalarField(volMesh, table, "e", 7, strandSize);
  apf::Field* nuField      = addScalarField(volMesh, table, "nu", 8, strandSize);
  printf("---- Adding Fields to Volume Mesh: Done. \n");

  printf("\n---- Printing Volume/Strand Mesh Stats. \n");
  printf("number of mesh verts: %zu\n", volMesh->count(0));
  printf("number of mesh regions(hexes): %zu\n", volMesh->count(3));
  printf("---- Printing Volume/Strand Mesh Stats: Done. \n");

  double constructionTime = pcu::Time();
  std::cout<<"TIMER: Finished converting capstone mesh to volume mesh "<<constructionTime-initialTime<<std::endl;

  //Get Size Field for Adapt
  apf::Field* lambdaMaxField = apf::createLagrangeField(volMesh,"lambdaMax",apf::SCALAR,1);
  apf::Field* finalSizeField = apf::createLagrangeField(mesh,"final_size",apf::SCALAR,1);

  //get current size field
  apf::Field* currentSize = samSz::isoSize(mesh);

  apf::MeshIterator* it;
  apf::MeshEntity* vert;

  //transfer rho and e field to surface mesh
  std::cout<<"Get surface fields from volume\n";
  apf::Field* surfaceRhoField = apf::createLagrangeField(mesh,"surfaceDensity",apf::SCALAR,1);
  apf::Field* surfaceEField = apf::createLagrangeField(mesh,"surfaceE",apf::SCALAR,1);

  //getSpeed
  apf::Field* speedField = apf::createLagrangeField(volMesh,"speed",apf::SCALAR,1);
  it = volMesh->begin(0);
  while( (vert = volMesh->iterate(it)) )
  {
    double rhoVal = apf::getScalar(rhoField,vert,0);
    apf::Vector3 rhoVelocity;
    apf::getVector(rho_uvwField,vert,0,rhoVelocity);
    double speed = rhoVelocity.getLength()/rhoVal;
    apf::setScalar(speedField,vert,0,speed);
  }
  volMesh->end(it);
  //End getSpeed

  apf::Field* surfaceSpeedField = apf::createLagrangeField(mesh,"surface_speed",apf::SCALAR,1);
  it = mesh->begin(0);
  apf::MeshIterator* itVol = volMesh->begin(0);
  while( (vert = mesh->iterate(it)) )
  {
    apf::MeshEntity* vertVol = volMesh->iterate(itVol);
    double rhoVal = apf::getScalar(rhoField,vertVol,0);
    double eVal = apf::getScalar(eField,vertVol,0);
    double speedVal = apf::getScalar(speedField,vertVol,0);
    apf::setScalar(surfaceRhoField,vert,0,rhoVal);
    apf::setScalar(surfaceEField,vert,0,eVal);
    apf::setScalar(surfaceSpeedField,vert,0,speedVal);
  }
  mesh->end(it);
  volMesh->end(itVol);
  std::cout<<"Finished surface fields from volume\n";

  double getSurfaceTime = pcu::Time();
  std::cout<<"TIMER: Finished computing speed and transferring fields to surface "<<getSurfaceTime-constructionTime <<std::endl;

  //get the true pressure field
  apf::Field* pressureField = apf::createLagrangeField(volMesh,"pressure",apf::SCALAR,1);
  itVol = volMesh->begin(0);
  while( (vert = volMesh->iterate(itVol)) )
  {
    double rho_e = apf::getScalar(eField,vert,0);
    double rho = apf::getScalar(rhoField,vert,0);
    double speedVal = apf::getScalar(speedField,vert,0);
    double pressure_gamma = rho_e-0.5*rho*speedVal*speedVal;
    apf::setScalar(pressureField,vert,0,pressure_gamma);
  }
  volMesh->end(itVol);
  apf::Field* gradPressureField = apf::recoverGradientByVolume(pressureField);
  apf::Field* hessianPressureField = apf::recoverGradientByVolume(gradPressureField);

  //get eigenvalues in the volume mesh
  getLambdaMax(volMesh,hessianPressureField,lambdaMaxField);


  FieldOfInterest eBased;
  eBased.lambdaMaxField = lambdaMaxField;
  eBased.lambdaStrandMax = apf::createLagrangeField(mesh,"surf_lambda_strandMax",apf::SCALAR,1);
  eBased.sizeField = apf::createLagrangeField(mesh,"surface_size_e",apf::SCALAR,1);

  getVolMaxPair(mesh,surfToStrandMap,eBased.lambdaMaxField,eBased.lambdaStrandMax,currentSize,eBased.lambda_max,eBased.h_lambdamax);

  //set size field

  setSizeField(mesh,eBased.lambdaStrandMax,eBased.sizeField,eBased.lambda_max,eBased.lambda_cutoff(),eBased.h_lambdamax,h_global,factor);

  double getPressureTime = pcu::Time();
  std::cout<<"TIMER: Finished pressure size field "<<getPressureTime-getSurfaceTime <<std::endl;

  //get surface shear stress for adaptivity
  std::cout<<"Getting surface shear\n";
  apf::Field* surfaceShearField = apf::createLagrangeField(mesh,"skin",apf::SCALAR,1);
  it = mesh->begin(0);
  int vIDcounter = 0; //better way to do this would be with a numbering system
  //apf::Numbering* nVID = volMesh->findNumbering; //this would be for a volume mesh

  while( (vert = mesh->iterate(it)) )
  {
    apf::Vector3 pointVol;
    apf::Vector3 pointSurf;
    apf::MeshEntity* volVert = surfToStrandMap[vIDcounter][1];
    double speed = apf::getScalar(speedField,volVert,0);
    mesh->getPoint(vert,0,pointSurf);
    volMesh->getPoint(volVert,0,pointVol);
    apf::Vector3 differenceVec = pointVol-pointSurf;
    double shearStress = speed/differenceVec.getLength();
    apf::setScalar(surfaceShearField,vert,0,shearStress);
    vIDcounter++;
  }
  mesh->end(it);

  FieldOfInterest shearBased;
  shearBased.lambdaMaxField = apf::createLagrangeField(mesh,"surflambdaMax",apf::SCALAR,1);
  shearBased.sizeField = apf::createLagrangeField(mesh,"surface_size",apf::SCALAR,1);

  //get surface shear gradient
  std::cout<<"Reached gradshearfield\n";
  apf::Field* gradShearField = apf::recoverGradientByVolume(surfaceShearField);
  std::cout<<"Got gradshearfield\n";

  it=mesh->begin(0);
  while( (vert = mesh->iterate(it)) )
  {
    apf::Vector3 gradShear;
    apf::getVector(gradShearField,vert,0,gradShear);
    apf::setScalar(shearBased.lambdaMaxField,vert,0,gradShear.getLength());
  }
  mesh->end(it);

  std::cout<<"Got surf lambda max\n";
  getSurfMaxPair(mesh,shearBased.lambdaMaxField,currentSize,shearBased.lambda_max,shearBased.h_lambdamax);

  std::cout<<"Got surf lambda max pair\n";
  //set size field
  setSizeField(mesh,shearBased.lambdaMaxField,shearBased.sizeField,shearBased.lambda_max,shearBased.lambda_cutoff(),shearBased.h_lambdamax,h_global,factor);

  std::cout<<"set size field\n";

  double getShearTime = pcu::Time();
  std::cout<<"TIMER: Finished skin friction size field "<<getShearTime-getPressureTime<<std::endl;

  //Mesh Intersection
  std::queue<apf::Field*> sizeFieldList;
  sizeFieldList.push(eBased.sizeField);
  sizeFieldList.push(shearBased.sizeField);
  apf::Field* finalChoiceField = apf::createLagrangeField(mesh,"finalChoice",apf::SCALAR,1);
  it = mesh->begin(0);
  while( (vert = mesh->iterate(it)) )
  {
    apf::setScalar(finalChoiceField,vert,0,0);
  }

  isotropicIntersect(mesh,sizeFieldList,bitFields,finalSizeField,finalChoiceField);

  double getIntersectionTime = pcu::Time();
  std::cout<<"TIMER: get mesh intersection "<<getIntersectionTime-getShearTime<<std::endl;

  //grade mesh
  gradeMesh(mesh,finalSizeField);
  std::cout<<"Exiting surface shear\n";

  double getGradationTime = pcu::Time();
  std::cout<<"TIMER: get graded size field "<<getGradationTime-getIntersectionTime<<std::endl;

  //adjust refinement level
  adjustRefinementLevel(mesh,finalSizeField,currentSize,maxLevel);

  //Save initial meshes
  apf::writeVtkFiles("initialSurfaceMesh", mesh);
  apf::writeVtkFiles("volumeMeshSizeField", volMesh);

  //adapt the mesh
  ma::Input* in;
  apf::Field* adaptSize  = apf::createFieldOn(mesh, "adapt_size", apf::SCALAR);
  apf::copyData(adaptSize, finalSizeField);
  in = ma::configure(mesh, adaptSize);
  ma::validateInput(in);
  in->shouldRunPreZoltan = true;
  in->shouldRunMidZoltan = true;
  in->shouldRunPostZoltan = true;
  in->maximumImbalance = 1.05;
  in->maximumIterations = log2(factor) + 2;
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  /* in->shouldTransferToClosestPoint = true; */
  in->shouldForceAdaptation = true;
  /* in->debugFolder = "debug"; */
  ma::adaptVerbose(in, false);

  double adaptTime = pcu::Time();
  std::cout<<"TIMER: adaptMesh "<<adaptTime-getGradationTime<<std::endl;


  // add size distribution based on area
  std::vector<int> binCount;
  std::vector<double> binArea;
  computeSizeDistribution(mesh, factor, binCount, binArea);

  for (std::size_t i = 0; i < binCount.size(); i++) {
    printf("bin %lu's count/area is %d/%g \n", i, binCount[i], binArea[i]);
  }


  // write the adapted mesh to vtk
  apf::writeVtkFiles("adaptedMesh", mesh);
  // write the adapted mesh to smb
  // Note: This is for statistic measurements using
  // measureIsoStats.cc
  // The last argument is the name of the field to be used
  // in measureIsoStats.cc
  writeMdsMesh(mesh, "adaptedMesh.smb", "adapt_size");
  // write the adapted mesh to cre
  writeCre(cs, "adaptedMesh.cre");

  // clean up and exit calls
  //destroy all fields...
  for(int i =0;i<mesh->countFields();i++)
  {
    apf::destroyField(mesh->getField(i));
  } 
  for(int i =0;i<volMesh->countFields();i++)
  {
    apf::destroyField(volMesh->getField(i));
  } 

  gmi_cap_stop();
  }
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


apf::Mesh2* createVolumeMesh(apf::Mesh2* m, const std::vector<row> &t, int s, std::vector<std::vector<apf::MeshEntity*> > &surfToStrandMap)
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

  // create the surface to strand map
  it = vMesh->begin(0);
  while ( (e = vMesh->iterate(it)) ) {
    int layer = apf::getNumber(nLayer, e, 0, 0) - 1;
    int vid   = apf::getNumber(nVID, e, 0, 0);
    PCU_ALWAYS_ASSERT((layer >= 0) && (layer < s));
    PCU_ALWAYS_ASSERT((vid   >= 0) && (vid   < (int)m->count(0)));
    surfToStrandMap[vid][layer] = e;
  }
  vMesh->end(it);

  // return the volume mesh
  return vMesh;
}

void writeCre(CapstoneModule& cs, const std::string& fileName)
{
  GeometryDatabaseInterface    *gdbi = cs.get_geometry();
  MeshDatabaseInterface        *mdbi = cs.get_mesh();
  AppContext                   *ctx = cs.get_context();

  // Get the VTK writer.
  Writer *creWriter = get_writer(ctx, "Create Native Writer");
  if (!creWriter)
          error(HERE, ERR_GENERIC, "Could not find the CRE writer");

  IdMapping idmapping;
  std::vector<M_MModel> mmodels;
  M_GModel gmodel;
  M_MModel mmodel;
  gdbi->get_current_model(gmodel);
  mdbi->get_current_model(mmodel);
  mmodels.clear();
  mmodels.push_back(mmodel);
  creWriter->write(ctx, gmodel, mmodels, fileName.c_str(), idmapping);
}

void writeMdsMesh(apf::Mesh2* m, const char* name, const char* fieldName)
{
  // this would be the smb mesh with a field associated with it
  apf::Mesh2* outMesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 2, false);

  // add numbering to the surface mash
  apf::Numbering* n = apf::createNumbering(m, "v_num", apf::getLagrange(1) , 1);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  int count = 0;
  while ( (e = m->iterate(it)) ) {
    apf::number(n, e, 0, 0, count);
    count++;
  }
  m->end(it);


  int numVerts = m->count(0);
  apf::MeshEntity* vMap[numVerts];

  apf::MeshEntity* v;
  it = m->begin(0);
  while( (v = m->iterate(it)) ) {
    apf::Vector3 p;
    m->getPoint(v, 0, p);
    apf::MeshEntity* nv = outMesh->createVertex(0, p, apf::Vector3(0, 0, 0));
    int id = apf::getNumber(n, v, 0, 0);
    vMap[id] = nv;
  }
  m->end(it);

  it = m->begin(2);
  while( (e = m->iterate(it)) ) {
    apf::MeshEntity* vs[3];
    m->getDownward(e, 0, vs);
    apf::MeshEntity* triVs[3];
    triVs[0] = vMap[apf::getNumber(n, vs[0], 0, 0)];
    triVs[1] = vMap[apf::getNumber(n, vs[1], 0, 0)];
    triVs[2] = vMap[apf::getNumber(n, vs[2], 0, 0)];
    apf::buildElement(outMesh, 0, apf::Mesh::TRIANGLE, triVs);
  }
  m->end(it);

  outMesh->acceptChanges();
  apf::deriveMdsModel(outMesh);
  outMesh->verify();
  apf::verify(outMesh);

  apf::destroyNumbering(n);

  // get the field on input mesh
  apf::Field* inputField = m->findField(fieldName);
  PCU_ALWAYS_ASSERT(inputField);

  // create the corresponding field for the out mesh
  apf::Field* outputField  = apf::createFieldOn(outMesh, "sizes", apf::SCALAR);

  it = m->begin(0);
  apf::MeshIterator* outIt = outMesh->begin(0);
  while( (v = m->iterate(it)) ) {
    double value = apf::getScalar(inputField, v, 0);
    apf::MeshEntity* outV = outMesh->iterate(outIt);
    apf::setScalar(outputField, outV, 0, value);
  }
  m->end(it);
  outMesh->end(outIt);
  outMesh->writeNative(name);

  // clean up outMesh
  outMesh->destroyNative();
  apf::destroyMesh(outMesh);
}


void computeSizeDistribution(apf::Mesh2* m, int factor,
    std::vector<int>& binCount, std::vector<double>& binArea)
{
  // find the min length edge
  apf::MeshEntity* e;
  apf::MeshIterator* it;

  double minLength = 1.0e12;
  it = m->begin(1);
  while ( (e = m->iterate(it)) ) {
    double l = apf::measure(m, e);
    if (l < minLength)
      minLength = l;
  }
  m->end(it);

  double minArea = minLength * minLength / 2.;

  // initialize the arrays
  for (int i = 1; i <= factor; i*=2) {
    binCount.push_back(0);
    binArea.push_back(0.0);
  }


  it = m->begin(2);
  while ( (e = m->iterate(it)) ) {
    double a = apf::measure(m, e);
    int count = 0;
    for (int i = 1; i <= factor; i*=2) {
      if (a >= i * minArea && a < 4 * i * minArea) {
      	binCount[count] += 1;
      	binArea[count] += a;
      	break;
      }
      count++;
    }
  }
  m->end(it);
}

void adjustRefinementLevel(apf::Mesh2* m, apf::Field* finalSize,
    apf::Field* currentSize, int maxLevel)
{
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);

  while( (v = m->iterate(it)) ) {
    double currnetS = apf::getScalar(currentSize, v, 0);
    double finalS   = apf::getScalar(finalSize, v, 0);
    if (currnetS < finalS) continue;
    if (log2(currnetS/finalS) > maxLevel)
      finalS = currnetS / pow(2, maxLevel);
    apf::setScalar(finalSize, v, 0, finalS);
  }
  m->end(it);
}
