/******************************************************************************

  (c) 2004-2010 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#include "LIIPBMod.h"
#include <PCU.h>
#include <apfMesh.h>
using namespace std;

int LIIPBMod::run(apf::Mesh* m)
{

  int ipart;
  int numParts = 1;

  //check the number of nodes on each parts and computes the average 
  //nodes number
  int **numRgnSend, **numRgnRecv;
  double **RatioSend, **RatioRecv;
  int *numNP = new int[numParts];
  int *numRgn = new int[numParts];
  RatioSend = new double* [numParts];
  RatioRecv = new double* [numParts];
  numRgnSend = new int* [numParts];
  numRgnRecv = new int* [numParts];

  double *NP_ratio = new double [numParts];

  int numNPTot, numRgnTot;
  int numNPTotonPart=0, numRgnTotonPart=0;
  int *NpA, *NpB;

  if(!PCU_Comm_Self()){
      NpA = new int[PCU_Comm_Peers()*numParts];
      NpB = new int[PCU_Comm_Peers()*numParts];
  }

  numNP[0] = m->count(0);
  numRgn[0] = m->count(m->getDimension());
  numNPTotonPart += numNP[0];
  numRgnTotonPart += numRgn[0];
  
  MPI_Allreduce(&numNPTotonPart, &numNPTot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&numRgnTotonPart, &numRgnTot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(PCU_Comm_Self()==0)
      printf("numnpTot before Boundary Modification: %d\n", numNPTot);

  MPI_Gather(numNP, numParts, MPI_INT, NpB, numParts, MPI_INT, 0, MPI_COMM_WORLD); 

  double numNPAve = numNPTot/PCU_Comm_Peers()/numParts;
//  double numRgnAve = numRgnTot/P_size()/numParts;
//  double numRgnMax = numRgnAve*1.010; //region could be 2% higher than average
  int tag;

  pMeshDataId POtoMoveTag;
  POtoMoveTag = MD_newMeshDataId("POtoMove");

  map <int, int> *Neigbors = new map<int,int>[numParts];
  int *numNeigbor = new int [numParts];
  pmModel::PEIter peiter;
//  pmEntity* pe;

/*  for(ipart=0;ipart<numParts;ipart++) {
      peiter=pmModel::Instance()->peBegin();
      numNeigbor[ipart] = 0;
      for(;peiter!=pmModel::Instance()->peEnd();++peiter)
      {
          pe=(*peiter);
          int isOnPart = 0;
          for (pmEntity::BPIter bpiter=pe->bpBegin();bpiter!=pe->bpEnd();++bpiter)
              if(*bpiter == P_pid()*numParts+ipart){
                  isOnPart = 1;
                  break;
              }
          if(isOnPart) 
              for (pmEntity::BPIter bpiter=pe->bpBegin();bpiter!=pe->bpEnd();++bpiter)                  
                  if(*bpiter!=P_pid()*numParts+ipart && Neigbors[ipart].find(*bpiter)==Neigbors[ipart].end())
                      Neigbors[ipart][*bpiter] = numNeigbor[ipart]++;
      }
      
      RatioRecv[ipart] = new double[numNeigbor[ipart]];
      numRgnRecv[ipart] = new int[numNeigbor[ipart]];

      for(int i=0;i<numNeigbor[ipart];i++){
          RatioRecv[ipart][i]=0;
          numRgnRecv[ipart][i]=0;
      }
      NP_ratio[ipart] = numNP[ipart]/numNPAve;
  }
*/  
  for(ipart=0;ipart<numParts;ipart++) {
      numNeigbor[ipart] = 0;
      VIter viter = M_vertexIter(meshes[ipart]);
      pVertex vertex;
      while(vertex=VIter_next(viter)){
          if(EN_duplicate(vertex)){
              std::vector<std::pair<pEntity,int> > remoteCopies;
              std::vector<std::pair<pEntity,int> >::iterator vecIter;
              EN_getCopies(vertex,remoteCopies);
              for(vecIter=remoteCopies.begin();vecIter!=remoteCopies.end();vecIter++)
                  if(vecIter->second!=PCU_Comm_Self()*numParts+ipart && Neigbors[ipart].find(vecIter->second)==Neigbors[ipart].end())
                      Neigbors[ipart][vecIter->second]= numNeigbor[ipart]++;
          }
      }
      VIter_delete(viter);
      RatioRecv[ipart] = new double[numNeigbor[ipart]];
      numRgnRecv[ipart] = new int[numNeigbor[ipart]];
      
      for(int i=0;i<numNeigbor[ipart];i++){
          RatioRecv[ipart][i]=0;
          numRgnRecv[ipart][i]=0;
      }
      NP_ratio[ipart] = numNP[ipart]/numNPAve;
  }

  liipbmod_commuInt(numRgn, numRgnRecv, Neigbors, numParts);
  liipbmod_commuDouble(NP_ratio, RatioRecv, Neigbors, numParts);

  zoltanCB cb;
  cb.setAlgorithm(pmZoltanCallbacks::PARMETIS);

  // if the current part_mesh have high nodes number, move some regions to 
  // its neighbor
  for(Iter=0; Iter<IterMax;Iter++){
      map<pEntity, int*> POtoMove;
      map<pEntity, int*>::iterator poIter;
      int needtoupdate=0, needtoupdateglobal;
      map<pEntity, int>::iterator mapIter;
          
      for(ipart=0;ipart<numParts;ipart++) {

//          double numNodesToMove = (NP_ratio[ipart]-tolerance1)*numNPAve;
          int numNodesMarked = 0;
 
          if(NP_ratio[ipart]>tolerance1){ //high nodes number part
              //loop over the boundary nodes, seach the ones has only small number of 
              //adjcent regions on the current part
              pVertex vertex;
              VIter vertices = M_vertexIter(meshes[ipart]);
              while(vertex=VIter_next(vertices)){ // loop over the boundary nodes
                  if(EN_duplicate(vertex)){ //looking for the ones on the boundary
                      int numV_R = 0;
                      pPList vRegions=V_regions(vertex);
                      void* tmp=0;
                      pRegion region;
                      numV_R = PList_size(vRegions); //number of adjcent regions of the
                      // current part
                      if(numV_R<=numVregionMax){ //small number of  adjacent regions 
                                        //on the current part
                          std::vector<std::pair<pEntity,int> > remoteCopies;
                          std::vector<std::pair<pEntity,int> >::iterator vecIter;
                          EN_getCopies(vertex,remoteCopies);
                          if(remoteCopies.size()==1){ //has only one remote
                              // copy. i.e. belong to two parts
                              vecIter=remoteCopies.begin();
                              int rank= vecIter->second;
                              int index = Neigbors[ipart][rank];
                              //       if((NP_ratio-RatioRecv[rank]>tolerance2||RatioRecv[rank]<tolerance3) && numRgnRecv[rank]<numRgnMax){
                              //   if((NP_ratio-RatioRecv[rank]>tolerance2) && numRgnRecv[rank]<numRgnMax){
                              if((NP_ratio[ipart]-RatioRecv[ipart][index]>tolerance2||RatioRecv[ipart][index]<tolerance3)){
                                  needtoupdate = 1;
                                  numNodesMarked++;                  
                                  int pidtomove=rank;
                                  while(region=(pRegion)PList_next(vRegions,&tmp))
                                      if(!EN_getDataInt((pEntity)region,POtoMoveTag,&tag)){ 
                                          int *pid = new int[3];
                                          pid[0] = ipart;
                                          pid[1] = pidtomove;
                                          pid[2] = pidtomove/numParts;
                                          POtoMove[(pEntity)region] = pid;
                                          EN_attachDataInt((pEntity)region,POtoMoveTag,1);
                                      }
                              }
                          }
                      }                 
                      PList_delete(vRegions);
                  }
//           if(numNodesToMove<numNodesMarked)
//               break;          
              }
              VIter_delete(vertices);    
          }
          
      }  // loop over all the parts on this process

      MPI_Allreduce(&needtoupdate, &needtoupdateglobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

      for(poIter=POtoMove.begin();poIter!=POtoMove.end();poIter++)
          EN_deleteData((pEntity)(*poIter).first,POtoMoveTag);

      if(needtoupdateglobal){                

          if(PCU_Comm_Self()==0)
              printf("[%2d]migrationMeshEntities....Iter%d\n",PCU_Comm_Self(),Iter);
          migrateMeshEntities(meshes, POtoMove, cb);          

          //update the inforamtion for next iteration
          for(ipart=0;ipart<numParts;ipart++){
              numNP[ipart] = M_numVertices(meshes[ipart]);
              numRgn[ipart] = M_numRegions(meshes[ipart]);
              NP_ratio[ipart] = numNP[ipart]/numNPAve;
              
              //updated the neibgorhood
              Neigbors[ipart].clear();
              numNeigbor[ipart] = 0;
/*              peiter=pmModel::Instance()->peBegin();
              pmEntity* pe;
               
              for(;peiter!=pmModel::Instance()->peEnd();++peiter)
              {

                  pe=(*peiter);
                  int isOnPart = 0;
                  for (pmEntity::BPIter bpiter=pe->bpBegin();bpiter!=pe->bpEnd();++bpiter)
                      if(*bpiter == P_pid()*numParts+ipart){
                          isOnPart = 1;
                          break;
                      }
                  if(isOnPart) 
                      for (pmEntity::BPIter bpiter=pe->bpBegin();bpiter!=pe->bpEnd();++bpiter)                  
                          if(*bpiter!=P_pid()*numParts+ipart && Neigbors[ipart].find(*bpiter)==Neigbors[ipart].end())
                              Neigbors[ipart][*bpiter] = numNeigbor[ipart]++;
              }
*/
              VIter viter = M_vertexIter(meshes[ipart]);
              pVertex vertex;
              while(vertex=VIter_next(viter)){
                  if(EN_duplicate(vertex)){
                      std::vector<std::pair<pEntity,int> > remoteCopies;
                      std::vector<std::pair<pEntity,int> >::iterator vecIter;
                      EN_getCopies(vertex,remoteCopies);
                      for(vecIter=remoteCopies.begin();vecIter!=remoteCopies.end();vecIter++)
                          if(vecIter->second!=PCU_Comm_Self()*numParts+ipart && Neigbors[ipart].find(vecIter->second)==Neigbors[ipart].end())
                              Neigbors[ipart][vecIter->second]= numNeigbor[ipart]++;
                  }
              }
              VIter_delete(viter);
              
              
              delete [] RatioRecv[ipart];
              delete [] numRgnRecv[ipart];
              RatioRecv[ipart] = new double[numNeigbor[ipart]];
              numRgnRecv[ipart] = new int [numNeigbor[ipart]];
              
              for(int i=0;i<numNeigbor[ipart];i++){
                  RatioRecv[ipart][i]=0;
                  numRgnRecv[ipart][i]=0;
              }
              
          }
          liipbmod_commuInt(numRgn, numRgnRecv, Neigbors, numParts);
          liipbmod_commuDouble(NP_ratio,RatioRecv, Neigbors, numParts);
              
//        sprintf(vtkfile, "vtkfile%d_",Iter);
//        M_writeVTKFile(mesh, vtkfile,NP_ratio);
      }
      else
          break;
  }

  numNPTotonPart = 0;
  for(ipart=0;ipart<numParts;ipart++) {
//      printf("[%2d] numnp after Boundary Modification: %d\n", P_pid()*numParts+ipart,numNP[ipart]);
//      printf("[%2d] numRgn after Boundary Modification: %d\n", P_pid()*numParts+ipart,numRgn[ipart]);
      numNPTotonPart += numNP[ipart];
  }
  
  MPI_Allreduce(&numNPTotonPart, &numNPTot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(PCU_Comm_Self()==0)
      printf("numnpTot after Boundary Modification: %d\n", numNPTot);  

  MPI_Gather(numNP, numParts, MPI_INT, NpA, numParts, MPI_INT, 0, MPI_COMM_WORLD); 

//  char sys[64];
//  sprintf(sys,"new_meshes");
//  mkdir(sys,"00755");

//  PM_write2(meshes,"geom_.sms");
  
  if(PCU_Comm_Self()==0){
      for(ipart=0;ipart<numParts*PCU_Comm_Peers();ipart++)
          printf("[%2d] numnp before Boundary Modification: %d\n",ipart,NpB[ipart]);
      for(ipart=0;ipart<numParts*PCU_Comm_Peers();ipart++)
          printf("[%2d] numnp after Boundary Modification: %d\n",ipart,NpA[ipart]);
  }

  delete [] numRgnRecv;
  delete [] RatioRecv;
  delete [] numRgnSend;
  delete [] RatioSend;
  delete [] NP_ratio;
  delete [] numRgn;
  delete [] numNP;
  delete [] Neigbors;
  delete [] numNeigbor;
  if(PCU_Comm_Self()==0){
      delete [] NpA;
      delete [] NpB;
  }

  return 1;
}


void liipbmod_commuInt(int* ns, int **nr, map<int,int> *Neigbors, int numParts)
{
    int numNeigbor=0, ipart;
    for(ipart=0;ipart<numParts;ipart++)
        numNeigbor  += Neigbors[ipart].size();

    MPI_Request* req = (MPI_Request *)malloc(sizeof(MPI_Request)*2*numNeigbor ); 
    MPI_Status* stat = (MPI_Status *)malloc(sizeof(MPI_Status)*2*numNeigbor );

    int m, tag = 0;
    map <int, int>::iterator neigb;

    for(ipart=0;ipart<numParts;ipart++){
        for (m=0, neigb=Neigbors[ipart].begin(); neigb != Neigbors[ipart].end(); neigb++)  {
            int neigbor = neigb->first;
            int index = neigb->second;
            int sender = neigbor/numParts;
            tag = neigbor - sender*numParts;
            MPI_Irecv(nr[ipart]+index, 1, MPI_INT, sender, tag, MPI_COMM_WORLD, &req[m++]);
            MPI_Isend(&ns[ipart], 1, MPI_INT, sender, ipart, MPI_COMM_WORLD, &req[m++]);        
        }
    }
    
    MPI_Waitall(m, req, stat);
    
/*      delete [] req; */
/*      delete [] stat;     */
    free(req);
    free(stat);
}

void liipbmod_commuDouble(double *ns, double **nr, map<int,int> *Neigbors, int numParts)
{
    int numNeigbor=0, ipart;
    for(ipart=0;ipart<numParts;ipart++)
        numNeigbor  += Neigbors[ipart].size();

    MPI_Request* req = (MPI_Request *)malloc(sizeof(MPI_Request)*2*numNeigbor ); 
    MPI_Status* stat = (MPI_Status *)malloc(sizeof(MPI_Status)*2*numNeigbor );
    
/*      MPI_Request* req = new MPI_Request[2*(PMU_size()-1)]; */
/*      MPI_Status*  stat= new MPI_Status[2*(PMU_size()-1)]; */
    int m, tag = 0;
    map <int, int>::iterator neigb;

    for(ipart=0;ipart<numParts;ipart++){ 
        for (m=0, neigb=Neigbors[ipart].begin(); neigb != Neigbors[ipart].end(); neigb++) {
            int neigbor = neigb->first;
            int index = neigb->second;      
            int sender = neigbor/numParts;
            tag = neigbor - sender*numParts;
            MPI_Irecv(nr[ipart]+index, 1, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, &req[m++]);
            MPI_Isend(&ns, 1, MPI_DOUBLE, sender, ipart, MPI_COMM_WORLD, &req[m++]);
        }
    }
    MPI_Waitall(m, req, stat);
    
/*      delete [] req; */
/*      delete [] stat;     */
    free(req);
    free(stat);
}

