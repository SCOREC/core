#include <PCU.h>
#include "phOutput.h"
#include "phIO.h"
#include "phiotimer.h"
#include <cstdio>
#include <sstream>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>
#include <string.h>
#include <assert.h>
#include "phRestart.h"
#include <apf.h>
#include <apfField.h>
#include "apfShape.h"


#ifdef HAVE_CGNS
//
#include <cgns_io.h>
#include <pcgnslib.h>
//
#endif
typedef int lcorp_t;
#define NCORP_MPI_T MPI_INTEGER
static cgsize_t nDbgCG=50;
static int nDbgI=50;

namespace  {

template<class T>
MPI_Datatype getMpiType(T) {
  MPI_Datatype mpitype;
  //determine the type based on what is being sent
  if( std::is_same<T, double>::value ) {
    mpitype = MPI_DOUBLE;
  } else if ( std::is_same<T, int64_t>::value ) {
    mpitype = MPI_INT64_T;
  } else if ( std::is_same<T, int32_t>::value ) {
    mpitype = MPI_INT32_T;
  } else {
    assert(false);
    fprintf(stderr, "Unknown type in %s... exiting\n", __func__);
    exit(EXIT_FAILURE);
  }
  return mpitype;
}

}

// https://www.geeksforgeeks.org/sorting-array-according-another-array-using-pair-stl/
// Sort an array according to
// other using pair in STL.  Modified to be real-int pair (for distance matching) and in a separate routine, two integers (for idx sort by surfID)  
#include <bits/stdc++.h>
using namespace std;
 
// Function to sort integer array b[]
// according to the order defined by a[]
void pairsortDI(double a[], int b[], int n)
{
       pair<double, int> *pairt = new pair<double, int>[n];  // when done    delete pairt;
 
    // Storing the respective array
    // elements in pairs.
    for (int i = 0; i < n; i++)
    {
        pairt[i].first = a[i];
        pairt[i].second = b[i];
    }
 
    // Sorting the pair array.
    sort(pairt, pairt + n);
     
    // Modifying original arrays
    for (int i = 0; i < n; i++)
    {
        a[i] = pairt[i].first;
        b[i] = pairt[i].second;
    }
    delete [] pairt;
}

// Function to sort integer array b[]
// according to the order defined by a[]
void pairsort(int a[], int b[], int n)
{
    pair<double, int> *pairt = new pair<double, int>[n];
 
    // Storing the respective array
    // elements in pairs.
    for (int i = 0; i < n; i++)
    {
        pairt[i].first = a[i];
        pairt[i].second = b[i];
    }
 
    // Sorting the pair array.
    sort(pairt, pairt + n);
     
    // Modifying original arrays
    for (int i = 0; i < n; i++)
    {
        a[i] = pairt[i].first;
        b[i] = pairt[i].second;
    }
    delete [] pairt;
}
void pairDeal6sort(int a[], int b[], int n)
{
    int c[6]={0}; 
    for (int i = 0; i < n; i++) c[a[i]-1]++;  // count number each type in a pre-scan
    int** p = new int*[6];
    for (int i = 0; i < 6; i++) p[i]=new int[c[i]];
    int** idx = new int*[6];
    for (int i = 0; i < 6; i++) idx[i]=new int[c[i]];
    for (int i = 0; i < 6; i++) c[i]=0;
    int isrfM1;
    for (int i = 0; i < n; i++)
    {
       isrfM1=a[i]-1;
       p[isrfM1][c[isrfM1]]=b[i];
       idx[isrfM1][c[isrfM1]]=a[i];
       c[isrfM1]++;
    }
    int igc=0;
    for (int j = 0; j < 6; j++){
      for (int i = 0; i < c[j]; i++) {
        b[igc] = p[j][i];
        a[igc] = idx[j][i];
        igc++;
      }
    }
    assert(igc==n);
    for (int i = 0; i < 6; i++) delete [] p[i];
    for (int i = 0; i < 6; i++) delete [] idx[i];
    delete [] idx;
    delete [] p;
}


namespace ph {

static lcorp_t count_owned(int* ilwork, int nlwork,cgsize_t* ncorp_tmp, int num_nodes);
static lcorp_t count_local(int* ilwork, int nlwork,cgsize_t* ncorp_tmp, int num_nodes);


void commuInt(Output& o, cgsize_t* global)
{ // translating a commuInt out from PHASTA to c
  int numtask=o.arrays.ilwork[0];
  int itkbeg=0;
  int maxseg=1;
  int numseg;
  for (int itask=0; itask<numtask; ++itask) {
    numseg = o.arrays.ilwork[itkbeg + 4];
    maxseg=std::max(numseg,maxseg);
    itkbeg+=4+2*numseg;
  }
         
  int itag, iacc, iother, isgbeg;
  MPI_Datatype sevsegtype[numtask];
//first do what ctypes does for setup
  int* isbegin;
  int* lenseg;
  int* ioffset;
  isbegin = (int*) malloc(sizeof(int) * maxseg);
  lenseg  = (int*) malloc(sizeof(int) * maxseg);
  ioffset = (int*) malloc(sizeof(int) * maxseg);
// no VLA        MPI_Request  req[numtask];
// no VLA        MPI_Status stat[numtask];
  int maxtask=1000;
  assert(maxtask>=numtask);
  MPI_Request  req[maxtask];
  MPI_Status stat[maxtask];
  int maxfront=0;
  int lfront;
  itkbeg=0;
  for (int itask=0; itask<numtask; ++itask) {
    iacc   = o.arrays.ilwork[itkbeg + 2];
    numseg = o.arrays.ilwork[itkbeg + 4];
    // ctypes.f decrements itkbeg+3 by one for rank 0-based.  do that where used below
    lfront=0; 
    for(int is=0; is<numseg; ++is){
      isbegin[is]= o.arrays.ilwork[itkbeg+3+2*(is+1)] -1 ; // ilwork was created for 1-based
      lenseg[is]= o.arrays.ilwork[itkbeg+4+2*(is+1)];
      lfront+=lenseg[is];
    }
    maxfront=std::max(maxfront,lfront);
    for ( int iseg=0; iseg<numseg; ++iseg) ioffset[iseg] = isbegin[iseg] - isbegin[0];
    auto type = getMpiType( cgsize_t() );
    MPI_Type_indexed (numseg, lenseg, ioffset,type, &sevsegtype[itask]);
    MPI_Type_commit (&sevsegtype[itask]);
    itkbeg+=4+2*numseg;
  }
  free(isbegin);
  free(lenseg);
  free(ioffset);

  int m = 0; 
  itkbeg=0;
  for (int itask=0; itask<numtask; ++itask) {
    itag   = o.arrays.ilwork[itkbeg + 1];
    iacc   = o.arrays.ilwork[itkbeg + 2];
    iother = o.arrays.ilwork[itkbeg + 3] - 1; // MPI is 0 based but this was prepped wrong
    numseg = o.arrays.ilwork[itkbeg + 4]; /// not used
    isgbeg = o.arrays.ilwork[itkbeg + 5] - 1;
    if (iacc==0){ 
      MPI_Irecv(&global[isgbeg], 1, sevsegtype[itask],iother, itag, MPI_COMM_WORLD, &req[m]);
    } else {
      MPI_Isend(&global[isgbeg], 1, sevsegtype[itask],iother, itag, MPI_COMM_WORLD, &req[m]);
    }
    itkbeg+=4+2*numseg;
    m      = m + 1; 
  }
  MPI_Waitall(m, req, stat);
}

void gen_ncorp(Output& o )
{
  apf::Mesh* m = o.mesh;
  int i;
  lcorp_t nilwork = o.nlwork;
  int num_nodes=m->count(0);
  o.arrays.ncorp = new cgsize_t[num_nodes]; 
  lcorp_t owned;
  lcorp_t local;
  lcorp_t* owner_counts;
  cgsize_t  local_start_id;
  cgsize_t  gid;

  const int num_parts = PCU_Comm_Peers();
  const int part = PCU_Comm_Self() ;

  for(int i=0; i < num_nodes; i++) o.arrays.ncorp[i]=0;
  owned = count_owned(o.arrays.ilwork, nilwork, o.arrays.ncorp, num_nodes);
  local = count_local(o.arrays.ilwork, nilwork, o.arrays.ncorp, num_nodes);
  o.iownnodes = owned+local;
#ifdef PRINT_EVERYTHING
  printf("%d: %d local only nodes\n", part, local);
  printf("%d: %d owned nodes\n", part, owned);
#endif
  assert( owned <= num_nodes );
  assert( owned+local <= num_nodes );

  owner_counts = (lcorp_t*) malloc(sizeof(lcorp_t)*num_parts);
  for(int i=0; i < num_parts; i++) owner_counts[i]=0;
  owner_counts[part] = owned+local;
#ifdef PRINT_EVERYTHING
  for(i=0;i<num_parts;i++)
    printf("%d,", owner_counts[i]);
  printf("\n");
#endif
  MPI_Allgather(MPI_IN_PLACE, 1, NCORP_MPI_T, owner_counts, 1, NCORP_MPI_T, MPI_COMM_WORLD);
#ifdef PRINT_EVERYTHING
  for(i=0;i<num_parts;i++)
    printf("%d,", owner_counts[i]);
  printf("\n");
#endif
  local_start_id=0;
  for(i=0;i<part;i++) //TODO: MPI_Exscan()?
    local_start_id += owner_counts[i];
  local_start_id++; //Fortran numbering
  o.local_start_id = local_start_id;
#ifdef PRINT_EVERYTHING
  printf("%d: %d\n", part, local_start_id);
#endif
  gid = local_start_id;
  if(gid<0) printf("part,gid, %d %ld",part,gid);
  assert(gid>=0);
  for(i=0;i<num_nodes;i++) //assign owned node's numbers
  { //if shared, owned 1 //if shared, slave -1 //if local only, 0
    if(o.arrays.ncorp[i] == 1)
    {
      o.arrays.ncorp[i]=gid;
      assert(o.arrays.ncorp[i]>=0);
      gid++;
      continue;
    }
    if(o.arrays.ncorp[i] == 0)
    {
      o.arrays.ncorp[i] = gid;
      assert(o.arrays.ncorp[i]>=0);
      gid++;
      continue;
    }
    if(o.arrays.ncorp[i] == -1)
      o.arrays.ncorp[i] = 0; //commu() adds, so zero slaves
  } //char code[] = "out"; //int ione = 1;

  if(num_parts > 1) 
    commuInt(o, o.arrays.ncorp);
if(0==1) {
  for (int ipart=0; ipart<num_parts; ++ipart){
    if(part==ipart) { // my turn
      for (int inod=0; inod<num_nodes; ++inod) printf("%ld ", o.arrays.ncorp[inod]);
        printf(" \n");
    }
    PCU_Barrier();
  }
}
free(owner_counts);
}

static lcorp_t count_local(int* ilwork, int nlwork,cgsize_t* ncorp_tmp, int num_nodes)
{
	int i;
	lcorp_t num_local = 0;
	for(i=0;i<num_nodes;i++)
	{
		if(ncorp_tmp[i] == 0)
			num_local++; //nodes away from part boundary
		assert(!(ncorp_tmp[i] < -1 || ncorp_tmp[i] > 1));
	}
	return(num_local);
}

static lcorp_t count_owned(int* ilwork, int nlwork,cgsize_t* ncorp_tmp, int num_nodes)
{
	int numtask = ilwork[0];
	int itkbeg = 0; //task offset
	int owned = 0;
	int i,j,k;
	for(i=0;i<numtask;i++)
	{
		int itag = ilwork[itkbeg+1]; //mpi tag
		int iacc = ilwork[itkbeg+2]; //0 for slave, 1 for master
		assert(iacc >= 0 && iacc <= 1);
		int iother = ilwork[itkbeg+3]-1; //other rank (see ctypes.f for off by one)
		int numseg = ilwork[itkbeg+4]; //number of segments
		for(j=0;j<numseg;j++)
		{
			int isgbeg = ilwork[itkbeg+5+(j*2)]; //first idx of seg
			int lenseg = ilwork[itkbeg+6+(j*2)]; //length of seg
			assert(iacc == 0 || iacc == 1);
			if(iacc)
			{
				for(k=0;k<lenseg;k++)
				{
					if(ncorp_tmp[isgbeg-1+k] == 0)
						owned++;
					//make sure we're not both master and slave
					assert(ncorp_tmp[isgbeg-1+k] != -1);
					ncorp_tmp[isgbeg-1+k] = 1;
					assert(isgbeg-1+k < num_nodes);
				}
				assert(owned <= num_nodes);
			}
			else
			{
				for(k=0;k<lenseg;k++)
				{
					ncorp_tmp[isgbeg-1+k] = -1;
					assert(isgbeg-1+k < num_nodes);
				}
			}
			//ncorp_tmp init'd to 0
			//if shared, owned 1
			//if shared, slave -1
			//if local only, 0

			assert(itkbeg+6+(j*2) < nlwork);
		}
		itkbeg+= 4+2*numseg;
	}
	return(owned);
}


// renamed, retained but not yet updated
static std::string buildCGNSFileName(std::string timestep_or_dat)
{
  std::stringstream ss;
  ss << "chefO." << timestep_or_dat;
  return ss.str();
}

// update is only a transpose to match CNGS.  
void getInteriorConnectivityCGNS(Output& o, int block, cgsize_t* c)
{
  int nelem = o.blocks.interior.nElements[block];
  int nvert = o.blocks.interior.keys[block].nElementVertices;
  size_t i = 0;
  if(nvert==4) { //prepped for PHASTA's negative volume tets so flip second and third vertex
    for (int elem = 0; elem < nelem; ++elem){
        c[i++] = o.arrays.ncorp[o.arrays.ien[block][elem][0]];
        c[i++] = o.arrays.ncorp[o.arrays.ien[block][elem][2]];
        c[i++] = o.arrays.ncorp[o.arrays.ien[block][elem][1]];
        c[i++] = o.arrays.ncorp[o.arrays.ien[block][elem][3]];
     }
  } else {
    for (int elem = 0; elem < nelem; ++elem)
      for (int vert = 0; vert < nvert; ++vert)
        c[i++] = o.arrays.ncorp[o.arrays.ien[block][elem][vert]]; // input is 0-based,  out is  1-based do drop the +1
  }
  PCU_ALWAYS_ASSERT(i == nelem*nvert);
}

// update is both a transpose to match CNGS and reduction to only filling the first number of vertices on the boundary whereas PHASTA wanted full volume
void getBoundaryConnectivityCGNS(Output& o, int block, cgsize_t* c, double* eCenx, double* eCeny, double* eCenz)
{
  int nelem = o.blocks.boundary.nElements[block];
  int nvertVol = o.blocks.boundary.keys[block].nElementVertices;
  int nvert = o.blocks.boundary.keys[block].nBoundaryFaceEdges;
  int num_nodes=o.mesh->count(0);
  size_t i = 0;
  size_t phGnod = 0;
  std::vector<int> lnode={0,1,2,3}; // Standard pattern of first 4 (or 3)
  // PHASTA's use of volume elements has an lnode array that maps the surface nodes from the volume numbering.  We need it here too
  //  see hierarchic.f but note that is fortran numbering
  if(nvertVol==4) lnode={0, 2, 1, -1};             // tet is first three but opposite normal of others to go with neg volume
  if(nvertVol==5 && nvert==3) lnode={0, 4, 1, -1}; // pyramid tri is a fortran map of 1 5 2 
  if(nvertVol==6 && nvert==4) lnode={0, 3, 4, 1};  // wedge quad is a fortran map of 1 4 5 2
  for (int elem = 0; elem < nelem; ++elem){
    eCenx[elem]=0;
    eCeny[elem]=0;
    eCenz[elem]=0;
    for (int vert = 0; vert < nvert; ++vert){
      phGnod=o.arrays.ienb[block][elem][lnode[vert]]; //actually it is on-rank Global 
      c[i++] = o.arrays.ncorp[phGnod]; // PETSc truely global
      eCenx[elem]+=o.arrays.coordinates[0*num_nodes+phGnod];
      eCeny[elem]+=o.arrays.coordinates[1*num_nodes+phGnod];
      eCenz[elem]+=o.arrays.coordinates[2*num_nodes+phGnod];
    }
    eCenx[elem]/=nvert; // only necessary if you really want to use this as a correct Centroid rather than comparison
    eCeny[elem]/=nvert; // only necessary if you really want to use this as a correct Centroid rather than comparison
    eCenz[elem]/=nvert; // only necessary if you really want to use this as a correct Centroid rather than comparison
   }
  PCU_ALWAYS_ASSERT(i == nelem*nvert);
}

void getInterfaceConnectivityCGNS // not extended yet other than transpose
(
  Output& o,
  int block,
  apf::DynamicArray<int>& c
)
{
  int nelem = o.blocks.interface.nElements[block];
  int nvert0 = o.blocks.interface.keys[block].nElementVertices;
  int nvert1 = o.blocks.interface.keys[block].nElementVertices1;
  c.setSize(nelem * (nvert0 + nvert1));
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert0; ++vert)
      c[i++] = o.arrays.ncorp[o.arrays.ienif0[block][elem][vert]];
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert1; ++vert)
      c[i++] = o.arrays.ncorp[o.arrays.ienif1[block][elem][vert]];
  PCU_ALWAYS_ASSERT(i == c.getSize());
}

// renamed stripped down to just give srfID
void getNaturalBCCodesCGNS(Output& o, int block, int* codes)
{
  int nelem = o.blocks.boundary.nElements[block];
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
      codes[i++] = o.arrays.ibcb[block][elem][1]; //srfID is the second number so 1
// if we wanted we could use PHASTA's bit in coding in the first number to us attributes to set
// arbitrary combinations of BCs but leaving that out for now
}

void topoSwitch(char* Ename, int nvert,int F,int B,int Z,int *E, cgsize_t e_startg,cgsize_t e_endg)
{
  int Ep;
  switch(nvert){
    case 4:
      snprintf(Ename, 4, "Tet");
      if (cgp_section_write(F, B, Z, Ename, CGNS_ENUMV(TETRA_4), e_startg, e_endg, 0, &Ep))
         cgp_error_exit();
      break;
    case 5:
      snprintf(Ename, 4, "Pyr");
      if (cgp_section_write(F, B, Z, Ename, CGNS_ENUMV(PYRA_5), e_startg, e_endg, 0, &Ep))
          cgp_error_exit();
      break;
    case 6:
      snprintf(Ename, 4, "Wdg");
      if (cgp_section_write(F, B, Z, Ename, CGNS_ENUMV(PENTA_6), e_startg, e_endg, 0, &Ep))
          cgp_error_exit();
      break;
    case 8:
      snprintf(Ename, 4, "Hex");
      if (cgp_section_write(F, B, Z, Ename, CGNS_ENUMV(HEXA_8), e_startg, e_endg, 0, &Ep))
          cgp_error_exit();
      break;
  }
if(0==1)  printf("%d %d %d %s %ld %ld %d\n",F,B,Z,Ename,e_startg,e_endg,Ep);
  *E=Ep;
}
void topoSwitchB(char* Ename, int nvert,int F,int B,int Z,int *E, cgsize_t e_startg,cgsize_t e_endg)
{
  int Ep;
  switch(nvert){
    case 3:
      snprintf(Ename, 4, "Tri");
      if (cgp_section_write(F, B, Z, Ename, CGNS_ENUMV(TRI_3), e_startg, e_endg, 0, &Ep))
            cgp_error_exit();
      break;
    case 4:
      snprintf(Ename, 5, "Quad");
      if (cgp_section_write(F, B, Z, Ename, CGNS_ENUMV(QUAD_4), e_startg, e_endg, 0, &Ep))
            cgp_error_exit();
      break;
  }
if(0==1)  printf("%d %d %d %s %ld %ld %d\n",F,B,Z,Ename,e_startg,e_endg,Ep);
  *E=Ep;
}

// renamed and calling the renamed functions above with output writes now to CGNS
void writeBlocksCGNSinteror(int F,int B,int Z, Output& o, cgsize_t *e_written)
{
  int E,S,Fs,Fs2,Fsb,Fsb2;
  cgsize_t e_owned, e_start,e_end;
  cgsize_t e_startg,e_endg;
  const int num_parts = PCU_Comm_Peers();
  const cgsize_t num_parts_cg=num_parts;
  const int part = PCU_Comm_Self() ;
  const cgsize_t part_cg=part;
  // create a centered solution 
  if (cg_sol_write(F, B, Z, "CellRank", CGNS_ENUMV(CellCenter), &S) ||
      cgp_field_write(F, B, Z, S, CGNS_ENUMV(Integer), "CellRank", &Fs))
      cgp_error_exit();
  int nblki= o.blocks.interior.getSize();
  int nvMap[4] = {4,5,6,8};
  int nvC,nvert,nvAll,invC,iblkC;
  for (int i = 0; i < 4; ++i) { // check all topologies
    nvAll=0;
    nvC=nvMap[i];
    for (int j = 0; j < nblki; ++j) { // check all blocks
      BlockKey& k = o.blocks.interior.keys[j];
      nvert = o.blocks.interior.keys[j].nElementVertices;
      if(nvC==nvert) {
        invC=1;
        iblkC=j;
        break;
      } else invC=0;
    }
    nvAll= PCU_Add_Int(invC); // add across all
    cgsize_t* e=NULL; // = (cgsize_t *)malloc(nvC * e_owned * sizeof(cgsize_t));
    if(nvAll!=0) { //nvC present on at least 1 rank
      if(invC!=0){  //nvC present on my rank
         e_owned = o.blocks.interior.nElements[iblkC];
         e = (cgsize_t *)malloc(nvC * e_owned * sizeof(cgsize_t));
        getInteriorConnectivityCGNS(o, iblkC, e);
      }
      else e_owned=0;
      long safeArg=e_owned; // e_owned is cgsize_t which could be an 32 or 64 bit int
      e_endg=*e_written + PCU_Add_Long(safeArg); // end for the elements of this topology
      e_startg=1+*e_written; // start for the elements of this topology
      char Ename[5];
      topoSwitch(Ename, nvC,F,B,Z,&E,e_startg,e_endg);
      e_start=0;
      auto type = getMpiType( cgsize_t() );
      MPI_Exscan(&e_owned, &e_start, 1, type , MPI_SUM, MPI_COMM_WORLD);
      e_start+=1+*e_written; // my parts global element start 1-based
      e_end=e_start+e_owned-1;  // my parts global element stop 1-based
      if (cgp_elements_write_data(F, B, Z, E, e_start, e_end, e))
          cgp_error_exit();
      *e_written=e_endg;
      if(invC!=0) free(e);
      //     create the field data for this process 
      int* d = NULL;
      if(invC!=0){  //nvC present on my rank
        d = (int *)malloc(e_owned * sizeof(int));
        for (int n = 0; n < e_owned; n++) 
          d[n] = part;
        //     write the solution field data in parallel 
      }
      if (cgp_field_write_data(F, B, Z, S, Fs, &e_start, &e_end, d))
          cgp_error_exit();
      if(invC!=0) free(d);
      char UserDataName[11];
      snprintf(UserDataName, 11, "n%sOnRank", Ename);
      // create Helper array for number of elements on part of a given topology 
      if ( cg_goto(F, B, "Zone_t", 1, NULL) ||
           cg_gorel(F, "User Data", 0, NULL) ||
           cgp_array_write(UserDataName, CGNS_ENUMV(Integer), 1, &num_parts_cg, &Fs2))
           cgp_error_exit();
        // create the field data for this process 
      int nIelVec=e_owned;
      cgsize_t  partP1=part+1;
if(0==1)      printf("Intr, %s,  %d, %d, %d, %d \n", UserDataName, nIelVec,part,Fs,Fs2);
      if ( cgp_array_write_data(Fs2, &partP1, &partP1, &nIelVec))
           cgp_error_exit();

if(0==1){
    printf("interior cnn %s %d %ld %ld \n", Ename,part, e_start, e_end);
//    for (int ne=0; ne<std::min(nDbgCG,e_owned); ++ne) {
//      printf("%d, %d ", part,(ne+1));
//      for(int nv=0; nv< nvC; ++nv) printf("%ld ", e[ne*nvC+nv]);
//      printf("\n");
//    }
}

    } // end if ANY rank has this topology
  } // end of loop over ALL topologies
  PCU_Barrier();
}
void writeBlocksCGNSboundary(int F,int B,int Z, Output& o, int* srfID, int* srfIDidx, double** srfIDCen1, double** srfIDCen2, int* srfID1OnBlk, int* srfID2OnBlk, int* startBelBlk, int* endBelBlk, cgsize_t *e_written, cgsize_t *totBel, int *nStackedOnRank, int nblkb)
{
    int E,Fsb,Fsb2, nvC,nvert,nvAll,invC;
    const int num_parts = PCU_Comm_Peers();
    const cgsize_t num_parts_cg=num_parts;
    const int part = PCU_Comm_Self() ;
    const cgsize_t part_cg=part;
    cgsize_t e_owned, e_start,e_end, e_startg,e_endg;
    cgsize_t eVolElm=*e_written;
    cgsize_t e_belWritten=0;
    int nvMap[2] = {3,4};
    int iblkC[2];
    int estart[2];
    for (int i = 0; i < 2; ++i) { // check all topologies
      nvAll=0;
      invC=0;
      nvC=nvMap[i];
      int icountB=0;
      for (int j = 0; j < nblkb; ++j) { // check all blocks
        BlockKey& k = o.blocks.boundary.keys[j];
        nvert = o.blocks.boundary.keys[j].nBoundaryFaceEdges;
        if(nvert==nvC) {
           invC=1;
           iblkC[icountB]=j; // mark the block numbers (could be more than one) that have current topology
           icountB++;
        } 
      }
      nvAll= PCU_Add_Int(invC); // add across all
      cgsize_t* e=NULL; double* eCenx=NULL; double* eCeny=NULL; double* eCenz=NULL; 
      if(nvAll!=0) { //nvC present on at least 1 rank
        e_owned=0;
        if(invC!=0){  //nvC present on my rank
           for (int j = 0; j < icountB; ++j) { // combine blocks 
             estart[j]=e_owned;
             e_owned += o.blocks.boundary.nElements[iblkC[j]];
           }
           e = (cgsize_t *)malloc(nvC * e_owned * sizeof(cgsize_t));
           eCenx = (double *)malloc( e_owned * sizeof(double));
           eCeny = (double *)malloc( e_owned * sizeof(double));
           eCenz = (double *)malloc( e_owned * sizeof(double));
           for (int j = 0; j < icountB; ++j) {// combine blocks 
             getBoundaryConnectivityCGNS(o, iblkC[j], &e[estart[j]], &eCenx[estart[j]], 
                                                  &eCeny[estart[j]], &eCenz[estart[j]]); // stack repeated topologies 
             getNaturalBCCodesCGNS(o, iblkC[j], &srfID[e_belWritten+estart[j]]); // note e_owned counts all same topo
           }
           (*nStackedOnRank)++;  // no longer have nblkb blocks so count them as you stack them
        }
        e_startg=1+*e_written; // start for the elements of this topology
        long safeArg=e_owned; // e_owned is cgsize_t which could be an 32 or 64 bit int
        cgsize_t  numBelTP = PCU_Add_Long(safeArg); // number of elements of this topology
        e_endg=*e_written + numBelTP; // end for the elements of this topology
        char Ename[6];
        topoSwitchB(Ename, nvC,F,B,Z,&E,e_startg,e_endg);
        e_start=0;
        auto type = getMpiType( cgsize_t() );
        MPI_Exscan(&e_owned, &e_start, 1, type , MPI_SUM, MPI_COMM_WORLD);
        e_start+=1+*e_written; // my parts global element start 1-based
        e_end=e_start+e_owned-1;  // my parts global element stop 1-based
        // write the element connectivity in parallel 
        if (cgp_elements_write_data(F, B, Z, E, e_start, e_end, e))
            cgp_error_exit();
if(0==1)        printf("boundary cnn %d, %ld, %ld \n", part, e_start, e_end);
if(0==1){
    for (int ne=0; ne<std::min(nDbgCG,e_owned); ++ne) { printf("%d, %d ", part,(ne+1)); for(int nv=0; nv< nvert; ++nv) printf("%ld ", e[ne*nvert+nv]); printf("\n"); }
}
        int idx =((*nStackedOnRank) - 1);

        int S, Fs;

//        if (cg_sol_write(F, B, Z, "BoundaryCellRank", CGNS_ENUMV(FaceCenter), &S) ||
//            cgp_field_write(F, B, Z, S, CGNS_ENUMV(Integer), "BoundaryCellRank", &Fs))
//            cgp_error_exit();
/*        int ec0=cg_sol_write(F, B, Z, "BoundaryCellRank2", CGNS_ENUMV(CellCenter), &S);
        int ec1=cg_sol_write(F, B, Z, "BoundaryCellRank", CGNS_ENUMV(FaceCenter), &S);
        int ec2= cgp_field_write(F, B, Z, S, CGNS_ENUMV(Integer), "BoundaryCellRank", &Fs);
      //     create the field data for this process 
        int* d = NULL;
        if(invC!=0){  //nvC present on my rank
          d = (int *)malloc(e_owned * sizeof(int));
          for (int n = 0; n < e_owned; n++) 
            d[n] = part;
        }
        if (cgp_field_write_data(F, B, Z, S, Fs, &e_start, &e_end, d))
          cgp_error_exit();
        if(invC!=0) free(d);
*/
// more tricky to put something into PV to visualize the above (approximately) through vertex field
        int* dv = (int *)malloc(o.iownnodes * sizeof(int));
        cgsize_t start=o.local_start_id;
        cgsize_t end=start+o.iownnodes-1;
        if (cg_sol_write(F, B, Z, "BoundaryVertexRank", CGNS_ENUMV(Vertex), &S) ||
            cgp_field_write(F, B, Z, S, CGNS_ENUMV(Integer), "BoundaryVertexRank", &Fs))
            cgp_error_exit();
      //     create the field data for this process 
        for (int inode = 0; inode < o.iownnodes; ++inode) dv[inode]= -1;
        int idxl, en;
        if(invC!=0) { 
           for (int ibel = 0; ibel < e_owned; ++ibel){ 
             for (int ilv=0; ilv < nvC; ilv++) {
               en=e[ibel*nvC+ilv]; 
               if(en>=start && en<=end) {
                 dv[en-start]= part; 
               }
             }
           }
         }
         if (cgp_field_write_data(F, B, Z, S, Fs, &start, &end, dv))
            cgp_error_exit();
// more tricky to put srfID on nodes to see in PV (approximately) through vertex field
        if (cg_sol_write(F, B, Z, "BoundaryVertexSrfID", CGNS_ENUMV(Vertex), &S) ||
            cgp_field_write(F, B, Z, S, CGNS_ENUMV(Integer), "BoundaryVertexSrfID", &Fs))
            cgp_error_exit();
      //     create the field data for this process 
        for (int inode = 0; inode < o.iownnodes; ++inode) dv[inode]= -1;
        if(invC!=0) { 
           for (int ibel = 0; ibel < e_owned; ++ibel){ 
             for (int ilv=0; ilv < nvC; ilv++) {
               en=e[ibel*nvC+ilv]; 
               if(en>=start && en<=end) {
                 dv[en-start]= srfID[ibel]; 
//                 printf("%d %d %d %d %d %d %d\n ", part,ibel, ilv, en, en-start, dv[en-start], srfID[ibel]);
               }
             }
           }
         }
         if (cgp_field_write_data(F, B, Z, S, Fs, &start, &end, dv))
            cgp_error_exit();
        free(dv);
        if(invC!=0) {
          free(e);
//moved above          getNaturalBCCodesCGNS(o, iblkC[, &srfID[e_belWritten]);
          int icnt1=0; int icnt2=0;
          for (int ne=0; ne<e_owned; ++ne){ //count srfID =1 and 2 on this part,block
             if(srfID[e_belWritten+ne]==1) icnt1++; 
             if(srfID[e_belWritten+ne]==2) icnt2++;
          } 
          srfIDCen1[idx]=new double[icnt1*3]; // icnt{1,2} will equal to zero OFTEN
          srfIDCen2[idx]=new double[icnt2*3];
          srfID1OnBlk[idx]=icnt1;
          srfID2OnBlk[idx]=icnt2;
          int j1=0; int j2=0;
          for (int ne=0; ne<e_owned; ++ne){
             if(srfID[e_belWritten+ne]==1){  // copy in the toplogies centroids
               srfIDCen1[idx][j1++]=eCenx[ne];
               srfIDCen1[idx][j1++]=eCeny[ne];
               srfIDCen1[idx][j1++]=eCenz[ne];
             }
             if(srfID[e_belWritten+ne]==2) {
               srfIDCen2[idx][j2++]=eCenx[ne];
               srfIDCen2[idx][j2++]=eCeny[ne];
               srfIDCen2[idx][j2++]=eCenz[ne];
             }
          } 
          free(eCenx); free(eCeny); free(eCenz);
if(0==1){      printf("CentroidCounts %d %d %d %d %d %d %d %d\n",part,icnt1, icnt2, j1, j2, e_owned, srfID1OnBlk[i],srfID2OnBlk[i]);}
          for (int j = 0; j < (int) e_owned; ++j) srfIDidx[e_belWritten+j]=e_start+j;
          startBelBlk[idx]=e_start-eVolElm; // provides start point for each block in srfID 
          endBelBlk[idx]=e_end-eVolElm; // provides end point for each block in srfID
        }
        *e_written=e_endg;
        e_belWritten+=e_owned; // this is tracking written by this rank as we unpack srfID later
        char UserDataName[12]; snprintf(UserDataName, 13, "n%sOnRank", Ename);
        if ( cg_goto(F, B, "Zone_t", 1, NULL) ||
             cg_gorel(F, "User Data", 0, NULL) ||
             cgp_array_write(UserDataName, CGNS_ENUMV(Integer), 1, &num_parts_cg, &Fsb2))
             cgp_error_exit();
if(0==1)        printf("Bndy %s, %ld, %d, %d \n", UserDataName, e_owned, part,Fsb2);
        cgsize_t partP1=part+1;
        if (cgp_array_write_data(Fsb2, &partP1, &partP1, &e_owned))
            cgp_error_exit();
      } // at least one part has this topo
    } // loop over both topos
    *totBel = *e_written-eVolElm;
}
void writeCGNS_UserData_srfID(int F,int B, Output& o, int* srfID,  int* startBelBlk, int *endBelBlk, cgsize_t *e_written, cgsize_t *totBel, cgsize_t *eVolElm, int nStackedOnRank)
{
    cgsize_t e_owned, e_start,e_end;
    int Fsb;
    // setup User Data for boundary faces 
    if ( cg_goto(F, B, "Zone_t", 1, NULL) ||
         cg_gorel(F, "User Data", 0, NULL) ||
         cgp_array_write("srfID", CGNS_ENUMV(Integer), 1,totBel, &Fsb)) 
         cgp_error_exit();
    // write the user data for this process 
    int nblkb = o.blocks.boundary.getSize(); 
    int nvMap[2] = {3,4};
    int iblkC[2];
    int estart[2];
    int nvC,nvert,nvAll,invC;
    for (int i = 0; i < 2; ++i) { // at most, 2 blocks (assumed as we have (untested) collapsed/stacked blocks withe same number of verts)
      nvAll=0;
      nvC=nvMap[i];
      invC=0;
      for (int j = 0; j < nblkb; ++j) { // check all blocks
        BlockKey& k = o.blocks.boundary.keys[j];
        nvert = o.blocks.boundary.keys[j].nBoundaryFaceEdges;
        if(nvert==nvC) {
           invC=1;
        } 
      }
      nvAll= PCU_Add_Int(invC); // add across all
      if(nvAll!=0) { //nvC present on at least 1 rank
        e_owned=0;
        int e_startB=0;
        if(i<nStackedOnRank) {
          e_owned=endBelBlk[i]-startBelBlk[i]+1;
          e_startB=startBelBlk[i]; // srfID is only for bel....matches linear order with eVolElm offset from 
        }                                 // bel# that starts from last volume element
        e_start=0;
        auto type = getMpiType( cgsize_t() );
        MPI_Exscan(&e_owned, &e_start, 1, type , MPI_SUM, MPI_COMM_WORLD);
        e_start+=1+*e_written; // my parts global element start 1-based
        e_end=e_start+e_owned-1;  // my parts global element stop 1-based
if(0==1)      printf("BndyUserData %s, %ld, %ld, %ld,  %d, %d %d \n", "srfID", e_start, e_end, e_owned, i, e_startB,*totBel);
        if(i<nStackedOnRank){ 
          if (cgp_array_write_data(Fsb, &e_start, &e_end, &srfID[e_startB]))
          cgp_error_exit();
        } else { 
          if (cgp_array_write_data(Fsb, &e_start, &e_end, NULL))
          cgp_error_exit();
        }
        long safeArg=e_owned; // is cgsize_t which could be an 32 or 64 bit int
        *e_written += PCU_Add_Long(safeArg); // number of elements of this topology
      }
    }

}
void sortID1andID2(double* srfID1GCen,double* srfID2GCen, int nmatchFace, int* imapD1, int*imapD2)
{
    int* imapD2v = (int *)malloc( nmatchFace * sizeof(int));
    double* srfID1distSq = (double *)malloc( nmatchFace * sizeof(double));
    double* srfID2distSq = (double *)malloc( nmatchFace * sizeof(double));
    const int part = PCU_Comm_Self() ;
    double xc=0.1; // true cubes with uniform meshes and other symmetries set up ties  (good for debugging/verifying that dumb search backup works)
    for (int i = 0; i < nmatchFace; ++i) {
      srfID1distSq[i]=(srfID1GCen[i*3+0]-xc)*(srfID1GCen[i*3+0]-xc) 
                   +srfID1GCen[i*3+1]*srfID1GCen[i*3+1] 
                   +srfID1GCen[i*3+2]*srfID1GCen[i*3+2]; // periodicity always in Z then could be omitted
      srfID2distSq[i]=(srfID2GCen[i*3+0]-xc)*(srfID2GCen[i*3+0]-xc) 
                   +srfID2GCen[i*3+1]*srfID2GCen[i*3+1] 
                   +srfID2GCen[i*3+2]*srfID2GCen[i*3+2]; // periodicity always in Z then could be omitted
      imapD1[i]=i;
      imapD2[i]=i;
    }
    if(0==1){ if(part==0) {
      printf(" srfID1dist GLOBAL B "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%f ", srfID1distSq[is]); printf("\n");
      printf(" imapD1 GLOBAL B     "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%d ", imapD1[is]); printf("\n"); }
    }
    if(0==1){ if(part==0) {
      printf(" srfID2dist GLOBAL B "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%f ", srfID2distSq[is]); printf("\n");
      printf(" imapD2 GLOBAL B     "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%d ", imapD2[is]); printf("\n"); }
    }
    pairsortDI(srfID1distSq,imapD1,nmatchFace); // imapD1 puts elements with srfID=1 in order of increasing distance from pt 0.1, 0 0 
    pairsortDI(srfID2distSq,imapD2,nmatchFace); // imapD1 puts elements with srfID=2 in order of increasing distance from pt 0.1, 0 0 

    if(0==1){ if(part==0) {
      printf(" srfID1dist GLOBAL "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%f ", srfID1distSq[is]); printf("\n");
      printf(" imapD1 GLOBAL     "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%d ", imapD1[is]); printf("\n"); }
    }
    if(0==1){ if(part==0) {
      printf(" srfID2dist GLOBAL "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%f ", srfID2distSq[is]); printf("\n");
      printf(" imapD2 GLOBAL     "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%d ", imapD2[is]); printf("\n"); }
    }

    double tol2=1.0e-14;
    int jclosest, iclose1, iclose2;
    double d1,d2,vDistSq,vDSmin;
    int DistFails=0;
    for (int i = 0; i < nmatchFace; ++i) {
        iclose1=imapD1[i];
        iclose2=imapD2[i];
        d1=srfID1GCen[(iclose1)*3+0]-srfID2GCen[(iclose2)*3+0];
        d2=srfID1GCen[(iclose1)*3+1]-srfID2GCen[(iclose2)*3+1];
        vDistSq= d1*d1+d2*d2;
        if(vDistSq < tol2) {
           imapD2v[i]=imapD2[i];
        } else {// Centroid for i did not match-> search list srfID=2 list to find true match
          vDSmin=vDistSq;
          DistFails++;
          for (int j = 0; j < nmatchFace; ++j) {   // if this turns out to be taken a lot then it could be narrowed e.g. j=max(0,i-50), j< i+min(matchFace,i+50),
            iclose2=imapD2[j];
            d1=srfID1GCen[(iclose1)*3+0]-srfID2GCen[(iclose2)*3+0];
            d2=srfID1GCen[(iclose1)*3+1]-srfID2GCen[(iclose2)*3+1];
            vDistSq= d1*d1+d2*d2;
            if(vDistSq<vDSmin) {
              vDSmin=vDistSq;
              jclosest=j;
              if(vDistSq<tol2) break;
            } 
          }
          assert(vDistSq<tol2);
          imapD2v[i]=imapD2[jclosest];
        } 
    } 
    for (int i = 0; i < nmatchFace; ++i) imapD2[i]=imapD2v[i];
    if(0==1&&part==0) {
      printf("Number of Distance Failures=%d\n ",DistFails);
      printf(" srfID1dist GLOBAL "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%f ", srfID1distSq[is]); printf("\n");
      printf(" imapD1 GLOBAL     "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%d ", imapD1[is]); printf("\n"); 
      printf(" srfID2dist GLOBAL "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%f ", srfID2distSq[is]); printf("\n");
      printf(" imapD2 GLOBAL     "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%d ", imapD2[is]); printf("\n"); }
    free(srfID1distSq); free(srfID2distSq);
    free(imapD2v);
}
void GatherCentroid(double** srfIDCen,int* srfIDOnBlk, double** srfIDGCen, int *nmatchFace, int nStackedOnRank)
{
// stack  connectivities on rank before gather (should preserve order)
    const int num_parts = PCU_Comm_Peers();
    const int part = PCU_Comm_Self() ;
    int* rcounts = (int *)malloc( num_parts * sizeof(int));
    int* displs = (int *)malloc( num_parts * sizeof(int));
    int numSurfIDOnRank=0;
    for (int i = 0; i < nStackedOnRank; ++i) numSurfIDOnRank+=srfIDOnBlk[i];
    double* srfIDCenAllBlocks = (double *)malloc(numSurfIDOnRank*3 * sizeof(double));
    int k1=0;
    for (int i = 0; i < nStackedOnRank; ++i) 
      for (int j = 0; j < srfIDOnBlk[i]*3; ++j) srfIDCenAllBlocks[k1++]=srfIDCen[i][j];
    int ncon=numSurfIDOnRank*3;
    auto type_i = getMpiType( int() );
    MPI_Gather(&ncon,1,type_i,rcounts,1,type_i,0,MPI_COMM_WORLD);
    displs[0]=0;
    for (int i = 1; i < num_parts; ++i) displs[i]=displs[i-1]+rcounts[i-1]; 
    int GsrfIDcnt=displs[num_parts-1]+rcounts[num_parts-1];
    *nmatchFace=GsrfIDcnt/3;
    if(part==0) *srfIDGCen = (double *)malloc( GsrfIDcnt * sizeof(double));
if(0==1){ printf("displs1 ");for(int ip=0; ip< num_parts; ++ip) printf("% ld ", displs[ip]); printf("\n"); }
    auto type_d = getMpiType( double() );
    MPI_Gatherv(srfIDCenAllBlocks,ncon,type_d,*srfIDGCen,rcounts,displs,type_d,0, MPI_COMM_WORLD);
    free(srfIDCenAllBlocks); free(rcounts); free(displs);
}
void Allgather2IntAndSort(int* srfID, int* srfIDidx,Output& o,int* srfIDG, int* srfIDGidx, int totOnRankBel)
{
    const int part = PCU_Comm_Self() ;
    const int num_parts = PCU_Comm_Peers();
    const cgsize_t part_cg=part;
    int* rcounts = (int *)malloc( num_parts * sizeof(int));
    int* displs = (int *)malloc( num_parts * sizeof(int));
    auto type_i = getMpiType( int() );
    MPI_Allgather(&totOnRankBel,1,type_i,rcounts,1,type_i,MPI_COMM_WORLD);
    displs[0]=0;
    for (int i = 1; i < num_parts; ++i) displs[i]=displs[i-1]+rcounts[i-1]; 
    int totBel=displs[num_parts-1]+rcounts[num_parts-1];
if(0==1){ for(int ip=0; ip< num_parts; ++ip) printf("%ld ", displs[ip]); printf("\n"); }
    MPI_Allgatherv(srfID,totOnRankBel,type_i,srfIDG,rcounts,displs,type_i,MPI_COMM_WORLD);
    MPI_Allgatherv(srfIDidx,totOnRankBel,type_i,srfIDGidx,rcounts,displs,type_i,MPI_COMM_WORLD);
    free(rcounts); free(displs);
if(0==1){ if(part==0) {
    printf(" srfID GLOBAL    "); for(int is=0; is< std::min(nDbgI,totBel); ++is)  printf("%d ", srfIDG[is]); printf("\n");
    printf(" srfIDidx GLOBAL "); for(int is=0; is< std::min(nDbgI,totBel); ++is)  printf("%d ", srfIDGidx[is]); printf("\n"); }
    printf("rank %d ",part); printf(" srfID on Part "); for(int is=0; is< std::min(nDbgI,totOnRankBel); ++is)  printf("%d ", srfID[is]); printf("\n");
    printf(" srfIDidx on Part "); for(int is=0; is< std::min(nDbgI,totOnRankBel); ++is)  printf("%d ", srfIDidx[is]); printf("\n"); }
//    pairsort(srfIDG,srfIDGidx,totBel);
    pairDeal6sort(srfIDG,srfIDGidx,totBel);
if(0==1){ if(part==0) {
    printf(" srfID GLOBAL    "); for(int is=0; is< std::min(nDbgI,totBel); ++is)  printf("%d ", srfIDG[is]); printf("\n");
    printf(" srfIDidx GLOBAL "); for(int is=0; is< std::min(nDbgI,totBel); ++is)  printf("%d ", srfIDGidx[is]); printf("\n"); } }
}
void writeCGNSboundary(int F,int B,int Z, Output& o, int* srfID, int* srfIDidx, double** srfIDCen1, double** srfIDCen2, int* srfID1OnBlk, int* srfID2OnBlk, int* startBelBlk, int *endBelBlk, cgsize_t *e_written, int totOnRankBel, cgsize_t *totBel, int nStackedOnRank)
{
// srfID is for ALL Boundary faces
    const int num_parts = PCU_Comm_Peers();
    const cgsize_t num_parts_cg=num_parts;
    const int part = PCU_Comm_Self() ;
    const cgsize_t part_cg=part;
    cgsize_t e_owned, e_start,e_end;
    int Fsb;
    cgsize_t  eVolElm = *e_written-*totBel;
    *e_written=0; //recycling  eVolElm holds 
    writeCGNS_UserData_srfID(F,B, o, srfID,  startBelBlk, endBelBlk, e_written, totBel, &eVolElm, nStackedOnRank );
    double* srfID1GCen; 
    double* srfID2GCen; 
    int nmatchFace1,nmatchFace;
    GatherCentroid(srfIDCen1,srfID1OnBlk,&srfID1GCen,&nmatchFace1, nStackedOnRank);
    GatherCentroid(srfIDCen2,srfID2OnBlk,&srfID2GCen,&nmatchFace, nStackedOnRank);
if(0==1)    if(part==0)  printf("matchface %d, %d\n", nmatchFace1, nmatchFace);
    if(part==0) assert(nmatchFace1==nmatchFace);
//   compute the translation while we still have ordered centroids data  Assuming Translation = donor minus periodic but documents unclear
    double  TranslationD[3];
    if (part==0){  TranslationD[0]=srfID2GCen[0]-srfID1GCen[0]; TranslationD[1]=srfID2GCen[1]-srfID1GCen[1];TranslationD[2]=srfID2GCen[2]-srfID1GCen[2];}
if(0==1){  printf("%d part srfID 1 xc ",part); for(int ip=0; ip< std::min(nDbgI,nmatchFace); ++ip) printf("%f ", srfID1GCen[ip*3+0]); printf("\n"); }
if(0==1){  printf("%d part srfID 1 yc ",part); for(int ip=0; ip< std::min(nDbgI,nmatchFace); ++ip) printf("%f ", srfID1GCen[ip*3+1]); printf("\n"); }
if(0==1){  printf("%d part srfID 1 zc ",part); for(int ip=0; ip< std::min(nDbgI,nmatchFace); ++ip) printf("%f ", srfID1GCen[ip*3+2]); printf("\n"); }
       PCU_Barrier();
if(0==1){  printf("%d part srfID 2 xc ",part); for(int ip=0; ip< std::min(nDbgI,nmatchFace); ++ip) printf("%f ", srfID2GCen[ip*3+0]); printf("\n"); }
if(0==1){  printf("%d part srfID 2 yc ",part); for(int ip=0; ip< std::min(nDbgI,nmatchFace); ++ip) printf("%f ", srfID2GCen[ip*3+1]); printf("\n"); }
if(0==1){  printf("%d part srfID 2 zc ",part); for(int ip=0; ip< std::min(nDbgI,nmatchFace); ++ip) printf("%f ", srfID2GCen[ip*3+2]); printf("\n"); }
    auto type_i = getMpiType( int() );
    MPI_Bcast(&nmatchFace,1,type_i,0, MPI_COMM_WORLD);
    int* imapD1 = (int *)malloc( nmatchFace * sizeof(int));
    int* imapD2 = (int *)malloc( nmatchFace * sizeof(int));
    if(part==0) sortID1andID2(srfID1GCen,srfID2GCen,nmatchFace, imapD1, imapD2);
    PCU_Barrier();
    MPI_Bcast(imapD1,nmatchFace,type_i,0, MPI_COMM_WORLD);
    MPI_Bcast(imapD2,nmatchFace,type_i,0, MPI_COMM_WORLD);
    auto type_d = getMpiType( double() );
    MPI_Bcast(TranslationD,3,type_d,0, MPI_COMM_WORLD);
    if(part==0) {free(srfID1GCen); free(srfID2GCen);}
// ZonalBC data 
    int* srfIDG = (int *)malloc( *totBel * sizeof(int));
    int* srfIDGidx = (int *)malloc( *totBel * sizeof(int));
    cgsize_t* donor2 = (cgsize_t *)malloc(nmatchFace * sizeof(cgsize_t));
    cgsize_t* periodic1 = (cgsize_t *)malloc(nmatchFace * sizeof(cgsize_t));
    Allgather2IntAndSort(srfID, srfIDidx,o,srfIDG, srfIDGidx,totOnRankBel);
    int BC_scan=0;
    cgsize_t* eBC = (cgsize_t *)malloc(*totBel * sizeof(cgsize_t));
    for (int BCid = 1; BCid < 7; BCid++) {
      int imatch=0;
      for (int ib = BC_scan; ib < *totBel; ib++) {
        if(srfIDG[ib]==BCid){
          eBC[imatch]=srfIDGidx[BC_scan];
          BC_scan++;
          imatch++;
        } else  break;
      }
//reorder SurfID = 1 and 2 using idmapD{1,2} based on distance to support periodicity 
      if(BCid==1) {
        for (int i = 0; i < nmatchFace; i++) periodic1[i]=eBC[imapD1[i]];
        for (int i = 0; i < nmatchFace; i++) eBC[i]=periodic1[i];
if(0==1&&part==1){ printf(" srfIDidx 1 "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%d ", eBC[is]); printf("\n"); }
      }       
      if(BCid==2) {
        for (int i = 0; i < nmatchFace; i++) donor2[i]=eBC[imapD2[i]];
        for (int i = 0; i < nmatchFace; i++) eBC[i]=donor2[i];
if(0==1&&part==1){ printf(" srfIDidx 2 "); for(int is=0; is< std::min(nDbgI,nmatchFace); ++is)  printf("%d ", eBC[is]); printf("\n"); }
      }       
if(0==1) {
      printf(" srfID =%d    ",BCid); for(int is=0; is< std::min(nDbgI,imatch); ++is)  printf("%d ", eBC[is]); printf("\n");
}
      int BC_index;
      char BC_name[33];
      snprintf(BC_name, 33, "SurfID_%d", BCid );
      if(cg_boco_write(F, B, Z, BC_name, CGNS_ENUMV(BCTypeUserDefined), CGNS_ENUMV(PointList), imatch, eBC,  &BC_index))
         cg_error_exit();
      if(cg_goto(F, B, "Zone_t", 1, "ZoneBC_t", 1, "BC_t", BC_index, "end")) cg_error_exit();;
      if(cg_gridlocation_write(CGNS_ENUMV(FaceCenter))) cg_error_exit();
    }
    int cgconn;
    if (cg_conn_write(F, B, Z, "Periodic Connectivity",
          CGNS_ENUMV(FaceCenter), CGNS_ENUMV(Abutting1to1),
          CGNS_ENUMV(PointList), nmatchFace, periodic1, "Zone",
          CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
          CGNS_ENUMV(Integer), nmatchFace, donor2, &cgconn)) cgp_error_exit();
    const float  RotationCenter[3]={0};
    const float  RotationAngle[3]={0};
    const float  Translation[3]={(float)TranslationD[0],(float)TranslationD[1],(float)TranslationD[2]};

    if (cg_conn_periodic_write(F, B, Z, cgconn, RotationCenter, RotationAngle,  Translation)) cgp_error_exit();
    free(imapD1); free(imapD2);
    free(eBC); free(srfIDG); free(srfIDGidx); free(donor2); free(periodic1);
} 
void CGNS_NodalSolution(int F,int B,int Z, Output& o)
{
  // create a nodal solution 
  char fieldName[12];
  snprintf(fieldName, 13, "solution");
if(0==1)  printf("solution=%s",fieldName);
  double* data;
  int size, S,Q;
  detachField(o.mesh, fieldName, data, size);
  assert(size==5);

//     create the field data for this process 
  double* p = (double *)malloc(o.iownnodes * sizeof(double));
  double* u = (double *)malloc(o.iownnodes * sizeof(double));
  double* v = (double *)malloc(o.iownnodes * sizeof(double));
  double* w = (double *)malloc(o.iownnodes * sizeof(double));
  double* T = (double *)malloc(o.iownnodes * sizeof(double));
  int icount=0;
  int num_nodes=o.mesh->count(0);
  cgsize_t gnod,start,end;
  start=o.local_start_id;
  end=start+o.iownnodes-1;
  for (int n = 0; n < num_nodes; n++) {
    gnod=o.arrays.ncorp[n];
    if(gnod >= start && gnod <= end) { // solution to write
         p[icount]= data[0*num_nodes+n];
         u[icount]= data[1*num_nodes+n];
         v[icount]= data[2*num_nodes+n];
         w[icount]= data[3*num_nodes+n];
         T[icount]= data[4*num_nodes+n];
         icount++;
    }
  }
//     write the solution field data in parallel 
  if (cg_sol_write(F, B, Z, "Solution", CGNS_ENUMV(Vertex), &S) ||
      cgp_field_write(F, B, Z, S, CGNS_ENUMV(RealDouble), "Pressure", &Q))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, S, Q, &start, &end, p))
      cgp_error_exit();
  if ( cgp_field_write(F, B, Z, S, CGNS_ENUMV(RealDouble), "VelocityX", &Q))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, S, Q, &start, &end, u))
      cgp_error_exit();
  if ( cgp_field_write(F, B, Z, S, CGNS_ENUMV(RealDouble), "VelocityY", &Q))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, S, Q, &start, &end, v))
      cgp_error_exit();
  if ( cgp_field_write(F, B, Z, S, CGNS_ENUMV(RealDouble), "VelocityZ", &Q))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, S, Q, &start, &end, w))
      cgp_error_exit();
  if ( cgp_field_write(F, B, Z, S, CGNS_ENUMV(RealDouble), "Temperature", &Q))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, S, Q, &start, &end, T))
      cgp_error_exit();
  free(p); free(u); free(v); free(w); free(T); free(data);
}
void CGNS_Coordinates(int F,int B,int Z,Output& o)
{
   int Cx,Cy,Cz;
  if (cgp_coord_write(F, B, Z, CGNS_ENUMV(RealDouble), "CoordinateX", &Cx) ||
      cgp_coord_write(F, B, Z, CGNS_ENUMV(RealDouble), "CoordinateY", &Cy) ||
      cgp_coord_write(F, B, Z, CGNS_ENUMV(RealDouble), "CoordinateZ", &Cz))
      cgp_error_exit();

// condense out vertices owned by another rank in a new array, x, whose slices are ready for CGNS.
  int num_nodes=o.mesh->count(0);
  cgsize_t gnod;
  cgsize_t start=o.local_start_id;
  cgsize_t end=start+o.iownnodes-1;
  double* x = (double *)malloc(o.iownnodes * sizeof(double));
  for (int j = 0; j < 3; ++j) {
    int icount=0;
    for (int inode = 0; inode < num_nodes; ++inode){
      gnod=o.arrays.ncorp[inode];
      if(gnod >= start && gnod <= end) { // coordinate to write
         x[icount]= o.arrays.coordinates[j*num_nodes+inode];
         icount++;
      }
    }
if(0==1) {
    printf("%ld, %ld \n", start, end);
    for (int ne=0; ne<std::min(icount,nDbgI); ++ne)
	printf("%d, %f \n", (ne+1), x[ne]);
}
    if(j==0) if(cgp_coord_write_data(F, B, Z, Cx, &start, &end, x)) cgp_error_exit();
    if(j==1) if(cgp_coord_write_data(F, B, Z, Cy, &start, &end, x)) cgp_error_exit();
    if(j==2) if(cgp_coord_write_data(F, B, Z, Cz, &start, &end, x)) cgp_error_exit();
  }
  free (x);
  int S2,Fs2;
  const int part = PCU_Comm_Self() ;
  if (cg_sol_write(F, B, Z, "VertexRank", CGNS_ENUMV(Vertex), &S2) ||
      cgp_field_write(F, B, Z, S2, CGNS_ENUMV(Integer), "VertexRank", &Fs2))
      cgp_error_exit();
  int* d = (int *)malloc(o.iownnodes * sizeof(int));
  for (int inode = 0; inode < o.iownnodes; ++inode) d[inode]= part;
  if (cgp_field_write_data(F, B, Z, S2, Fs2, &start, &end, d))
          cgp_error_exit();
  free(d);
}
void writeCGNS(Output& o, std::string path)
{
  double t0 = PCU_Time();
  apf::Mesh* m = o.mesh;
  const int num_parts = PCU_Comm_Peers();
  const int part = PCU_Comm_Self() ;
  const cgsize_t  num_parts_cg=num_parts;
  std::string timestep_or_dat;
  static char outfile[] = "chefOut.cgns";
  int  F, B, Z, E, S, Fs, Fs2, A, Cx, Cy, Cz;
  cgsize_t sizes[3],*e, start, end;
  int num_nodes=m->count(0);
if(0==1){  // ilwork debugging
    for (int ipart=0; ipart<num_parts; ++ipart){
        if(part==ipart) { // my turn
           printf("ilwork %d, %d, %d \n", part, o.nlwork,o.arrays.ilwork[0]);
           int ist=0;
           for (int itask=0; itask<o.arrays.ilwork[0]; ++itask) {
              printf("%d  ",itask);
              for (int itt=1; itt<5; ++itt)  printf("%d ", o.arrays.ilwork[ist+itt]);
              printf(" \n");
              int pnumseg=o.arrays.ilwork[ist+4];
              for (int is=0; is<pnumseg; ++is) { 
                 printf("%d, %d, %d \n",is,o.arrays.ilwork[ist+5+2*is],o.arrays.ilwork[ist+6+2*is]);
              } 
           }
       }
       PCU_Barrier();
     }
}
if(0==1){
  for (int ipart=0; ipart<num_parts; ++ipart){
    if(part==part) { // my turn    
    printf("xyz %d, %d \n", part, num_nodes);
    for (int inode = 0; inode < std::min(nDbgI,num_nodes); ++inode){
      printf("%d ",inode+1);
      for (int j=0; j<3; ++j) printf("%f ", o.arrays.coordinates[j*num_nodes+inode]);
      printf(" \n");
   }
       }
       PCU_Barrier();
     }
}
// copied gen_ncorp from PHASTA to help map on-rank numbering to CGNS/PETSC friendly global numbering
  gen_ncorp( o );
//  o carries
//     o.arrays.ncorp[on-rank-node-number(0-based)] => PETSc global node number (1-based)
//     o.iownnodes => nodes owned by this rank
//     o.local_start_id => this rank's first node number (1-based and also which must be a long long int)
  long safeArg=o.iownnodes; // cgsize_t could be an int
  sizes[0]=PCU_Add_Long(safeArg);
  int ncells=m->count(m->getDimension()); // this ranks number of elements
  safeArg=ncells; // cgsize_t could be an int
  sizes[1]=PCU_Add_Long(safeArg);
  sizes[2]=0;
  if(cgp_mpi_comm(MPI_COMM_WORLD)) cgp_error_exit;
  if ( cgp_open(outfile, CG_MODE_WRITE, &F) ||
       cg_base_write(F, "Base", 3, 3, &B) ||
       cg_zone_write(F, B, "Zone", sizes, CGNS_ENUMV(Unstructured), &Z))
       cgp_error_exit();
    // create data nodes for coordinates 
  cg_set_file_type(CG_FILE_HDF5);
  CGNS_Coordinates(F,B,Z,o);
  CGNS_NodalSolution(F,B,Z,o);
  // create Helper array for number of elements on rank 
  if ( cg_goto(F, B, "Zone_t", 1, NULL) ||
       cg_user_data_write("User Data") ||
       cg_gorel(F, "User Data", 0, NULL) ||
       cgp_array_write("nCoordsOnRank", CGNS_ENUMV(Integer), 1, &num_parts_cg, &Fs2))
       cgp_error_exit();
  // create the field data for this process 
  int nCoordVec=o.iownnodes;
  cgsize_t partP1=part+1;
if(0==1)  printf("Coor %d, %d, %d, \n", nCoordVec,part,Fs2);
  if ( cgp_array_write_data(Fs2, &partP1, &partP1, &nCoordVec))
       cgp_error_exit();
  cgsize_t e_written=0; 
  if(o.writeCGNSFiles > 2) 
    writeBlocksCGNSinteror(F,B,Z,o,&e_written);
  if(o.writeCGNSFiles > 2) {
    cgsize_t totBel;
    int nblkb = o.blocks.boundary.getSize(); 
    double** srfIDCen1 = new double*[nblkb]; // might not all be used
    double** srfIDCen2 = new double*[nblkb];
    int totOnRankBel=0;
    for (int i = 0; i < nblkb; ++i) 
      totOnRankBel += o.blocks.boundary.nElements[i];
    int* srfID = (int *)malloc( totOnRankBel * sizeof(int));
    int* srfID1OnBlk = (int *)malloc( nblkb * sizeof(int));
    int* srfID2OnBlk = (int *)malloc( nblkb * sizeof(int));
    int* startBelBlk = (int *)malloc( nblkb * sizeof(int));
    int* endBelBlk = (int *)malloc( nblkb * sizeof(int));
    int* srfIDidx = (int *)malloc( totOnRankBel * sizeof(int));
    int nStackedOnRank=0;
    writeBlocksCGNSboundary(F,B,Z,o, srfID, srfIDidx, srfIDCen1, srfIDCen2, srfID1OnBlk, srfID2OnBlk, startBelBlk, endBelBlk, &e_written, &totBel, &nStackedOnRank, nblkb);
    writeCGNSboundary      (F,B,Z,o, srfID, srfIDidx, srfIDCen1, srfIDCen2, srfID1OnBlk, srfID2OnBlk, startBelBlk, endBelBlk, &e_written, totOnRankBel, &totBel,  nStackedOnRank);
    free(srfID); free(srfIDidx);
    free(srfID1OnBlk); free(srfID2OnBlk);
    free(startBelBlk); free(endBelBlk);
    for (int i = 0; i < nStackedOnRank; ++i) delete [] srfIDCen1[i];
    for (int i = 0; i < nStackedOnRank; ++i) delete [] srfIDCen2[i];
    delete [] srfIDCen1; delete [] srfIDCen2;
    if(cgp_close(F)) cgp_error_exit();
    double t1 = PCU_Time();
    if (!PCU_Comm_Self())
      lion_oprint(1,"CGNS file written in %f seconds\n", t1 - t0);
  }
}
} // namespace
