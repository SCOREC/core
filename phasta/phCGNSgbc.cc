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
// other using pair in STL.
#include <bits/stdc++.h>
using namespace std;
 
// Function to sort integer array b[]
// according to the order defined by a[]
void pairsortDI(double a[], int b[], int n)
{
    pair<double, int> pairt[n];
 
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
}

// Function to sort integer array b[]
// according to the order defined by a[]
void pairsort(int a[], int b[], int n)
{
    pair<int, int> pairt[n];
 
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
    delete idx;
    delete p;
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
  o.arrays.ncorp = (cgsize_t *)malloc(num_nodes * sizeof(cgsize_t)); //FIXME where to deallocate 
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
if(1==0) {
  for (int ipart=0; ipart<num_parts; ++ipart){
    if(part==ipart) { // my turn
      for (int inod=0; inod<num_nodes; ++inod) printf("%ld ", o.arrays.ncorp[inod]);
        printf(" \n");
    }
    PCU_Barrier();
  }
}
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

enum {
  MAX_PARAMS = 12
};

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
void getBoundaryConnectivityCGNS(Output& o, int block, cgsize_t* c, double* ecenx, double* eceny, double* ecenz)
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
    ecenx[elem]=0;
    eceny[elem]=0;
    ecenz[elem]=0;
    for (int vert = 0; vert < nvert; ++vert){
      phGnod=o.arrays.ienb[block][elem][lnode[vert]]; //actually it is on-rank Global 
      c[i++] = o.arrays.ncorp[phGnod]; // PETSc truely global
      ecenx[elem]+=o.arrays.coordinates[0*num_nodes+phGnod];
      eceny[elem]+=o.arrays.coordinates[1*num_nodes+phGnod];
      ecenz[elem]+=o.arrays.coordinates[2*num_nodes+phGnod];
    }
    ecenx[elem]/=nvert; // only necessary if you really want to use this as a correct centroid rather than comparison
    eceny[elem]/=nvert; // only necessary if you really want to use this as a correct centroid rather than comparison
    ecenz[elem]/=nvert; // only necessary if you really want to use this as a correct centroid rather than comparison
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

// renamed and calling the renamed functions above with output writes now to CGNS

void writeBlocksCGNS(int F,int B,int Z, Output& o)
{
  int params[MAX_PARAMS];
  int E,S,Fs,Fs2,Fsb,Fsb2;
  cgsize_t e_owned, e_start,e_end;
  cgsize_t e_startg,e_endg;
  cgsize_t e_written=0;
  const int num_parts = PCU_Comm_Peers();
  const cgsize_t num_parts_cg=num_parts;
  const int part = PCU_Comm_Self() ;
  const cgsize_t part_cg=part;
  // create a centered solution 
  if (cg_sol_write(F, B, Z, "RankOfWriter", CG_CellCenter, &S) ||
      cgp_field_write(F, B, Z, S, CG_Integer, "RankOfWriter", &Fs))
      cgp_error_exit();
  for (int i = 0; i < o.blocks.interior.getSize(); ++i) {
    BlockKey& k = o.blocks.interior.keys[i];
    std::string phrase = getBlockKeyPhrase(k, "connectivity interior ");
    params[0] = o.blocks.interior.nElements[i];
    e_owned = o.blocks.interior.nElements[i];
    int nvert = o.blocks.interior.keys[i].nElementVertices;
    cgsize_t* e = (cgsize_t *)malloc(nvert * e_owned * sizeof(cgsize_t));
    getInteriorConnectivityCGNS(o, i, e);
    // create data node for elements 
    e_startg=1+e_written; // start for the elements of this topology
    long safeArg=e_owned; // e_owned is cgsize_t which could be an 32 or 64 bit int
    e_endg=e_written + PCU_Add_Long(safeArg); // end for the elements of this topology
    char Ename[5];
    switch(nvert){
      case 4:
        snprintf(Ename, 4, "Tet");
        if (cgp_section_write(F, B, Z, Ename, CG_TETRA_4, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
      case 5:
        snprintf(Ename, 4, "Pyr");
        if (cgp_section_write(F, B, Z, Ename, CG_PYRA_5, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
      case 6:
        snprintf(Ename, 4, "Wdg");
        if (cgp_section_write(F, B, Z, Ename, CG_PENTA_6, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
      case 8:
        snprintf(Ename, 4, "Hex");
        if (cgp_section_write(F, B, Z, Ename, CG_HEXA_8, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
    }
    e_start=0;
    auto type = getMpiType( cgsize_t() );
    MPI_Exscan(&e_owned, &e_start, 1, type , MPI_SUM, MPI_COMM_WORLD);
    e_start+=1+e_written; // my parts global element start 1-based
    e_end=e_start+e_owned-1;  // my parts global element stop 1-based
    // write the element connectivity in parallel 
    if (cgp_elements_write_data(F, B, Z, E, e_start, e_end, e))
        cgp_error_exit();
    e_written=e_endg; // update count of elements written

if(1==0){
    printf("interior cnn %d, %ld, %ld \n", part, e_start, e_end);
    for (int ne=0; ne<e_owned; ++ne) {
      printf("%d, %d ", part,(ne+1));
      for(int nv=0; nv< nvert; ++nv) printf("%ld ", e[ne*nvert+nv]);
      printf("\n");
    }
}
    free(e);

//     create the field data for this process 
    int* d = (int *)malloc(e_owned * sizeof(int));
    for (int n = 0; n < e_owned; n++) 
            d[n] = part;
//     write the solution field data in parallel 
    if (cgp_field_write_data(F, B, Z, S, Fs, &e_start, &e_end, d))
        cgp_error_exit();
    free(d);

    char UserDataName[11];
        snprintf(UserDataName, 11, "n%sOnRank", Ename);
        // create Helper array for number of elements on part of a given topology 
    if ( cg_goto(F, B, "Zone_t", 1, NULL) ||
         cg_gorel(F, "User Data", 0, NULL) ||
         cgp_array_write(UserDataName, CG_Integer, 1, &num_parts_cg, &Fs2))
         cgp_error_exit();
    // create the field data for this process 
    int nIelVec=e_owned;
    cgsize_t  partP1=part+1;
    printf("Intr, %s,  %d, %d, %d, %d \n", UserDataName, nIelVec,part,Fs,Fs2);
    if ( cgp_array_write_data(Fs2, &partP1, &partP1, &nIelVec))
        cgp_error_exit();
  } // end of loop over interior blocks


  if(o.writeCGNSFiles > 2) {
    cgsize_t eVolElm=e_written;
    cgsize_t e_belWritten=0;
    int totOnRankBel=0;
    int triCount=0;
    int quadCount=0;
    int nblkb = o.blocks.boundary.getSize(); 
    for (int i = 0; i < nblkb; ++i) 
      totOnRankBel += o.blocks.boundary.nElements[i];
    int* srfID = (int *)malloc( totOnRankBel * sizeof(int));
    int* srfIDidx = (int *)malloc( totOnRankBel * sizeof(int));
    double** srfIDCen1 = new double*[nblkb];
    double** srfIDCen2 = new double*[nblkb];
    int* srfID1OnBlk = (int *)malloc( nblkb * sizeof(int));
    int* srfID2OnBlk = (int *)malloc( nblkb * sizeof(int));

    int* startBelBlk = (int *)malloc( nblkb * sizeof(int));
    int* endBelBlk = (int *)malloc( nblkb * sizeof(int));
    for (int i = 0; i < o.blocks.boundary.getSize(); ++i) {
      BlockKey& k = o.blocks.boundary.keys[i];
      params[0] = o.blocks.boundary.nElements[i];
      e_owned = params[0];
      int nvert = o.blocks.boundary.keys[i].nBoundaryFaceEdges;
      cgsize_t* e = (cgsize_t *)malloc(nvert * e_owned * sizeof(cgsize_t));
      double* ecenx = (double *)malloc( e_owned * sizeof(double));
      double* eceny = (double *)malloc( e_owned * sizeof(double));
      double* ecenz = (double *)malloc( e_owned * sizeof(double));
      getBoundaryConnectivityCGNS(o, i, e,ecenx,eceny,ecenz);
      e_startg=1+e_written; // start for the elements of this topology
      long safeArg=e_owned; // e_owned is cgsize_t which could be an 32 or 64 bit int
      cgsize_t  numBelTP = PCU_Add_Long(safeArg); // number of elements of this topology
      e_endg=e_written + numBelTP; // end for the elements of this topology
      if(nvert==3) triCount++;
      if(nvert==4) quadCount++;
      char Ename[7];
      switch(nvert){
        case 3:
          snprintf(Ename, 5, "Tri%d",triCount);
          if (cgp_section_write(F, B, Z, Ename, CG_TRI_3, e_startg, e_endg, 0, &E))
            cgp_error_exit();
          break;
        case 4:
          snprintf(Ename, 6, "Quad%d",quadCount);
          if (cgp_section_write(F, B, Z, Ename, CG_QUAD_4, e_startg, e_endg, 0, &E))
            cgp_error_exit();
          break;
      }
      e_start=0;
      auto type = getMpiType( cgsize_t() );
      MPI_Exscan(&e_owned, &e_start, 1, type , MPI_SUM, MPI_COMM_WORLD);
      e_start+=1+e_written; // my parts global element start 1-based
      e_end=e_start+e_owned-1;  // my parts global element stop 1-based
      // write the element connectivity in parallel 
      if (cgp_elements_write_data(F, B, Z, E, e_start, e_end, e))
          cgp_error_exit();
      printf("boundary cnn %d, %ld, %ld \n", part, e_start, e_end);
if(1==0){
    for (int ne=0; ne<e_owned; ++ne) { printf("%d, %d ", part,(ne+1)); for(int nv=0; nv< nvert; ++nv) printf("%ld ", e[ne*nvert+nv]); printf("\n"); }
}
      getNaturalBCCodesCGNS(o, i, &srfID[e_belWritten]);
      int icnt1=0;
      int icnt2=0;
      for (int ne=0; ne<e_owned; ++ne){ //count srfID =1 and 2 on this part,block
         if(srfID[e_belWritten+ne]==1) icnt1++; 
         if(srfID[e_belWritten+ne]==2) icnt2++;
      } 
      srfIDCen1[i]=new double[icnt1*3];
      srfIDCen2[i]=new double[icnt2*3];
      srfID1OnBlk[i]=icnt1;
      srfID2OnBlk[i]=icnt2;
      int j1=0;
      int j2=0;
      for (int ne=0; ne<e_owned; ++ne){
         if(srfID[e_belWritten+ne]==1){ 
           srfIDCen1[i][j1++]=ecenx[ne];
           srfIDCen1[i][j1++]=eceny[ne];
           srfIDCen1[i][j1++]=ecenz[ne];
         }
         if(srfID[e_belWritten+ne]==2) {
           srfIDCen2[i][j2++]=ecenx[ne];
           srfIDCen2[i][j2++]=eceny[ne];
           srfIDCen2[i][j2++]=ecenz[ne];
         }
      } 
      free(ecenx); free(eceny); free(ecenz);
if(1==0){      printf("CentroidCounts %d %d %d %d %d %d %d %d\n",part,icnt1, icnt2, j1, j2, e_owned, srfID1OnBlk[i],srfID2OnBlk[i]);}
      for (int j = 0; j < (int) e_owned; ++j) srfIDidx[e_belWritten+j]=e_start+j;
      startBelBlk[i]=e_start; // provides start point for each block in srfID
      endBelBlk[i]=e_end; // provides end point for each block in srfID
      e_written=e_endg;
      e_belWritten+=e_owned; // this is tracking written by this rank as we unpack srfID later
      char UserDataName[12]; snprintf(UserDataName, 13, "n%sOnRank", Ename);
      if ( cg_goto(F, B, "Zone_t", 1, NULL) ||
           cg_gorel(F, "User Data", 0, NULL) ||
           cgp_array_write(UserDataName, CG_Integer, 1, &num_parts_cg, &Fsb2))
           cgp_error_exit();
      printf("Bndy %s, %ld, %d, %d \n", UserDataName, e_owned, part,Fsb2);
      cgsize_t partP1=part+1;
      if (cgp_array_write_data(Fsb2, &partP1, &partP1, &e_owned))
          cgp_error_exit();
    }
// srfID is for ALL Boundary faces
    cgsize_t  totBel = e_written-eVolElm;
    // setup User Data for boundary faces 
    if ( cg_goto(F, B, "Zone_t", 1, NULL) ||
         cg_gorel(F, "User Data", 0, NULL) ||
         cgp_array_write("srfID", CG_Integer, 1,&totBel, &Fsb)) 
         cgp_error_exit();
    // write the user data for this process 
    e_written=0; //recycling  eVolElm holds 
    for (int i = 0; i < nblkb; ++i) {
      int e_startB=startBelBlk[i]-eVolElm; // srfID is only for bel....matches linear order with eVolElm offset from 
                                       // bel# that starts from last volume element
      e_owned=endBelBlk[i]-startBelBlk[i]+1;
      e_start=0;
      auto type = getMpiType( cgsize_t() );
      MPI_Exscan(&e_owned, &e_start, 1, type , MPI_SUM, MPI_COMM_WORLD);
      e_start+=1+e_written; // my parts global element start 1-based
      e_end=e_start+e_owned-1;  // my parts global element stop 1-based
      printf("Bndy %s, %ld, %ld, %ld, %d, %d, %d \n", "srfID", e_start, e_end, e_owned, i, part,Fsb);
      if (cgp_array_write_data(Fsb, &e_start, &e_end, &srfID[e_startB]))
        cgp_error_exit();
      long safeArg=e_owned; // is cgsize_t which could be an 32 or 64 bit int
      e_written += PCU_Add_Long(safeArg); // number of elements of this topology
    }
// stack  connectivities on rank before gather (should preserve order)
    int* rcounts = (int *)malloc( num_parts * sizeof(int));
    int* displs = (int *)malloc( num_parts * sizeof(int));
    int numsurfID1onRank=0;
    int numsurfID2onRank=0;
    for (int i = 0; i < nblkb; ++i) numsurfID1onRank+=srfID1OnBlk[i];
    for (int i = 0; i < nblkb; ++i) numsurfID2onRank+=srfID2OnBlk[i];
    double* srfIDCen1AllBlocks = (double *)malloc(numsurfID1onRank*3 * sizeof(double));
    double* srfIDCen2AllBlocks = (double *)malloc(numsurfID2onRank*3 * sizeof(double));
    int k1=0;
    int k2=0;
    for (int i = 0; i < nblkb; ++i) {
      for (int j = 0; j < srfID1OnBlk[i]*3; ++j) srfIDCen1AllBlocks[k1++]=srfIDCen1[i][j];
      for (int j = 0; j < srfID2OnBlk[i]*3; ++j) srfIDCen2AllBlocks[k2++]=srfIDCen2[i][j];
    }
    free(srfID1OnBlk); free(srfID2OnBlk);
    delete srfIDCen1; delete srfIDCen2;
    int ncon=numsurfID1onRank*3;
    auto type_i = getMpiType( int() );
    MPI_Allgather(&ncon,1,type_i,rcounts,1,type_i,MPI_COMM_WORLD);
    displs[0]=0;
    for (int i = 1; i < num_parts; ++i) displs[i]=displs[i-1]+rcounts[i-1]; 
if(1==0){ printf("displs1 %d ",part);for(int ip=0; ip< num_parts; ++ip) printf("% ld ", displs[ip]); printf("\n"); }
    int GsrfID1cnt=displs[num_parts-1]+rcounts[num_parts-1];
if(1==0){    printf("Stack1 %d %d, %d, %d, %d, %d\n",part, GsrfID1cnt, ncon, nblkb, numsurfID1onRank, numsurfID2onRank);}
    double* srfID1Gcen = (double *)malloc( GsrfID1cnt * sizeof(double));
    auto type_d = getMpiType( double() );
    MPI_Allgatherv(srfIDCen1AllBlocks,ncon,type_d,srfID1Gcen,rcounts,displs,type_d,MPI_COMM_WORLD);
// srfID=2 repeats
    ncon=numsurfID2onRank*3;
    MPI_Allgather(&ncon,1,type_i,rcounts,1,type_i,MPI_COMM_WORLD);
    displs[0]=0;
    for (int i = 1; i < num_parts; ++i) displs[i]=displs[i-1]+rcounts[i-1]; 
if(1==0){ printf("displs2 %d ",part);for(int ip=0; ip< num_parts; ++ip) printf("% ld ", displs[ip]); printf("\n"); }
    int GsrfID2cnt=displs[num_parts-1]+rcounts[num_parts-1];
if(1==0){     printf("Stack2 %d %d, %d, %d, %d, %d\n",part, GsrfID2cnt, ncon, nblkb, numsurfID1onRank, numsurfID2onRank);}
    assert(GsrfID1cnt==GsrfID2cnt);
    int nmatchFace=GsrfID1cnt/3;
    double* srfID2Gcen = (double *)malloc( GsrfID2cnt * sizeof(double));
    MPI_Allgatherv(srfIDCen2AllBlocks,ncon,type_d,srfID2Gcen,rcounts,displs,type_d,MPI_COMM_WORLD);
    const  float  Lz=abs(srfID2Gcen[2]-srfID1Gcen[2]);
if(1==0){  printf("%d part srfID 1 xc ",part); for(int ip=0; ip< nmatchFace; ++ip) printf("%f ", srfID1Gcen[ip*3+0]); printf("\n"); }
if(1==0){  printf("%d part srfID 1 yc ",part); for(int ip=0; ip< nmatchFace; ++ip) printf("%f ", srfID1Gcen[ip*3+1]); printf("\n"); }
if(1==0){  printf("%d part srfID 1 zc ",part); for(int ip=0; ip< nmatchFace; ++ip) printf("%f ", srfID1Gcen[ip*3+2]); printf("\n"); }
       PCU_Barrier();
if(1==0){  printf("%d part srfID 2 xc ",part); for(int ip=0; ip< nmatchFace; ++ip) printf("%f ", srfID2Gcen[ip*3+0]); printf("\n"); }
if(1==0){  printf("%d part srfID 2 yc ",part); for(int ip=0; ip< nmatchFace; ++ip) printf("%f ", srfID2Gcen[ip*3+1]); printf("\n"); }
if(1==0){  printf("%d part srfID 2 zc ",part); for(int ip=0; ip< nmatchFace; ++ip) printf("%f ", srfID2Gcen[ip*3+2]); printf("\n"); }
    free(srfIDCen1AllBlocks); free(srfIDCen2AllBlocks);
    double* srfID1distSq = (double *)malloc( nmatchFace * sizeof(double));
    double* srfID2distSq = (double *)malloc( nmatchFace * sizeof(double));
    int* imapD1 = (int *)malloc( nmatchFace * sizeof(int));
    int* imapD2 = (int *)malloc( nmatchFace * sizeof(int));
    int* imapD2v = (int *)malloc( nmatchFace * sizeof(int));
    double xc=10; // true cubes with uniform meshes set up ties  (good for debugging/verifying that dumb search backup works)
    for (int i = 0; i < nmatchFace; ++i) {
      srfID1distSq[i]=(srfID1Gcen[i*3+0]-xc)*(srfID1Gcen[i*3+0]-xc) 
                   +srfID1Gcen[i*3+1]*srfID1Gcen[i*3+1] 
                   +srfID1Gcen[i*3+2]*srfID1Gcen[i*3+2]; 
      srfID2distSq[i]=(srfID2Gcen[i*3+0]-xc)*(srfID2Gcen[i*3+0]-xc) 
                   +srfID2Gcen[i*3+1]*srfID2Gcen[i*3+1] 
                   +srfID2Gcen[i*3+2]*srfID2Gcen[i*3+2]; 
      imapD1[i]=i;
      imapD2[i]=i;
    }
    if(1==0){ if(part==0) {
      printf(" srfID1dist GLOBAL B "); for(int is=0; is< nmatchFace; ++is)  printf("%f ", srfID1distSq[is]); printf("\n");
      printf(" imapD1 GLOBAL B     "); for(int is=0; is< nmatchFace; ++is)  printf("%d ", imapD1[is]); printf("\n"); }
    }
    if(1==0){ if(part==0) {
      printf(" srfID2dist GLOBAL B "); for(int is=0; is< nmatchFace; ++is)  printf("%f ", srfID2distSq[is]); printf("\n");
      printf(" imapD2 GLOBAL B     "); for(int is=0; is< nmatchFace; ++is)  printf("%d ", imapD2[is]); printf("\n"); }
    }
    pairsortDI(srfID1distSq,imapD1,nmatchFace);
    pairsortDI(srfID2distSq,imapD2,nmatchFace);

    if(1==0){ if(part==0) {
      printf(" srfID1dist GLOBAL "); for(int is=0; is< nmatchFace; ++is)  printf("%f ", srfID1distSq[is]); printf("\n");
      printf(" imapD1 GLOBAL     "); for(int is=0; is< nmatchFace; ++is)  printf("%d ", imapD1[is]); printf("\n"); }
    }
    if(1==0){ if(part==0) {
      printf(" srfID2dist GLOBAL "); for(int is=0; is< nmatchFace; ++is)  printf("%f ", srfID2distSq[is]); printf("\n");
      printf(" imapD2 GLOBAL     "); for(int is=0; is< nmatchFace; ++is)  printf("%d ", imapD2[is]); printf("\n"); }
    }

    double tol=1.0e-12;
    double tol2=1.0e-14;
    int jclosest, iclose1, iclose2;
    double d1,d2,vDistSq,vDSmin;
    int DistFails=0;
    for (int i = 0; i < nmatchFace; ++i) {
        iclose1=imapD1[i];
        iclose2=imapD2[i];
        d1=srfID1Gcen[(iclose1)*3+0]-srfID2Gcen[(iclose2)*3+0];
        d2=srfID1Gcen[(iclose1)*3+1]-srfID2Gcen[(iclose2)*3+1];
        vDistSq= d1*d1+d2*d2;
        if(vDistSq < tol2) {
           imapD2v[i]=imapD2[i];
        } else {// centroid for i-1 did not match-> search list srfID=2 list to find true match
          vDSmin=vDistSq;
          DistFails++;
          for (int j = 0; j < nmatchFace; ++j) {   // if this turns out to be taken a lot then it could be narrowed e.g. j=max(0,i-50), j< i+min(matchFace,i+50),
            iclose2=imapD2[j];
            d1=srfID1Gcen[(iclose1)*3+0]-srfID2Gcen[(iclose2)*3+0];
            d2=srfID1Gcen[(iclose1)*3+1]-srfID2Gcen[(iclose2)*3+1];
            vDistSq= d1*d1+d2*d2;
            if(vDistSq<vDSmin) {
              vDSmin=vDistSq;
              jclosest=j;
            } 
          }
          imapD2v[i]=imapD2[jclosest];
        } 
    } 
    for (int i = 0; i < nmatchFace; ++i) imapD2[i]=imapD2v[i];
    if(1==1&&part==0) {
      printf("Number of Distance Failures=%d\n ",DistFails);
      printf(" srfID1dist GLOBAL "); for(int is=0; is< nmatchFace; ++is)  printf("%f ", srfID1distSq[is]); printf("\n");
      printf(" imapD1 GLOBAL     "); for(int is=0; is< nmatchFace; ++is)  printf("%d ", imapD1[is]); printf("\n"); 
      printf(" srfID2dist GLOBAL "); for(int is=0; is< nmatchFace; ++is)  printf("%f ", srfID2distSq[is]); printf("\n");
      printf(" imapD2 GLOBAL     "); for(int is=0; is< nmatchFace; ++is)  printf("%d ", imapD2[is]); printf("\n"); }
    free(srfID1Gcen); free(srfID2Gcen);
    free(srfID1distSq); free(srfID2distSq);
    free(imapD2v);

// ZonalBC data 
    int* srfIDG = (int *)malloc( totBel * sizeof(int));
    int* srfIDGidx = (int *)malloc( totBel * sizeof(int));
    cgsize_t* donor2 = (cgsize_t *)malloc(nmatchFace * sizeof(cgsize_t));
    cgsize_t* periodic1 = (cgsize_t *)malloc(nmatchFace * sizeof(cgsize_t));
    auto type_cg = getMpiType( cgsize_t() );
    MPI_Allgather(&totOnRankBel,1,type_i,rcounts,1,type_i,MPI_COMM_WORLD);
    displs[0]=0;
    for (int i = 1; i < num_parts; ++i) displs[i]=displs[i-1]+rcounts[i-1]; 
if(0==1){ for(int ip=0; ip< num_parts; ++ip) printf("%ld ", displs[ip]); printf("\n"); }
    MPI_Allgatherv(srfID,totOnRankBel,type_i,srfIDG,rcounts,displs,type_i,MPI_COMM_WORLD);
    MPI_Allgatherv(srfIDidx,totOnRankBel,type_i,srfIDGidx,rcounts,displs,type_i,MPI_COMM_WORLD);
    free(rcounts); free(displs);
if(0==1){ if(part==0) {
    printf(" srfID GLOBAL    "); for(int is=0; is< totBel; ++is)  printf("%d ", srfIDG[is]); printf("\n");
    printf(" srfIDidx GLOBAL "); for(int is=0; is< totBel; ++is)  printf("%d ", srfIDGidx[is]); printf("\n"); }
    printf("rank %d ",part); printf(" srfID on Part "); for(int is=0; is< totOnRankBel; ++is)  printf("%d ", srfID[is]); printf("\n");
    printf(" srfIDidx on Part "); for(int is=0; is< totOnRankBel; ++is)  printf("%d ", srfIDidx[is]); printf("\n"); }
    free(srfID); free(srfIDidx); 
//    pairsort(srfIDG,srfIDGidx,totBel);
    pairDeal6sort(srfIDG,srfIDGidx,totBel);
if(1==0){ if(part==0) {
    printf(" srfID GLOBAL    "); for(int is=0; is< totBel; ++is)  printf("%d ", srfIDG[is]); printf("\n");
    printf(" srfIDidx GLOBAL "); for(int is=0; is< totBel; ++is)  printf("%d ", srfIDGidx[is]); printf("\n"); }
}
    int BC_scan=0;
    int imatch1;
    cgsize_t* eBC = (cgsize_t *)malloc(totBel * sizeof(cgsize_t));
    for (int BCid = 1; BCid < 7; BCid++) {
      int imatch=0;
      while (srfIDG[BC_scan]==BCid) {
        eBC[imatch]=srfIDGidx[BC_scan];
        BC_scan++;
        imatch++;
      }
//reorder SurfID = 1 and 2 using idmapD{1,2} based on distance to support periodicity 
      if(BCid==1) {
        for (int i = 0; i < nmatchFace; i++) periodic1[i]=eBC[imapD1[i]];
        for (int i = 0; i < nmatchFace; i++) eBC[i]=periodic1[i];
if(1==1&&part==1){ printf(" srfIDidx 1 "); for(int is=0; is< nmatchFace; ++is)  printf("%d ", eBC[is]); printf("\n"); }
      }       
      if(BCid==2) {
        for (int i = 0; i < nmatchFace; i++) donor2[i]=eBC[imapD2[i]];
        for (int i = 0; i < nmatchFace; i++) eBC[i]=donor2[i];
if(1==1&&part==1){ printf(" srfIDidx 2 "); for(int is=0; is< nmatchFace; ++is)  printf("%d ", eBC[is]); printf("\n"); }
      }       
if(0==1) {
      printf(" srfID =%d    ",BCid); for(int is=0; is< imatch; ++is)  printf("%d ", eBC[is]); printf("\n");
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
    const float  Translation[3]={0,0,-Lz};
    if (cg_conn_periodic_write(F, B, Z, cgconn, RotationCenter, RotationAngle, Translation)) cgp_error_exit();
    free(imapD1); free(imapD2);
    free(eBC); free(srfIDG); free(srfIDGidx);
  } // processing boundary elments
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
  int Fp, Fu, Fv, Fw, FT;
  int Sp, Su, Sv, Sw, ST;
  cgsize_t sizes[3],*e, start, end;

  int num_nodes=m->count(0);

if(1==0){  // ilwork debugging
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
if(1==0){
  for (int ipart=0; ipart<num_parts; ++ipart){
    if(part==part) { // my turn    
    printf("xyz %d, %d \n", part, num_nodes);
    for (int inode = 0; inode < num_nodes; ++inode){
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
       cg_zone_write(F, B, "Zone", sizes, CG_Unstructured, &Z))
       cgp_error_exit();
    // create data nodes for coordinates 
  cg_set_file_type(CG_FILE_HDF5);

  if (cgp_coord_write(F, B, Z, CG_RealDouble, "CoordinateX", &Cx) ||
      cgp_coord_write(F, B, Z, CG_RealDouble, "CoordinateY", &Cy) ||
      cgp_coord_write(F, B, Z, CG_RealDouble, "CoordinateZ", &Cz))
      cgp_error_exit();

// condense out vertices owned by another rank in a new array, x, whose slices are ready for CGNS.
  cgsize_t gnod;
  start=o.local_start_id;
  end=start+o.iownnodes-1;
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
    for (int ne=0; ne<num_nodes; ++ne)
	printf("%d, %f \n", (ne+1), x[ne]);
}
    if(j==0) if(cgp_coord_write_data(F, B, Z, Cx, &start, &end, x)) cgp_error_exit();
    if(j==1) if(cgp_coord_write_data(F, B, Z, Cy, &start, &end, x)) cgp_error_exit();
    if(j==2) if(cgp_coord_write_data(F, B, Z, Cz, &start, &end, x)) cgp_error_exit();
  }
  free (x);
/* abort....matcched mesh was needed but this breaks our approach to building ncorp
  cgsize_t* gizmin = (cgsize_t *)malloc(num_nodes * sizeof(cgsize_t));
  cgsize_t* gizmax = (cgsize_t *)malloc(num_nodes * sizeof(cgsize_t));
  double zmin=-0.5;
  double zmax=65;
  double eps=1e-3;
  double z;
  for (int inode = 0; inode < num_nodes; ++inode){
      gnod=o.arrays.ncorp[inode];
      z= o.arrays.coordinates[2*num_nodes+inode];
      if( abs(z-zmin) < eps) {
         gizmin[inode]=gnod;
         gizmax[inode]=-1*(part+1);
      } else if(abs(zmax-z)<eps){
         gizmax[inode]=gnod;
         gizmin[inode]=-1*(part+1);
      } else {
         gizmin[inode]=-1*(part+1);
         gizmax[inode]=-1*(part+1);
      }
  }
if(1==1) {
      printf(" gizmin bc on part %d ",part); for(int is=0; is< num_nodes; ++is)  printf("%d ", gizmin[is]); printf("\n");
      printf(" gizmax bc on part %d ",part); for(int is=0; is< num_nodes; ++is)  printf("%d ", gizmax[is]); printf("\n");
}
  commuInt(o, gizmax);
  commuInt(o, gizmin);
if(1==1) {
      printf(" gizmin ac on part %d ",part); for(int is=0; is< num_nodes; ++is)  printf("%d ", gizmin[is]); printf("\n");
      printf(" gizmax ac on part %d ",part); for(int is=0; is< num_nodes; ++is)  printf("%d ", gizmax[is]); printf("\n");
}
 */ 
  // create a nodal solution 
  char fieldName[12];
  snprintf(fieldName, 13, "solution");
  printf("solution=%s",fieldName);
  double* data;
  int size;
  detachField(o.mesh, fieldName, data, size);
  assert(size==5);

//     create the field data for this process 
  double* p = (double *)malloc(o.iownnodes * sizeof(double));
  double* u = (double *)malloc(o.iownnodes * sizeof(double));
  double* v = (double *)malloc(o.iownnodes * sizeof(double));
  double* w = (double *)malloc(o.iownnodes * sizeof(double));
  double* T = (double *)malloc(o.iownnodes * sizeof(double));
  int icount=0;
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
  if (cg_sol_write(F, B, Z, "Solution", CG_Vertex, &Sp) ||
      cgp_field_write(F, B, Z, Sp, CG_RealDouble, "Pressure", &Fp))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, Sp, Fp, &start, &end, p))
      cgp_error_exit();
  if ( cgp_field_write(F, B, Z, Sp, CG_RealDouble, "VelocityX", &Fu))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, Sp, Fu, &start, &end, u))
      cgp_error_exit();
  if ( cgp_field_write(F, B, Z, Sp, CG_RealDouble, "VelocityY", &Fv))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, Sp, Fv, &start, &end, v))
      cgp_error_exit();
  if ( cgp_field_write(F, B, Z, Sp, CG_RealDouble, "VelocityZ", &Fw))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, Sp, Fw, &start, &end, w))
      cgp_error_exit();
  if ( cgp_field_write(F, B, Z, Sp, CG_RealDouble, "Temperature", &FT))
      cgp_error_exit();
  if (cgp_field_write_data(F, B, Z, Sp, FT, &start, &end, T))
      cgp_error_exit();
  free(p);
  free(u);
  free(v);
  free(w);
  free(T);
  free(data);
  // create Helper array for number of elements on rank 
  if ( cg_goto(F, B, "Zone_t", 1, NULL) ||
       cg_user_data_write("User Data") ||
       cg_gorel(F, "User Data", 0, NULL) ||
       cgp_array_write("nCoordsOnRank", CG_Integer, 1, &num_parts_cg, &Fs2))
       cgp_error_exit();
  // create the field data for this process 
  int nCoordVec=o.iownnodes;
  cgsize_t partP1=part+1;
  printf("Coor %d, %d, %d, \n", nCoordVec,part,Fs2);
  if ( cgp_array_write_data(Fs2, &partP1, &partP1, &nCoordVec))
       cgp_error_exit();
  if(o.writeCGNSFiles > 1) 
  writeBlocksCGNS(F,B,Z, o);
  if(cgp_close(F)) cgp_error_exit();
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    lion_oprint(1,"CGNS file written in %f seconds\n", t1 - t0);
}
} // namespace
