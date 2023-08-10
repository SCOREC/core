#include <PCU.h>
#include "phOutput.h"
#include "phIO.h"
#include "phiotimer.h"
#include <sstream>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>
#include <string.h>
#include <assert.h>
#ifdef HAVE_CGNS
//
#include <cgns_io.h>
#include <pcgnslib.h>
//
#endif
typedef int lcorp_t;
#define NCORP_MPI_T MPI_INTEGER

namespace ph {


static lcorp_t count_owned(int* ilwork, int nlwork,cgsize_t* ncorp_tmp, int num_nodes);
static lcorp_t count_local(int* ilwork, int nlwork,cgsize_t* ncorp_tmp, int num_nodes);


void gen_ncorp(Output& o )
{
        apf::Mesh* m = o.mesh;
	int part;
	int num_parts;
	int i;
	lcorp_t nilwork = o.nlwork;
        int num_nodes=m->count(0);
	o.arrays.ncorp = new cgsize_t[num_nodes];
	lcorp_t owned;
	lcorp_t local;
	lcorp_t* owner_counts;
	cgsize_t  local_start_id;
	cgsize_t  gid;

	MPI_Comm_rank(MPI_COMM_WORLD, &part);
	MPI_Comm_size(MPI_COMM_WORLD, &num_parts);

	memset(o.arrays.ncorp, 0, sizeof(cgsize_t)*(num_nodes));
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
	memset(owner_counts, 0, sizeof(lcorp_t)*num_parts);
	owner_counts[part] = owned+local;
#ifdef PRINT_EVERYTHING
	for(i=0;i<num_parts;i++)
	{
		printf("%d,", owner_counts[i]);
	}
	printf("\n");
#endif
	MPI_Allgather(MPI_IN_PLACE, 1, NCORP_MPI_T, owner_counts,
		       	1, NCORP_MPI_T, MPI_COMM_WORLD);
#ifdef PRINT_EVERYTHING
	for(i=0;i<num_parts;i++)
	{
		printf("%d,", owner_counts[i]);
	}
	printf("\n");
#endif
	local_start_id=0;
	for(i=0;i<part;i++) //TODO: MPI_Exscan()?
	{
// global so needs long long
		local_start_id += owner_counts[i];
	}
	local_start_id++; //Fortran numbering
        o.local_start_id = local_start_id;

// also get the global number of nodes
	o.numGlobalNodes=0;
	for(i=0;i<num_parts;i++)
	   o.numGlobalNodes += owner_counts[i];

#ifdef PRINT_EVERYTHING
	printf("%d: %d\n", part, local_start_id);
#endif
// global so needs long long
	gid = local_start_id;
        if(gid<0) printf("part,gid, %d %ld",part,gid);
        assert(gid>=0);
	for(i=0;i<num_nodes;i++) //assign owned node's numbers
	{
		//if shared, owned 1
			//if shared, slave -1
			//if local only, 0
		if(o.arrays.ncorp[i] == 1)
		{
// global so needs long long
			o.arrays.ncorp[i]=gid;
                        assert(o.arrays.ncorp[i]>=0);

// global so needs long long
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
		{
			o.arrays.ncorp[i] = 0; //commu() adds, so zero slaves
		}

	}
	//char code[] = "out";
	//int ione = 1;
     cgsize_t* ncorp = new cgsize_t[num_nodes];

     if(num_parts > 1) {
// translating a commuInt out from PHASTA to c
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
        int isbegin[maxseg];
        int lenseg[maxseg];
        int ioffset[maxseg];
        MPI_Request  req[numtask];
        MPI_Status stat[numtask];
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
          MPI_Type_indexed (numseg, lenseg, ioffset,MPI_LONG_LONG_INT, &sevsegtype[itask]);
          MPI_Type_commit (&sevsegtype[itask]);
          itkbeg+=4+2*numseg;
        }

        int m = 0; 
        itkbeg=0;
        for (int itask=0; itask<numtask; ++itask) {
          itag   = o.arrays.ilwork[itkbeg + 1];
          iacc   = o.arrays.ilwork[itkbeg + 2];
          iother = o.arrays.ilwork[itkbeg + 3] - 1; // MPI is 0 based but this was prepped wrong
          numseg = o.arrays.ilwork[itkbeg + 4]; /// not used
          isgbeg = o.arrays.ilwork[itkbeg + 5] - 1;
          if (iacc==0){ 
             MPI_Irecv(&o.arrays.ncorp[isgbeg], 1, sevsegtype[itask],iother, itag, MPI_COMM_WORLD, &req[m]);
          } else {
             MPI_Isend(&o.arrays.ncorp[isgbeg], 1, sevsegtype[itask],iother, itag, MPI_COMM_WORLD, &req[m]);
          }
          itkbeg+=4+2*numseg;
          m      = m + 1; 
        }
        MPI_Waitall(m, req, stat);
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
  int rank = PCU_Comm_Self() + 1;
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
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert; ++vert)
      c[i++] = o.arrays.ncorp[o.arrays.ien[block][elem][vert]]; // input is 0-based,  out is  1-based do drop the +1
  PCU_ALWAYS_ASSERT(i == nelem*nvert);
}

// update is both a transpose to match CNGS and reduction to only filling the first number of vertices on the boundary whereas PHASTA wanted full volume
void getBoundaryConnectivityCGNS(Output& o, int block, cgsize_t* c)
{
  int nelem = o.blocks.boundary.nElements[block];
// CGNS wants surface elements  int nvert = o.blocks.boundary.keys[block].nElementVertices;
  int nvert = o.blocks.boundary.keys[block].nBoundaryFaceEdges;
  //c.setSize(nelem * nvert);
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert; ++vert)
      c[i++] = o.arrays.ncorp[o.arrays.ienb[block][elem][vert]];
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
      c[i++] = o.arrays.ncorp[o.arrays.ienif0[block][elem][vert]-1];
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert1; ++vert)
      c[i++] = o.arrays.ncorp[o.arrays.ienif1[block][elem][vert]-1];
  PCU_ALWAYS_ASSERT(i == c.getSize());
}

// renamed but not updated yet
void getNaturalBCCodesCGNS(Output& o, int block, int* codes)
{
  int nelem = o.blocks.boundary.nElements[block];
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
      codes[i++] = o.arrays.ibcb[block][elem][1]; //srfID is the second number so 1
// if we wanted we could use PHASTA's bit in coding in the first number to us attributes to set
// arbitrary combinations of BCs but leaving that out for now
}

// renamed and calling the renamed functions above with output writes commented as they are PHASTA file style
void writeBlocksCGNS(int F,int B,int Z, Output& o)
{
  int params[MAX_PARAMS];
  int E;
  cgsize_t e_owned, e_start,e_end;
  cgsize_t e_startg,e_endg;
  cgsize_t e_written=0;
  for (int i = 0; i < o.blocks.interior.getSize(); ++i) {
    BlockKey& k = o.blocks.interior.keys[i];
    std::string phrase = getBlockKeyPhrase(k, "connectivity interior ");
    params[0] = o.blocks.interior.nElements[i];
    e_owned = o.blocks.interior.nElements[i];
    int nvert = o.blocks.interior.keys[i].nElementVertices;
    cgsize_t* e = (cgsize_t *)malloc(nvert * e_owned * sizeof(cgsize_t));
    getInteriorConnectivityCGNS(o, i, e);
    /* create data node for elements */
    e_startg=1+e_written; // start for the elements of this topology
    e_endg=e_written + PCU_Add_Long(e_owned); // end for the elements of this topology
    switch(nvert){
      case 4:
        if (cgp_section_write(F, B, Z, "Tet", CG_TETRA_4, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
      case 5:
        if (cgp_section_write(F, B, Z, "Pyr", CG_PYRA_5, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
      case 6:
        if (cgp_section_write(F, B, Z, "Wdg", CG_PENTA_6, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
      case 8:
        if (cgp_section_write(F, B, Z, "Hex", CG_HEXA_8, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
    }
    e_start=0;
    MPI_Exscan(&e_owned, &e_start, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    e_start+=1+e_written; // my ranks global element start 1-based
    e_end=e_start+e_owned-1;  // my ranks global element stop 1-based
    /* write the element connectivity in parallel */
    if (cgp_elements_write_data(F, B, Z, E, e_start, e_end, e))
        cgp_error_exit();
    e_written=e_endg; // update count of elements written
    free(e);
  }
  if(o.writeCGNSFiles > 2) {
  for (int i = 0; i < o.blocks.boundary.getSize(); ++i) {
    BlockKey& k = o.blocks.boundary.keys[i];
    params[0] = o.blocks.boundary.nElements[i];
    e_owned = params[0];
    int nvert = o.blocks.boundary.keys[i].nBoundaryFaceEdges;
    cgsize_t* e = (cgsize_t *)malloc(nvert * e_owned * sizeof(cgsize_t));
    getBoundaryConnectivityCGNS(o, i, e);
    e_startg=1+e_written; // start for the elements of this topology
    e_endg=e_written + PCU_Add_Long(e_owned); // end for the elements of this topology
    switch(nvert){
      case 3:
        if (cgp_section_write(F, B, Z, "Tri", CG_TETRA_4, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
      case 4:
        if (cgp_section_write(F, B, Z, "Quad", CG_QUAD_4, e_startg, e_endg, 0, &E))
           cgp_error_exit();
        break;
    }
    e_start=0;
    MPI_Exscan(&e_owned, &e_start, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    e_start+=1+e_written; // my ranks global element start 1-based
    e_end=e_start+e_owned-1;  // my ranks global element stop 1-based
    /* write the element connectivity in parallel */
    if (cgp_elements_write_data(F, B, Z, E, e_start, e_end, e))
        cgp_error_exit();
    free(e);
    int* srfID = (int *)malloc(nvert * e_owned * sizeof(int));
    getNaturalBCCodesCGNS(o, i, srfID);
//  I am not sure if you want to put the code here to generate the face BC "node" but srfID has
//  a number from 1 to 6 for the same numbered surfaces as we use in the box

  }
 }
}

void writeCGNS(Output& o, std::string path)
{
  double t0 = PCU_Time();
  apf::Mesh* m = o.mesh;
  int rank = PCU_Comm_Self() + 0;

  std::string timestep_or_dat;
  static char outfile[] = "chefOut.cgns";
  int  F, B, Z, E, S, Fs, A, Cx, Cy, Cz;
  cgsize_t sizes[3],*e, start, end, ncells;

    int num_nodes=m->count(0);


// copied gen_ncorp from PHASTA to help map on-rank numbering to CGNS/PETSC friendly global numbering
    gen_ncorp( o );
//  o carries
//     o.arrays.ncorp[on-rank-node-number(0-based)] => PETSc global node number (1-based)
//     o.iownnodes => nodes owned by this rank
//     o.local_start_id => this rank's first node number (1-based and also which must be a long long int)
//     o.numGlobalNodes
    ncells=m->count(m->getDimension());
    ncells=PCU_Add_Long(ncells);
    o.numGlobalVolumeElements = ncells;

    sizes[0]=o.numGlobalNodes;
    sizes[1]=ncells;
    sizes[2]=0;
    if(cgp_mpi_comm(MPI_COMM_WORLD)) cgp_error_exit;
    if ( cgp_open(outfile, CG_MODE_WRITE, &F) ||
        cg_base_write(F, "Base", 3, 3, &B) ||
        cg_zone_write(F, B, "Zone", sizes, CG_Unstructured, &Z))
        cgp_error_exit();
    /* create data nodes for coordinates */
    cg_set_file_type(CG_FILE_HDF5);

    if (cgp_coord_write(F, B, Z, CG_RealDouble, "CoordinateX", &Cx) ||
        cgp_coord_write(F, B, Z, CG_RealDouble, "CoordinateY", &Cy) ||
        cgp_coord_write(F, B, Z, CG_RealDouble, "CoordinateZ", &Cz))
        cgp_error_exit();

// condense out vertices owned by another rank in a new array, x, whose slices are ready for CGNS.  Seeing now PETSc CGNS writer did one coordinate at a time which is probably better....feel free to rewrite.
  cgsize_t gnod;
  start=o.local_start_id;
  end=start+o.iownnodes-1;
  double* x = new double[o.iownnodes];
  for (int j = 0; j < 3; ++j) {
    int icount=0;
    for (int inode = 0; inode < num_nodes; ++inode){
      gnod=o.arrays.ncorp[inode];
      if(gnod >= start && gnod <= end) { // coordinate to write
         x[icount]= o.arrays.coordinates[j*num_nodes+inode];
         icount++;
      }
    }
    if(j==0) if(cgp_coord_write_data(F, B, Z, Cx, &start, &end, x)) cgp_error_exit();
    if(j==1) if(cgp_coord_write_data(F, B, Z, Cy, &start, &end, x)) cgp_error_exit();
    if(j==2) if(cgp_coord_write_data(F, B, Z, Cz, &start, &end, x)) cgp_error_exit();
  }
  if(o.writeCGNSFiles > 1) 
  writeBlocksCGNS(F,B,Z, o);
  if(cgp_close(F)) cgp_error_exit();
}
} // namespace
