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
typedef long long int gcorp_t;

namespace ph {


static lcorp_t count_owned(int* ilwork, int nlwork,gcorp_t* ncorp_tmp, int num_nodes);
static lcorp_t count_local(int* ilwork, int nlwork,gcorp_t* ncorp_tmp, int num_nodes);


void gen_ncorp(Output& o )
{
        apf::Mesh* m = o.mesh;
	int part;
	int num_parts;
	int i;
	lcorp_t nilwork = o.nlwork;
        int num_nodes=m->count(0);
	o.arrays.ncorp = new gcorp_t[num_nodes];
	lcorp_t owned;
	lcorp_t local;
	lcorp_t* owner_counts;
	gcorp_t  local_start_id;
	gcorp_t  gid;

	MPI_Comm_rank(MPI_COMM_WORLD, &part);
	MPI_Comm_size(MPI_COMM_WORLD, &num_parts);

	memset(o.arrays.ncorp, 0, sizeof(gcorp_t)*(num_nodes));
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

}

static lcorp_t count_local(int* ilwork, int nlwork,gcorp_t* ncorp_tmp, int num_nodes)
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
static lcorp_t count_owned(int* ilwork, int nlwork,gcorp_t* ncorp_tmp, int num_nodes)
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
//  ss << "geombc." << timestep_or_dat << "." << rank;
  ss << "chefO." << timestep_or_dat;
  return ss.str();
}

enum {
  MAX_PARAMS = 12
};

// renamed, update is only a transpose to match CNGS.  Parallel will require mapping here or later to global numbering
void getInteriorConnectivityCGNS(Output& o, int block, apf::DynamicArray<int>& c)
{
  int nelem = o.blocks.interior.nElements[block];
  int nvert = o.blocks.interior.keys[block].nElementVertices;
  c.setSize(nelem * nvert);
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert; ++vert)
      c[i++] = o.arrays.ncorp[o.arrays.ien[block][elem][vert]-1]; // input is 0-based,  out is  1-based do drop the +1
  PCU_ALWAYS_ASSERT(i == c.getSize());
}

//renamed, update is both a transpose to match CNGS and reduction to only filling the first number of vertices on the boundary whereas PHAST wanted full volume
void getBoundaryConnectivityCGNS(Output& o, int block, apf::DynamicArray<int>& c)
{
  int nelem = o.blocks.boundary.nElements[block];
// CGNS wants surface elements  int nvert = o.blocks.boundary.keys[block].nElementVertices;
  int nvert = o.blocks.boundary.keys[block].nBoundaryFaceEdges;
  c.setSize(nelem * nvert);
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert; ++vert)
      c[i++] = o.arrays.ncorp[o.arrays.ienb[block][elem][vert]-1]; 
  PCU_ALWAYS_ASSERT(i == c.getSize());
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
void getNaturalBCCodesCGNS(Output& o, int block, apf::DynamicArray<int>& codes)
{
  int nelem = o.blocks.boundary.nElements[block];
  codes.setSize(nelem * 2); 
  size_t i = 0;
  for (int j = 0; j < 2; ++j)
    for (int elem = 0; elem < nelem; ++elem)
      codes[i++] = o.arrays.ibcb[block][elem][j];
  PCU_ALWAYS_ASSERT(i == codes.getSize());
}

// renamed and calling the renamed functions above with output writes commented as they are PHASTA file style
void writeBlocksCGNS(FILE* f, Output& o)
{
  apf::DynamicArray<int> c;
  int params[MAX_PARAMS];
  for (int i = 0; i < o.blocks.interior.getSize(); ++i) {
    BlockKey& k = o.blocks.interior.keys[i];
    std::string phrase = getBlockKeyPhrase(k, "connectivity interior ");
    params[0] = o.blocks.interior.nElements[i];
//    fillBlockKeyParams(params, k);
    getInteriorConnectivityCGNS(o, i, c);
//    ph_write_ints(f, phrase.c_str(), &c[0], c.getSize(), 7, params);
  }
  for (int i = 0; i < o.blocks.boundary.getSize(); ++i) {
    BlockKey& k = o.blocks.boundary.keys[i];
    std::string phrase = getBlockKeyPhrase(k, "connectivity boundary ");
    params[0] = o.blocks.boundary.nElements[i];
//    fillBlockKeyParams(params, k);
    getBoundaryConnectivityCGNS(o, i, c);
//    ph_write_ints(f, phrase.c_str(), &c[0], c.getSize(), 8, params);
// this is probably the easiest path to getting the list that tells us the face (through surfID of smd) that each boundary element face is on
    phrase = getBlockKeyPhrase(k, "nbc codes ");
    apf::DynamicArray<int> codes;
    getNaturalBCCodesCGNS(o, i, codes);
//    ph_write_ints(f, phrase.c_str(), &codes[0], codes.getSize(), 8, params);
  }

}



// WIP
void writeCGNS(Output& o, std::string path)
{
  double t0 = PCU_Time();
  apf::Mesh* m = o.mesh;
  int         cgid = -1;

  std::string timestep_or_dat;
//  if (! timestep)
    timestep_or_dat = "cgns";
//  else {
//    tss << timestep;   
//    timestep_or_dat = tss.str();
//  }
//  cgp_mpi_comm();
//  cgp_open('chefOut.cgns', CG_MODE_WRITE, &cgid);
//static std::string buildCGNSFileName(std::string timestep_or_dat)
//  path += buildCGNSFileName(timestep_or_dat);
  static char *outfile = "chefOut.cgns";
 // if (!PCU_Comm_Self())
    cgp_mpi_comm(MPI_COMM_WORLD);
    cgp_open(outfile, CG_MODE_READ, &cgid);
  
//FAILED    cgp_open('chefO.cgns', CG_MODE_READ, &cgid);
//    PetscCheck(cgid > 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "cg_open(\"%s\",...) did not return a valid file ID", filename);
     
// copied gen_ncorp from PHASTA to help map on-rank numbering to CGNS/PETSC friendly global numbering
  gen_ncorp( o );
//  o carries
//     o.arrays.ncorp[on-rank-node-number(0-based)] => PETSc global node number (1-based)
//     o.iownnodes => nodes owned by this rank
//     o.local_start_id => this rank's first node number (1-based and also which must be a long long int)
//     o.numGlobalNodes
       int numel=m->count(m->getDimension());
       PCU_Add_Ints(&numel,1);
       o.numGlobalVolumeElements = numel;


// condense out vertices owned by another rank in a new array, x, whose slices are ready for CGNS.  Seeing now PETSc CGNS writer did one coordinate at a time which is probably better....feel free to rewrite. 
  int num_nodes=m->count(0);
  int icount=0;
  gcorp_t gnod; 
  double* x = new double[o.iownnodes * 3];
  for (int inode = 0; inode < num_nodes; ++inode){
    gnod=o.arrays.ncorp[inode];
    if(gnod >= o.local_start_id && gnod <= o.local_start_id + o.iownnodes -1) { // coordinate to write 
       for (int j = 0; j < 3; ++j) 
         x[j*o.iownnodes+icount]= o.arrays.coordinates[j*num_nodes+inode];
       icount++;
    }
  }

  
//  path += buildCGNSFileName(timestep_or_dat);
//  phastaio_setfile(GEOMBC_WRITE);
//  FILE* f = o.openfile_write(o, path.c_str());
//  if (!f) {
//    lion_eprint(1,"failed to open \"%s\"!\n", path.c_str());
//    abort();
//  }
//  ph_write_preamble(f);
  int params[MAX_PARAMS];
  
/* all of these strings are looked for by the other programs
   reading this format, so don't fix spelling errors or
   other silliness, it has already been set in stone */
/*
  writeInt(f, "number of nodes", m->count(0));
  writeInt(f, "number of modes", o.nOverlapNodes);
  writeInt(f, "number of shapefunctions soved on processor", 0);
  writeInt(f, "number of global modes", 0);
  writeInt(f, "number of interior elements", m->count(m->getDimension()));
  writeInt(f, "number of boundary elements", o.nBoundaryElements);
  writeInt(f, "maximum number of element nodes", o.nMaxElementNodes);
  writeInt(f, "number of interior tpblocks", o.blocks.interior.getSize());
  writeInt(f, "number of boundary tpblocks", o.blocks.boundary.getSize());
  writeInt(f, "number of nodes with Dirichlet BCs", o.nEssentialBCNodes);

  params[0] = m->count(0);
  params[1] = 3;
  ph_write_doubles(f, "co-ordinates", o.arrays.coordinates,
      params[0] * params[1], 2, params);
  writeInt(f, "number of processors", PCU_Comm_Peers());
  writeInt(f, "size of ilwork array", o.nlwork);
  params[0] = m->count(0);
  writeInts(f, " mode number map from partition to global",
      o.arrays.globalNodeNumbers, m->count(0));
  writeBlocksCGNS(f, o);
  writeInts(f, "bc mapping array", o.arrays.nbc, m->count(0));
  writeInts(f, "bc codes array", o.arrays.ibc, o.nEssentialBCNodes);
  apf::DynamicArray<double> bc;
  PHASTAIO_CLOSETIME(fclose(f);)
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    lion_oprint(1,"geombc file written in %f seconds\n", t1 - t0);
*/
}
}
