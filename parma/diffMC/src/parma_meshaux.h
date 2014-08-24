#ifndef PARMA_MESHAUX_H_
#define PARMA_MESHAUX_H_

#include "apf.h"
#include "apfMesh.h"

typedef apf::DynamicArray<apf::MeshEntity*> dynEntArr;

inline void clearTag(apf::Mesh*& m, apf::MeshTag* t) {
   apf::MeshEntity* e;
   for(int d=0; d <= m->getDimension(); d++) {
      apf::MeshIterator* itr = m->begin(d);
      while( (e = m->iterate(itr)) ) 
	 if( m->hasTag(e, t) ) 
	    m->removeTag(e, t);
      m->end(itr);
   }
}

inline int getEntDim(apf::Mesh* m, apf::MeshEntity* e) {
   assert(e);
   assert(m);
   return apf::Mesh::typeDimension[m->getType(e)];
}

inline apf::MeshEntity* getUpEnt(apf::Mesh* m, apf::MeshEntity* e) {
   const int upDim = getEntDim(m, e) + 1;
   dynEntArr adjEnt;
   m->getAdjacent(e, upDim, adjEnt);
   assert( NULL != adjEnt[0] );
   return adjEnt[0];
}

inline void getEdgeAdjVtx(apf::Mesh* m, apf::MeshEntity* vtx, dynEntArr& adjVtx) {
   assert( 0 == getEntDim(m, vtx) );   
   dynEntArr adjE;
   apf::Downward adjDown;
  
   m->getAdjacent(vtx, 1, adjE);
   APF_ITERATE(dynEntArr, adjE, eit) {
      const int na = m->getDownward(*eit, 0, adjDown);
      for(int i=0; i<na; i++)
         if ( adjDown[i] != vtx ) 
            adjVtx.append(adjDown[i]);
   }
   assert( 0 != adjVtx.getSize() );
}

inline void getDwn2ndAdj(apf::Mesh* m, apf::MeshEntity* ent, dynEntArr& adj) {
   const int dim = getEntDim(m, ent);
   assert(dim >= 1);
   apf::Downward adjDown;
   const int na = m->getDownward(ent, dim-1, adjDown);
   for(int i=0; i<na; i++) {
      dynEntArr adjUp;  
      m->getAdjacent(adjDown[i], dim, adjUp);
      APF_ITERATE(dynEntArr, adjUp, upIt) 
         if ( *upIt != ent ) 
            adj.append(*upIt);
   }
}

inline int getNumTaggedEnts(apf::Mesh* m, apf::MeshTag* vtag, int visited, apf::Downward& ents, const int numEnts) {
   int taggedEnts = 0;
   int tval = -1;
   for(int i=0; i<numEnts; i++) {
      if( m->hasTag(ents[i], vtag) ) {
         m->getIntTag(ents[i], vtag, &tval);
         if ( tval == visited )
	    taggedEnts++; 
      }
   }
   return taggedEnts;
}

inline int getNumFaceOnPb(apf::Mesh* m, const int destPid, apf::Downward& adjF, const int numFaces) {
   int numPbFaces = 0;
 
   apf::Parts pbPid; 
   pbPid.insert(m->getId());
   pbPid.insert(destPid);

   //APF_ITERATE(dynEntArr, adjF, fIt) {
   for(int i=0; i<numFaces; i++) {
      assert(m->getDimension()-1 == getEntDim(m, adjF[i]));
      apf::Parts resPid; 
      m->getResidence(adjF[i], resPid);
      if ( resPid == pbPid ) 
         numPbFaces++;
   }
   return numPbFaces;
}

#endif
