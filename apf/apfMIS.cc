#include <time.h>
#include <cstdlib>
#include "apfMIS.h"
#include "apf.h"
#include <lionPrint.h>
#include <stdio.h>
namespace apf {    
  MIS::MIS(Mesh* mesh, int vtx_dim_, int edge_dim_)
    : m(mesh), vtx_dim(vtx_dim_), edge_dim(edge_dim_),
      ents(NULL), n(0),color(0) {
    srand(time(NULL));
  }
  MIS* initializeMIS(Mesh* mesh, int vtx_dim, int edge_dim) {
    MeshTag* coloring = mesh->findTag("coloring");
    if (coloring!=NULL) {
      lion_oprint(1,"[ERROR] Only one MIS can exist at a time.\n");
      return NULL;
    }
    MIS* mis = new MIS(mesh,vtx_dim,edge_dim);
    coloring = mesh->createIntTag("coloring",1);
    MeshIterator* vitr = mesh->begin(vtx_dim);
    MeshEntity* ent;
    int uncolored = -2;
    while ((ent = mesh->iterate(vitr))) {
      mesh->setIntTag(ent,coloring,&uncolored);
    }
    mesh->end(vitr);
    mesh->createIntTag("degrees",1);
    return mis;
  }
  void finalizeMIS(MIS* mis) {
    MeshTag* coloring = mis->m->findTag("coloring");
    MeshTag* degrees = mis->m->findTag("degrees");
    mis->m->destroyTag(coloring);
    mis->m->destroyTag(degrees);
    if (mis->ents) {
      delete [] mis->ents;
      mis->ents = NULL;
    }
    delete mis;
  }
  double getRandomNumber()
  {
    return (double)rand() / (double)RAND_MAX ;
  }
  void prepareForMIS(Mesh* m, MeshTag* coloring, int vtx_dim,
                     int uncolored,int neighbor);
  void selectVertices(Mesh* m, MeshTag* coloring, MeshTag* degrees,
                      int vtx_dim, int edge_dim, int temp_color, int neighbor);
  void trimColoring(Mesh* m,MeshTag* coloring, MeshTag* degrees,
                    int vtx_dim, int edge_dim, int uncolored,
                    int temp_color, int neighbor);
  int setColor(Mesh* m, MeshTag* coloring, int vtx_dim,
               int edge_dim, int temp_color, int neighbor, int color);
  void constructMIS(Mesh* m, MeshTag* coloring, int vtx_dim, int color,
                    MeshEntity** ents);
  bool getIndependentSet(MIS* mis) {
    Mesh* m = mis->m;
    int vtx_dim = mis->vtx_dim;
    int edge_dim = mis->edge_dim;
    int n=1;
    int uncolored = -2;
    int temp_color = -1;
    int neighbor = 0;
    mis->color++;
    MeshTag* coloring = m->findTag("coloring");
    MeshTag* degrees = m->findTag("degrees");
    if (mis->ents) {
      //Cleanup previous MIS if one exists
      delete [] mis->ents;
      mis->ents = NULL;
    }
    
    prepareForMIS(m,coloring,vtx_dim,uncolored,neighbor);
    
    mis->n=0;    
    while (n>0) {
      //Choose the vertices to try to color
      selectVertices(m, coloring, degrees, vtx_dim,
                     edge_dim, temp_color, neighbor);

      //Remove vertices that are adjacent with the new color
      trimColoring(m, coloring, degrees, vtx_dim, edge_dim,
                   uncolored, temp_color, neighbor);
      
      //Finalize vertices to be colored and set neighbors
      n = setColor(m, coloring, vtx_dim, edge_dim,
                   temp_color, neighbor, mis->color);
      mis->n+=n;
    }

    if (mis->n==0)
      return false;    

    mis->ents = new MeshEntity*[mis->n];

    constructMIS(m, coloring, vtx_dim, mis->color, mis->ents);
    return true;
  }

  void prepareForMIS(Mesh* m, MeshTag* coloring, int vtx_dim, int uncolored,
                     int neighbor) {
    MeshIterator* vitr = m->begin(vtx_dim);
    MeshEntity* ent;
    while ((ent = m->iterate(vitr))) {
      int tag;
      m->getIntTag(ent,coloring,&tag);
      if (tag<=neighbor)
        m->setIntTag(ent,coloring,&uncolored);
    }
    m->end(vitr);
  }
  void selectVertices(Mesh* m, MeshTag* coloring, MeshTag* degrees,
                      int vtx_dim, int edge_dim, int temp_color, int neighbor) {
    //Thread Parallize This:
    //Choose random set of vertices (temporarily color them to temp_color)
    MeshIterator* vitr = m->begin(vtx_dim);
    MeshEntity* ent;
    while ((ent = m->iterate(vitr))) {
      int tag;
      m->getIntTag(ent,coloring,&tag);
      if (tag>=neighbor)
        continue;
      //use set to uniquely retrieve neighbors
      std::set<MeshEntity*> neighbors;
      //get edge_dim adjacent of ent-> edge_adj
      Adjacent edge_adj;
      m->getAdjacent(ent,edge_dim,edge_adj);
      APF_ITERATE(Adjacent,edge_adj,edge_itr) {
        //get vtx_dim adjacent of edge_adj -> vtx_adj
        Adjacent vtx_adj;
        m->getAdjacent(*edge_itr,vtx_dim,vtx_adj);
        APF_ITERATE(Adjacent,vtx_adj,vtx_itr) {
          //count unique number of vtx_adj that are not colored -> count
          MeshEntity* vtx = *vtx_itr;
          if (vtx == ent)
            continue;
          int tag;
          m->getIntTag(vtx, coloring, &tag);
          if (tag>temp_color)
            continue;
          neighbors.insert(vtx);
        }
      }
      int count = neighbors.size();
      m->setIntTag(ent,degrees,&count);
      //set degree tag of this entity to be count
      if (count <= 1) {
        m->setIntTag(ent,coloring,&temp_color);
      }
      else {
        //Get random number in [0,1]
        double rand = getRandomNumber();
        if (rand > 1.0 / count) {
          m->setIntTag(ent, coloring, &temp_color);
        }
      }
    }
    m->end(vitr);
  }

  void trimColoring(Mesh* m,MeshTag* coloring, MeshTag* degrees,
                    int vtx_dim, int edge_dim, int uncolored,
                    int temp_color, int neighbor) {
    //Thread Parallize this:
    MeshIterator* eitr = m->begin(edge_dim);
    MeshEntity* ent;
    while ((ent = m->iterate(eitr))) {
      //Get vtx_dim adjacent of ent -> vtx_adj
      Adjacent vtx_adj;
      m->getAdjacent(ent,vtx_dim,vtx_adj);
      //Find vtx of max degree colored with temp_color
      int max_deg = -1;
      int max_index = -1;
      for (unsigned int i = 0;i < vtx_adj.size(); i++) {
        MeshEntity* vtx = vtx_adj[i];
        int color;
        m->getIntTag(vtx, coloring, &color);
        if (color != temp_color)
          continue;
        int deg;
        m->getIntTag(vtx,degrees,&deg);
        if (deg > max_deg) {
          max_deg = deg;
          max_index = i;
        }
      }
      if (max_index==-1)
        continue;
      //Set colors of all other vtxs to be a neighbor(0)
      for (unsigned int i = 0; i < vtx_adj.size(); i++) {
        if (i==(unsigned int)max_index)
          continue;
        MeshEntity* vtx = vtx_adj[i];
        int color;
        m->getIntTag(vtx, coloring, &color);
        if (color >= neighbor)
          continue;
        m->setIntTag(vtx,coloring,&uncolored);
      }
    }
    m->end(eitr);
  }
  int setColor(Mesh* m, MeshTag* coloring, int vtx_dim,
               int edge_dim, int temp_color, int neighbor, int color) {
    int n=0;
    //Set all temp_color to color and neighbors of temp_color to neighbor
    MeshIterator* vitr = m->begin(vtx_dim);
    MeshEntity* ent;
    while ((ent = m->iterate(vitr))) {
      int tag;
      m->getIntTag(ent,coloring,&tag);
      if (tag!=temp_color)
        continue;
      m->setIntTag(ent,coloring,&color);
      n++;
      //get edge_dim adjacent of ent-> edge_adj
      Adjacent edge_adj;
      m->getAdjacent(ent,edge_dim,edge_adj);
      APF_ITERATE(Adjacent,edge_adj,edge_itr) {
        //get vtx_dim adjacent of edge_adj -> vtx_adj
        Adjacent vtx_adj;
        m->getAdjacent(*edge_itr,vtx_dim,vtx_adj);
        APF_ITERATE(Adjacent,vtx_adj,vtx_itr) {
          MeshEntity* vtx = *vtx_itr;
          if (vtx == ent)
            continue;
          int tag;
          m->getIntTag(vtx, coloring, &tag);
          if (tag>temp_color)
            continue;
          m->setIntTag(vtx,coloring,&neighbor);
        }
      }
    }
    m->end(vitr);
    return n;
  }
  void constructMIS(Mesh* m, MeshTag* coloring, int vtx_dim, int color,
                    MeshEntity** ents) {
    //Thread Parallize this: (use an atomic to make index for ents)
    int n=0;
    MeshIterator* vitr = m->begin(vtx_dim);
    MeshEntity* ent;
    while ((ent = m->iterate(vitr))) {
      int c;
      m->getIntTag(ent, coloring, &c);
      if (c == color)
        ents[n++] = ent;
    }
    m->end(vitr);
  }
    
}
