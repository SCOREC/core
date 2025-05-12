#include <PCU.h>
#include <apf.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <set>
#include <map>
#include <vector>
#include <ma.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <apfConvert.h>
#include "sizeFields.h"

#define TOLERANCE 1e-9f

int vectorEqual(const apf::Vector3& v1, const apf::Vector3& v2) {
	bool a, b, c;
	a = std::fabs(v1[0] - v2[0]) < TOLERANCE;
	b = std::fabs(v1[1] - v2[1]) < TOLERANCE;
	c = std::fabs(v1[2] - v2[2]) < TOLERANCE;
	return a && b && c;
}

//Return the desired component of the vector normal of v1 and v2
double getVec3Normal(const apf::Vector3& v1, const apf::Vector3& v2, int dir) {
    double x, y, z;
    x = v1[1]*v2[2] - v1[2]*v2[1];
    y = v1[2]*v2[0] - v1[0]*v2[2];
    z = v1[0]*v2[1] - v1[1]*v2[0];
    apf::Vector3 out(x, y, z);
    return out[dir];
  }

int main(int argc, char* argv[]) {
  if(argc < 3) {
    printf("USAGE: ./a.out <dmg> <smb>\n");
    return 1;
  }

  MPI_Init(&argc, &argv);
  pcu::PCU* PCUObj = new pcu::PCU(MPI_COMM_WORLD);
  gmi_register_mesh();
  gmi_register_null();

  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2], PCUObj);

  //Adapt the mesh
  bool adaptOn = true;
  ma::AnisotropicFunction* sf = new Shock(m);
  apf::Field* frameField = nullptr;
  apf::Field* scaleField = nullptr;
  ma::Input *in = nullptr;
  if(adaptOn) {
    frameField = apf::createFieldOn(m, "adapt_frames", apf::MATRIX);
    scaleField = apf::createFieldOn(m, "adapt_scales", apf::VECTOR);

    ma::Entity *v;
    ma::Iterator* it = m->begin(0);
    while( (v = m->iterate(it)) ) {
      ma::Vector s;
      ma::Matrix f;
      sf->getValue(v, f, s);
      apf::setVector(scaleField, v, 0, s);
      apf::setMatrix(frameField, v, 0, f);
    }
    m->end(it);
    in = ma::makeAdvanced(ma::configure(m, scaleField, frameField));
    in->maximumIterations = 10;
    ma::adapt(in);
  }

  //First, go through the mesh faces on the cube. Use them to determine the normal axis for each cube face (X, Y, or Z)
  int posFaces, negFaces;
  std::map<int, int> cubeFaces;

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  std::map<int, int> dirsByFace;
  std::map<int, std::pair<int, std::vector<apf::MeshEntity*>>> seen;
  std::set<int> done;
  std::map<int, int> pos;
  std::map<int, int> neg;
  while((e = m->iterate(it))) {
    if(m->getModelType(m->toModel(e)) != 2) {
      continue;
    }
    int id = m->getModelTag(m->toModel(e));
    if(done.find(id) == done.end()) {
      //Identify and store two different mesh faces per cube face
      if(seen.find(id) == seen.end()) {
        seen[id] = std::pair<int, std::vector<apf::MeshEntity*>>(0, std::vector<apf::MeshEntity*>());
      }
      seen[id].first++;
      seen[id].second.push_back(e);
      
      //Once two mesh faces have been stored 
      if(seen[id].first == 2) {
        std::vector<apf::Vector3> one, two;
        apf::Downward down1, down2;
        m->getDownward(seen[id].second[0], 0, down1);
        m->getDownward(seen[id].second[1], 0, down2);

        //Use the stored vertices of the two mesh faces on this cube face to identify normal axis
        for(int i = 0; i < 3; i++) {
            one.push_back(ma::getPosition(m, down1[i]));
            two.push_back(ma::getPosition(m, down2[i]));
        }
        if(std::abs(one[2][0] - one[1][0]) < TOLERANCE && std::abs(two[2][0] - two[1][0]) < TOLERANCE) {
            dirsByFace[id] = 0;
        }
        else if(std::abs(one[2][1] - one[1][1]) < TOLERANCE && std::abs(two[2][1] - two[1][1]) < TOLERANCE) {
            dirsByFace[id] = 1;
        }
        else if(std::abs(one[2][2] - one[1][2]) < TOLERANCE && std::abs(two[2][2] - two[1][2]) < TOLERANCE) {
            dirsByFace[id] = 2;
        }
        //Catch all for the internal faces, we only care at boundary layer
        else {
            dirsByFace[id] = -1;
        }
        //This cube face is finished and does not need to be inspected again
        done.insert(id);
      } 
    }
  }

  FILE* writeFile;
  if(argc > 3) {
    writeFile = fopen(argv[3], "w");
  }
  else {
    writeFile = fopen("cubeData.txt", "w");
  }

  //Loop through each mesh face, determining whether its normal is positive or negative on the axis of its cube face
  int numMeshFaces = 0;
  int goodId = -1, badId = -1;
  int badDir, goodDir;
  apf::MeshEntity* badE, *goodE;
  it = m->begin(2);
  while((e = m->iterate(it))) {
    if(m->getModelType(m->toModel(e)) != 2) {
      continue;
    }
    int id = m->getModelTag(m->toModel(e));
    if(dirsByFace[id] == -1) {
        continue;
    }
    if(cubeFaces.find(id) == cubeFaces.end()) {
      cubeFaces[id] = 0;
      pos[id] = 0;
      neg[id] = 0;
    }
    cubeFaces[id] += 1;
    fprintf(writeFile, "Face %d Model face %d\n", reinterpret_cast<char*>(e) - ((char*)1), id);
    apf::Downward verts;
    m->getDownward(e, 0, verts);
    std::vector<apf::Vector3> all;
    for(int i = 0; i < 3; i++) {
      apf::Vector3 pos = ma::getPosition(m, verts[i]);
      fprintf(writeFile, "\tV%d: %e %e %e\n", (reinterpret_cast<char*>(verts[i]) - ((char*)1)) / 8, pos[0], pos[1], pos[2]);
      all.push_back(pos);
    }
    int direction = getVec3Normal(all[2] - all[1], all[1] - all[0], dirsByFace[id]) > 0 ? 1 : -1;
    fprintf(writeFile, "Direction: %d\n", direction);
    if(direction > 0) {
        pos[id]++;
        if(neg[id] > 0 && badId == -1) {
          badId = id;
          badE = e;
          badDir = direction;
        }
        else if(id == badId && goodId == -1){
          goodId = id;
          goodE = e;
          goodDir = direction;
        }
    }
    else {
        neg[id]++;
        if(pos[id] > 0 && badId == -1) {
          badId = id;
          badE = e;
          badDir = direction;
        }
        else if(id == badId && goodId == -1) {
          goodId = id;
          goodE = e;
          goodDir = direction;
        }
    }
    numMeshFaces++;
  }
  fclose(writeFile);

  //Report data by cube face
  for(std::map<int, int>::const_iterator i = cubeFaces.begin(); i != cubeFaces.end(); i++) {
    printf("\nCube Face %d\nNum Mesh Faces %d\n", i->first, i->second);
    printf("Direction +: %d -: %d\n", pos[i->first], neg[i->first]);
    printf("Normal Axis: %d\n", dirsByFace[i->first]);
  }

  apf::writeVtkFiles("cubeVtk", m);

  if(adaptOn) {
    m->writeNative("adapted.smb");
  }


  delete PCUObj;
  MPI_Finalize();

}