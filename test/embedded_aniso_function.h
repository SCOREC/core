#include <unistd.h>
#include <cstdlib>
#include <PCU.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfConvert.h>
#include <ma.h>
#include <float.h>

#include <CapstoneModule.h>
#include "cap_analytic_closest_point.h"
 
class EmbeddedShockFunction : public ma::AnisotropicFunction {
    public:

    EmbeddedShockFunction(ma::Mesh* m, gmi_model* g, std::list<gmi_ent*> surfs, double nsize,
        double AR, double h0, double thickness) : mesh(m),
        norm_size(nsize), init_size(h0), ref(g), shock_surfaces(surfs), thickness(thickness) {
        //shock_surface = reinterpret_cast<apf::ModelEntity*>(toGmiEntity(shock));
        thickness_tol = thickness * thickness / 4;
        tan_size = norm_size * AR;
    }
    EmbeddedShockFunction(ma::Mesh* m, std::list<gmi_ent*> surfs, double nsize,
        double AR, double h0, double thickness) : 
        EmbeddedShockFunction(m, m->getModel(), surfs, nsize, AR, h0, thickness) {};

    void getValue(ma::Entity* vtx, ma::Matrix& frame, ma::Vector& scale);
    double getMaxEdgeLengthAcrossShock();

    protected:

    double getZoneIsoSize(ma::Entity* vtx, apf::Vector3 closestPt, bool inShockBand, bool& inTipRef);
    double h_global = 0.2113125;
    ma::Mesh* mesh;
    gmi_model* ref;
    double thickness_tol, thickness, norm_size, init_size, tan_size;
    std::list<gmi_ent*> shock_surfaces;
    int nInShockBand;

}; // class EmbeddedSizeFunction

double EmbeddedShockFunction::getZoneIsoSize(ma::Entity* vtx, apf::Vector3 closestPt, bool inShockBand, bool& inTipRef) {
    apf::Vector3 pos;

    //double h_tip = h_global/4; // h_global/4
    double h_tip = norm_size;
    double h_tip_min = h_global/32;
    double h_upstream = 4 * h_global;

    mesh->getPoint(vtx, 0, pos);
    double sphere_size = 4*h_global;
    //apf::Vector3 sphere_cent(-h_global,0,0);
    apf::Vector3 sphere_cent(0,0,0);

    apf::Vector3 dist = pos - sphere_cent;
    apf::Vector3 vecToPos = pos - closestPt;

    double sphere_dist_sqr = std::abs(dist * dist);    
    if (sphere_dist_sqr < sphere_size*sphere_size) {
        inTipRef = true;
        return h_tip_min + (h_tip-h_tip_min)*(std::sqrt(sphere_dist_sqr)/sphere_size);
    }

    // (h_norm-h_global)*exp(-abs(testx)/smooth_dist) + h_global;
    #define EXP_SMOOTH(initial, final, x, distance) ( initial - final )*exp(-abs(x)/ distance ) + final

    double sphere_smooth_pos = std::sqrt(sphere_dist_sqr)-sphere_size;
    double sphere_smooth_dist = 4*h_global;
    double sphere_smooth_size = EXP_SMOOTH(h_tip, h_global, sphere_smooth_pos, sphere_smooth_dist);

    double sphere_fs_smooth_pos = std::sqrt(sphere_dist_sqr)-sphere_size;
    double sphere_fs_smooth_dist = 6*h_global;
    double sphere_fs_smooth_size = EXP_SMOOTH(h_tip, h_upstream, sphere_smooth_pos, sphere_smooth_dist);

    double fs_smooth_pos = std::sqrt(std::abs(vecToPos * vecToPos))-0.5*thickness;
    double fs_smooth_dist = 6*h_global;
    double fs_smooth_size = EXP_SMOOTH(sphere_smooth_size, sphere_fs_smooth_size, fs_smooth_pos, fs_smooth_dist);
    if (!inShockBand && vecToPos.x() > -1e-3) { // slight negative tolerance for outer outlet edge
        return fs_smooth_size;
    }
    return sphere_smooth_size;
}

void EmbeddedShockFunction::getValue(ma::Entity* vtx, ma::Matrix& frame, ma::Vector& scale) {
    apf::Vector3 pos;
    mesh->getPoint(vtx, 0, pos);
    double posArr[3];
    pos.toArray(posArr);

    double clsArr[3], clsParArr[2];
    double nrmArr[3];
    apf::Vector3 clsVec;
    apf::Vector3 norm;
    gmi_ent* closestSurf;
    PCU_ALWAYS_ASSERT(gmi_can_get_closest_point(ref));
    doubleConeClosestPointAnalytic(posArr, clsArr, nrmArr, h_global);
    clsVec.fromArray(clsArr);
    norm.fromArray(nrmArr);
    double shockDistSquare = std::abs((pos-clsVec)*(pos-clsVec));

    // Negate largest component to get tangent.
    apf::Vector3 trial(norm[2], norm[1], norm[0]);
    int largest = trial[0] > trial[1] && trial[0] > trial[2] ? 0 : (trial[1] > trial[0] && trial[1] > trial[2] ? 1 : 2);
    trial[largest] *= -1;
    apf::Vector3 tan1 = apf::cross(norm, trial).normalize();
    apf::Vector3 tan2 = apf::cross(norm, tan1);

    frame[0][0] = nrmArr[0]; frame[0][1] = tan1[0]; frame[0][2] = tan2[0];
    frame[1][0] = nrmArr[1]; frame[1][1] = tan1[1]; frame[1][2] = tan2[1];
    frame[2][0] = nrmArr[2]; frame[2][1] = tan1[2]; frame[2][2] = tan2[2];

    double shockDist = std::sqrt(shockDistSquare);
    //#define EXP_SMOOTH2(hi, hf, x, d, dAR) (hi-hf)*exp(-x*x/(-d*d/std::log((dAR-1)*hf/(hi-hf))))+hf
    //#define TANH_SMOOTH(min, max, x, slope) min+(max-min)*std::tanh(slope*x)
    //TANH_SMOOTH(norm_size, h_global, shockDist, 0.333d);
    #define LINE(mins, maxs, x, slope) std::min(std::max(slope*x,mins),maxs)
    bool inTipRef = false;
    double zoneIsoSize = getZoneIsoSize(vtx, clsVec, shockDistSquare < thickness_tol, inTipRef);
    scale[0] = inTipRef ? zoneIsoSize : LINE(norm_size, h_global, shockDist, 1);
    //scale[0] = zoneIsoSize;
    scale[1] = std::max(zoneIsoSize, scale[0]);
    scale[2] = std::max(zoneIsoSize, scale[0]);

    if(inTipRef) {
        frame[0][0] = 1; frame[0][1] = 0; frame[0][2] = 0;
        frame[1][0] = 0; frame[1][1] = 1; frame[1][2] = 0;
        frame[2][0] = 0; frame[2][1] = 0; frame[2][2] = 1;
    }

    if(lion_get_verbosity() >= 1 && nInShockBand % 500 == 0){
        std::cout << posArr[0] << " " << posArr[1] << " " << posArr[2] << " ";
        std::cout << clsArr[0] << " " << clsArr[1] << " " << clsArr[2] << " ";
        std::cout << nrmArr[0] << " " << nrmArr[1] << " " << nrmArr[2] << " ";
        std::cout << tan1[0] << " " << tan1[1] << " " << tan1[2] << " ";
        std::cout << tan2[0] << " " << tan2[1] << " " << tan2[2] << std::endl;
    }
}

double EmbeddedShockFunction::getMaxEdgeLengthAcrossShock() {
    double max_length = -1.0;
    apf::MeshIterator* it = mesh->begin(1);
    apf::MeshEntity* e;
    while ((e = mesh->iterate(it))) {
        for (gmi_ent* surf : shock_surfaces) {
        apf::MeshEntity* adj_pts[2];
        mesh->getDownward(e, 0, adj_pts);
        apf::Vector3 pointA, pointB, closest;
        mesh->getPoint(adj_pts[0], 0, pointA);
        mesh->getPoint(adj_pts[1], 0, pointB);

        double posArr[3], clsArr[3], clsParArr[2];
        pointA.toArray(posArr);
        gmi_closest_point(ref, surf, posArr, clsArr, clsParArr);
        closest.fromArray(clsArr);

        apf::Vector3 vecA = pointA - closest;
        apf::Vector3 vecB = pointB - closest;
        if (vecA * vecB < 0) {
            max_length = std::max(max_length, (vecB-vecA).getLength());
            break;
        }
        }
    }
    mesh->end(it);
    return max_length;
}