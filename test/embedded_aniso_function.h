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

    EmbeddedShockFunction(ma::Mesh* m, gmi_model* g, std::list<gmi_ent*> surfs, double anisosize, double t) : 
        mesh(m), ref(g), shock_surfaces(surfs) {
        if(anisosize > 0) {
            norm_size = anisosize;
            std::cout << "Overriding normal size with valid value " << norm_size << std::endl;
        }
        if(t > 0) {
            thickness = t;
            std::cout << "Overriding thickness with valid value " << thickness << std::endl;
        }
        thickness_tol = thickness * thickness / 4;
    }
    EmbeddedShockFunction(ma::Mesh* m, std::list<gmi_ent*> surfs, double anisosize, double t) : 
        EmbeddedShockFunction(m, m->getModel(), surfs, anisosize, t) {};
    EmbeddedShockFunction(ma::Mesh* m, std::list<gmi_ent*> surfs) : 
        EmbeddedShockFunction(m, m->getModel(), surfs, -1, -1) {};

    void getValue(ma::Entity* vtx, ma::Matrix& frame, ma::Vector& scale);
    void getValue(apf::Vector3& pos, ma::Matrix& frame, ma::Vector& scale);
    double getMaxEdgeLengthAcrossShock();

    protected:

    double getZoneIsoSize(apf::Vector3 vtx, apf::Vector3 closest_pt, bool in_shock_band, bool& in_tip_ref);

    double getClosestPointAndNormal(double pos_arr[3], double cls_arr[3], double nrm_arr[3]);

    double h_global = 0.2113125;
    double norm_size = h_global/16;
    double thickness = 0.721796;

    ma::Mesh* mesh;
    gmi_model* ref;
    double thickness_tol;
    std::list<gmi_ent*> shock_surfaces;
    int n_evals;

}; // class EmbeddedSizeFunction

double EmbeddedShockFunction::getZoneIsoSize(apf::Vector3 pos, apf::Vector3 closest_pt, bool in_shock_band, bool& in_tip_ref) {

    double h_tip = h_global/4;
    //double h_tip = norm_size;
    double h_tip_min = h_global/16;
    double h_upstream = 4 * h_global;

    double sphere_size = 4*h_global;
    //apf::Vector3 sphere_cent(-h_global,0,0);
    apf::Vector3 sphere_cent(0,0,0);

    apf::Vector3 dist = pos - sphere_cent;
    apf::Vector3 vecToPos = pos - closest_pt;

    double sphere_dist_sqr = std::abs(dist * dist);
    if (sphere_dist_sqr < sphere_size*sphere_size) {
        in_tip_ref = true;
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
    if (!in_shock_band && vecToPos.x() > -1e-3) { // slight negative tolerance for outer outlet edge
        return fs_smooth_size;
    }
    return sphere_smooth_size;
}

void EmbeddedShockFunction::getValue(ma::Entity* vtx, ma::Matrix& frame, ma::Vector& scale) {
    apf::Vector3 pos;
    mesh->getPoint(vtx, 0, pos);
    getValue(pos, frame, scale);
}

void EmbeddedShockFunction::getValue(apf::Vector3& pos, ma::Matrix& frame, ma::Vector& scale) {
    double pos_arr[3];
    pos.toArray(pos_arr);

    double cls_arr[3];
    double nrm_arr[3];
    apf::Vector3 cls_vec;
    apf::Vector3 norm;
    double shock_dist_square = getClosestPointAndNormal(pos_arr, cls_arr, nrm_arr);
    cls_vec.fromArray(cls_arr);
    norm.fromArray(nrm_arr);

    // Setting scale 
    double shockDist = std::sqrt(shock_dist_square);
    bool in_tip_ref = false;
    #define LINE(mins, maxs, x, slope) std::min(std::max(slope*x,mins),maxs)
    double zoneIsoSize = getZoneIsoSize(pos, cls_vec, shock_dist_square < thickness_tol, in_tip_ref);
    scale[0] = in_tip_ref ? zoneIsoSize : LINE(norm_size, h_global, shockDist, 1);
    //scale[0] = zoneIsoSize;
    scale[1] = std::max(zoneIsoSize, scale[0]);
    scale[2] = std::max(zoneIsoSize, scale[0]);

    if(in_tip_ref) {
        frame[0][0] = 1; frame[0][1] = 0; frame[0][2] = 0;
        frame[1][0] = 0; frame[1][1] = 1; frame[1][2] = 0;
        frame[2][0] = 0; frame[2][1] = 0; frame[2][2] = 1;
    } else {
        // Negate largest component to get tangent.
        apf::Vector3 trial(norm[2], norm[1], norm[0]);
        int largest = trial[0] > trial[1] && trial[0] > trial[2] ? 0 : (trial[1] > trial[0] && trial[1] > trial[2] ? 1 : 2);
        trial[largest] *= -1;
        apf::Vector3 tan1 = apf::cross(norm, trial).normalize();
        apf::Vector3 tan2 = apf::cross(norm, tan1);

        frame[0][0] = nrm_arr[0]; frame[0][1] = tan1[0]; frame[0][2] = tan2[0];
        frame[1][0] = nrm_arr[1]; frame[1][1] = tan1[1]; frame[1][2] = tan2[1];
        frame[2][0] = nrm_arr[2]; frame[2][1] = tan1[2]; frame[2][2] = tan2[2];
    }

    n_evals++;
    if(lion_get_verbosity() >= 1 && n_evals % 500 == 0){
        std::cout << pos_arr[0] << " " << pos_arr[1] << " " << pos_arr[2] << " ";
        std::cout << cls_arr[0] << " " << cls_arr[1] << " " << cls_arr[2] << " ";
        std::cout << nrm_arr[0] << " " << nrm_arr[1] << " " << nrm_arr[2] << " ";
        //std::cout << tan1[0] << " " << tan1[1] << " " << tan1[2] << " ";
        //std::cout << tan2[0] << " " << tan2[1] << " " << tan2[2] << std::endl;
        std::cout << frame[0][1] << " " << frame[1][1] << " " << frame[2][1] << " ";
        std::cout << frame[0][2] << " " << frame[1][2] << " " << frame[2][2] << std::endl;
    }
}

double EmbeddedShockFunction::getClosestPointAndNormal(double pos_arr[3], double cls_arr[3], double nrm_arr[3]) {
    if (shock_surfaces.size() == 0) {
        doubleConeClosestPointAnalytic(pos_arr, cls_arr, nrm_arr, h_global);
        //return std::abs((pos-cls_vec)*(pos-cls_vec));
        return (pos_arr[0]-cls_arr[0])*(pos_arr[0]-cls_arr[0])+
            (pos_arr[1]-cls_arr[1])*(pos_arr[1]-cls_arr[1])+
            (pos_arr[2]-cls_arr[2])*(pos_arr[2]-cls_arr[2]);
    }
    
    double shock_dist_square = DBL_MAX;
    double cls_par_arr[2];
    PCU_ALWAYS_ASSERT(gmi_can_get_closest_point(ref));
    gmi_ent* closest_surf;
    for(gmi_ent* surf : shock_surfaces) {
        double cur_cls_arr[3];
        double cur_cls_par_arr[2];
        gmi_closest_point(ref, surf, pos_arr, cur_cls_arr, cur_cls_par_arr);
        double cur_shock_dist_square = (pos_arr[0]-cur_cls_arr[0])*(pos_arr[0]-cur_cls_arr[0])+
            (pos_arr[1]-cur_cls_arr[1])*(pos_arr[1]-cur_cls_arr[1])+
            (pos_arr[2]-cur_cls_arr[2])*(pos_arr[2]-cur_cls_arr[2]);
        if (cur_shock_dist_square < shock_dist_square) {
            shock_dist_square = cur_shock_dist_square;
            closest_surf = surf;
            cls_arr[0] = cur_cls_arr[0]; cls_arr[1] = cur_cls_arr[1]; cls_arr[2] = cur_cls_arr[2];
            cls_par_arr[0] = cur_cls_par_arr[0]; cls_par_arr[1] = cur_cls_par_arr[1];
        }
    }
    PCU_ALWAYS_ASSERT(gmi_has_normal(ref));
    gmi_normal(ref, closest_surf, cls_par_arr, nrm_arr);

    return shock_dist_square;
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

        double pos_arr[3], cls_arr[3], cls_par_arr[2];
        pointA.toArray(pos_arr);
        gmi_closest_point(ref, surf, pos_arr, cls_arr, cls_par_arr);
        closest.fromArray(cls_arr);

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