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

    #define SIZING_PARAMS iso, global, anisosize, t
    #define SIZING_DEFAULTS bool iso = false, double global = -1, double anisosize = -1, double t = -1
    EmbeddedShockFunction(ma::Mesh* m, gmi_model* g, std::list<gmi_ent*> surfs, SIZING_DEFAULTS) 
        : mesh(m), ref(g), shock_surfaces(surfs), test_iso(iso) {
        if (global > 0) {
            h_global = global;
            std::cout << "Overriding h_global size with valid value " << global << std::endl;
        }
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
    EmbeddedShockFunction(ma::Mesh* m, std::list<gmi_ent*> surfs, SIZING_DEFAULTS) : 
        EmbeddedShockFunction(m, m->getModel(), surfs, SIZING_PARAMS) {};
    EmbeddedShockFunction(ma::Mesh* m, SIZING_DEFAULTS) : 
        EmbeddedShockFunction(m, m->getModel(), {}, SIZING_PARAMS) {};
    #undef SIZING_PARAMS
    #undef SIZING_DEFAULTS

    void getValue(ma::Entity* vtx, ma::Matrix& frame, ma::Vector& scale);
    void getValue(apf::Vector3& pos, ma::Matrix& frame, ma::Vector& scale);
    double getMaxEdgeLengthAcrossShock();

    protected:

    double getZoneIsoSize(apf::Vector3 vtx, apf::Vector3 closest_pt, bool in_shock_band, bool& in_tip_ref);

    double getClosestPointAndNormal(double pos_arr[3], double cls_arr[3], double nrm_arr[3]);

    double h_global = 0.2113125;
    double norm_size = h_global/16;
    double thickness = norm_size*2;

    //double h_tip = h_global/4;
    double h_tip = norm_size;
    double h_upstream = 4 * h_global;

    bool test_iso;
    ma::Mesh* mesh;
    gmi_model* ref;
    double thickness_tol;
    std::list<gmi_ent*> shock_surfaces;
    int n_evals;

    private:

    double sizeLerp(double initial, double final, double x, double dist);

    bool inSphere(double c_x, double c_y, double c_z, double radius, apf::Vector3 pos, double& sphere_dist_sqr);

    bool inSphere(double c_x, double c_y, double c_z, double radius, apf::Vector3 pos);

}; // class EmbeddedSizeFunction

double EmbeddedShockFunction::sizeLerp(double initial, double final, double x, double dist) {
    return std::min(std::max(initial, x*(final-initial)/dist), final);
}

bool EmbeddedShockFunction::inSphere(double c_x, double c_y, double c_z, double radius, apf::Vector3 pos, double& sphere_dist_sqr) {
    sphere_dist_sqr = std::pow(pos.x()-c_x, 2) + std::pow(pos.y()-c_y, 2) + std::pow(pos.z()-c_z, 2);
    return sphere_dist_sqr < radius * radius;
}

bool EmbeddedShockFunction::inSphere(double c_x, double c_y, double c_z, double radius, apf::Vector3 pos) {
    double dummy = 0;
    return inSphere(c_x, c_y, c_z, radius, pos, dummy);
}

double EmbeddedShockFunction::getZoneIsoSize(apf::Vector3 pos, apf::Vector3 closest_pt, bool in_shock_band, bool& in_tip_ref) {
    double sphere_size = 8*h_tip;

    apf::Vector3 vecToPos = pos - closest_pt;

    double sphere_dist_sqr;
    in_tip_ref = inSphere(0, 0, 0, sphere_size, pos, sphere_dist_sqr);
    double h_tip_min = test_iso ? norm_size/2 : norm_size;
    if (test_iso || in_tip_ref) {
        return std::min(h_tip_min + (h_tip-h_tip_min)*(std::sqrt(sphere_dist_sqr)/sphere_size), 4*h_global);
    }

    // (h_norm-h_global)*exp(-abs(testx)/smooth_dist) + h_global;
    //#define EXP_SMOOTH(initial, final, x, distance) ( initial - final )*exp(-abs(x)/ distance ) + final

    double sphere_smooth_pos = std::sqrt(sphere_dist_sqr)-sphere_size;
    double sphere_smooth_dist = 8*h_global;
    //double sphere_smooth_size = EXP_SMOOTH(h_tip, h_global, sphere_smooth_pos, sphere_smooth_dist);
    double sphere_smooth_size = sizeLerp(h_tip, h_global, sphere_smooth_pos, sphere_smooth_dist);

    double sphere_fs_smooth_pos = std::sqrt(sphere_dist_sqr)-sphere_size;
    double sphere_fs_smooth_dist = 12*h_global;
    //double sphere_fs_smooth_size = EXP_SMOOTH(h_tip, h_upstream, sphere_smooth_pos, sphere_smooth_dist);
    double sphere_fs_smooth_size = sizeLerp(h_tip, h_upstream, sphere_smooth_pos, sphere_smooth_dist);

    double fs_smooth_pos = std::sqrt(std::abs(vecToPos * vecToPos))-0.5*thickness;
    double fs_smooth_dist = 12*h_global;
    //double fs_smooth_size = EXP_SMOOTH(sphere_smooth_size, sphere_fs_smooth_size, fs_smooth_pos, fs_smooth_dist);
    double fs_smooth_size = sizeLerp(sphere_smooth_size, sphere_fs_smooth_size, fs_smooth_pos, fs_smooth_dist);
    if (!in_shock_band && vecToPos.x() > -1e-3 * h_global) { // slight negative tolerance for outer outlet edge
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
    double zoneIsoSize = getZoneIsoSize(pos, cls_vec, shock_dist_square < thickness_tol, in_tip_ref);

    double h_n_max = h_global;
    if (inSphere(-1.25*h_global, 0, 0, h_global*2, pos)) {
        h_n_max = 2*h_tip;
    }
    if (inSphere(-2.75*h_global, 0, 0, h_global*4, pos)) {
        h_n_max = 4*h_tip;
    }

    scale[0] = test_iso || in_tip_ref ? zoneIsoSize : std::min(std::max(norm_size, shockDist - thickness/2), h_n_max);
    scale[1] = std::max(zoneIsoSize, scale[0]);
    scale[2] = std::max(zoneIsoSize, scale[0]);

    if(test_iso || in_tip_ref) {
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