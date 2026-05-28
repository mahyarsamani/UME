/*
  Copyright (c) 2023, Triad National Security, LLC. All rights reserved.

  This is open source software; you can redistribute it and/or modify it under
  the terms of the BSD-3 License. If software is modified to produce derivative
  works, such modified software should be clearly marked, so as not to confuse
  it with the version available from LANL. Full text of the BSD-3 License can be
  found in the LICENSE.md file, and the full assertion of copyright in the
  NOTICE.md file.
*/

/*!
   \file Ume/hov_context.hh
*/
#ifndef UME_HOV_CONTEXT_HH
#define UME_HOV_CONTEXT_HH 1

#ifdef HOV
#include "hov.h"
#include "Ume/SOA_Idx_Mesh.hh"

namespace Ume {

struct GradzatpHOVContext {
    DS_Types::DBLV_T point_volume;
    DS_Types::DBLV_T pg_x;
    DS_Types::DBLV_T pg_y;
    DS_Types::DBLV_T pg_z;
    
    hov_pair_t pair_zf;
    hov_pair_t pair_pv;
    hov_pair_t pair_pg_x;
    hov_pair_t pair_pg_y;
    hov_pair_t pair_pg_z;

    GradzatpHOVContext() = default;

    void init(Ume::SOA_Idx::Mesh &mesh, DS_Types::DBLV_T const &zone_field) {
        int const pll = mesh.points.size();
        int const cl = mesh.corners.local_size();
        
        point_volume.assign(pll, 0.0);
        pg_x.assign(pll, 0.0);
        pg_y.assign(pll, 0.0);
        pg_z.assign(pll, 0.0);

        auto const &c_to_p_map = mesh.ds->caccess_intv("m:c>p");
        auto const &c_to_z_map = mesh.ds->caccess_intv("m:c>z");

        pair_zf = hov_create_pair((void*)zone_field.data(), (void*)c_to_z_map.data(), sizeof(double), sizeof(int), cl);
        pair_pv = hov_create_pair((void*)point_volume.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), cl);
        pair_pg_x = hov_create_pair((void*)pg_x.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), cl);
        pair_pg_y = hov_create_pair((void*)pg_y.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), cl);
        pair_pg_z = hov_create_pair((void*)pg_z.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), cl);
    }

    void destroy() {
        hov_destroy_pair(&pair_zf);
        hov_destroy_pair(&pair_pv);
        hov_destroy_pair(&pair_pg_x);
        hov_destroy_pair(&pair_pg_y);
        hov_destroy_pair(&pair_pg_z);
    }
};

struct GradzatzHOVContext {
    DS_Types::DBLV_T zone_volume;
    DS_Types::DBLV_T zg_x;
    DS_Types::DBLV_T zg_y;
    DS_Types::DBLV_T zg_z;
    DS_Types::DBLV_T pg_x;
    DS_Types::DBLV_T pg_y;
    DS_Types::DBLV_T pg_z;

    hov_pair_t pair_zv;
    hov_pair_t pair_zg_x;
    hov_pair_t pair_zg_y;
    hov_pair_t pair_zg_z;
    hov_pair_t pair_pg_x;
    hov_pair_t pair_pg_y;
    hov_pair_t pair_pg_z;

    GradzatzHOVContext() = default;

    void init(Ume::SOA_Idx::Mesh &mesh) {
        int const num_local_corners = mesh.corners.local_size();
        
        zone_volume.assign(mesh.zones.size(), 0.0);
        zg_x.assign(mesh.zones.size(), 0.0);
        zg_y.assign(mesh.zones.size(), 0.0);
        zg_z.assign(mesh.zones.size(), 0.0);
        pg_x.assign(mesh.points.size(), 0.0);
        pg_y.assign(mesh.points.size(), 0.0);
        pg_z.assign(mesh.points.size(), 0.0);

        auto const &c_to_z_map = mesh.ds->caccess_intv("m:c>z");
        auto const &c_to_p_map = mesh.ds->caccess_intv("m:c>p");

        pair_zv = hov_create_pair((void*)zone_volume.data(), (void*)c_to_z_map.data(), sizeof(double), sizeof(int), num_local_corners);
        pair_zg_x = hov_create_pair((void*)zg_x.data(), (void*)c_to_z_map.data(), sizeof(double), sizeof(int), num_local_corners);
        pair_zg_y = hov_create_pair((void*)zg_y.data(), (void*)c_to_z_map.data(), sizeof(double), sizeof(int), num_local_corners);
        pair_zg_z = hov_create_pair((void*)zg_z.data(), (void*)c_to_z_map.data(), sizeof(double), sizeof(int), num_local_corners);
        
        pair_pg_x = hov_create_pair((void*)pg_x.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), num_local_corners);
        pair_pg_y = hov_create_pair((void*)pg_y.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), num_local_corners);
        pair_pg_z = hov_create_pair((void*)pg_z.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), num_local_corners);
    }

    void destroy() {
        hov_destroy_pair(&pair_zv);
        hov_destroy_pair(&pair_zg_x);
        hov_destroy_pair(&pair_zg_y);
        hov_destroy_pair(&pair_zg_z);
        hov_destroy_pair(&pair_pg_x);
        hov_destroy_pair(&pair_pg_y);
        hov_destroy_pair(&pair_pg_z);
    }
};

struct GradzatpInvertHOVContext {
    hov_pair_t pair_zf;

    GradzatpInvertHOVContext() = default;

    void init(Ume::SOA_Idx::Mesh &mesh, DS_Types::DBLV_T const &zone_field) {
        auto const &c_to_z_map = mesh.ds->caccess_intv("m:c>z");
        pair_zf = hov_create_pair((void*)zone_field.data(), (void*)c_to_z_map.data(), sizeof(double), sizeof(int), c_to_z_map.size());
    }

    void destroy() {
        hov_destroy_pair(&pair_zf);
    }
};

struct GradzatzInvertHOVContext {
    DS_Types::DBLV_T pg_x;
    DS_Types::DBLV_T pg_y;
    DS_Types::DBLV_T pg_z;
    
    hov_pair_t pair_pg_x;
    hov_pair_t pair_pg_y;
    hov_pair_t pair_pg_z;

    GradzatzInvertHOVContext() = default;

    void init(Ume::SOA_Idx::Mesh &mesh) {
        pg_x.assign(mesh.points.size(), 0.0);
        pg_y.assign(mesh.points.size(), 0.0);
        pg_z.assign(mesh.points.size(), 0.0);

        auto const &c_to_p_map = mesh.ds->caccess_intv("m:c>p");
        
        pair_pg_x = hov_create_pair((void*)pg_x.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), c_to_p_map.size());
        pair_pg_y = hov_create_pair((void*)pg_y.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), c_to_p_map.size());
        pair_pg_z = hov_create_pair((void*)pg_z.data(), (void*)c_to_p_map.data(), sizeof(double), sizeof(int), c_to_p_map.size());
    }

    void destroy() {
        hov_destroy_pair(&pair_pg_x);
        hov_destroy_pair(&pair_pg_y);
        hov_destroy_pair(&pair_pg_z);
    }
};

struct FaceAreaHOVContext {
    DS_Types::INTV_T side_tag;

    hov_pair_t pair_fct;
    hov_pair_t pair_fa;
    hov_pair_t pair_st;

    FaceAreaHOVContext() = default;

    void init(Ume::SOA_Idx::Mesh &mesh, DS_Types::DBLV_T &face_area) {
        int const sll = mesh.sides.size();
        int const sl = mesh.sides.local_size();

        side_tag.assign(sll, 0);

        auto const &face_comm_type = mesh.faces.comm_type;
        auto const &s_to_f_map = mesh.ds->caccess_intv("m:s>f");
        auto const &s_to_s2_map = mesh.ds->caccess_intv("m:s>s2");

        pair_fct = hov_create_pair((void*)face_comm_type.data(), (void*)s_to_f_map.data(), sizeof(int), sizeof(int), sl);
        pair_fa = hov_create_pair((void*)face_area.data(), (void*)s_to_f_map.data(), sizeof(double), sizeof(int), sl);
        pair_st = hov_create_pair((void*)side_tag.data(), (void*)s_to_s2_map.data(), sizeof(int), sizeof(int), sl);
    }

    void destroy() {
        hov_destroy_pair(&pair_fct);
        hov_destroy_pair(&pair_fa);
        hov_destroy_pair(&pair_st);
    }
};

} // namespace Ume
#endif

#endif
