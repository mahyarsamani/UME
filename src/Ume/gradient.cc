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
  \file Ume/gradient.cc
*/
#include "Ume/gradient.hh"
#include "Ume/Comm_MPI.hh"

#ifdef HOV
#include "hov.h"
#endif

#ifdef ANNOTATE
#include <annotate.h>
#endif

#include <iostream>

namespace Ume {

using DBLV_T = DS_Types::DBLV_T;
using VEC3V_T = DS_Types::VEC3V_T;
using VEC3_T = DS_Types::VEC3_T;

void gradzatp(Ume::SOA_Idx::Mesh &mesh, DBLV_T const &zone_field, VEC3V_T &point_gradient, bool measure) {
  auto const &csurf = mesh.ds->caccess_vec3v("corner_csurf");
  auto const &corner_volume = mesh.ds->caccess_dblv("corner_vol");
  auto const &point_normal = mesh.ds->caccess_vec3v("point_norm");
  auto const &c_to_p_map = mesh.ds->caccess_intv("m:c>p");
  auto const &c_to_z_map = mesh.ds->caccess_intv("m:c>z");
  auto const &corner_type = mesh.corners.mask;
  auto const &point_type = mesh.points.mask;

  int const pll = mesh.points.size();
  int const pl = mesh.points.local_size();
  int const cl = mesh.corners.local_size();

  DBLV_T point_volume(pll, 0.0);
  point_gradient.assign(pll, VEC3_T(0.0));

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(1);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int c = 0; c < cl; ++c) {
    if (corner_type[c] < 1)
      continue; // Only operate on interior corners
    int const z = c_to_z_map[c];
    int const p = c_to_p_map[c];
    point_volume[p] += corner_volume[c];
    point_gradient[p] += csurf[c] * zone_field[z];
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  mesh.points.gathscat(Ume::Comm::Op::SUM, point_volume);
  mesh.points.gathscat(Ume::Comm::Op::SUM, point_gradient);

  for (int p = 0; p < pl; ++p) {
    if (point_type[p] > 0) {
      // Internal point
      point_gradient[p] /= point_volume[p];
    } else if (point_type[p] == -1) {
      // Mesh boundary point
      double const ppdot = dotprod(point_gradient[p], point_normal[p]);
      point_gradient[p] =
        (point_gradient[p] - point_normal[p] * ppdot) / point_volume[p];
    }
  }
  mesh.points.scatter(point_gradient);
}

void gradzatz(Ume::SOA_Idx::Mesh &mesh, DBLV_T const &zone_field, VEC3V_T &zone_gradient, VEC3V_T &point_gradient, bool measure) {
  auto const &c_to_z_map = mesh.ds->caccess_intv("m:c>z");
  auto const &c_to_p_map = mesh.ds->caccess_intv("m:c>p");
  int const num_local_corners = mesh.corners.local_size();
  auto const &corner_type = mesh.corners.mask;
  auto const &corner_volume = mesh.ds->caccess_dblv("corner_vol");

  // Get the field gradient at each mesh point.
  gradzatp(mesh, zone_field, point_gradient, measure);

  DBLV_T zone_volume(mesh.zones.size(), 0.0);
  for (int corner_idx = 0; corner_idx < num_local_corners; ++corner_idx) {
    if (corner_type[corner_idx] < 1)
      continue; // Only operate on interior corners
    int const zone_idx = c_to_z_map[corner_idx];
    zone_volume[zone_idx] += corner_volume[corner_idx];
  }

  zone_gradient.assign(mesh.zones.size(), VEC3_T(0.0));

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(1);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int corner_idx = 0; corner_idx < num_local_corners; ++corner_idx) {
    if (corner_type[corner_idx] < 1)
      continue; // Only operate on interior corners
    int const zone_idx = c_to_z_map[corner_idx];
    int const point_idx = c_to_p_map[corner_idx];
    double const c_z_vol_ratio =
        corner_volume[corner_idx] / zone_volume[zone_idx];
    zone_gradient[zone_idx] += point_gradient[point_idx] * c_z_vol_ratio;
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  mesh.zones.scatter(zone_gradient);
}

void gradzatp_invert(Ume::SOA_Idx::Mesh &mesh, DBLV_T const &zone_field, VEC3V_T &point_gradient, bool measure) {
  auto const &csurf = mesh.ds->caccess_vec3v("corner_csurf");
  auto const &corner_volume = mesh.ds->caccess_dblv("corner_vol");
  auto const &point_normal = mesh.ds->caccess_vec3v("point_norm");
  auto const &p_to_c_map = mesh.ds->caccess_intrr("m:p>rc");
  auto const &c_to_z_map = mesh.ds->caccess_intv("m:c>z");
  auto const &point_type = mesh.points.mask;

  int const num_points = mesh.points.size();
  int const num_local_points = mesh.points.local_size();

  DBLV_T point_volume(num_points, 0.0);
  point_gradient.assign(num_points, VEC3_T(0.0));

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ_INVERT)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(1);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int point_idx = 0; point_idx < num_local_points; ++point_idx) {
    for (int const &corner_idx : p_to_c_map[point_idx]) {
      int const zone_idx = c_to_z_map[corner_idx];
      point_volume[point_idx] += corner_volume[corner_idx];
      point_gradient[point_idx] += csurf[corner_idx] * zone_field[zone_idx];
    }
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ_INVERT)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  mesh.points.gathscat(Ume::Comm::Op::SUM, point_volume);
  mesh.points.gathscat(Ume::Comm::Op::SUM, point_gradient);

  for (int point_idx = 0; point_idx < num_local_points; ++point_idx) {
    if (point_type[point_idx] > 0) {
      // Internal point
      point_gradient[point_idx] /= point_volume[point_idx];
    } else if (point_type[point_idx] == -1) {
      // Mesh boundary point
      double const ppdot =
          dotprod(point_gradient[point_idx], point_normal[point_idx]);
      point_gradient[point_idx] =
          (point_gradient[point_idx] - point_normal[point_idx] * ppdot) /
          point_volume[point_idx];
    }
  }
  mesh.points.scatter(point_gradient);
}

void gradzatz_invert(Ume::SOA_Idx::Mesh &mesh, DBLV_T const &zone_field, VEC3V_T &zone_gradient, VEC3V_T &point_gradient, bool measure) {
  auto const &z_to_c_map = mesh.ds->caccess_intrr("m:z>c");
  auto const &c_to_p_map = mesh.ds->caccess_intv("m:c>p");
  int const num_local_zones = mesh.zones.local_size();
  auto const &zone_type = mesh.zones.mask;
  auto const &corner_volume = mesh.ds->caccess_dblv("corner_vol");

  // Get the field gradient at each mesh point.
  gradzatp_invert(mesh, zone_field, point_gradient, measure);

  zone_gradient.assign(mesh.zones.size(), VEC3_T(0.0));

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ_INVERT)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(1);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int zone_idx = 0; zone_idx < num_local_zones; ++zone_idx) {
    if (zone_type[zone_idx] < 1)
      continue; // Only operate on local interior zones

    // Accumulate the (local) zone volume
    double zone_volume{0.0}; // Only need a local volume
    for (int const &corner_idx : z_to_c_map[zone_idx]) {
      zone_volume += corner_volume[corner_idx];
    }

    for (int const &corner_idx : z_to_c_map[zone_idx]) {
      int const point_idx = c_to_p_map[corner_idx];
      double const c_z_vol_ratio = corner_volume[corner_idx] / zone_volume;
      zone_gradient[zone_idx] += point_gradient[point_idx] * c_z_vol_ratio;
    }
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ_INVERT)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  mesh.zones.scatter(zone_gradient);
}

} // namespace Ume

#ifdef HOV
namespace Ume {

void gradzatp_hov(Ume::SOA_Idx::Mesh &mesh, DS_Types::DBLV_T const &zone_field, DS_Types::VEC3V_T &point_gradient, GradzatpHOVContext &ctx, bool measure) {
  auto const &csurf = mesh.ds->caccess_vec3v("corner_csurf");
  auto const &corner_volume = mesh.ds->caccess_dblv("corner_vol");
  auto const &point_normal = mesh.ds->caccess_vec3v("point_norm");
  auto const &corner_type = mesh.corners.mask;
  auto const &point_type = mesh.points.mask;

  int const pll = mesh.points.size();
  int const pl = mesh.points.local_size();
  int const cl = mesh.corners.local_size();

  ctx.point_volume.assign(pll, 0.0);
  ctx.pg_x.assign(pll, 0.0);
  ctx.pg_y.assign(pll, 0.0);
  ctx.pg_z.assign(pll, 0.0);

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(1);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int c = 0; c < cl; ++c) {
    if (corner_type[c] < 1)
      continue;
    
    hov_result_f64_u32_t pv_res, zf_res, pgx_res, pgy_res, pgz_res;
    hov_gather_f64_u32(&pv_res, &ctx.pair_pv, c);
    pv_res.data_val += corner_volume[c];
    hov_scatter_f64_u32(pv_res.data_val, &ctx.pair_pv, c);

    hov_gather_f64_u32(&zf_res, &ctx.pair_zf, c);

    hov_gather_f64_u32(&pgx_res, &ctx.pair_pg_x, c);
    pgx_res.data_val += csurf[c][0] * zf_res.data_val;
    hov_scatter_f64_u32(pgx_res.data_val, &ctx.pair_pg_x, c);

    hov_gather_f64_u32(&pgy_res, &ctx.pair_pg_y, c);
    pgy_res.data_val += csurf[c][1] * zf_res.data_val;
    hov_scatter_f64_u32(pgy_res.data_val, &ctx.pair_pg_y, c);

    hov_gather_f64_u32(&pgz_res, &ctx.pair_pg_z, c);
    pgz_res.data_val += csurf[c][2] * zf_res.data_val;
    hov_scatter_f64_u32(pgz_res.data_val, &ctx.pair_pg_z, c);
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  point_gradient.assign(pll, VEC3_T(0.0));
  for (int p = 0; p < pl; ++p) {
    point_gradient[p] = VEC3_T(ctx.pg_x[p], ctx.pg_y[p], ctx.pg_z[p]);
  }

  mesh.points.gathscat(Ume::Comm::Op::SUM, ctx.point_volume);
  mesh.points.gathscat(Ume::Comm::Op::SUM, point_gradient);

  for (int p = 0; p < pl; ++p) {
    if (point_type[p] > 0) {
      point_gradient[p] /= ctx.point_volume[p];
    } else if (point_type[p] == -1) {
      double const ppdot = dotprod(point_gradient[p], point_normal[p]);
      point_gradient[p] =
        (point_gradient[p] - point_normal[p] * ppdot) / ctx.point_volume[p];
    }
  }
  mesh.points.scatter(point_gradient);
}

void gradzatz_hov(Ume::SOA_Idx::Mesh &mesh, DBLV_T const &zone_field, VEC3V_T &zone_gradient, VEC3V_T &point_gradient, GradzatpHOVContext &p_ctx, GradzatzHOVContext &z_ctx, bool measure) {
  auto const &c_to_z_map = mesh.ds->caccess_intv("m:c>z");
  int const num_local_corners = mesh.corners.local_size();
  auto const &corner_type = mesh.corners.mask;
  auto const &corner_volume = mesh.ds->caccess_dblv("corner_vol");

  gradzatp_hov(mesh, zone_field, point_gradient, p_ctx, measure);

  z_ctx.zone_volume.assign(mesh.zones.size(), 0.0);
  for (int corner_idx = 0; corner_idx < num_local_corners; ++corner_idx) {
    if (corner_type[corner_idx] < 1)
      continue;
    int const zone_idx = c_to_z_map[corner_idx];
    z_ctx.zone_volume[zone_idx] += corner_volume[corner_idx];
  }

  z_ctx.zg_x.assign(mesh.zones.size(), 0.0);
  z_ctx.zg_y.assign(mesh.zones.size(), 0.0);
  z_ctx.zg_z.assign(mesh.zones.size(), 0.0);
  z_ctx.pg_x.assign(mesh.points.size(), 0.0);
  z_ctx.pg_y.assign(mesh.points.size(), 0.0);
  z_ctx.pg_z.assign(mesh.points.size(), 0.0);
  for (size_t i=0; i<mesh.points.size(); ++i) {
      z_ctx.pg_x[i] = point_gradient[i][0];
      z_ctx.pg_y[i] = point_gradient[i][1];
      z_ctx.pg_z[i] = point_gradient[i][2];
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(1);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int corner_idx = 0; corner_idx < num_local_corners; ++corner_idx) {
    if (corner_type[corner_idx] < 1)
      continue;
    
    hov_result_f64_u32_t zv_res;
    hov_gather_f64_u32(&zv_res, &z_ctx.pair_zv, corner_idx);
    double const c_z_vol_ratio = corner_volume[corner_idx] / zv_res.data_val;

    hov_result_f64_u32_t pgx_res, pgy_res, pgz_res;
    hov_gather_f64_u32(&pgx_res, &z_ctx.pair_pg_x, corner_idx);
    hov_gather_f64_u32(&pgy_res, &z_ctx.pair_pg_y, corner_idx);
    hov_gather_f64_u32(&pgz_res, &z_ctx.pair_pg_z, corner_idx);

    hov_result_f64_u32_t zgx_res, zgy_res, zgz_res;
    hov_gather_f64_u32(&zgx_res, &z_ctx.pair_zg_x, corner_idx);
    zgx_res.data_val += pgx_res.data_val * c_z_vol_ratio;
    hov_scatter_f64_u32(zgx_res.data_val, &z_ctx.pair_zg_x, corner_idx);

    hov_gather_f64_u32(&zgy_res, &z_ctx.pair_zg_y, corner_idx);
    zgy_res.data_val += pgy_res.data_val * c_z_vol_ratio;
    hov_scatter_f64_u32(zgy_res.data_val, &z_ctx.pair_zg_y, corner_idx);

    hov_gather_f64_u32(&zgz_res, &z_ctx.pair_zg_z, corner_idx);
    zgz_res.data_val += pgz_res.data_val * c_z_vol_ratio;
    hov_scatter_f64_u32(zgz_res.data_val, &z_ctx.pair_zg_z, corner_idx);
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  zone_gradient.assign(mesh.zones.size(), VEC3_T(0.0));
  for (size_t i=0; i<mesh.zones.size(); ++i) {
      zone_gradient[i] = VEC3_T(z_ctx.zg_x[i], z_ctx.zg_y[i], z_ctx.zg_z[i]);
  }
  mesh.zones.scatter(zone_gradient);
}


void gradzatp_invert_hov(Ume::SOA_Idx::Mesh &mesh, DBLV_T const &zone_field, VEC3V_T &point_gradient, GradzatpInvertHOVContext &ctx, bool measure) {
  auto const &csurf = mesh.ds->caccess_vec3v("corner_csurf");
  auto const &corner_volume = mesh.ds->caccess_dblv("corner_vol");
  auto const &point_normal = mesh.ds->caccess_vec3v("point_norm");
  auto const &p_to_c_map = mesh.ds->caccess_intrr("m:p>rc");
  auto const &point_type = mesh.points.mask;

  int const num_points = mesh.points.size();
  int const num_local_points = mesh.points.local_size();

  DBLV_T point_volume(num_points, 0.0);
  point_gradient.assign(num_points, VEC3_T(0.0));

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ_INVERT)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(1);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int point_idx = 0; point_idx < num_local_points; ++point_idx) {
    for (int const &corner_idx : p_to_c_map[point_idx]) {
      hov_result_f64_u32_t zf_res;
      hov_gather_f64_u32(&zf_res, &ctx.pair_zf, corner_idx);
      point_volume[point_idx] += corner_volume[corner_idx];
      point_gradient[point_idx] += csurf[corner_idx] * zf_res.data_val;
    }
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ_INVERT)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  mesh.points.gathscat(Ume::Comm::Op::SUM, point_volume);
  mesh.points.gathscat(Ume::Comm::Op::SUM, point_gradient);

  for (int point_idx = 0; point_idx < num_local_points; ++point_idx) {
    if (point_type[point_idx] > 0) {
      point_gradient[point_idx] /= point_volume[point_idx];
    } else if (point_type[point_idx] == -1) {
      double const ppdot =
          dotprod(point_gradient[point_idx], point_normal[point_idx]);
      point_gradient[point_idx] =
          (point_gradient[point_idx] - point_normal[point_idx] * ppdot) /
          point_volume[point_idx];
    }
  }
  mesh.points.scatter(point_gradient);
}


void gradzatz_invert_hov(Ume::SOA_Idx::Mesh &mesh, DBLV_T const &zone_field, VEC3V_T &zone_gradient, VEC3V_T &point_gradient, GradzatpInvertHOVContext &p_ctx, GradzatzInvertHOVContext &z_ctx, bool measure) {
  auto const &z_to_c_map = mesh.ds->caccess_intrr("m:z>c");
  int const num_local_zones = mesh.zones.local_size();
  auto const &zone_type = mesh.zones.mask;
  auto const &corner_volume = mesh.ds->caccess_dblv("corner_vol");

  gradzatp_invert_hov(mesh, zone_field, point_gradient, p_ctx, measure);

  z_ctx.pg_x.assign(mesh.points.size(), 0.0);
  z_ctx.pg_y.assign(mesh.points.size(), 0.0);
  z_ctx.pg_z.assign(mesh.points.size(), 0.0);
  for (size_t i=0; i<mesh.points.size(); ++i) {
      z_ctx.pg_x[i] = point_gradient[i][0];
      z_ctx.pg_y[i] = point_gradient[i][1];
      z_ctx.pg_z[i] = point_gradient[i][2];
  }

  zone_gradient.assign(mesh.zones.size(), VEC3_T(0.0));

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ_INVERT)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(1);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int zone_idx = 0; zone_idx < num_local_zones; ++zone_idx) {
    if (zone_type[zone_idx] < 1)
      continue; // Only operate on local interior zones

    double zone_volume{0.0};
    for (int const &corner_idx : z_to_c_map[zone_idx]) {
      zone_volume += corner_volume[corner_idx];
    }

    for (int const &corner_idx : z_to_c_map[zone_idx]) {
      double const c_z_vol_ratio = corner_volume[corner_idx] / zone_volume;
      hov_result_f64_u32_t pgx_res, pgy_res, pgz_res;
      hov_gather_f64_u32(&pgx_res, &z_ctx.pair_pg_x, corner_idx);
      hov_gather_f64_u32(&pgy_res, &z_ctx.pair_pg_y, corner_idx);
      hov_gather_f64_u32(&pgz_res, &z_ctx.pair_pg_z, corner_idx);
      
      zone_gradient[zone_idx][0] += pgx_res.data_val * c_z_vol_ratio;
      zone_gradient[zone_idx][1] += pgy_res.data_val * c_z_vol_ratio;
      zone_gradient[zone_idx][2] += pgz_res.data_val * c_z_vol_ratio;
    }
  }

  if (measure) {
#if defined(ANNOTATE) && defined(ROI_GRADZATZ_INVERT)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  mesh.zones.scatter(zone_gradient);
}

} // namespace Ume
#endif
