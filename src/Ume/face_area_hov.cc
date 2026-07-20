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
  \file Ume/face_area_hov.cc

  HOV (Hardware Observability Vector) variant of face area computation.
  This file is compiled per-executable so that KERNEL_* macros are
  available for controlling ROI annotations.
*/
#include "Ume/face_area.hh"
#include "hov.h"

#ifdef ANNOTATE
extern "C" {
#include <annotate.h>
}
#endif

#include <iostream>

namespace Ume {

using Mesh = SOA_Idx::Mesh;
using DBLV_T = DS_Types::DBLV_T;
using INTV_T = DS_Types::INTV_T;
using VEC3V_T = DS_Types::VEC3V_T;

void calc_face_area_hov(SOA_Idx::Mesh &mesh, DBLV_T &face_area, FaceAreaHOVContext &ctx, bool measure) {
  auto const &side_type = mesh.sides.mask;
  auto const &surz = mesh.ds->caccess_vec3v("side_surz");

  int const sll = mesh.sides.size();
  int const sl = mesh.sides.local_size();

  std::fill(face_area.begin(), face_area.end(), 0.0);
  ctx.side_tag.assign(sll, 0);

  if (measure) {
#if defined(ANNOTATE) && defined(KERNEL_FACE_AREA)
    roi_begin_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(2);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  for (int s = 0; s < sl; ++s) {
    if (side_type[s] < 1)
      continue;
    if (ctx.side_tag[s] == 1)
      continue;

    hov_result_u32_u32_t fct_res;
    hov_gather_u32_u32(&fct_res, &ctx.pair_fct, s);

    if ((int)fct_res.data_val < 3) {
      double const side_area = vectormag(surz[s]);

      hov_result_f64_u32_t fa_res;
      hov_gather_f64_u32(&fa_res, &ctx.pair_fa, s);

      fa_res.data_val += side_area;
      hov_scatter_f64_u32(fa_res.data_val, &ctx.pair_fa, s);

      hov_scatter_u32_u32(1, &ctx.pair_st, s);
    }
  }

  if (measure) {
#if defined(ANNOTATE) && defined(KERNEL_FACE_AREA)
    roi_end_();
#ifdef SYNC_ON_ROI
    annotate_synchronize_(3);
#endif // SYNC_ON_ROI
#endif // ANNOTATE
  }

  mesh.faces.scatter(face_area);
}

} // namespace Ume
