/*
  This file is part of nbody.

  Copyright (C) 2025 Fredy W. Aquino

  nbody is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  nbody is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with nbody.  If not, see <https://www.gnu.org/licenses/>.

  Description : Penetration of a lossless dielectric sphere by a plane wave
                Using FDTD (Finite-Difference Time-Domain)
                to solve Maxwell's equations in time domain
                using modified Yee's algorithm.
                Computing |Ey| / [Einc| (ey1) vs j and |Ez| / [Einc| (ez1) vs j
                within lossless dielectric sphere.
  Reference   : Computational Electromagnetics with Matlab by Matthew N.O. Sadiku
  Date        : 02-23-26
*/

#include <cmath>
#include <numbers>
#include <string>
#include <vector>
#include "fdtd.h"

using namespace std;

const constexpr static double Pi = std::numbers::pi;

namespace fdtd_auxiliar {
int mod(int x, int y) {
  if (y == 0) {
    return x;
  }
  // MATLAB mod definition: m = x - floor(x./y) .* y
  return x - y * std::floor((double)x / (double)y);
}

template<typename T>
void ndgrid_3d(
  const size_t imax,
  const size_t jmax,
  const size_t kmax,
  std::vector<T>& X_grid,
  std::vector<T>& Y_grid,
  std::vector<T>& Z_grid) {
  const size_t nx = imax + 2;
  const size_t ny = jmax + 2;
  const size_t nz = kmax + 2;
  X_grid.resize(nx * ny * nz);
  Y_grid.resize(nx * ny * nz);
  Z_grid.resize(nx * ny * nz);

  // Fill X_grid
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      for (size_t k = 0; k < nz; ++k) {
        X_grid[i * ny * nz + j * nz + k] = (T)i;
      }
    }
  }
  // Fill Y_grid
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      for (size_t k = 0; k < nz; ++k) {
        Y_grid[i * ny * nz + j * nz + k] = (T)j;
      }
    }
  }
  // Fill Z_grid
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      for (size_t k = 0; k < nz; ++k) {
        Z_grid[i * ny * nz + j * nz + k] = (T)k;
      }
    }
  }
}
} // end namespace

template<typename T>
void fdtd::fdtd<T>::compute_media_parameters() {
  const T e0 = (1e-9) / (36 * Pi);
  const T u0 = (1e-7) * 4 * Pi;
  const T dt = delta_ / (2 * cl_);
  const T r = dt / e0;
  const T dt2 = dt * dt;
  const T delta2 = delta_ * delta_;
  const T ra = dt2 / (u0 * e0 * delta2);
  rb_ = dt / (u0 * delta_);
  tpifdt_ = 2.0 * Pi * f_ * dt;
  for (int i = 0; i != (int)sig_.size(); ++i) {
    ca_[i] = 1 - r * sig_[i] / er_[i];
    cb_[i] = ra / er_[i];
    cbmrb_[i] = cb_[i] / rb_;
  }
}

template<typename T>
void fdtd::fdtd<T>::compute_media_arrays() {
  std::vector<T> mx, my, mz;
  fdtd_auxiliar::ndgrid_3d(imax_, jmax_, kmax_, mx, my, mz);
  for (size_t i = 0; i != (size_t)nx_; ++i) {
    for (size_t j = 0; j != (size_t)ny_; ++j) {
      for (size_t k = 0; k != (size_t)nz_; ++k) {
        const int ind = i * ny_ * nz_ + j * nz_ + k;
        const T mx_oi = mx[ind] - oi_;
        const T my_oj = my[ind] - oj_;
        const T mz_ok = mz[ind] - ok_;
        const T mxo = mx_oi * mx_oi;
        const T mxd = (mx_oi + 0.5) * (mx_oi + 0.5);
        const T myo = my_oj * my_oj;
        const T myd = (my_oj + 0.5) * (my_oj + 0.5);
        const T mzo = mz_ok * mz_ok;
        const T mzd = (mz_ok + 0.5) * (mz_ok + 0.5);
        ixmed_[i * ny_ * nz_ + j * nz_ + k] = (sqrt(mxd + myo + mzo) <= radius_);
        iymed_[i * ny_ * nz_ + j * nz_ + k] = (sqrt(mxo + myd + mzo) <= radius_);
        izmed_[i * ny_ * nz_ + j * nz_ + k] = (sqrt(mxo + myo + mzd) <= radius_);
      }
    }
  }
}

template<typename T>
void fdtd::fdtd<T>::fdtd_solver() {
  int ncur = 3;
  int npr1 = 2;
  int npr2 = 1;
  for (int nn = 1; nn <= nnmax_; ++nn) { // time loop
    npr2 = npr1;
    npr1 = ncur;
    ncur = fdtd_auxiliar::mod(ncur, 3) + 1;
    const int shift1 = nmax_ + 1;
    const int shift2 = (kmax_ + 2) * shift1;
    const int shift3 = (jmax_ + 2) * shift2;
    for (int k = 0; k <= kmax_; ++k) { // z loop
      for (int j = 0; j <= jmax_; ++j) { // y loop
        for (int i = 0; i <= imax_; ++i) { // x loop
          // (ii)-apply soft lattice truncation conditions
          //---x=delta/2
          if (i == 0) {
            const int ind1 = (0 + 1) * shift3 + (j + 1) * shift2 + (k + 1) * shift1 + ncur;
            const int ind2 = (1 + 1) * shift3 + (j + 1) * shift2 +  k      * shift1 + npr2;
            const int ind3 = (1 + 1) * shift3 + (j + 1) * shift2 + (k + 1) * shift1 + npr2;
            const int ind4 = (1 + 1) * shift3 + (j + 1) * shift2 + (k + 2) * shift1 + npr2;
            if ((k != kmax_) && (k != 0)) {
              hy_[ind1] = (hy_[ind2] + hy_[ind3] + hy_[ind4]) / 3;
              hz_[ind1] = (hz_[ind2] + hz_[ind3] + hz_[ind4]) / 3;
            } else {
              if (k == kmax_) {
                const int ind5 = (0 + 1) * shift3 + (j + 1) * shift2 + (kmax_ + 1) * shift1 + ncur;
                const int ind6 = (1 + 1) * shift3 + (j + 1) * shift2 +  kmax_      * shift1 + npr2;
                const int ind7 = (1 + 1) * shift3 + (j + 1) * shift2 + (kmax_ + 1) * shift1 + npr2;
                hy_[ind5] = (hy_[ind6] + hy_[ind7]) / 2;
                hz_[ind1] = (hz_[ind2] + hz_[ind3]) / 2;
              } else {
                const int ind8  = (1 + 1) * shift3 + (j + 1) * shift2 +     shift1 + npr2;
                const int ind9  = (1 + 1) * shift3 + (j + 1) * shift2 + 2 * shift1 + npr2;
                const int ind10 = (0 + 1) * shift3 + (j + 1) * shift2 +     shift1 + npr2;
                hy_[ind1] = (hy_[ind3] + hy_[ind4]) / 2;
                hz_[ind10] = (hz_[ind8] + hz_[ind9]) / 2;
              }
            }
          }
          // ---y=0
          if (j == 0) {
            const int ind11 = (i + 1) * shift3 + (0 + 1) * shift2 + (k + 1) * shift1 + ncur;
            const int ind12 = (i + 1) * shift3 + (1 + 1) * shift2 + (k + 1) * shift1 + npr2;
            ex_[ind11] = ex_[ind12];
            ez_[ind11] = ez_[ind12];
          } else {
            //---y=ymax
            if (j == jmax_) {
              const int ind13 = (i + 1) * shift3 + (jmax_ + 1) * shift2 + (k + 1) * shift1 + ncur;
              const int ind14 = (i + 1) * shift3 +  jmax_      * shift2 + (k + 1) * shift1 + npr2;
              ex_[ind13] = ex_[ind14];
              ez_[ind13] = ez_[ind14];
            }
          }
          //---z=0
          if(k == 0) {
            const int ind15 = (i + 1) * shift3 + (j + 1) * shift2 +     shift1 + ncur;
            const int ind16 =  i      * shift3 + (j + 1) * shift2 + 2 * shift1 + npr2;
            const int ind17 = (i + 1) * shift3 + (j + 1) * shift2 + 2 * shift1 + npr2;
            const int ind18 = (i + 2) * shift3 + (j + 1) * shift2 + 2 * shift1 + npr2;
            if ((i != 0) && (i != imax_)) {
              ex_[ind15] = (ex_[ind16] + ex_[ind17] + ex_[ind18]) / 3;
              ey_[ind15] = (ey_[ind16] + ey_[ind17] + ey_[ind18]) / 3;
            } else {
              if (i == 0) {
                const int ind19 =     shift3 + (j + 1) * shift2 +     shift1 + ncur;
                const int ind20 =     shift3 + (j + 1) * shift2 + 2 * shift1 + npr2;
                const int ind21 = 2 * shift3 + (j + 1) * shift2 + 2 * shift1 + npr2;
                ex_[ind19] = (ex_[ind20] + ex_[ind21]) / 2;
                ey_[ind15] = (ey_[ind17] + ey_[ind18]) / 2;
              } else {
                ex_[ind15] = (ex_[ind16] + ex_[ind17]) / 2;
                ey_[ind15] = (ey_[ind16] + ey_[ind17]) / 2;
              }
            }
          }
          //  (iii)  apply  fd/td algorithm
          const int ind27 = (i + 1) * shift3 + (j + 1) * shift2 + (k + 1) * shift1 + ncur;
          const int ind28 = (i + 1) * shift3 + (j + 1) * shift2 + (k + 1) * shift1 + npr1;
          const int ind29 = (i + 1) * shift3 + (j + 1) * shift2 + (k + 2) * shift1 + npr1;
          const int ind30 = (i + 2) * shift3 + (j + 1) * shift2 + (k + 1) * shift1 + npr1;
          const int ind31 = (i + 1) * shift3 + (j + 2) * shift2 + (k + 1) * shift1 + npr1;
          const int ind32 = (i + 1) * shift3 +  j      * shift2 + (k + 1) * shift1 + ncur;
          const int ind33 =  i      * shift3 + (j + 1) * shift2 + (k + 1) * shift1 + ncur;
          const int ind34 = (i + 1) * shift3 + (j + 1) * shift2 +  k      * shift1 + ncur;
          //-----a. hx  generation:
          hx_[ind27] = hx_[ind28] + rb_ * (ey_[ind29] - ey_[ind28] + ez_[ind28] - ez_[ind31]);
          //-----b. hy  generation:
          hy_[ind27] = hy_[ind28] + rb_ * (ez_[ind30] - ez_[ind28] + ex_[ind28] - ex_[ind29]);
          //-----c. hz  generation:
          hz_[ind27] = hz_[ind28] + rb_ * (ex_[ind31] - ex_[ind28] + ey_[ind28] - ey_[ind30]);
          //---k=kmax ! symmetry
          if (k == kmax_) {
            const int ind25 = (i + 1) * shift3 + (j + 1) * shift2 +  kmax_      * shift1 + ncur;
            const int ind26 = (i + 1) * shift3 + (j + 1) * shift2 + (kmax_ + 1) * shift1 + ncur;
            hx_[ind26] = hx_[ind25];
            hy_[ind26] = hy_[ind25];
          }
          // -----d. ex  generation:
          const int ind_med = i * ny_ * nz_ + j * nz_ + k;
          if ((j != 0) && (j != jmax_) && (k != 0)) {
            const int m = ixmed_[ind_med];
            ex_[ind27] = ca_[m] * ex_[ind28] + cbmrb_[m] * (hz_[ind27] - hz_[ind32] + hy_[ind34] - hy_[ind27]);
          }
          //-----e. ey  generation:
          if(k != 0) {
            const int m = iymed_[ind_med];
            if (i != 0) {
              ey_[ind27] = ca_[m] * ey_[ind28] + cbmrb_[m] * (hx_[ind27] - hx_[ind34] + hz_[ind33] - hz_[ind27]);
            } else {
              ey_[ind27] = ca_[m] * ey_[ind28] + cbmrb_[m] * (hx_[ind27] - hx_[ind34] +   0        - hz_[ind27]);
            }
          }
          //-----f. ez  generation:
          if ((j != 0) && (j != jmax_)) {
            const int m = izmed_[ind_med];
            //  sig(ext)=0 for ez only from taflove[14]
            const T cam = (m == 0) ? 1.0 : ca_[m];
            if (i != 0) {
              ez_[ind27] = cam * ez_[ind28] + cbmrb_[m] * (hy_[ind27] - hy_[ind33] + hx_[ind32] - hx_[ind27]);
            } else {
              ez_[ind27] = cam * ez_[ind28] + cbmrb_[m] * (hy_[ind27] -     0      + hx_[ind32] - hx_[ind27]);
            }

            // (iv)  apply the plane-wave source
            if (j == js_) {
              const int ind37 = (i + 1) * shift3 + (js_ + 1) * shift2 + (k + 1) * shift1 + ncur;
              ez_[ind37] += sin(tpifdt_ * (T)nn);
            }
          }
          //---i=imax+1/2  ! symmetry
          if(i == imax_) {
            const int ind22 = (imax_ + 2) * shift3 + (j + 1) * shift2 + (k + 1) * shift1 + ncur;
            const int ind23 = (imax_ + 1) * shift3 + (j + 1) * shift2 + (k + 1) * shift1 + ncur;
            ey_[ind22] = ey_[ind23];
            ez_[ind22] = ez_[ind23];
          }
          //---k=kmax
          if(k == kmax_) {
            const int ind24 = (i + 1) * shift3 + (j + 1) * shift2 + (kmax_ + 2) * shift1 + ncur;
            const int ind25 = (i + 1) * shift3 + (j + 1) * shift2 +  kmax_      * shift1 + ncur;
            ex_[ind24] = ex_[ind25];
            ey_[ind24] = ey_[ind25];
          }
        } // x loop
// ********************************************************
//  step # 4 - retain the maximum absolute values during
//              the last half-wave
// ********************************************************
        if ((k == kmax_) && (nn > (nnmax_ - nhw_))) {
          const int ind35 = (imax_ + 1) * shift3 + (j + 1) * shift2 + kmax_ * shift1 + ncur;
          const T temp = fabs(ey_[ind35]);
          if (temp > ey1_[j]) {
            ey1_[j] = temp;
          }
          const int ind36 = (imax_ + 1) * shift3 + (j + 1) * shift2 + (kmax_ + 1) * shift1 + ncur;
          const T temp1 = fabs(ez_[ind36]);
          if (temp1 > ez1_[j]) {
            ez1_[j] = temp1;
          }
        }
      } // y loop
    } // z loop
  } // time loop
}

template<typename T>
fdtd::fdtd<T>::fdtd(const nlohmann::json input_pars0) :
  input_pars_(input_pars0), imax_(input_pars0["imax"]),
  jmax_(input_pars0["jmax"]), kmax_(input_pars0["kmax"]),
  nx_(imax_ + 2), ny_(jmax_ + 2), nz_(kmax_ + 2),
  nmax_(input_pars0["nmax"]), nnmax_(input_pars0["nnmax"]),
  nhw_(input_pars0["nhw"]), med_(input_pars0["med"]),
  js_(input_pars0["js"]), delta_(input_pars0["delta"]),
  cl_(input_pars0["cl"]), f_(input_pars0["f"]),
  oi_(input_pars0["oi"]), oj_(input_pars0["oj"]),
  ok_(input_pars0["ok"]), radius_(input_pars0["radius"]),
  er_(input_pars0["er"]), sig_(input_pars0["sig"]),
  ca_(2, 0), cb_(2, 0), cbmrb_(2, 0),
  ixmed_(nx_ * ny_ * nz_, 0), iymed_(nx_ * ny_ * nz_, 0),
  izmed_(nx_ * ny_ * nz_, 0),
  ey1_(jmax_ + 3, 0.0), ez1_(jmax_ + 3, 0),
  ex_((imax_ + 3) * (jmax_ + 3) * (kmax_ + 3) * (nmax_ + 1), 0),
  ey_((imax_ + 3) * (jmax_ + 3) * (kmax_ + 3) * (nmax_ + 1), 0),
  ez_((imax_ + 3) * (jmax_ + 3) * (kmax_ + 3) * (nmax_ + 1), 0),
  hx_((imax_ + 3) * (jmax_ + 3) * (kmax_ + 3) * (nmax_ + 1), 0),
  hy_((imax_ + 3) * (jmax_ + 3) * (kmax_ + 3) * (nmax_ + 1), 0),
  hz_((imax_ + 3) * (jmax_ + 3) * (kmax_ + 3) * (nmax_ + 1), 0)
{
  compute_media_parameters();
  compute_media_arrays();
  fdtd_solver();
}

template class fdtd::fdtd<double>;  // INSTANTIATION
