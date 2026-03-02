/*
  This file is part of helmholtz.

  Copyright (C) 2025 Fredy W. Aquino

  helmholtz is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  helmholtz is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with helmholtz.  If not, see <https://www.gnu.org/licenses/>.

  Description : Solving homogeneous (e.g. g = 0), scalar Helmholtz's equation
                \nabla^2 \Phi + k^2 \Phi = g
                using Finite Element (FE) for a rectangular region, using
                a grid mesh with triangular elements.
                Ref: Sakiku, M.N.O., Computation Electromagnetics with Matlab
  Date        : 03-02-26
*/

#include <cmath>
#include <numbers>
#include <algorithm>
#include <vector>
#include <array>
#include <numeric>
#include <iterator>
#include <mkl.h>
#include "helmholtz.h"

using namespace std;

const constexpr static double Pi = std::numbers::pi;

template<typename T>
vector<int> helmholtz::Helmholtz<T>::get_ndp(
  const array<int, 3>& bottom,
  const array<int, 3>& right,
  const array<int, 3>& top,
  const array<int, 3>& left) {
  vector<int> out;
  auto update_out = [](vector<int>& out, const array<int, 3>& side) {
    if (side[1] < 0) {
      for (int i = side[0]; i >= side[2]; i += side[1]) {
        out.push_back(i);
      }
    } else {
      for (int i = side[0]; i <= side[2]; i += side[1]) {
        out.push_back(i);
      }
    }
  };
  update_out(out, bottom);
  update_out(out, right);
  update_out(out, top);
  update_out(out, left);
  return out;
}

template<typename T>
vector<T> helmholtz::Helmholtz<T>::get_mat_for_mesh(
  vector<T>& seed,
  const int n,
  const bool transp) {
  const int n1 = seed.size();
  vector<T> out;

  if (transp) {
    for (int j = 0; j != n1; ++j) {
      for (int i = 0; i != n; ++i) {
        out.push_back(seed[j]);
      }
    }
  } else {
    for (int j = 0; j != n; ++j) {
      for (int i = 0; i != n1; ++i) {
        out.push_back(seed[i]);
      }
    }
  }
  return out;
}

template<typename T>
void helmholtz::Helmholtz<T>::update_vect(
  vector<int>& vec,
  vector<int>& update,
  const int inc_to_update,
  const int ini,
  const int inc,
  const int col) {
  for (int i = 0, ind = ini; i != (int)update.size(); ind += inc, ++i) {
    vec[ind * 3 + col] = update[i] + inc_to_update;
  }
}

template<typename T>
void helmholtz::Helmholtz<T>::gridmesh(
  int& ne,
  int& nd,
  int& np,
  vector<int>& nl,
  vector<T>& x,
  vector<T>& y,
  vector<int>& ndp,
  const int nx,
  const int ny,
  const vector<T>& dx,
  const vector<T>& dy)
{

// Generation of a rectangular mesh using triangular elements (nx by ny).
// ne = no. of elements in the mesh
// nd = no. of nodes in the mesh
// np = no. of boundary (prescribed) nodes
// x(i)  & y(i)  are global coordinates  of node i
// dx(i) & dy(i) are distances between nodes along x & y axes
// nl(i,j)       is the list of nodes for element i, j= 1,2,3 are local numbers
// ndp(i) = list of prescribed nodes i
//
// Ref:  Reddy, NJ, “An introduction to the finite element
//                   method.”, 1984, p. 436

  ne  = 2 * nx * ny;      // number of elements
  np  = 2 * (nx + ny);    // number of prescribed nodes (edges)
  const int nx1 = nx + 1; // number of nodes in the x direction
  const int ny1 = ny + 1; // number of nodes in the y direction
  nd  = nx1 * ny1;        // number of nodes in mesh

  vector<int> nl1;
  for (int i = 1; i <= nx1 * ny; ++i) {
    if (i % nx1 != 0) {
      nl1.push_back(i);
    }
  }

  vector<int> nl3;
  for (int i = nx + 2; i <= nx + nx1 * ny + 1; ++i) {
    if (i % nx1 != 0) {
      nl3.push_back(i);
    }
  }

  nl.resize(ne * 3);
  update_vect(nl, nl1, 0, 0, 2, 0);
  update_vect(nl, nl1, 0, 1, 2, 0);
  update_vect(nl, nl3, 1, 0, 2, 1);
  update_vect(nl, nl1, 1, 1, 2, 1);
  update_vect(nl, nl3, 0, 0, 2, 2);
  update_vect(nl, nl3, 1, 1, 2, 2);

  vector<T> dx_ac;
  dx_ac.push_back(0.0);
  std::partial_sum(dx.begin(), dx.end(), std::back_inserter(dx_ac));
  vector<T> dy_ac;
  dy_ac.push_back(0.0);
  std::partial_sum(dy.begin(), dy.end(), std::back_inserter(dy_ac));
  x = get_mat_for_mesh(dx_ac, ny1, false);
  y = get_mat_for_mesh(dy_ac, nx1, true);
  const array<int, 3> bottom = {1, 1, nx1};
  const array<int, 3> right  = {2 * nx1, nx1, nd};
  const array<int, 3> top    = {nd - 1, -1, nd - nx};
  const array<int, 3> left   = {nd - 2 * nx - 1, -nx1, nx + 2};
  ndp = get_ndp(bottom, right, top, left);
}

template<typename T>
void helmholtz::Helmholtz<T>::helmholtz_solver() {
  const int ny = (choice_ == 0) ? nx_ : (choice_ == 1) ? 2 * nx_: -1;
  const T aa = 1.0;
  const T bb = (choice_ == 0) ? 1.0 : (choice_ == 1) ? 2.0 : -1;
  if (!(choice_ == 0 || choice_ == 1)) {
    throw std::runtime_error("ny or bb = -1, wrong choice, choice ne 0 or 1");
  }

  const T deltax = aa / nx_;
  const T deltay = bb / ny;
  vector<T> dx(nx_, deltax);
  vector<T> dy(ny, deltay);

  int nd, np;
  vector<int> nl, ndp;
  vector<T> x, y;
  gridmesh(ne_, nd, np, nl, x, y, ndp, nx_, ny, dx, dy);
  vector<T> c(nd * nd, 0);
  vector<T> t(nd * nd, 0);
  array<T, 9> t1 = {2,1,1,
                    1,2,1,
                    1,1,2};
  for (int i = 0; i < ne_; ++i) {
    array<T, 3> xl;
    array<T, 3> yl;
    for (int k = 0; k != 3; ++k) {
      const int ind = nl[3 * i + k] - 1;
      xl[k] = x[ind];
      yl[k] = y[ind];
    }
    array<T, 3> p = {yl[1] - yl[2], yl[2] - yl[0], yl[0] - yl[1]};
    array<T, 3> q = {xl[2] - xl[1], xl[0] - xl[2], xl[1] - xl[0]};
    const T area = 0.5 * fabs(p[1] * q[2] - q[1] * p[2]);
    for (int j = 0; j != 3; ++j) {
      const int ir = nl[3 * i + j] - 1;
      for (int l = 0; l != 3; ++l) {
        const int ic = nl[3 * i + l] - 1;
        c[ir * nd + ic] += (p[j] * p[l] + q[j] * q[l]) / (4 * area);
        t[ir * nd + ic] += area / 12.0 * t1[j * 3 + l];
      }
    }
  }

  vector<int> lf;
  for (int i = 1; i <= nd; ++i) {
    auto it = std::find(ndp.begin(), ndp.end(), i);
    if (it == ndp.end()) {
      lf.push_back(i - 1);
    }
  }

  vector<T> c_ff, t_ff;
  for (int& i : lf) {
    for (int& j : lf) {
      c_ff.push_back(c[i * nd + j]);
      t_ff.push_back(t[i * nd + j]);
    }
  }

  const int n_lf = lf.size();
  const int lda = n_lf;
  const int ldb = n_lf;
  vector<int> ipiv(n_lf);
  vector<T> a(n_lf * n_lf);
  copy_n(c_ff.data(), n_lf * n_lf, a.data());
  int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n_lf, n_lf, t_ff.data(), lda, ipiv.data(), a.data(), ldb);
  if (info < 0) {
    std::ostringstream tag;
    tag << "In helmholtz_solver:: Failing linear solver, dgesv info = " << info;
    throw std::runtime_error(tag.str());
  }
  vector<T> alam_r(n_lf), alam_i(n_lf);
  T vl[1], vr[1];
  MKL_INT ldvl = 1, ldvr = n_lf;
  MKL_INT n_lf1 = (MKL_INT)n_lf;
  // Using MKL dgeev to compute only eigenvalues (alam_r) we discard imaginary part (alam_i) which is zero.
  const int info1 = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', n_lf1, a.data(), n_lf1, alam_r.data(), alam_i.data(), vl, ldvl, vr, ldvr);
  if (info1 < 0) {
    std::ostringstream tag;
    tag << "In helmholtz_solver:: Failing finding eigenvalues, dgeev info1 = " << info1;
    throw std::runtime_error(tag.str());
  }

  typename std::vector<T>::iterator min_itr = std::min_element(alam_r.begin(), alam_r.end());
  const T alam_min = *min_itr;
  kcalc_ = sqrt(alam_min);
  const T kc = Pi * sqrt(1.0 / (aa * aa) + 1.0 / (bb * bb));
  error_ = 100.0 * (kcalc_ - kc) / kc;
}

template<typename T>
helmholtz::Helmholtz<T>::Helmholtz(const int nx, const int choice) :
  nx_(nx), choice_(choice)
{
  helmholtz_solver();
}

template class helmholtz::Helmholtz<double>;  // INSTANTIATION
