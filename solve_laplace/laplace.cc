/*
  This file is part of laplace.

  Copyright (C) 2025 Fredy W. Aquino

  laplace is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  laplace is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with laplace.  If not, see <https://www.gnu.org/licenses/>.

  Description : Solving Laplace's equation (\nabla^2 V = 0) using Finite Element (FE)               
                Ref: Sakiku, M.N.O., Computation Electromagnetics with Matlab
  Date        : 03-03-26
*/

#include <cmath>
#include <array>
#include <algorithm>
#include "laplace.h"

using namespace std;

template<typename T>
vector<int> laplace::Laplace<T>::get_lf() {
  vector<int> lf;
  for (int i = 1; i <= nd_; ++i) {
    auto it = std::find(ndp_.begin(), ndp_.end(), i);
    if (it == ndp_.end()) {
      lf.push_back(i - 1);
    }
  }
  return lf;
}

template<typename T>
void laplace::Laplace<T>::get_matrix_c() {
  for (int i = 0; i < ne_; ++i) {
    array<T, 3> xl;
    array<T, 3> yl;
    for (int k = 0; k != 3; ++k) {
      const int ind = nl_[3 * i + k] - 1;
      xl[k] = x_[ind];
      yl[k] = y_[ind];
    }
    array<T, 3> p = {yl[1] - yl[2], yl[2] - yl[0], yl[0] - yl[1]};
    array<T, 3> q = {xl[2] - xl[1], xl[0] - xl[2], xl[1] - xl[0]};
    const T area = 0.5 * fabs(p[1] * q[2] - q[1] * p[2]);
    for (int j = 0; j != 3; ++j) {
      const int ir = nl_[3 * i + j] - 1;
      for (int l = 0; l != 3; ++l) {
        const int ic = nl_[3 * i + l] - 1;
        c_[ir * nd_ + ic] += (p[j] * p[l] + q[j] * q[l]) / (4 * area);
      }
    }
  }
}

template<typename T>
void laplace::Laplace<T>::laplace_solver() {
  const vector<int> lf = get_lf();
  for (int k = 0; k != (int)ndp_.size(); ++k) {
    v_[ndp_[k] - 1] = (T)val_[k];
  }

  for (int n = 1; n <= ni_; ++n) {
    for (int ii = 0; ii != (int)lf.size(); ++ii) {
      vector<T> rowc(nd_);
      for (int k = 0; k != nd_; ++k) {
        rowc[k] = c_[lf[ii] * nd_ + k];
      }
      rowc[lf[ii]] = 0;
      T ac = 0;
      for (int l = 0; l != (int)rowc.size(); ++l) {
        ac += rowc[l] * v_[l];
      }
      v_[lf[ii]] = -1.0 / c_[lf[ii] * nd_ + lf[ii]] * ac;
    }
  }
}

template<typename T>
laplace::Laplace<T>::Laplace(const nlohmann::json input_pars0) :
  input_pars_(input_pars0),
  x_(input_pars0.at("x")), y_(input_pars0.at("y")), nd_(x_.size()), v_(nd_, 0),
  np_(input_pars0.at("np")), ni_(input_pars0.at("ni")), ne_(input_pars0.at("ne")),
  ndp_(input_pars0.at("ndp")), val_(input_pars0.at("val")), nl_(input_pars0.at("nl")),
  c_(nd_ * nd_, 0)
{
  get_matrix_c();
  laplace_solver();
}

template class laplace::Laplace<double>;  // INSTANTIATION
