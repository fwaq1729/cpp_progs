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

#ifndef HELMHOLTZ_H_
#define HELMHOLTZ_H_
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

namespace helmholtz {
template<typename T>
class Helmholtz {
  public:
    explicit Helmholtz(const int nx, const int choice);
    void helmholtz_solver();
    int get_ne() const { return ne_; }
    T get_kcalc() const { return kcalc_; }
    T get_error() const { return error_; }
  private:
    int nx_;
    int choice_;
    int ne_;
    T kcalc_;
    T error_;
    nlohmann::json input_pars_;

    std::vector<T> get_mat_for_mesh(
      std::vector<T>& seed,
      const int n,
      const bool transp);
    std::vector<int> get_ndp(
      const std::array<int, 3>& bottom,
      const std::array<int, 3>& right,
      const std::array<int, 3>& top,
      const std::array<int, 3>& left);
    void update_vect(
      std::vector<int>& vec,
      std::vector<int>& update,
      const int inc_to_update,
      const int ini,
      const int inc,
      const int col);
    void gridmesh(
      int& ne,
      int& nd,
      int& np,
      std::vector<int>& nl,
      std::vector<T>& x,
      std::vector<T>& y,
      std::vector<int>& ndp,
      const int nx,
      const int ny,
      const std::vector<T>& dx,
      const std::vector<T>& dy);
  protected:
};
}
#endif // HELMHOLTZ_H_
