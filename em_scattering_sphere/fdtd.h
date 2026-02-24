/*
  This file is part of fdtd.

  Copyright (C) 2025 Fredy W. Aquino

  fdtd is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  fdtd is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with fdtd.  If not, see <https://www.gnu.org/licenses/>.

  Description : Penetration of a lossless dielectric sphere by a plane wave
                Using FDTD (Finite-Difference Time-Domain)
                to solve Maxwell's equations in time domain
                using modified Yee's algorithm.
                Computing |Ey| / [Einc| (ey1) vs j and |Ez| / [Einc| (ez1) vs j
                within lossless dielectric sphere.
  Reference   : Computational Electromagnetics with Matlab by Matthew N.O. Sadiku
  Date        : 02-23-26
*/
#ifndef FDTD_H_
#define FDTD_H_
#include <vector>
#include <nlohmann/json.hpp>

namespace fdtd {
template<typename T>
class fdtd {
  public:
    explicit fdtd(const nlohmann::json input_pars);
    void fdtd_solver();
    const T get_rb() { return rb_; }
    const T get_tpifdt() { return tpifdt_; }
    const std::vector<T>& get_ca() { return ca_; }
    const std::vector<T>& get_cb() { return cb_; }
    const std::vector<T>& get_cbmrb() { return cbmrb_; }
    const std::vector<T>& get_ey1() { return ey1_; }
    const std::vector<T>& get_ez1() { return ez1_; }
  private:
    nlohmann::json input_pars_;
    int imax_;
    int jmax_;
    int kmax_;
    int nx_;
    int ny_;
    int nz_;
    int nmax_;
    int nnmax_;
    int nhw_;
    int med_;
    int js_;
    T delta_;
    T cl_;
    T f_;
    T oi_;
    T oj_;
    T ok_;
    T radius_;
    std::vector<T> er_;
    std::vector<T> sig_;
    std::vector<T> ca_;
    std::vector<T> cb_;
    std::vector<T> cbmrb_;
    std::vector<T> ixmed_;
    std::vector<T> iymed_;
    std::vector<T> izmed_;
    std::vector<T> ey1_;
    std::vector<T> ez1_;
    std::vector<T> ex_;
    std::vector<T> ey_;
    std::vector<T> ez_;
    std::vector<T> hx_;
    std::vector<T> hy_;
    std::vector<T> hz_;
    T rb_;
    T tpifdt_;
    void compute_media_parameters();
    void compute_media_arrays();
  protected:
};
} // namespace fdtd

#endif // FDTD_H_
