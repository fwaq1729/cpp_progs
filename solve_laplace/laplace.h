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

#ifndef LAPLACE_H_
#define LAPLACE_H_
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

namespace laplace {
template<typename T>
class Laplace {
  public:
    explicit Laplace(const nlohmann::json input_pars);
    void laplace_solver();
    std::vector<T> get_x() const { return x_; }
    std::vector<T> get_y() const { return y_; }
    std::vector<T> get_v() const { return v_; }
  private:
    nlohmann::json input_pars_;
    std::vector<T> x_;
    std::vector<T> y_;
    int nd_;
    std::vector<T> v_;
    int np_;
    int ni_;
    int ne_;
    std::vector<int> ndp_;
    std::vector<int> val_;
    std::vector<int> nl_;
    std::vector<T> c_;
    std::vector<int> get_lf();
    void get_matrix_c();

  protected:
};
}
#endif // LAPLACE_H_

