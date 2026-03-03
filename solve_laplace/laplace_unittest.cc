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

#include "laplace.h"
#include <vector>
#include <fstream>
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

using namespace std;

using json = nlohmann::json;

namespace {
TEST(test_laplace, test1) {
  const double tolerance = 1e-9;
  const nlohmann::json input_pars =
    nlohmann::json::parse(std::ifstream("laplace_data_test.json"));
  laplace::Laplace solver = laplace::Laplace<double>(input_pars);
  const vector<double> x = solver.get_x();
  const vector<double> y = solver.get_y();
  const vector<double> v = solver.get_v();
  const vector<double> x_ref = input_pars.at("x").get<std::vector<double>>();
  const vector<double> y_ref = input_pars.at("y").get<std::vector<double>>();
  const vector<double> v_ref = input_pars.at("v").get<std::vector<double>>();
  const int ndata = x_ref.size();
  for (int i = 0; i != ndata; ++i) {
    EXPECT_NEAR(x[i], x_ref[i], tolerance);
    EXPECT_NEAR(y[i], y_ref[i], tolerance);
    EXPECT_NEAR(v[i], v_ref[i], tolerance);
  }
}
} // end-namespace
