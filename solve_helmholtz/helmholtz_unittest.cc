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

#include "helmholtz.h"
#include <vector>
#include <fstream>
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace {
TEST(test_helmholtz, test1) {
  const double tolerance = 1e-9;
  const nlohmann::json data = nlohmann::json::parse(std::ifstream("helmholtz_data_test.json"));
  const int ncases = 2;
  for (int i = 0; i != ncases; ++i) {
    const std::string tag = "case_" + std::to_string(i + 1);
    const nlohmann::json input_pars = data.at(tag);
    const std::vector<int> nx_ref = input_pars.at("nx").get<std::vector<int>>();
    const std::vector<int> ne_ref = input_pars.at("ne").get<std::vector<int>>();
    const std::vector<double> kcalc_ref = input_pars.at("kcalc").get<std::vector<double>>();
    const std::vector<double> error_ref = input_pars.at("error").get<std::vector<double>>();
    for (int j = 0; j != (int)nx_ref.size(); ++j) {
      helmholtz::Helmholtz solver = helmholtz::Helmholtz<double>(nx_ref[j], i);
      EXPECT_EQ(solver.get_ne(), ne_ref[j]);
      EXPECT_NEAR(solver.get_kcalc(), kcalc_ref[j], tolerance);
      EXPECT_NEAR(solver.get_error(), error_ref[j], tolerance);
    }
  }
}
} // end namespace
