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
#include "fdtd.h"
#include <string>
#include <vector>
#include <fstream>
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

using namespace std;

using json = nlohmann::json;

namespace {
TEST(test_fdtd, test1) {
  const double tolerance = 1e-6;
  const nlohmann::json input_pars =
    nlohmann::json::parse(std::ifstream("fdtd_input.json"));
  const nlohmann::json data1 =
    nlohmann::json::parse(std::ifstream("fdtd_data_for_test.json"));
  fdtd::fdtd solver = fdtd::fdtd<double>(input_pars);
  EXPECT_NEAR(solver.get_rb(), data1.at("rb"), tolerance);
  EXPECT_NEAR(solver.get_tpifdt(), data1.at("tpifdt"), tolerance);

  const vector<double> ca_ref = data1.at("ca").get<std::vector<double>>();
  const vector<double> ca = solver.get_ca();
  for (int i = 0; i != (int)ca_ref.size(); ++i) {
    EXPECT_NEAR(ca[i], ca_ref[i], tolerance);
  }

  const vector<double> cb_ref = data1.at("cb").get<std::vector<double>>();
  const vector<double> cb = solver.get_cb();
  for (int i = 0; i != (int)cb_ref.size(); ++i) {
    EXPECT_NEAR(cb[i], cb_ref[i], tolerance);
  }

  const vector<double> cbmrb_ref = data1.at("cbmrb").get<std::vector<double>>();
  const vector<double> cbmrb = solver.get_cbmrb();
  for (int i = 0; i != (int)cbmrb_ref.size(); ++i) {
    EXPECT_NEAR(cbmrb[i], cbmrb_ref[i], tolerance);
  }

  const vector<double> ey1_ref = data1.at("ey1").get<std::vector<double>>();
  const vector<double> ey1 = solver.get_ey1();
  for (int i = 0; i != (int)ey1_ref.size(); ++i) {
    EXPECT_NEAR(ey1[i], ey1_ref[i], tolerance);
  }

  const vector<double> ez1_ref = data1.at("ez1").get<std::vector<double>>();
  const vector<double> ez1 = solver.get_ez1();
  for (int i = 0; i != (int)ez1_ref.size(); ++i) {
    EXPECT_NEAR(ez1[i], ez1_ref[i], tolerance);
  }
}
} // namespace
