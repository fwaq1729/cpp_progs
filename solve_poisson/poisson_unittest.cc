#include "poisson.h"
#include <vector>
#include <fstream>
#include <limits.h>
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>

using namespace std;

using json = nlohmann::json;

namespace {
TEST(test_poisson, test1) {
  const double tolerance = 1e-2;
  const double tolerance_tight = 1e-10;
  const nlohmann::json data =
    nlohmann::json::parse(std::ifstream("poisson_input.json"));
  const nlohmann::json data1 =
    nlohmann::json::parse(std::ifstream("poisson_data_for_test.json"));
  const int ncases = 3;
  for (int i = 0; i != ncases; ++i) {
    const std::string tag = "case_" + std::to_string(i + 1);
    const nlohmann::json input_pars = data.at(tag);
    poisson::Poisson solver = poisson::Poisson<double>(input_pars);

    const double h = solver.get_h();
    const double w = solver.get_relaxation_factor();
    const int n_iter_conv = solver.get_niter_conv();
    const double h_ref = data1.at(tag).at("h").get<double>();
    const double w_ref = data1.at(tag).at("relaxation_factor").get<double>();
    const int n_iter_conv_ref = data1.at(tag).at("niter_conv").get<int>();
    EXPECT_NEAR(h, h_ref, 1e-6);
    EXPECT_NEAR(w, w_ref, 1e-6);
    EXPECT_EQ(n_iter_conv, n_iter_conv_ref);
    const std::vector<double> v_num = solver.get_v_num();
    const std::vector<double> vnum_x = solver.extract_data_for_test(v_num); 
    const std::vector<double> vnum_x_ref =
      data1.at(tag).at("v_num_x").get<std::vector<double>>();
    for (int j = 0; j != (int)vnum_x.size(); ++j) {
      EXPECT_NEAR(vnum_x[j], vnum_x_ref[j], tolerance);
    }
    const std::vector<double> vnum_x1_ref =
      data1.at(tag).at("v_num_x1").get<std::vector<double>>();
    for (int j = 0; j != (int)vnum_x.size(); ++j) {
      EXPECT_NEAR(vnum_x[j], vnum_x1_ref[j], tolerance_tight);
    }
    
    solver.compute_exact_solution_case1();
    const std::vector<double> v_exact = solver.get_v_exact();
    const std::vector<double> vexact_x = solver.extract_data_for_test(v_exact);
    const std::vector<double> vexact_x_ref =
      data1.at(tag).at("v_exact_x").get<std::vector<double>>();
    for (int j = 0; j != (int)vexact_x.size(); ++j) {
      EXPECT_NEAR(vexact_x[j], vexact_x_ref[j], tolerance);
    }
    const std::vector<double> vexact_x1_ref =
      data1.at(tag).at("v_exact_x1").get<std::vector<double>>();
    for (int j = 0; j != (int)vexact_x.size(); ++j) {
      EXPECT_NEAR(vexact_x[j], vexact_x1_ref[j], tolerance_tight);
    }
  }
}
} // end-namespace
