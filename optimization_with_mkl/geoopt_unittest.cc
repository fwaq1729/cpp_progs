#include "geoopt.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <limits.h>
#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
using namespace std;

using json = nlohmann::json;

namespace {
#define EXPECT_VECTORS_NEAR(expected, actual, abs_error) \
    ASSERT_EQ(expected.size(), actual.size()) << "Vector sizes differ."; \
    for (size_t i = 0; i < expected.size(); ++i) { \
        EXPECT_NEAR(expected[i], actual[i], abs_error) << "at index " << i; \
    }
TEST(test_get_proj, test1)
{
  const double tolerance = 1e-15;
  const string filename = "data_Hproj.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<int> size_B = data["size_B"].get<vector<int>>();
  const int nrow_B = size_B[0];
  const int ncol_B = size_B[1];
  const vector<double> B = data["B_5"].get<vector<double>>();
  const vector<int> size_Binv = data["size_Binv"].get<vector<int>>();
  const int ncol_Binv = size_Binv[1];
  const vector<double> B_inv = data["Binv_5"].get<vector<double>>();
  const vector<double> proj_ref = data["proj_5"].get<vector<double>>();
  const int ndata = nrow_B;
  const vector<double> proj = geoopt_routines_with_mkl::get_proj(ndata, nrow_B, ncol_B, ncol_Binv, B, B_inv);
  EXPECT_VECTORS_NEAR(proj, proj_ref, tolerance);
}

TEST(test_get_g_new, test2)
{
  const double tolerance = 1e-15;
  const string filename = "data_Hproj.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<double> proj_ref = data["proj_5"].get<vector<double>>();
  const vector<double> g_interp = data["interpolated_g_5"].get<vector<double>>();
  const vector<double> g_new_ref = data["gnew_5"].get<vector<double>>();
  const vector<double> g_new = geoopt_routines_with_mkl::get_g_new(proj_ref, g_interp);
  EXPECT_VECTORS_NEAR(g_new, g_new_ref, tolerance);
}

TEST(test_get_Hproj, test3)
{
  const double tolerance = 1e-12;
  const string filename = "data_Hproj.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<double> proj_ref = data["proj_5"].get<vector<double>>();
  const vector<double> H = data["H_5"].get<vector<double>>();
  const vector<double> Hproj_ref = data["Hproj_5"].get<vector<double>>();
  const int ndata = (int)sqrt(Hproj_ref.size());
  const vector<double> Hproj = geoopt_routines_with_mkl::get_Hproj(ndata, proj_ref, H);
  EXPECT_VECTORS_NEAR(Hproj, Hproj_ref, tolerance);
}

TEST(test_do_quadratic_step, test4)
{
  const double tolerance = 1e-10;
  const string filename = "data_Hproj.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<double> g_new_ref = data["gnew_5"].get<vector<double>>();
  const vector<double> Hproj_ref = data["Hproj_5"].get<vector<double>>();
  const vector<double> dq_ref = data["dq_5"].get<vector<double>>();
  const vector<double> E3 = data["E_5"].get<vector<double>>();
  const double dE_ref = E3[0];
  double dE;
  vector<double> dq;
  bool on_sphere;
  double trust = 0.3;
  tie (dE, dq, on_sphere) = geoopt_routines_with_mkl::do_quadratic_step(trust, g_new_ref, Hproj_ref);
  EXPECT_EQ(on_sphere, false);
  EXPECT_NEAR(dE, dE_ref, tolerance);
  EXPECT_VECTORS_NEAR(dq, dq_ref, tolerance);
}

TEST(test_update_Hessian_BFGS, test5)
{
  const double tolerance = 1e-10;
  const string filename = "data_update_hessian.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<double> q = data["q_5"].get<vector<double>>();
  const vector<double> best_q = data["best_q_5"].get<vector<double>>();
  const vector<double> g = data["g_5"].get<vector<double>>();
  const vector<double> best_g = data["best_g_5"].get<vector<double>>();
  const vector<double> H = data["H_5"].get<vector<double>>();
  const vector<double> H_updated_ref = data["H_updated_5"].get<vector<double>>();
  const vector<double> H_updated = geoopt_routines_with_mkl::update_Hessian_BFGS(q, best_q, g, best_g, H);  
  EXPECT_VECTORS_NEAR(H_updated, H_updated_ref, tolerance);
}
}
