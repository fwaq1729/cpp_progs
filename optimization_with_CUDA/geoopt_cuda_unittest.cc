#include "geoopt_cuda.h"
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
  const string filename = "../optimization_with_mkl/data_Hproj.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<int> size_B = data["size_B"].get<vector<int>>();
  const int nrow_B = size_B[0];
  const int ncol_B = size_B[1];
  const vector<double> B = data["B_5"].get<vector<double>>();
  vector<double> B1 = geoopt_routines_with_cuda::get_colwise_matrix(B, nrow_B, ncol_B);
  const vector<int> size_Binv = data["size_Binv"].get<vector<int>>();
  const int ncol_Binv = size_Binv[1];
  const vector<double> B_inv = data["Binv_5"].get<vector<double>>();
  vector<double> B_inv1 = geoopt_routines_with_cuda::get_colwise_matrix(B_inv, size_Binv[0], ncol_Binv);
  const vector<double> proj_ref = data["proj_5"].get<vector<double>>();
  const int ndata = nrow_B;
  vector<double> proj_ref1 = geoopt_routines_with_cuda::get_colwise_matrix(proj_ref, ndata, ndata);

  geoopt_routines_with_cuda::Geoopt_cuda opt_cuda = geoopt_routines_with_cuda::Geoopt_cuda();
  vector<double> proj = opt_cuda.get_proj(ndata, nrow_B, ncol_B, ncol_Binv, B1, B_inv1);
  EXPECT_VECTORS_NEAR(proj, proj_ref1, tolerance);
}

TEST(test_get_g_new, test2)
{
  const double tolerance = 1e-15;
  const string filename = "../optimization_with_mkl/data_Hproj.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<double> proj_ref = data["proj_5"].get<vector<double>>();
  const vector<double> g_interp = data["interpolated_g_5"].get<vector<double>>();
  const int ndata = g_interp.size();
  vector<double> proj_ref1 = geoopt_routines_with_cuda::get_colwise_matrix(proj_ref, ndata, ndata);

  const vector<double> g_new_ref = data["gnew_5"].get<vector<double>>();
  geoopt_routines_with_cuda::Geoopt_cuda opt_cuda = geoopt_routines_with_cuda::Geoopt_cuda();
  const vector<double> g_new = opt_cuda.get_g_new(proj_ref1, g_interp);
  EXPECT_VECTORS_NEAR(g_new, g_new_ref, tolerance);
}

TEST(test_get_Hproj, test3)
{
  const double tolerance = 1e-12;
  const string filename = "../optimization_with_mkl/data_Hproj.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<double> proj_ref = data["proj_5"].get<vector<double>>();
  const vector<double> H = data["H_5"].get<vector<double>>();
  const vector<double> Hproj_ref = data["Hproj_5"].get<vector<double>>();
  const int ndata = (int)sqrt(Hproj_ref.size());
  vector<double> proj_ref1 = geoopt_routines_with_cuda::get_colwise_matrix(proj_ref, ndata, ndata);
  vector<double> H1 = geoopt_routines_with_cuda::get_colwise_matrix(H, ndata, ndata);
  vector<double> Hproj_ref1 = geoopt_routines_with_cuda::get_colwise_matrix(Hproj_ref, ndata, ndata);
  geoopt_routines_with_cuda::Geoopt_cuda opt_cuda = geoopt_routines_with_cuda::Geoopt_cuda();
  const vector<double> Hproj = opt_cuda.get_Hproj(ndata, proj_ref1, H1);
  EXPECT_VECTORS_NEAR(Hproj, Hproj_ref1, tolerance);
}

TEST(test_update_Hessian_BFGS, test5)
{
  const double tolerance = 1e-10;
  const string filename = "../optimization_with_mkl/data_update_hessian.json";
  const json data = json::parse(std::ifstream(filename));
  const vector<double> q = data["q_5"].get<vector<double>>();
  const vector<double> best_q = data["best_q_5"].get<vector<double>>();
  const vector<double> g = data["g_5"].get<vector<double>>();
  const vector<double> best_g = data["best_g_5"].get<vector<double>>();
  const vector<double> H = data["H_5"].get<vector<double>>();
  const int ndata = q.size();
  vector<double> H1 = geoopt_routines_with_cuda::get_colwise_matrix(H, ndata, ndata);
  const vector<double> H_updated_ref = data["H_updated_5"].get<vector<double>>();
  geoopt_routines_with_cuda::Geoopt_cuda opt_cuda = geoopt_routines_with_cuda::Geoopt_cuda();
  const vector<double> H_updated = opt_cuda.update_Hessian_BFGS(q, best_q, g, best_g, H1);
  EXPECT_VECTORS_NEAR(H_updated, H_updated_ref, tolerance);
}

}
