#include "nbody.h"
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
TEST(test_nbody_simulator, test1)
{
  const double tolerance = 1e-6;
  const nlohmann::json data = nlohmann::json::parse(std::ifstream("nbody_input.json"));
  const nlohmann::json data1 = nlohmann::json::parse(std::ifstream("nbody_data_for_test.json"));
  const int ncases = 7;
  for (int i = 0; i != ncases; ++i)
  {
    const string tag = "case_" + to_string(i + 1);
    const nlohmann::json input_pars = data.at(tag);
    nbody::Nbody simulator = nbody::Nbody(input_pars);
    simulator.nbody_simulator();
    for (double en_t : simulator.get_energy())
    {
      EXPECT_NEAR(en_t, data1.at(tag).at("energy"), tolerance);
    }
    const vector<double> ang_mom_ref = data1.at(tag).at("angular_momentum").get<vector<double>>();
    const vector<double> ang_mom_t = simulator.get_angular_momentum();
    for (int j = 0, count = 0; j != simulator.get_ndata(); ++j)
    {
      for (int k = 0; k != 3; ++count, ++k)
      {
        EXPECT_NEAR(ang_mom_ref[k], ang_mom_t[count], tolerance);
      }
    }
    const vector<double> ang_mom_i_ref = data1.at(tag).at("angular_momentum_i").get<vector<double>>();
    const vector<double> ang_mom_i_t = simulator.get_angular_momentum_i();
    for (int j = 0, count = 0; j != simulator.get_ndata(); ++j)
    {
      for (int k = 0; k != 3 * simulator.get_nbody(); ++count, ++k)
      {
        EXPECT_NEAR(ang_mom_i_ref[k], ang_mom_i_t[count], tolerance);
      }
    }

    for (double val : simulator.get_c_alpha())
    {
      EXPECT_NEAR(val, data1.at(tag).at("c_alpha"), tolerance);
    }
  }
}
}
