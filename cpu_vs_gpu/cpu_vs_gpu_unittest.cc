#include "integration_cpu.h"
#include "integration_gpu.h"
#include "gtest/gtest.h"

using namespace std;

const int steps = 1000000;
const int terms = 1000;

namespace {
TEST(test_cpu, test1)
{
  double sum;
  double duration;
  tie(sum, duration) = integration_cpu_routines::cpu_integration(steps, terms);
  const double tolerance = 20;
  EXPECT_NEAR(sum, 2.0, 1e-8);
  EXPECT_NEAR(duration, 870.0, tolerance);
}

TEST(test_cpu_omp, test2)
{
  const int threads = 4;
  double sum;
  double duration;
  tie(sum, duration) = integration_cpu_routines::cpu_omp_integration(steps, terms, threads);
  const double tolerance = 50;
  EXPECT_NEAR(sum, 2.0, 1e-8);
  EXPECT_NEAR(duration, 220.0, tolerance);
}

TEST(test_cpu_omp1, test3)
{
  const int threads = 8;
  double sum;
  double duration;
  tie(sum, duration) = integration_cpu_routines::cpu_omp_integration(steps, terms, threads);
  const double tolerance = 50;
  EXPECT_NEAR(sum, 2.0, 1e-8);
  EXPECT_NEAR(duration, 110.0, tolerance);
}

TEST(test_gpu, test4)
{
  double sum;
  double duration;
  tie(sum, duration) = integration_gpu_routines::gpu_integration(steps, terms);
  const double tolerance = 5;
  EXPECT_NEAR(sum, 2.0, 3e-7);
  EXPECT_NEAR(duration, 1.4, tolerance);
}
}
