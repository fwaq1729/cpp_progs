#include <cmath>
#include <omp.h>
#include <chrono>

#include "integration_cpu.h"

namespace
{
  const static double my_pi = std::acos(-1.0);
  inline float sinsum(float x, int terms) { // sin(x) = x - x^3/3! + x^5/5! ...
  float x2 = x * x;
  float term = x;
  float sum = term;
  for(int n = 1; n < terms; ++n) {
    term *= -x2 / (2 * n * (2 * n + 1));
    sum += term;
  }
  return sum;
}

void trapezoidal_integration_correction(double &sum, int terms, double step_size) {
  sum -= 0.5f * (sinsum(0.0, terms) + sinsum(my_pi, terms));
  sum *= step_size;
}
} // namespace

std::tuple<double, double> integration_cpu_routines::cpu_integration(int steps, int terms) {
  double step_size = my_pi / (steps - 1);
  auto start_time = std::chrono::high_resolution_clock::now();
  double cpu_sum = 0.0;
  for (int step = 0; step < steps; ++step) {
    float x = (float)(step_size * step);
    cpu_sum += sinsum(x, terms);
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);
  trapezoidal_integration_correction(cpu_sum, terms, step_size);
  return std::make_tuple(cpu_sum, (double)duration.count());
}

std::tuple<double, double> integration_cpu_routines::cpu_omp_integration(int steps, int terms, int threads) {
  double step_size = my_pi / (steps - 1);
  auto start_time = std::chrono::high_resolution_clock::now();
  double omp_sum = 0.0;
  omp_set_num_threads(threads);
  #pragma omp parallel for reduction (+:omp_sum)
  for(int step = 0; step < steps; ++step) {
    float x = (float)(step_size * step);
    omp_sum += sinsum(x, terms);
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);
  trapezoidal_integration_correction(omp_sum, terms, step_size);
  return std::make_tuple(omp_sum, (double)duration.count());
}
