#include <cmath>
#include <chrono>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#include "integration_gpu.h"

namespace
{
  const static double my_pi = std::acos(-1.0);
__host__ __device__ inline float sinsum(float x, int terms) { // sin(x) = x - x^3/3! + x^5/5! ...
  float x2 = x * x;
  float term = x;
  float sum = term;
  for(int n = 1; n < terms; ++n) {
    term *= -x2 / (2 * n * (2 * n + 1));
    sum += term;
  }
  return sum;
}

__global__ void gpu_sin(float *sums, int steps, int terms, float step_size) {
  int step = blockIdx.x * blockDim.x + threadIdx.x;
  if (step < steps) {
    float x = step_size * step;
    sums[step] = sinsum(x, terms);
  }
}

void trapezoidal_integration_correction(double &sum, int terms, double step_size) {
  sum -= 0.5f * (sinsum(0.0, terms) + sinsum(my_pi, terms));
  sum *= step_size;
}
} // namespace

std::tuple<double, double> integration_gpu_routines::gpu_integration(int steps, int terms) {
  double step_size = my_pi / (steps - 1);

  int threads = 256;
  int blocks = (steps + threads - 1) / threads;

  thrust::device_vector<float> dsums(steps);
  float *dptr = thrust::raw_pointer_cast(&dsums[0]);

  auto start_time = std::chrono::high_resolution_clock::now();
  gpu_sin<<<blocks,threads>>>(dptr, steps, terms, (float)step_size);
  double gpu_sum = thrust::reduce(dsums.begin(),dsums.end());
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);
  trapezoidal_integration_correction(gpu_sum, terms, step_size);
  return std::make_tuple(gpu_sum, (double)duration.count());
}
