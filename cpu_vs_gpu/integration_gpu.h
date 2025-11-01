#ifndef __INTEGRATION_GPU_H
#define __INTEGRATION_GPU_H
#include <tuple>
namespace integration_gpu_routines
{
std::tuple<double, double> gpu_integration(int steps, int terms);
} // end namespace
#endif

