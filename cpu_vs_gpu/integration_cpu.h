#ifndef __INTEGRATION_CPU_H
#define __INTEGRATION_CPU_H
#include <tuple>
namespace integration_cpu_routines
{
std::tuple<double, double> cpu_integration(int steps, int terms);
std::tuple<double, double> cpu_omp_integration(int steps, int terms, int threads);
} // end namespace
#endif

