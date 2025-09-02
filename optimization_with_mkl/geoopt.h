#ifndef __SRC_GEOOPT_H
#define __SRC_GEOOPT_H
#include <vector>
#include <tuple>

namespace geoopt_routines_with_mkl
{
  std::vector<double> get_proj(
  const int ndata,
  const int nrow_B,
  const int ncol_B,
  const int ncol_Binv,
  const std::vector<double>& B,
  const std::vector<double>& B_inv);

  std::vector<double> get_g_new(
  const std::vector<double>& proj,
  const std::vector<double>& g_interp);

  std::vector<double> get_Hproj(
  const int ndata,
  const std::vector<double>& proj,
  const std::vector<double>& H);

  std::tuple<double, std::vector<double>, bool> do_quadratic_step(
  const double trust,
  const std::vector<double>& g_new,
  const std::vector<double>& Hproj);

  std::vector<double> update_Hessian_BFGS(
  const std::vector<double>& q,
  const std::vector<double>& best_q,
  const std::vector<double>& g,
  const std::vector<double>& best_g,
  const std::vector<double>& H);
}

#endif
