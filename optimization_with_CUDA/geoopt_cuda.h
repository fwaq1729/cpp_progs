#ifndef __SRC_GEOOPT_CUDA_H
#define __SRC_GEOOPT_CUDA_H
#include <string>
#include <vector>
#include <tuple>

namespace geoopt_routines_with_cuda
{
  void fill_mat(std::vector<double>& A, const int nrow, const int ncol);

  void prn_mat(const std::string tag, const std::vector<double>& A, const int nrow, const int ncol);

  std::vector<double> get_colwise_matrix(const std::vector<double>& in, const int nrow, const int ncol);

  class Geoopt_cuda {
  private:
  public:
  Geoopt_cuda() {}

  void matmul_cublas(
  const int m,
  const int k,
  const int n,
  const double* h_A, const double* h_B, double* h_C);

  std::vector<double> linear_solver(const std::vector<double>& h_A, const std::vector<double>& h_B);

  void daxpy_cublas(
  std::vector<double>& h_y,
  const std::vector<double>& h_x,
  const double alpha);

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

  std::vector<double> update_Hessian_BFGS(
  const std::vector<double>& q,
  const std::vector<double>& best_q,
  const std::vector<double>& g,
  const std::vector<double>& best_g,
  const std::vector<double>& H);

  };
}

#endif
