#include <cmath>
#include <algorithm>
#include <mkl.h>
#include "geoopt.h"

// Documentation:
// Original routines located in: https://github.com/jhrmnn/pyberny
// written in python using numpy
// Adapted to C++ using MKL routines
// List of MKL functions used:
// cblas_dgemm, matrix products
// cblas_ddot, dot product
// cblas_dscal, scaling a vector
// cblas_daxpy, ax_p_py routine
// LAPACKE_dsyev, solving a linear system: A x = b

using namespace std;

// Computing: proj = dot(B, B_inv)
vector<double> geoopt_routines_with_mkl::get_proj(
  const int ndata,
  const int nrow_B,
  const int ncol_B,
  const int ncol_Binv,
  const vector<double>& B,     // Wilson matrix
  const vector<double>& B_inv) // inverse of Wilson matrix
{
  vector<double> proj(ndata * ndata);
  // A(m,k) B(k,n) -> C(m, n)
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
              nrow_B,       // m      
              ncol_Binv,    // n
              ncol_B,       // k
              1.0,          // alpha
              B.data(),     // A
              nrow_B,       // m
              B_inv.data(), // B
              ncol_B,       // k
              0.0,          // beta  
              proj.data(),  // C
              nrow_B);      // m

  return proj;
}

// Computing: g_new = dot(proj, s.interpolated.g)
vector<double> geoopt_routines_with_mkl::get_g_new(
  const vector<double>& proj,
  const vector<double>& g_interp)
{
  const int ndata = g_interp.size();
  vector<double> g_new(ndata);            
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
              ndata,           // m      
              1,               // n
              ndata,           // k
              1.0,             // alpha
              proj.data(),     // A
              ndata,           // m
              g_interp.data(), // B
              ndata,           // k
              0.0,             // beta  
              g_new.data(),    // C
              ndata);          // m      

  return g_new;
}

// Computing: H_proj = proj.dot(s.H).dot(proj) + 1000 * (eye(len(s.coords)) - proj)
vector<double> geoopt_routines_with_mkl::get_Hproj(
  const int ndata,
  const vector<double>& proj,
  const vector<double>& H)
{
  vector<double> proj1(ndata * ndata);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              ndata,        // rows in A              
              ndata,        // cols in B
              ndata,        // cols in A
              1.0,          // alpha
              proj.data(),  // A
              ndata,        // cols in A
              H.data(),     // B
              ndata,        // cols in B
              0.0,          // beta  
              proj1.data(), // C
              ndata);       // cols in C
  vector<double> proj2(ndata * ndata);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              ndata,        // rows in A              
              ndata,        // cols in B
              ndata,        // cols in A
              1.0,          // alpha
              proj1.data(), // A
              ndata,        // cols in A
              proj.data(),  // B
              ndata,        // cols in B
              0.0,          // beta  
              proj2.data(), // C
              ndata);       // cols in C
  
  vector<double> v(ndata * ndata);
  copy_n(proj.data(), ndata * ndata, v.data());
  cblas_dscal(ndata * ndata, -1.0, v.data(), 1);
  for (int i = 0; i != ndata; ++i)
  {
    v[i + ndata * i] = 1.0 + v[i + ndata * i];
  }

  vector<double> Hproj(ndata * ndata);
  copy_n(proj2.data(), ndata * ndata, Hproj.data());
  cblas_daxpy(ndata * ndata,  1000.0, v.data(), 1, Hproj.data(), 1);

  return Hproj;
}

// Computing:
// dq, dE, on_sphere = quadratic_step(dot(proj, s.interpolated.g), H_proj, s.weights, s.trust, log=log)
tuple<double, vector<double>, bool> geoopt_routines_with_mkl::do_quadratic_step(
  const double trust,
  const vector<double>& g_new,
  const vector<double>& Hproj)
{
  const int ndata = g_new.size();
  // Get H + H^t:
  vector<double> Hproj_hermitian(ndata * ndata);
  for (int j = 0; j != ndata; ++j)
  {
    Hproj_hermitian[j + ndata * j] = Hproj[j + ndata * j];
    for (int i = j + 1; i != ndata; ++i)
    {
      Hproj_hermitian[i + ndata * j] = 0.5 * (Hproj[i + ndata * j] + Hproj[j + ndata * i]);
      Hproj_hermitian[j + ndata * i] = Hproj_hermitian[i + ndata * j];
    }
  }
  vector<double> eigenvals(ndata);
  vector<double> eigenvects(ndata * ndata);
  // Warning: Hproj_hermitian, will be lost: now it contains eigenvectors
  // dsyev, diagonalize a real-symmetric matrix
  int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', ndata, Hproj_hermitian.data(), ndata,
                           eigenvals.data());
  const int ndata1 = ndata + 1;
  vector<double> rfo_symm(ndata1 * ndata1);
  for (int j = 0; j != ndata; ++j)
  {
    rfo_symm[ndata1 - 1 + ndata1 * j] = g_new[j];
    rfo_symm[j + ndata1 * (ndata1 - 1)] = g_new[j];
    rfo_symm[j + ndata1 * j] = Hproj[j + ndata * j];
    for (int i = j + 1; i != ndata; ++i)
    {
      rfo_symm[i + ndata1 * j] = 0.5 * (Hproj[i + ndata * j] + Hproj[j + ndata * i]);
      rfo_symm[j + ndata1 * i] = rfo_symm[i + ndata1 * j];
    }
  }
  vector<double> eigenvals1(ndata1);
  int info1 = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', ndata1, rfo_symm.data(), ndata1,
                            eigenvals1.data());
  vector<double> dq(ndata);
  for (int j = 0; j != ndata; j++)
  {
    dq[j] = rfo_symm[0 + ndata1 * j] / rfo_symm[0 + ndata1 * ndata];
  }
  vector<double> Hdq(ndata * 1);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              ndata,        // rows in A              
              1,            // cols in B
              ndata,        // cols in A
              1.0,          // alpha
              Hproj.data(), // A
              ndata,        // cols in A
              dq.data(),    // B
              1,            // cols in B
              0.0,          // beta  
              Hdq.data(),   // C
              1);           // cols in C

  const double dE = cblas_ddot(ndata, dq.data(), 1, g_new.data(), 1) +
                    0.5 * cblas_ddot(ndata, Hdq.data(), 1, dq.data(), 1);
  const double norm_dq = sqrt(cblas_ddot(ndata, dq.data(), 1, dq.data(), 1));
  const bool on_sphere = (norm_dq <= trust) ? false : true;
  return make_tuple(dE, dq, on_sphere);
}

vector<double> geoopt_routines_with_mkl::update_Hessian_BFGS(
  const vector<double>& q,
  const vector<double>& best_q,
  const vector<double>& g,
  const vector<double>& best_g,
  const vector<double>& H)
{
  const int ndata = q.size();
  vector<double> q_scratch(ndata);
  copy_n(q.data(), ndata, q_scratch.data());
  vector<double> g_scratch(ndata);
  copy_n(g.data(), ndata, g_scratch.data());
  cblas_daxpy(ndata, -1.0, best_q.data(), 1, q_scratch.data(), 1);
  cblas_daxpy(ndata, -1.0, best_g.data(), 1, g_scratch.data(), 1);
  // outer(dg,dg):
  vector<double> dH1(ndata * ndata);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              ndata,      // rows in A              
              ndata,      // cols in B
              1,          // cols in A
              1.0,        // alpha
              g_scratch.data(),   // A
              1,          // cols in A
              g_scratch.data(),   // B
              ndata,      // cols in B
              0.0,        // beta  
              dH1.data(), // C
              ndata);     // cols in C
  const double dqdg = 1.0 / cblas_ddot(ndata, g_scratch.data(), 1, q_scratch.data(), 1);
  cblas_dscal(ndata * ndata, dqdg, dH1.data(), 1);
  vector<double> dH2_0(ndata * ndata);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              ndata,      // rows in A              
              ndata,      // cols in B
              1,          // cols in A
              1.0,        // alpha
              q_scratch.data(),   // A
              1,          // cols in A
              q_scratch.data(),   // B
              ndata,      // cols in B
              0.0,        // beta  
              dH2_0.data(), // C
              ndata);     // cols in C
  vector<double> HdH2_0(ndata * ndata);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              ndata,        // rows in A              
              ndata,        // cols in B
              ndata,        // cols in A
              1.0,          // alpha
              H.data(),     // A
              ndata,        // cols in A
              dH2_0.data(), // B
              ndata,        // cols in B
              0.0,          // beta  
              HdH2_0.data(),// C
              ndata);       // cols in C
  vector<double> HdH2_0H(ndata * ndata);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              ndata,        // rows in A              
              ndata,        // cols in B
              ndata,        // cols in A
              1.0,          // alpha
              HdH2_0.data(),// A
              ndata,        // cols in A
              H.data(),     // B
              ndata,        // cols in B
              0.0,          // beta  
              HdH2_0H.data(),// C
              ndata);       // cols in C
  vector<double> dqH(ndata);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              1,            // rows in A              
              ndata,        // cols in B
              ndata,        // cols in A
              1.0,          // alpha
              q_scratch.data(),     // A
              ndata,        // cols in A
              H.data(),     // B
              ndata,        // cols in B
              0.0,          // beta  
              dqH.data(),   // C
              ndata);       // cols in C
  const double dqHdq = 1.0 / cblas_ddot(ndata, dqH.data(), 1, q_scratch.data(), 1);
  cblas_dscal(ndata * ndata, dqHdq, HdH2_0H.data(), 1);
  cblas_daxpy(ndata * ndata, -1.0, HdH2_0H.data(), 1, dH1.data(), 1);
  vector<double> H_updated(ndata * ndata);
  copy_n(H.data(), ndata * ndata, H_updated.data());
  cblas_daxpy(ndata * ndata,  1.0, dH1.data(), 1, H_updated.data(), 1);
  return H_updated;
}
