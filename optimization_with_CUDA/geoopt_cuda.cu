#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cublas_v2.h>

#include "geoopt_cuda.h"

using namespace std;

// Error checking macro for CUDA calls
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error at " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(err) << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

// Error checking macro for cuSOLVER calls
#define CUSOLVER_CHECK(call) \
    do { \
        cusolverStatus_t status = call; \
        if (status != CUSOLVER_STATUS_SUCCESS) { \
            std::cerr << "cuSOLVER Error at " << __FILE__ << ":" << __LINE__ << ": " << status << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

void geoopt_routines_with_cuda::fill_mat(vector<double>& A, const int nrow, const int ncol)
{
  for (int i = 0, count = 0; i != nrow; ++i)
  {
    for (int j = 0; j != ncol; ++count, ++j)
    {
       A[count] = rand() / static_cast<double>(RAND_MAX);
    }
  }
}

vector<double> geoopt_routines_with_cuda::get_colwise_matrix(const vector<double>& in, const int nrow, const int ncol)
{
  vector<double> ou(nrow * ncol);
  for (int i = 0, count = 0; i != nrow; ++i)
  {
    for (int j = 0; j != ncol; ++count, ++j)
    {
      ou[count] = in[i + nrow * j];
    }
  }

  return ou;
}

void geoopt_routines_with_cuda::prn_mat(const string tag, const vector<double>& A, const int nrow, const int ncol)
{
  cout << tag << " = np.array([\n";
  for (int i = 0; i != nrow; ++i)
  {
    cout << "  [";
    for (int j = 0; j != ncol; ++j)
    {
       if (j != ncol - 1)
       {
         cout << setw(15) << setprecision(8) << scientific << A[j + ncol * i] << ",";
       } else if (i != nrow - 1)
       {
         cout << setw(15) << setprecision(8) << scientific << A[j + ncol * i] << "],\n";
       } else
       {
         cout << setw(15) << setprecision(8) << scientific << A[j + ncol * i] << "]]);\n";
       }
    }
  }
}

void geoopt_routines_with_cuda::Geoopt_cuda::matmul_cublas(
  const int m,
  const int k,
  const int n,
  const double* h_A, const double* h_B, double* h_C)
{
  cublasHandle_t handle;
  cublasCreate(&handle);
  double *d_A, *d_B, *d_C;
  cudaMalloc(&d_A, m * k * sizeof(double));
  cudaMalloc(&d_B, k * n * sizeof(double));
  cudaMalloc(&d_C, m * n * sizeof(double));
  cudaMemcpy(d_A, h_A, m * k * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, k * n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_C, h_C, m * n * sizeof(double), cudaMemcpyHostToDevice);
  // Perform matrix multiplication C = alpha * A * B + beta * C
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
              n, m, k, &alpha, d_B, n, d_A, k, &beta, d_C, n);
  cudaMemcpy(h_C, d_C, m * n * sizeof(double), cudaMemcpyDeviceToHost);
  cublasDestroy(handle);
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
}

vector<double> geoopt_routines_with_cuda::Geoopt_cuda::linear_solver(const vector<double>& h_A, const vector<double>& h_B)
{
  const int n = h_B.size();
  vector<double> h_X(n);
  // Device-side pointers
  double *d_A, *d_B;
  int *d_P;
  double *d_Workspace; // Workspace for cuSOLVER routines
  int *d_Info;    // For error codes from cuSOLVER

  // Allocate device memory
  CUDA_CHECK(cudaMalloc(&d_A, n * n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_B, n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_P, n * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&d_Info, sizeof(int)));

  // Copy data from host to device
  CUDA_CHECK(cudaMemcpy(d_A, h_A.data(), n * n * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_B, h_B.data(), n * sizeof(double), cudaMemcpyHostToDevice));

  // Initialize cuSOLVER handle
  cusolverDnHandle_t cusolverH = NULL;
  CUSOLVER_CHECK(cusolverDnCreate(&cusolverH));

  // Query workspace size for LU factorization
  int lwork = 0;
  CUSOLVER_CHECK(cusolverDnDgetrf_bufferSize(cusolverH, n, n, d_A, n, &lwork));
  CUDA_CHECK(cudaMalloc(&d_Workspace, lwork * sizeof(double)));

  // Perform LU factorization
  CUSOLVER_CHECK(cusolverDnDgetrf(cusolverH, n, n, d_A, n, d_Workspace, d_P, d_Info));

  // Check for singularity (if *d_Info != 0, matrix is singular)
  int info;
  CUDA_CHECK(cudaMemcpy(&info, d_Info, sizeof(int), cudaMemcpyDeviceToHost));
  if (info != 0) {
      std::cerr << "Matrix A is singular. Cannot solve the system." << std::endl;
      // Handle singularity appropriately (e.g., return, use a different solver)
  } else {
      // Solve the linear system A*X = B using the LU factors
      CUSOLVER_CHECK(cusolverDnDgetrs(cusolverH, CUBLAS_OP_N, n, 1, d_A, n, d_P, d_B, n, d_Info));
      // Copy solution from device to host
      CUDA_CHECK(cudaMemcpy(h_X.data(), d_B, n * sizeof(double), cudaMemcpyDeviceToHost));
    }
  // Clean up
  CUDA_CHECK(cudaFree(d_A));
  CUDA_CHECK(cudaFree(d_B));
  CUDA_CHECK(cudaFree(d_P));
  CUDA_CHECK(cudaFree(d_Workspace));
  CUDA_CHECK(cudaFree(d_Info));
  CUSOLVER_CHECK(cusolverDnDestroy(cusolverH));

  return h_X;
}

tuple<vector<double>, vector<double>> geoopt_routines_with_cuda::Geoopt_cuda::eigenvalue_solver(
  const int ndata,
  const vector<double>& h_A)
{
    const int ndata2 = ndata * ndata;
    // Pointers for device memory
    double *d_A, *d_W, *d_work;
    int *d_info;

    // cuSOLVER variables
    cusolverDnHandle_t cusolverH = NULL;
    int lwork = 0;

    // 1. Initialize cuSOLVER handle
    CUSOLVER_CHECK(cusolverDnCreate(&cusolverH));

    // 2. Allocate device memory
    CUDA_CHECK(cudaMalloc((void**)&d_A, ndata2 * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_W, ndata * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_info, sizeof(int)));

    // 3. Copy host matrix to device
    CUDA_CHECK(cudaMemcpy(d_A, h_A.data(), ndata2 * sizeof(double), cudaMemcpyHostToDevice));

    // 4.1. Query for the optimal work buffer size
    CUSOLVER_CHECK(cusolverDnDsyevd_bufferSize(
        cusolverH,
        CUSOLVER_EIG_MODE_VECTOR, // CUSOLVER_EIG_MODE_VECTOR for both eigenvalues and eigenvectors
        CUBLAS_FILL_MODE_LOWER,   // Lower triangular part of the matrix
        ndata,                    // Matrix dimension
        d_A,                      // Device pointer to the matrix
        ndata,                    // Leading dimension of d_A
        d_W,                      // Device pointer for eigenvalues
        &lwork                    // Pointer to workspace size
    ));

    // 4.2. Allocate device work buffer
    CUDA_CHECK(cudaMalloc((void**)&d_work, lwork * sizeof(double)));

    // 4.3. Compute eigenvalues and eigenvectors
    // The eigenvectors are stored in the input matrix d_A
    CUSOLVER_CHECK(cusolverDnDsyevd(
        cusolverH,
        CUSOLVER_EIG_MODE_VECTOR,
        CUBLAS_FILL_MODE_LOWER,
        ndata,
        d_A,
        ndata,
        d_W,
        d_work,
        lwork,
        d_info
    ));

    // 5. Copy results back to host
    std::vector<double> h_eigenvalues(ndata);
    std::vector<double> h_eigenvectors(ndata2);

    CUDA_CHECK(cudaMemcpy(h_eigenvalues.data(), d_W, ndata * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_eigenvectors.data(), d_A, ndata2 * sizeof(double), cudaMemcpyDeviceToHost));

    int h_info;
    CUDA_CHECK(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));

    if (h_info != 0) {
        const string tag = "cuSOLVER failed to converge. Info = " + to_string(h_info);
        throw runtime_error(tag);
    }

    // 6. Clean up
    CUSOLVER_CHECK(cusolverDnDestroy(cusolverH));
    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_work));
    CUDA_CHECK(cudaFree(d_info));

    return make_tuple(h_eigenvalues, h_eigenvectors);
}

void geoopt_routines_with_cuda::Geoopt_cuda::daxpy_cublas(
  vector<double>& h_y,
  const vector<double>& h_x,
  const double alpha)
{
    // Device-side pointers
    double *d_x, *d_y;

    const int n = h_y.size();
    // Allocate memory on the device
    cudaMalloc(&d_x, n * sizeof(double));
    cudaMalloc(&d_y, n * sizeof(double));

    // Copy host data to device
    cudaMemcpy(d_x, h_x.data(), n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, h_y.data(), n * sizeof(double), cudaMemcpyHostToDevice);

    // Initialize cuBLAS handle
    cublasHandle_t handle;
    cublasCreate(&handle);

    // Perform the AXPY operation: d_y = alpha * d_x + d_y
    // Parameters: handle, n, alpha, d_x, incx, d_y, incy
    // incx and incy are strides (1 for contiguous elements)
    cublasDaxpy(handle, n, &alpha, d_x, 1, d_y, 1);

    // Copy results back to host
    cudaMemcpy(h_y.data(), d_y, n * sizeof(double), cudaMemcpyDeviceToHost);

    // Clean up
    cudaFree(d_x);
    cudaFree(d_y);
    cublasDestroy(handle);
}

// Computing: proj = dot(B, B_inv)
vector<double> geoopt_routines_with_cuda::Geoopt_cuda::get_proj(
  const int ndata,
  const int nrow_B,
  const int ncol_B,
  const int ncol_Binv,
  const vector<double>& B,     // Wilson matrix
  const vector<double>& B_inv) // inverse of Wilson matrix
{
  vector<double> proj(ndata * ndata);
  matmul_cublas(nrow_B, ncol_B, ncol_Binv, B.data(), B_inv.data(), proj.data());

  return proj;
}

// Computing: g_new = dot(proj, s.interpolated.g)
vector<double> geoopt_routines_with_cuda::Geoopt_cuda::get_g_new(
  const vector<double>& proj,
  const vector<double>& g_interp)
{
  const int ndata = g_interp.size();
  vector<double> g_new(ndata);
  matmul_cublas(ndata, ndata, 1, proj.data(), g_interp.data(), g_new.data());

  return g_new;
}

// Computing: H_proj = proj.dot(s.H).dot(proj) + 1000 * (eye(len(s.coords)) - proj)
vector<double> geoopt_routines_with_cuda::Geoopt_cuda::get_Hproj(
  const int ndata,
  const vector<double>& proj,
  const vector<double>& H)
{
  cublasHandle_t handle;
  cublasCreate(&handle);
  double *d_A, *d_B, *d_C;
  const int ndata2 = ndata * ndata;
  cudaMalloc(&d_A, ndata2 * sizeof(double));
  cudaMalloc(&d_B, ndata2 * sizeof(double));
  cudaMalloc(&d_C, ndata2 * sizeof(double));
  cudaMemcpy(d_A, proj.data(), ndata2 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, H.data(), ndata2 * sizeof(double), cudaMemcpyHostToDevice);
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
              ndata, ndata, ndata, &alpha, d_B, ndata, d_A, ndata, &beta, d_C, ndata);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
              ndata, ndata, ndata, &alpha, d_A, ndata, d_C, ndata, &beta, d_C, ndata);
  double alpha1 = -1000.0;
  cublasDaxpy(handle, ndata2, &alpha1, d_A, 1, d_C, 1);
  vector<double> Hproj(ndata2);
  cudaMemcpy(Hproj.data(), d_C, ndata2 * sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
  cublasDestroy(handle);
  for (int i = 0; i != ndata; ++i)
  {
    Hproj[i + ndata * i] -= alpha1;
  }

  return Hproj;
}

// Computing:
// dq, dE, on_sphere = quadratic_step(dot(proj, s.interpolated.g), H_proj, s.weights, s.trust, log=log)
tuple<double, vector<double>, bool> geoopt_routines_with_cuda::Geoopt_cuda::do_quadratic_step(
  const double trust,
  const vector<double>& g_new,
  const vector<double>& Hproj)
{
  const int ndata = g_new.size();
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
  vector<double> eigenvectors1(ndata1 * ndata1);
  tie(eigenvals1, eigenvectors1) = eigenvalue_solver(ndata1, rfo_symm);

  vector<double> dq(ndata);
  for (int j = 0; j != ndata; j++)
  {
    dq[j] = eigenvectors1[j + ndata1 * 0] / eigenvectors1[ndata + ndata1 * 0];
  }
  vector<double> Hdq(ndata);
  matmul_cublas(ndata, ndata, 1, Hproj.data(), dq.data(), Hdq.data());
  double dE1, dE2, norm_dq;
  double *g_new_d, *dq_d, *Hdq_d;
  cublasHandle_t handle;
  cublasCreate(&handle);
  cudaMalloc(&g_new_d, ndata * sizeof(double));
  cudaMalloc(&dq_d, ndata * sizeof(double));
  cudaMalloc(&Hdq_d, ndata * sizeof(double));
  cudaMemcpy(g_new_d, g_new.data(), ndata * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dq_d, dq.data(), ndata * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Hdq_d, Hdq.data(), ndata * sizeof(double), cudaMemcpyHostToDevice);
  cublasDdot(handle, ndata, g_new_d, 1, dq_d, 1, &dE1);
  cublasDdot(handle, ndata, Hdq_d, 1, dq_d, 1, &dE2);
  const double dE = dE1 + 0.5 * dE2;
  cublasDdot(handle, ndata, dq.data(), 1, dq.data(), 1, &norm_dq);
  cublasDestroy(handle);
  cudaFree(g_new_d);
  cudaFree(dq_d);
  cudaFree(Hdq_d);
  norm_dq = sqrt(norm_dq);
  const bool on_sphere = (norm_dq <= trust) ? false : true;

  return make_tuple(dE, dq, on_sphere);
}

vector<double> geoopt_routines_with_cuda::Geoopt_cuda::update_Hessian_BFGS(
  const vector<double>& q,
  const vector<double>& best_q,
  const vector<double>& g,
  const vector<double>& best_g,
  const vector<double>& H)
{
  const int ndata = q.size();
  const int ndata2 = ndata * ndata;
  double *q_h, *g_h, *q_h1, *g_h1, *dH1, *dH2_0, *H_h, *HdH2_0, *HdH2_0H, *dqH;
  cublasHandle_t handle;
  cublasCreate(&handle);
  cudaMalloc(&q_h, ndata * sizeof(double));
  cudaMalloc(&g_h, ndata * sizeof(double));
  cudaMalloc(&q_h1, ndata * sizeof(double));
  cudaMalloc(&g_h1, ndata * sizeof(double));
  cudaMalloc(&dH1, ndata2 * sizeof(double));
  cudaMalloc(&dH2_0, ndata2 * sizeof(double));
  cudaMalloc(&H_h, ndata2 * sizeof(double));
  cudaMalloc(&HdH2_0, ndata2 * sizeof(double));
  cudaMalloc(&HdH2_0H, ndata2 * sizeof(double));
  cudaMalloc(&dqH, ndata * sizeof(double));
  cudaMemcpy(q_h, q.data(), ndata * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(g_h, g.data(), ndata * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(q_h1, best_q.data(), ndata * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(g_h1, best_g.data(), ndata * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(H_h, H.data(), ndata2 * sizeof(double), cudaMemcpyHostToDevice);
  double alpha0 = -1.0;
  cublasDaxpy(handle, ndata, &alpha0, q_h1, 1, q_h, 1);
  cublasDaxpy(handle, ndata, &alpha0, g_h1, 1, g_h, 1);

  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
              ndata, ndata, 1, &alpha, g_h, ndata, g_h, 1, &beta, dH1, ndata);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
              ndata, ndata, 1, &alpha, q_h, ndata, q_h, 1, &beta, dH2_0, ndata);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
              ndata, ndata, ndata, &alpha, dH2_0, ndata, H_h, ndata, &beta, HdH2_0, ndata);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
              ndata, ndata, ndata, &alpha, H_h, ndata, HdH2_0, ndata, &beta, HdH2_0H, ndata);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
              1, ndata, ndata, &alpha, q_h, 1, H_h, ndata, &beta, dqH, 1);

  double qd_dot, dqHdq_dot;
  cublasDdot(handle, ndata, g_h, 1, q_h, 1, &qd_dot);
  cublasDdot(handle, ndata, dqH, 1, q_h, 1, &dqHdq_dot);
  const double dqdg = 1.0 / qd_dot;
  const double dqHdq = 1.0 / dqHdq_dot;
  cublasDscal(handle, ndata2, &dqdg, dH1, 1);
  cublasDscal(handle, ndata2, &dqHdq, HdH2_0H, 1);
  double alpha1 = -1.0;
  double alpha2 =  1.0;
  cublasDaxpy(handle, ndata2, &alpha1, HdH2_0H, 1, dH1, 1);
  cublasDaxpy(handle, ndata2, &alpha2, dH1, 1, H_h, 1);
  vector<double> H_updated(ndata2);
  cudaMemcpy(H_updated.data(), H_h, ndata2 * sizeof(double), cudaMemcpyDeviceToHost);
  cublasDestroy(handle);
  cudaFree(q_h);
  cudaFree(g_h);
  cudaFree(q_h1);
  cudaFree(g_h1);
  cudaFree(dH1);
  cudaFree(dH2_0);
  cudaFree(H_h);
  cudaFree(HdH2_0);
  cudaFree(HdH2_0H);
  cudaFree(dqH);

  return H_updated;
}
