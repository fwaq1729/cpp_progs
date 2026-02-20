#include <cmath>
#include <numbers>
#include "poisson.h"

#include<iostream>
#include<iomanip>
using namespace std;

const constexpr static double Pi = std::numbers::pi;

template<typename T>
void poisson::Poisson<T>::prn_mat(const vector<T>& v, const int nx, const int ny, const string tag) {
  cout << tag << ":\n";
  for (int i = 0; i != nx; ++i) {
    for (int j = 0; j != ny; ++j) {
      if (j != ny - 1) {
        cout << setw(8) << setprecision(3) << fixed << v[i * ny + j] << ",";
      } else if (i != nx - 1) {
        cout << setw(8) << setprecision(3) << fixed << v[i * ny + j] << ",\n";
      } else {
        cout << setw(8) << setprecision(3) << fixed << v[i * ny + j] << "\n";
      }
    }
  }
}

template<typename T>
vector<T> poisson::Poisson<T>::extract_data_for_test(const vector<T>& v) {
  const T data_ini = xdata_[0];
  const T data_step = xdata_[1];
  const T data_end = xdata_[2];
  vector<T> v_x;
  for (T x = data_ini; x <= data_end; x+= data_step) {
    const int i = (int)(x / h_);
    for (T y = data_ini; y <= data_end; y+= data_step) {
      const int j = (int)(y / h_);
      v_x.push_back(v[i * ny1_ + j]);
    }
  }
  return v_x;
}

template<typename T>
void poisson::Poisson<T>::set_potentials() {
  const T v1 = v1234_[0];
  const T v2 = v1234_[1];
  const T v3 = v1234_[2];
  const T v4 = v1234_[3];
  for (int i = 1; i <= nx_ - 1; ++i) {
    for (int j = 1; j <= ny_ - 1; ++j) {
      v_num_[i * ny1_ + j] = (v1 + v2 + v3 + v4) / 4;
    }
  }
  // Set potentials at fixed nodes
  for (int i = 1; i <= nx_ - 1; ++i) {
    v_num_[i * ny1_ + 0] = v1;
    v_num_[i * ny1_ + ny_] = v3; 
  }
  for (int j = 1; j <= ny_ - 1; ++j) {
      v_num_[0 * ny1_ + j] = v4;
      v_num_[nx_ * ny1_ + j] = v2;
  }
  v_num_[0 * ny1_ + 0] = (v1 + v4) / 2;
  v_num_[nx_ * ny1_ + 0] = (v1 + v2) / 2;
  v_num_[0 * ny1_ + ny_] = (v3 + v4) / 2;
  v_num_[nx_ * ny1_ + ny_] = (v2 + v3) / 2;
}

template<typename T>
void poisson::Poisson<T>::compute_relaxation_factor() {
  const T t = cos(Pi / (T)nx_) + cos(Pi / (T)ny_);
  const T t2 = t * t;
  w_ = 4 * (2 - sqrt(4 - t2)) / t2;
  w4_ = w_ / 4;
}


template<typename T>
void poisson::Poisson<T>::poisson_solver() {
  int ncount = 0;
  int loop = 1;
  while (loop == 1) {
    T rmin = 0.0;
    for (int i = 1; i <= nx_ - 1; ++i) {
      const T x = h_ * i;
      for (int j = 1; j <= ny_ - 1; ++j) {
        const T y = h_ * j;
        const T g = -36 * Pi * x * (y - 1);
        const T r = w4_ * (v_num_[(i + 1) * ny1_ + j] + v_num_[(i - 1) * ny1_ + j] + v_num_[i * ny1_ + j + 1] + v_num_[i * ny1_ + j - 1] - 4 * v_num_[i * ny1_ + j] - g * h2_);
        rmin += fabs(r);
        v_num_[i * ny1_ + j] += r;
      }
    }
    rmin /= (T)(nx_ * ny_);
    if (rmin >= tol_) {
      ++ncount;
      if (ncount > niter_max_) {
        loop = 0;
        std::ostringstream os;
        os << "solution does not converge in " << niter_max_ << " iterations";
        throw std::runtime_error(os.str());
      }
    } else {
      loop = 0;
      niter_conv_ = ncount;
      //cout << "Solution converges in " << ncount << " iterations " << endl;
      //cout << "h = " << setw(15) << setprecision(8) << fixed << h_ << endl;
    }
  } // end-while
}

template<typename T>
void poisson::Poisson<T>::compute_exact_solution_case1() {
  const int n_series = 10; // number of terms in series expansion method
  for (int i = 1; i <= nx_ - 1; ++i) {
    const T x = h_ * (T)i;
    for (int j = 1; j <= ny_ - 1; ++j) {
      const T y = h_ * (T)j;
      double sum = 0;
      for (int m = 1; m <= n_series; ++m) {
        const T fm = (T)m;
        for (int n = 1; n <= n_series; ++n) {
          const T fn = (T)n;
          const T fma = fm * Pi / a_;
          const T fnb = fn * Pi / b_;
          const T factor1 = fma * fma + fnb * fnb;
          const T factor2 = pow(-1.0, (T)(m + n)) * 144.0 * a_ * b_ / (Pi * fm * fn);
          const T factor3 = 1.0 - (1.0 - pow(-1.0, (T)n)) / b_;
          const T factor = factor2 * factor3 / factor1;
          sum += factor * sin (fma * x) * sin(fnb * y);
        }
      }
      const T vh = sum;
      const T v1 = v1234_[0];
      const T v2 = v1234_[1];
      const T v3 = v1234_[2];
      const T v4 = v1234_[3];
      const T c1 = 4 * v1 / Pi;
      const T c2 = 4 * v2 / Pi;
      const T c3 = 4 * v3 / Pi;
      const T c4 = 4 * v4 / Pi;
      sum = 0.0;
      for (int k = 1; k <= n_series; ++k) {
        const T n = (T)(2 * k - 1);
        const T a1 = sin(n * Pi * x / b_);
        const T a2 = sinh(n * Pi * (a_ - y) / b_);
        const T a3 = n * sinh(n * Pi * a_ / b_);
        const T term1 = c1 * a1 * a2 / a3;
        const T b1 = sinh(n * Pi * x / a_);
        const T b2 = sin(n * Pi * y / a_);
        const T b3 = n * sinh(n * Pi * b_ / a_);
        const T term2 = c2 * b1 * b2 / b3;
        const T d1 = sin(n * Pi * x / b_);
        const T d2 = sinh(n * Pi * y / b_);
        const T d3 = n * sinh(n * Pi * a_ / b_);
        const T term3 = c3 * d1 * d2 / d3;
        const T e1 = sinh(n * Pi * (b_ - x) / a_);
        const T e2 = sin(n * Pi * y / a_);
        const T e3 = n * sinh(n * Pi * b_ / a_);
        const T term4 = c4 * e1 * e2 / e3;
        const T term = term1 + term2 + term3 + term4;
        sum += term;
      } // end-loop-k
      const T vi = sum;
      v_exact_[i * ny1_ + j] = vh + vi;
    }
  }
}

template<typename T>
poisson::Poisson<T>::Poisson(const nlohmann::json input_pars0) :
  input_pars_(input_pars0), nx_(input_pars0["nx"]), ny_(input_pars0["ny"]), ny1_(ny_ + 1),
  a_(input_pars0["a"]), b_(input_pars0["b"]), v1234_(input_pars0["v1234"]),
  tol_(input_pars0["tolerance"]), niter_max_(input_pars0["niter_max"]),
  h_(a_ / nx_), h2_(h_ * h_),
  v_num_((nx_ + 1) * (ny_ + 1), 0.0),
  v_exact_((nx_ + 1) * (ny_ + 1), 0.0),
  xdata_(input_pars0["xdata"]),
  niter_conv_(0)
  {
    set_potentials();
    compute_relaxation_factor();
    poisson_solver();
    //prn_mat(v_num_, nx_ + 1, ny_ + 1, "v_num");
    //prn_mat(v_exact_, nx_ + 1, ny_ + 1, "v_exact");
  }

template class poisson::Poisson<double>;  // INSTANTIATION
