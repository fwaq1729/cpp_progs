#ifndef POISSON_H_
#define POISSON_H_
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

namespace poisson {
template<typename T>
class Poisson {
  public:
    explicit Poisson(const nlohmann::json input_pars);
    void poisson_solver();
    void compute_exact_solution_case1();
    std::vector<T> extract_data_for_test(const std::vector<T>& v);
    T get_relaxation_factor() const { return w_; }
    T get_h() const { return h_; }
    int get_niter_conv() const { return niter_conv_; }
    std::vector<T> get_v_num() const { return v_num_; }
    std::vector<T> get_v_exact() const { return v_exact_;}
  private:
    nlohmann::json input_pars_;
    int nx_;
    int ny_;
    int ny1_;
    T a_;
    T b_;
    std::vector<T> v1234_;
    T tol_;
    int niter_max_;
    T h_;
    T h2_;
    std::vector<T> v_num_;
    std::vector<T> v_exact_;
    T w_;
    T w4_;
    std::vector<T> xdata_;
    int niter_conv_;
    
    void prn_mat(const std::vector<T>& v, const int nx, const int ny, const std::string tag);
    void set_potentials();
    void compute_relaxation_factor();

  protected:
};
}
#endif // POISSON_H_

