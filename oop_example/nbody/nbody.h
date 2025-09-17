#ifndef __SRC_NBODY_H
#define __SRC_NBODY_H
#include <string>
#include <array>
#include <vector>
#include <tuple>
#include <sstream>
#include <fstream>
#include <nlohmann/json.hpp>

namespace nbody{
class Nbody {
  public:
    Nbody(const nlohmann::json input_pars);
    void nbody_simulator();
    int get_nbody() { return nbody_; }
    int get_ndata() { return ndata_; }
    const std::vector<double>& get_c_alpha() { return c_alpha_; }
    const std::vector<double>& get_energy() { return energy_; }
    const std::vector<double>& get_angular_momentum() { return angular_momentum_; }
    const std::vector<double>& get_angular_momentum_i() { return angular_momentum_i_; }

  private:
    nlohmann::json input_pars_;
    int nbody_;
    int ndata_;
    std::vector<double> c_alpha_;
    std::vector<double> pos_;
    std::vector<double> vel_;
    std::vector<double> mass_;
    std::vector<double> energy_;
    std::vector<double> angular_momentum_;
    std::vector<double> angular_momentum_i_;

    void get_ini_posvel(
      const nlohmann::json& input_pars);
    std::tuple<double, std::vector<double>> linear_configuration(
      const int config,
      const std::vector<double>& r_ini,
      const std::vector<double>& mass);

    std::vector<double> gen_vertices(
      const double radius,
      const int nbody);

    double fderiv(
      const std::array<double, 6>& a,
      const double xi,
      const int m);

    double newton_raphson(const std::array<double, 6>& c);

    void positions_wrt_cm(
      std::vector<double>& pos,
      const std::vector<double>& mass);

    std::vector<double> get_alpha_for_3bodies(
      const std::vector<double>& mass);

    std::vector<double> get_alpha_for_linear_configuration(
      const int config,
      const double x,
      const std::vector<double>& mass);

    std::vector<double> get_alpha_for_4bodies(
      const std::vector<double>& pos,
      const std::vector<double>& mass);

    std::vector<double> get_alpha_for_regular_polygon(const std::vector<double>& mass);

    std::tuple<std::vector<double>, std::vector<double>> gen_vertices_for_4bodies(
      const double rbase,
      const double kappa,
      const double mass3,
      const int selec_nu);

    std::vector<double> get_forces(
      const std::vector<double>& mass,
      const std::vector<double>& pos);

    std::vector<double> get_velocity_i(
      const double lambda,
      const double thetap,
      const std::vector<double>& pos);

    void rungekutta_for_nbody(
      std::vector<double>& r,
      std::vector<double>& v,
      const std::vector<double>& mass,
      const double dt,
      const int rk_order);

    void prn_vector(
      std::ofstream& ou,
      const std::string tag,
      const std::vector<double>& vec);

    double get_excentricity(
      const double lambda,
      const double c_alpha,
      const double theta_p);

    std::tuple<std::vector<double>, std::vector<double>> check_angular_momentum(
      const std::vector<double>& pos,
      const std::vector<double>& vel,
      const std::vector<double>& mass);

    double check_total_energy(
      const std::vector<double>& pos,
      const std::vector<double>& vel,
      const std::vector<double>& mass);
  protected:
};
}

#endif

