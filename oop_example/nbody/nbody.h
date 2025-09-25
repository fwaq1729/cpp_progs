/*
  This file is part of nbody.
 
  Copyright (C) 2025 Fredy W. Aquino
 
  nbody is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  nbody is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.
 
  You should have received a copy of the GNU Lesser General Public License
  along with nbody.  If not, see <https://www.gnu.org/licenses/>.

  Description : Analytical solutions to Gravitational N-Body problem
  Date        : 09-16-25
*/

#ifndef OOP_EXAMPLE_NBODY_NBODY_H_
#define OOP_EXAMPLE_NBODY_NBODY_H_
#include <string>
#include <array>
#include <vector>
#include <tuple>
#include <sstream>
#include <fstream>
#include <nlohmann/json.hpp>

namespace nbody {
template<typename T>
class Nbody {
 public:
  explicit Nbody(const nlohmann::json input_pars);
  void nbody_simulator();
  int get_nbody() { return nbody_; }
  int get_ndata() { return ndata_; }
  const std::vector<T>& get_c_alpha() { return c_alpha_; }
  const std::vector<T>& get_energy() { return energy_; }
  const std::vector<T>& get_angular_momentum()
    { return angular_momentum_; }
  const std::vector<T>& get_angular_momentum_i()
    { return angular_momentum_i_; }

 private:
  nlohmann::json input_pars_;
  int nbody_;
  int ndata_;
  std::vector<T> c_alpha_;
  std::vector<T> pos_;
  std::vector<T> vel_;
  std::vector<T> mass_;
  std::vector<T> energy_;
  std::vector<T> angular_momentum_;
  std::vector<T> angular_momentum_i_;

  void get_ini_posvel(
    const nlohmann::json& input_pars);
  std::tuple<T, std::vector<T>> linear_configuration(
    const int config,
    const std::vector<T>& r_ini,
    const std::vector<T>& mass);

  std::vector<T> gen_vertices(
    const T radius,
    const int nbody);

  T fderiv(
    const std::array<T, 6>& a,
    const T xi,
    const int m);

  T newton_raphson(const std::array<T, 6>& c);

  void positions_wrt_cm(
    std::vector<T>& pos,
    const std::vector<T>& mass);

  std::vector<T> get_alpha_for_3bodies(
    const std::vector<T>& mass);

  std::vector<T> get_alpha_for_linear_configuration(
    const int config,
    const T x,
    const std::vector<T>& mass);

  std::vector<T> get_alpha_for_4bodies(
    const std::vector<T>& pos,
    const std::vector<T>& mass);

  std::vector<T>
    get_alpha_for_regular_polygon(const std::vector<T>& mass);

  std::tuple<std::vector<T>, std::vector<T>>
    gen_vertices_for_4bodies(
    const T rbase,
    const T kappa,
    const T mass3,
    const int selec_nu);

  std::vector<T> get_forces(
    const std::vector<T>& mass,
    const std::vector<T>& pos);

  std::vector<T> get_velocity_i(
    const T lambda,
    const T thetap,
    const std::vector<T>& pos);

  void rungekutta_for_nbody(
    std::vector<T>& r,
    std::vector<T>& v,
    const std::vector<T>& mass,
    const T dt,
    const int rk_order);

  void prn_vector(
    std::ofstream& ou,
    const std::string tag,
    const std::vector<T>& vec);

  T get_excentricity(
    const T lambda,
    const T c_alpha,
    const T theta_p);

  std::tuple<std::vector<T>, std::vector<T>>
    check_angular_momentum(
    const std::vector<T>& pos,
    const std::vector<T>& vel,
    const std::vector<T>& mass);

  T check_total_energy(
    const std::vector<T>& pos,
    const std::vector<T>& vel,
    const std::vector<T>& mass);

 protected:
};
}  // namespace nbody

#endif  // OOP_EXAMPLE_NBODY_NBODY_H_

