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

#include <cmath>
#include <numbers>
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <string>
#include <vector>
#include "nbody.h"

const constexpr static double Gconst = 1.0;
const constexpr static double Pi = std::numbers::pi;

double nbody::Nbody::fderiv(
  const std::array<double, 6>& a,
  const double xi,
  const int m) {
// description:
//   find m-th derivative at x = xi of p(x)
//   multiplied by c(m) : c(m) * d^m / dx^m p(x)  (x = xi)
//   where c(m) = 1 / m!
//   p(x)= a[0] x^n + a[1] x^(n-1) +...+ a[n-1] x + a[n]
  const int n = a.size() - 1;
  const int m1 = m + 1;
  const int n1 = n + 1;
  std::vector<double> l(n1 * n1);

  for (int k = 1; k <= n; ++k)
    l[0 * n1 + k] = a[k];

  for (int k = 1; k <= m1; ++k) {
    l[k * n1 + 0] = (k == 1) ? a[0] :  l[(k - 1) * n1 + 0];
    for (int j = 1; j <= n - k + 1; ++j) {
      l[k * n1 + j] = l[(k - 1) * n1 + j] + xi * l[k * n1 + j - 1];
    }
  }

  return l[m1 * n1 + n - m1 + 1];
}

double nbody::Nbody::newton_raphson(const std::array<double, 6>& c) {
// Description:
//   Finds root x of P(x) = 0 using Newton-Raphson method
//   P(X) = c[0] X^N + c[1] X^(N-1) + ...+ c[N-1] X + c[N] = 0
//   c, polynomial coefficients
  const double x0 = 0.0;
  const int niterations = 300;
  double x = x0;
  const double thresh = 1e-12;
  for (int i = 0; i != niterations; ++i) {
    const double dx = fderiv(c, x, 0) / fderiv(c, x, 1);
    if (fabs(dx) < thresh)
      break;
    x -= dx;
  }

  return x;
}

void nbody::Nbody::positions_wrt_cm(
  std::vector<double>& pos,
  const std::vector<double>& mass) {
  const int nbody = pos.size() / 3;
  double total_mass = 0.0;
  for (int i = 0; i != nbody; ++i)
    total_mass += mass[i];
  std::array<double, 3> rcm{0.0, 0.0, 0.0};  // rcm, center-of-mass
  for (int k = 0; k != 3; ++k) {
    for (int i = 0; i != nbody; ++i)
      rcm[k] += mass[i] * pos[i + nbody * k];
    rcm[k] /= total_mass;
  }
  for (int i = 0; i != nbody; ++i)
    for (int k = 0; k != 3; ++k)
      pos[i + nbody * k] -= rcm[k];
}

std::tuple<double, std::vector<double>> nbody::Nbody::linear_configuration(
  const int config,
  const std::vector<double>& r_ini,
  const std::vector<double>& mass) {
// Purpose: generate three body linear configuration
// Note.- Bodies are located in X-axis
// config = 1 (321), 2 (231), 3 (213)
  std::array<int, 3> k{(config == 1) ? (3 - 1) : 2 - 1,
                       (config == 1) ? (2 - 1) : (config == 2) ? 3 - 1 : 1 - 1,
                       (config == 1 || config == 2) ? (1 - 1) : 3 - 1};
  const std::array<double, 6> a{mass[k[2]] + mass[k[1]],
                                3.0 * mass[k[2]] + 2.0 * mass[k[1]],
                                3.0 * mass[k[2]] + mass[k[1]],
                                -(mass[k[1]] + 3.0 * mass[k[0]]),
                                -(2.0 * mass[k[1]] + 3.0 * mass[k[0]]),
                                -(mass[k[1]] + mass[k[0]])};
  const double x = newton_raphson(a);  // Find x for P(x) = 0
  // condition: pos[k[1]] < pos[k[2]] < pos[k[3}],  for x-coordinate
  const int nbody = 3;
  std::vector<double> pos(3 * nbody);
  pos[0] = r_ini[0];
  pos[1] = r_ini[1];
  // Note.- only work out initial position along X-axis
  if (config == 1) pos[2] = pos[1] - x * (pos[0] - pos[1]);
  if (config == 2) pos[2] = (pos[1] + x * pos[0]) / ( 1.0 + x);
  if (config == 3) pos[2] = pos[0] + 1.0 / x * (pos[0] - pos[1]);

  return make_tuple(x, pos);
}

std::vector<double> nbody::Nbody::gen_vertices(
  const double radius,
  const int nbody) {
  const double dbeta = 2.0 * Pi / static_cast<double>(nbody);
  std::vector<double> pos(3 * nbody);
  double beta_i = 0.0;
  for (int i = 0; i != nbody; ++i) {
    pos[i + nbody * 0] = radius * cos(beta_i);
    pos[i + nbody * 1] = radius * sin(beta_i);
    pos[i + nbody * 2] = 0.0;
    beta_i += dbeta;
  }

  return pos;
}

std::tuple<std::vector<double>, std::vector<double>>
  nbody::Nbody::gen_vertices_for_4bodies(
  const double rbase,    // = r_{12}
  const double kappa,    // kappa >= 1
  const double mass3,    // third mass value
  const int selec_nu) {  // = 1, 2
  const double k3 = kappa * kappa * kappa;
  const double ksq = 1.0 + kappa * kappa;
  const double kspq = exp(1.5 * log(fabs(ksq)));
  const double f = (k3 * kspq + 1.0) / (k3 + kspq);
  const double radic = f * f - 1.0;
  if (radic < 0) {
    std::ostringstream tag;
    tag << "trapezoidal configuration:: f^2 < 1, f = "
        << std::setw(6) << std::setprecision(2) << std::fixed << f;
    throw std::runtime_error(tag.str());
  }
  const double coef = (selec_nu == 1) ? -1.0 : (selec_nu == 2) ? 1.0 : 0.0;
  if (coef == 0.0) {
    const std::string tag = "error selec_nu ne 1 or 2, selec_nu = " +
                            std::to_string(selec_nu);
    throw std::runtime_error(tag);
  }
  const double nu = f + coef * sqrt(radic);
  const double rnu_1 = exp(-0.666666667 * log(fabs(nu)));
  const double r34 = rbase * rnu_1;
  const double rnu_2 = exp(-0.333333333 * log(fabs(nu)));
  const double radicl = exp(0.5 * log(fabs(ksq)));
  const double r13 = rbase * radicl * rnu_2;
  // const double r23 = rbase * kappa * rnu_2;
  const double r12 = rbase;

  const double tq = (r12 - r34) / 2.0;
  const double rt = sqrt(r13 * r13 - (r12 - tq) * (r12 - tq));
  const std::vector<double> pos{0.0, r12, r12 - tq, tq,    // X
                                0.0, 0.0, rt      , rt,    // Y
                                0.0, 0.0, 0.0     , 0.0};  // Z
  const double g = (nu - 1.0 / kspq) / (1.0 - nu / kspq) / rnu_2;
  const std::vector<double> mass{mass3 * g, mass3 * g, mass3, mass3};

  return make_tuple(mass, pos);
}

std::vector<double> nbody::Nbody::get_alpha_for_3bodies(
  const std::vector<double>& mass) {
  const double mt = mass[0] + mass[1] + mass[2];
  const double mt2 = mt * mt;
  std::vector<double> mv(3);
  const double m2_1 = mass[0] * mass[0];
  const double m2_2 = mass[1] * mass[1];
  const double m2_3 = mass[2] * mass[2];
  mv[0] = m2_2 + m2_3 + mass[1] * mass[2];
  mv[1] = m2_1 + m2_3 + mass[0] * mass[2];
  mv[2] = m2_1 + m2_2 + mass[0] * mass[1];
  mv[0] = Gconst * mv[0] * sqrt(mv[0]) / mt2;
  mv[1] = Gconst * mv[1] * sqrt(mv[1]) / mt2;
  mv[2] = Gconst * mv[2] * sqrt(mv[2]) / mt2;

  return mv;
}

std::vector<double> nbody::Nbody::get_alpha_for_linear_configuration(
  const int config,
  const double x,
  const std::vector<double>& mass) {
// config = 1 (321), 2 (231), 3 (213)
  std::array<int, 3> k{(config == 1) ? (3 - 1) : 2 - 1,
                       (config == 1) ? (2 - 1) : (config == 2) ? 3 - 1 : 1 - 1,
                       (config == 1 || config == 2) ? (1 - 1) : 3 - 1};
  const double total_mass = mass[0] + mass[1] + mass[2];
  std::array<double, 3> v{1 + x, -x, x / (1.0 + x)};

  std::vector<double> mv(3);
  for (int i = 0; i != 3; ++i) {
    const int i1 = (i == 0) ? k[1] : k[2];
    const int i2 = (i == 2) ? k[1] : k[0];
    double s = (mass[i1] + mass[i2] * v[i]) / total_mass;
    s *= s;
    const double ss = (i == 1) ? -1.0 : 1.0;
    const double c1 = ((i == 1) && (mass[k[0]] < mass[k[2]])) ? -1.0 : 1.0;
    mv[k[2 - i]] = Gconst * s * c1 * (ss * mass[i1] + mass[i2] / (v[i] * v[i]));
  }

  return mv;
}

std::vector<double> nbody::Nbody::get_alpha_for_4bodies(
  const std::vector<double>& pos,
  const std::vector<double>& mass) {
  const int nbody = pos.size() / 3;
  const std::vector<double> acel = get_forces(mass, pos);
  std::vector<double> mv(nbody);
  for (int i = 0; i != nbody; ++i) {
    double acum = 0.0;
    for (int j = 0; j != 3; ++j) {
      acum += pos[i + nbody * j] * acel[i + nbody * j];
    }
    const int ind_x = i + nbody * 0;
    const int ind_y = i + nbody * 1;
    const double dist_i = sqrt(pos[ind_x] * pos[ind_x] +
                               pos[ind_y] * pos[ind_y]);
    mv[i] = -Gconst * dist_i * acum;
  }

  return mv;
}

std::vector<double> nbody::Nbody::get_alpha_for_regular_polygon(
  const std::vector<double>& mass) {
  const int nbody = mass.size();
  const double nc = nbody % 2 ? nbody : nbody + 1;
  const int nn = static_cast<int>(nc / 2.0 - 1);
  const double dtheta = Pi / static_cast<double>(nbody);
  double s = 0.0;
  double theta_j = 0.0;
  for (int j = 1; j <= nn; ++j) {
    theta_j += dtheta;
    s += 1.0 / sin(theta_j);
  }
  const double s1 = nbody % 2 ? 1.0 : 0.0;
  const double beta = Gconst * (s + s1 / 2.0) / 2.0;

  std::vector<double> mv(nbody);
  for (int i = 0; i != nbody; ++i)
    mv[i] = mass[i] * beta;

  return mv;
}

void nbody::Nbody::get_ini_posvel(
  const nlohmann::json& input_pars) {
// config = 0, equilateral triangle
//        = 1,2,3 en linear configuration types 321,231,213 respectively
//        = 4 four bodies in isosceles-trapezoid
//        = 5 bodies on vertices of regular polygon
// radius, regular polygon radius
// lambda, relative angle between vector position and vector velocity
//       (same angle for all bodies)
// thetap, angular velocity (same for all bodies)
  const int config = input_pars["config"];
  const int nbody = input_pars["nbody"];
  std::vector<double> mv(nbody);
  switch (config) {
    case 0:  // equilateral triangle configuration
      {
        mass_ = input_pars["mass"].get<std::vector<double>>();
        const double radius = input_pars["radius"];
        pos_ = gen_vertices(radius, nbody);
        mv = get_alpha_for_3bodies(mass_);
      }
      break;
    case 1: case 2: case 3:  // linear configuration (types: 321,231,213)
      {
        double x;
        mass_ = input_pars["mass"].get<std::vector<double>>();
        const std::vector<double> x_ini = input_pars["x_ini"];
        tie(x, pos_) = linear_configuration(config, x_ini, mass_);
        mv = get_alpha_for_linear_configuration(config, x, mass_);
      }
      break;
    case 4:  // isosceles trapezoid configuration (4 bodies)
      {
        const double rbase = input_pars["r12"];
        const double mass_3 = input_pars["mass_3"];
        const double kappa = input_pars["kappa"];
        const int nu = input_pars["nu"];
        tie(mass_, pos_) = gen_vertices_for_4bodies(rbase, kappa, mass_3, nu);
      }
      break;
    case 5:  // regular polygon configuration
      {
        const double mass_i = input_pars["mass_i"];
        const double radius = input_pars["radius"];
        for (int i = 0; i != nbody; ++i) mass_[i] = mass_i;
        pos_ = gen_vertices(radius, nbody);
        mv = get_alpha_for_regular_polygon(mass_);
      }
      break;
  }

  positions_wrt_cm(pos_, mass_);

  if (config == 4) {
    mv = get_alpha_for_4bodies(pos_, mass_);
  }

  const double lambda = input_pars["lambda"];
  for (int i = 0; i != nbody; ++i) {
    const int ind_x = i + nbody * 0;
    const int ind_y = i + nbody * 1;
    const double md = sqrt(pos_[ind_x] * pos_[ind_x] +
                           pos_[ind_y] * pos_[ind_y]);
    const double to_c_alpha = sqrt(mv[i] / md) / md * sin(lambda * Pi /180.0);
    // Note.- c_alpha[i] (i = 0, nbody - 1) should be equal
    c_alpha_[i] = to_c_alpha;
  }
  const double dtheta0 = sqrt(2.0) * c_alpha_[0];
  const double thetap = input_pars["angular_velocity"];
  const double thetap_min = -dtheta0;
  const double thetap_max =  dtheta0;
  if (!(thetap >= thetap_min && thetap <= thetap_max)) {
    std::ostringstream os;
    os << "initial angular velocity, thetap not in range: ["
       << std::setw(8) << std::setprecision(3)
       << std::fixed << thetap_min << ","
       << std::setw(8) << std::setprecision(3)
       << std::fixed << thetap_max << "]";
    throw std::runtime_error(os.str());
  }
  vel_ = get_velocity_i(lambda, thetap, pos_);
}

std::vector<double> nbody::Nbody::get_velocity_i(
  const double lambda,
  const double thetap,
  const std::vector<double>& pos) {
  const double ang = lambda * Pi / 180.0;
  const int nbody = pos.size() / 3;
  std::vector<double> v(nbody * 3);
  for (int i = 0; i < nbody; ++i) {
    const int ind_x = i + nbody * 0;
    const int ind_y = i + nbody * 1;
    const double ri = sqrt(pos[ind_x] * pos[ind_x] + pos[ind_y] * pos[ind_y]);
    const std::array<double, 2> u{pos[ind_x] / ri, pos[ind_y] / ri};
    const double mod_vel = ri * thetap / sin(ang);
    v[ind_x] = mod_vel * (u[0] * cos(ang) - u[1] * sin(ang));
    v[ind_y] = mod_vel * (u[1] * cos(ang) + u[0] * sin(ang));
  }

  return v;
}

void nbody::Nbody::rungekutta_for_nbody(
  std::vector<double>& r,
  std::vector<double>& v,
  const std::vector<double>& mass,
  const double dt,
  const int rk_order) {
  const int nbody = r.size() / 3;
  std::array<double, 4> c{0.0, 0.0, 0.0, 0.0};
  std::array<double, 4> d{0.0, 0.0, 0.0, 0.0};
  switch (rk_order) {
    case 1:  // rk_order=1; Euler's method
      d[0] = 0.0;
      c[0] = 1.0;
      break;
    case 2:  // rk_order=2; Improved Euler's method
      d[0] = 0.0;
      d[1] = 1.0 / 2.0;
      c[0] = 0.0;
      c[1] = 1.0;
      break;
    case 3:  // rk_order = 3
      d[0] = 0.0;
      d[1] = 1.0 / 2.0;
      d[2] = 3.0 / 4.0;
      c[0] = 2.0 / 9.0;
      c[1] = 3.0 / 9.0;
      c[2] = 4.0 / 9.0;
      break;
    case 4:  // rk_order = 4
      d[0] = 0.0;
      d[1] = 1.0 / 2.0;
      d[2] = 1.0 / 2.0;
      d[3] = 1.0;
      c[0] = 1.0 / 6.0;
      c[1] = 2.0 / 6.0;
      c[2] = 2.0 / 6.0;
      c[3] = 1.0 / 6.0;
      break;
  }

  std::vector<double> r0(3 * nbody);
  std::vector<double> v0(3 * nbody);
  std::copy_n(r.data(), r.size(), r0.data());
  std::copy_n(v.data(), v.size(), v0.data());
  std::vector<double> dr(3 * nbody);
  std::vector<double> dv(3 * nbody);
  std::vector<double> k1(3 * nbody);
  std::vector<double> k2(3 * nbody);
  for (int k = 0; k < rk_order; ++k) {
    for (int i = 0; i < nbody; ++i) {
      for (int j = 0; j < 3; ++j) {
        const int ind = i + nbody * j;
        k1[ind] = r0[ind] + dr[ind] * d[k];
        k2[ind] = v0[ind] + dv[ind] * d[k];
      }
    }
    const std::vector<double> acel = get_forces(mass, k1);
    for (int i = 0; i < nbody; ++i) {
      for (int j = 0; j < 3; ++j) {
        const int ind = i + nbody * j;
        dr[ind] = dt * k2[ind];
        dv[ind] = dt * acel[ind];
        r[ind] += dr[ind] * c[k];
        v[ind] += dv[ind] * c[k];
      }
    }
  }
}

std::vector<double> nbody::Nbody::get_forces(
  const std::vector<double>& mass,
  const std::vector<double>& pos) {
  const int nbody = pos.size() / 3;
  std::vector<double> forces(3 * nbody);
  for (int i = 0; i != nbody; ++i) {
    std::array<double, 3> ac{0.0, 0.0, 0.0};
    for (int j = 0; j != nbody; ++j) {
      if (i != j) {
        const std::array<double, 3> rij{
          pos[0 * nbody + j] - pos[0 * nbody + i],
          pos[1 * nbody + j] - pos[1 * nbody + i],
          pos[2 * nbody + j] - pos[2 * nbody + i]};
        const double rij_mod0 = rij[0] * rij[0] +
                                rij[1] * rij[1] +
                                rij[2] * rij[2];
        const double rij3 = sqrt(rij_mod0 * rij_mod0 * rij_mod0);
        for (int k = 0; k != 3; ++k) {
          ac[k] += mass[j] * rij[k] / rij3;
        }
      }
    }
    for (int k = 0; k != 3; ++k) {
      forces[i + k * nbody] = Gconst * ac[k];
    }
  }

  return forces;
}

double nbody::Nbody::get_excentricity(
  const double lambda,
  const double c_alpha,
  const double theta_p) {
  const double c1 = theta_p / sin(Pi / 180.0 * lambda) / 2.0;
  const double c3 = theta_p / c_alpha;

  return sqrt(1.0 + 2.0 * (c1 - c_alpha) * c3 * c3);
}

std::tuple<std::vector<double>, std::vector<double>>
  nbody::Nbody::check_angular_momentum(
  const std::vector<double>& pos,
  const std::vector<double>& vel,
  const std::vector<double>& mass) {
  const int nbody = pos.size() / 3;
  std::vector<double> ang(3 * nbody);
  std::vector<double> total_ang(3);
  for (int i = 0; i != nbody; ++i) {
    ang[i + nbody * 0] = mass[i] * (pos[i + nbody * 1] * vel[i +nbody * 2] -
                                    pos[i + nbody * 2] * vel[i +nbody * 1]);
    ang[i + nbody * 1] = mass[i] * (pos[i + nbody * 2] * vel[i +nbody * 0] -
                                    pos[i + nbody * 0] * vel[i +nbody * 2]);
    ang[i + nbody * 2] = mass[i] * (pos[i + nbody * 0] * vel[i +nbody * 1] -
                                    pos[i + nbody * 1] * vel[i +nbody * 0]);
    for (int k = 0; k != 3; ++k) {
      total_ang[k] += ang[i + nbody * k];
    }
  }

  return make_tuple(total_ang, ang);
}

double nbody::Nbody::check_total_energy(
  const std::vector<double>& pos,
  const std::vector<double>& vel,
  const std::vector<double>& mass) {
  double energy = 0.0;
  const int nbody = pos.size() / 3;
  for (int i = 0; i != nbody; ++i) {
    const double vel_i2 = vel[i + nbody * 0] * vel[i + nbody * 0] +
                          vel[i + nbody * 1] * vel[i + nbody * 1] +
                          vel[i + nbody * 2] * vel[i + nbody * 2];
    const double kinetic_i = mass[i] / 2.0 * vel_i2;
    double potential_energy_i = 0.0;
    for (int j = 0; j != i; ++j) {
      const std::array<double, 3> rij{
        pos[0 * nbody + j] - pos[0 * nbody + i],
        pos[1 * nbody + j] - pos[1 * nbody + i],
        pos[2 * nbody + j] - pos[2 * nbody + i]};
      const double rij_mod0 = sqrt(rij[0] * rij[0] +
                                   rij[1] * rij[1] +
                                   rij[2] * rij[2]);
      potential_energy_i -= mass[j] / rij_mod0;
    }
    potential_energy_i *= Gconst * mass[i];
    energy += kinetic_i + potential_energy_i;
  }
  return energy;
}

void nbody::Nbody::prn_vector(
  std::ofstream& ou,
  const std::string tag,
  const std::vector<double>& vec) {
  ou << "  \"" + tag + "\" : [";
  for (int j = 0; j != static_cast<int>(vec.size()); ++j) {
    if (j != static_cast<int>(vec.size() - 1)) {
      ou << std::setw(15) << std::setprecision(8)
         << std::scientific << vec[j] << ",";
    } else {
      ou << std::setw(15) << std::setprecision(8)
         << std::scientific << vec[j] << "],\n";
    }
  }
}

void nbody::Nbody::nbody_simulator() {
  get_ini_posvel(input_pars_);
  const double dt = input_pars_["dt"];
  const double tfinal = input_pars_["tfinal"];
  const int nsteps = tfinal / dt;
  std::ofstream ou("nbody_output.json");
  if (!ou.is_open()) {
    throw std::runtime_error("file nbody_output.json could not be opened");
  }
  ou << "{\n";
  const int ndata = input_pars_["ndata"];
  if (ndata > nsteps) {
    const std::string tag = "ndata > nsteps, (ndata,nsteps) = " +
                            std::to_string(ndata) + "," +
                            std::to_string(nsteps) + ")";
    throw std::runtime_error(tag);
  }
  double t = 0.0;
  for (int i = 0; i != nsteps; ++i) {
    const std::string tag = "\"step_" + std::to_string(i) + "\" : {\n";
    rungekutta_for_nbody(pos_, vel_, mass_, dt, 4);
    std::vector<double> total_ang;
    std::vector<double> ang;
    const double energy = check_total_energy(pos_, vel_, mass_);
    tie(total_ang, ang) = check_angular_momentum(pos_, vel_, mass_);
    if (i < ndata) {
      energy_[i] = energy;
      std::copy_n(total_ang.data(), total_ang.size(),
                  angular_momentum_.data() + i * total_ang.size());
      std::copy_n(ang.data(), ang.size(),
                  angular_momentum_i_.data() + i * ang.size());
    }
    // output to a json file:
    ou << tag;
    prn_vector(ou, "positions", pos_);
    prn_vector(ou, "velocities", vel_);
    prn_vector(ou, "angular_momentum", total_ang);
    prn_vector(ou, "angular_momentum_i", ang);
    ou << "  \"energy\" : "
       << std::setw(15) << std::setprecision(8)
       << std::scientific << energy << "\n";
    if (i != nsteps - 1)
      ou << "},\n";
    else
      ou << "}\n";
    t += dt;
  }
  ou << "}\n";
  ou.close();
}

nbody::Nbody::Nbody(const nlohmann::json input_pars0) :
  input_pars_(input_pars0), nbody_(input_pars0["nbody"]),
  ndata_(input_pars0["ndata"]),
  c_alpha_(nbody_), pos_(3 * nbody_), vel_(3 * nbody_), mass_(nbody_),
  energy_(ndata_), angular_momentum_(3 * ndata_),
  angular_momentum_i_(3 * ndata_ * nbody_) {
  nbody_simulator();
}
