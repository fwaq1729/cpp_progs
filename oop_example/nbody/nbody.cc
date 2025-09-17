
#include <cmath>
#include <numbers>
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
#include "nbody.h"

using namespace std;

const constexpr static double Gconst = 1.0;
const constexpr static double Pi = std::numbers::pi;

double nbody::Nbody::fderiv(
  const array<double, 6>& a,
  const double xi,
  const int m)
{
// description:
//   find m-th derivative at x = xi of p(x)
//   multiplied by c(m) : c(m) * d^m / dx^m p(x)  (x = xi)
//   where c(m) = 1 / m!
//   p(x)= a[0] x^n + a[1] x^(n-1) +...+ a[n-1] x + a[n]
  const int n = a.size() - 1;
  const int m1 = m + 1;
  const int n1 = n + 1;
  vector<double> l(n1 * n1);

  for (int k = 1; k <= n; ++k)
    l[0 * n1 + k] = a[k];

  for (int k = 1; k <= m1; ++k)
  {
    l[k * n1 + 0] = (k == 1) ? a[0] :  l[(k - 1) * n1 + 0];
    for (int j = 1; j <= n - k + 1; ++j)
    {
      l[k * n1 + j] = l[(k - 1) * n1 + j] + xi * l[k * n1 + j - 1];
    }
  }

  return l[m1 * n1 + n - m1 + 1];
}

double nbody::Nbody::newton_raphson(const array<double, 6>& c)
{
// Description:
//   Finds root x of P(x) = 0 using Newton-Raphson method
//   P(X) = c[0] X^N + c[1] X^(N-1) + ...+ c[N-1] X + c[N] = 0
//   c, polynomial coefficients
  const double x0 = 0.0;
  const int niterations = 300;
  double x = x0;
  for (int i = 0; i != niterations; ++i)
  {
    x -= fderiv(c, x, 0) / fderiv(c, x, 1);
    // Note.- here I need to put a thresh so that it skips loop
  }

  return x;
}

void nbody::Nbody::positions_wrt_cm(
  vector<double>& pos,
  const vector<double>& mass)
{
  const int nbody = pos.size() / 3;
  double total_mass = 0.0;
  for (int i = 0; i != nbody; ++i)
    total_mass += mass[i];
  array<double, 3> rcm{0.0, 0.0, 0.0}; // rcm, center-of-mass
  for (int k = 0; k != 3; ++k)
  {
    for (int i = 0; i != nbody; ++i)
      rcm[k] += mass[i] * pos[i + nbody * k];
    rcm[k] /= total_mass;
  }
  for (int i = 0; i != nbody; ++i)
    for (int k = 0; k != 3; ++k)
      pos[i + nbody * k] -= rcm[k];
}

tuple<double, vector<double>> nbody::Nbody::linear_configuration(
  const int config,
  const vector<double>& r_ini,
  const vector<double>& mass)
{
// Purpose: generate three body linear configuration
// Note.- Bodies are located in X-axis
  array<int, 3> k;
  if (config == 1) { k[0] = 3 - 1; k[1] = 2 - 1; k[2] = 1 - 1;} // config. 321
  if (config == 2) { k[0] = 2 - 1; k[1] = 3 - 1; k[2] = 1 - 1;} // config. 231
  if (config == 3) { k[0] = 2 - 1; k[1] = 1 - 1; k[2] = 3 - 1;} // config. 213
  const array<double, 6> a{mass[k[2]] + mass[k[1]],
                           3.0 * mass[k[2]] + 2.0 * mass[k[1]],
                           3.0 * mass[k[2]] + mass[k[1]],
                           -(mass[k[1]] + 3.0 * mass[k[0]]),
                           -(2.0 * mass[k[1]] + 3.0 * mass[k[0]]),
                           -(mass[k[1]] + mass[k[0]])};
  const double x = newton_raphson(a); // Find x for P(x) = 0
  // condition: pos[k[1]] < pos[k[2]] < pos[k[3}],  for x-coordinate
  const int nbody = 3;
  vector<double> pos(3 * nbody);
  pos[0] = r_ini[0];
  pos[1] = r_ini[1];
  // Note.- only work out initial position along X-axis
  if (config == 1) pos[2] = pos[1] - x * (pos[0] - pos[1]);
  if (config == 2) pos[2] = (pos[1] + x * pos[0]) / ( 1.0 + x);
  if (config == 3) pos[2] = pos[0] + 1.0 / x * (pos[0] - pos[1]);

  return make_tuple(x, pos);
}

vector<double> nbody::Nbody::gen_vertices(
  const double radius,
  const int nbody)
{
  const double dbeta = 2.0 * Pi / (double)nbody;
  vector<double> pos(3 * nbody);
  for (int i = 0; i != nbody; ++i)
  {
    const double beta_i = (double)i * dbeta;
    pos[i + nbody * 0] = radius * cos(beta_i);
    pos[i + nbody * 1] = radius * sin(beta_i);
    pos[i + nbody * 2] = 0.0;
  }

  return pos;
}

tuple<vector<double>, vector<double>> nbody::Nbody::gen_vertices_for_4bodies(
  const double rbase, // = r_{12}
  const double kappa, // kappa >= 1
  const double mass3, // third mass value
  const int selec_nu) // = 1, 2
{
  const double k3 = kappa * kappa * kappa;
  const double ksq = 1.0 + kappa * kappa;
  const double kspq = exp(1.5 * log(fabs(ksq)));
  const double f = (k3 * kspq + 1.0) / (k3 + kspq);
  const double radic = f * f - 1.0;
  if (radic < 0)
  {
    ostringstream tag;
    tag << "trapezoidal configuration:: f^2 < 1, f = "
        << setw(6) << setprecision(2) << fixed << f;
    throw runtime_error(tag.str());
  }
  const double coef = (selec_nu == 1) ? -1.0 : (selec_nu == 2) ? 1.0 : 0.0;
  if (coef == 0.0)
  {
    const string tag = "error selec_nu ne 1 or 2, selec_nu = " + to_string(selec_nu);
    throw runtime_error(tag);
  }
  const double nu = f + coef * sqrt(radic);
  const double rnu_1 = exp(-0.666666667 * log(fabs(nu)));
  const double r34 = rbase * rnu_1;
  const double rnu_2 = exp(-0.333333333 * log(fabs(nu)));
  const double radicl = exp(0.5 * log(fabs(ksq)));
  const double r13 = rbase * radicl * rnu_2;
  //const double r23 = rbase * kappa * rnu_2;
  const double r12 = rbase;

  const double tq = (r12 - r34) / 2.0;
  const double rt = sqrt(r13 * r13 - (r12 - tq) * (r12 - tq));
  const vector<double> pos{0.0, r12, r12 - tq, tq,   // X
                           0.0, 0.0, rt      , rt,   // Y
                           0.0, 0.0, 0.0     , 0.0}; // Z
  const double g = (nu - 1.0 / kspq) / (1.0 - nu / kspq) / rnu_2;
  const vector<double> mass{mass3 * g, mass3 * g, mass3, mass3};

  return make_tuple(mass, pos);
}

vector<double> nbody::Nbody::get_alpha_for_3bodies(
  const vector<double>& mass)
{
  const double mt = mass[0] + mass[1] + mass[2];
  const double mt2 = mt * mt;
  vector<double> mv(3);
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

vector<double> nbody::Nbody::get_alpha_for_linear_configuration(
  const int config,
  const double x,
  const vector<double>& mass)
{
  array<int, 3> k = {0, 0, 0};
  if (config == 1) { k[0] = 3 - 1; k[1] = 2 - 1; k[2] = 1 - 1;} // config. 321
  if (config == 2) { k[0] = 2 - 1; k[1] = 3 - 1; k[2] = 1 - 1;} // config. 231
  if (config == 3) { k[0] = 2 - 1; k[1] = 1 - 1; k[2] = 3 - 1;} // config. 213
  const double total_mass = mass[0] + mass[1] + mass[2];
  array<double, 3> v{1 + x, -x, x / (1.0 + x)};

  vector<double> mv(3);
  for (int i = 0; i != 3; ++i)
  {
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

vector<double> nbody::Nbody::get_alpha_for_4bodies(
  const vector<double>& pos,
  const vector<double>& mass)
{
  const int nbody = pos.size() / 3;
  const vector<double> acel = get_forces(mass, pos);
  vector<double> mv(nbody);
  for (int i = 0; i != nbody; ++i)
  {
    double acum = 0.0;
    for (int j = 0; j != 3; ++j)
    {
      acum += pos[i + nbody * j] * acel[i + nbody * j];
    }
    const int ind_x = i + nbody * 0;
    const int ind_y = i + nbody * 1;
    const double dist_i = sqrt(pos[ind_x] * pos[ind_x] + pos[ind_y] * pos[ind_y]);
    mv[i] = -Gconst * dist_i * acum;
  }

  return mv;
}

vector<double> nbody::Nbody::get_alpha_for_regular_polygon(
  const vector<double>& mass)
{
  const int nbody = mass.size();
  const double nc = nbody % 2 ? nbody : nbody + 1;
  const int nn = (int)(nc / 2.0 - 1);
  const double dtheta = Pi / (double)nbody;
  double s = 0.0;
  for (int j = 1; j <= nn; ++j)
  {
    s += 1.0 / sin((double)j * dtheta);
  }
  const double s1 = nbody % 2 ? 1.0 : 0.0;
  const double beta = Gconst * (s + s1 / 2.0) / 2.0;

  vector<double> mv(nbody);
  for (int i = 0; i != nbody; ++i)
    mv[i] = Gconst * mass[i] * beta;

  return mv;
}

void nbody::Nbody::get_ini_posvel(
  const nlohmann::json& input_pars)
// config = 0, equilateral triangle
//        = 1,2,3 en linear configuration types 321,231,213 respectively
//        = 4 four bodies in isosceles-trapezoid
//        = 5 bodies on vertices of regular polygon
// radius, regular polygon radius
// lambda, relative angle between vector position and vector velocity
//       (same angle for all bodies)
// thetap, angular velocity (same for all bodies)
{
  const int config = input_pars["config"];
  const int nbody = input_pars["nbody"];
  vector<double> mv(nbody);
  switch(config)
  {
    case 0: // equilateral triangle configuration
      {
        mass_ = input_pars["mass"].get<vector<double>>();
        const double radius = input_pars["radius"];
        pos_ = gen_vertices(radius, nbody);
        mv = get_alpha_for_3bodies(mass_);
      }
      break;
    case 1: case 2: case 3: // linear configuration (types: 321,231,213)
      {
        double x;
        mass_ = input_pars["mass"].get<vector<double>>();
        const vector<double> x_ini = input_pars["x_ini"];
        tie(x, pos_) = linear_configuration(config, x_ini, mass_);
        mv = get_alpha_for_linear_configuration(config, x, mass_);
      }
      break;
    case 4: // isosceles trapezoid configuration (4 bodies)
      {
        const double rbase = input_pars["r12"];
        const double mass_3 = input_pars["mass_3"];
        const double kappa = input_pars["kappa"];
        const int nu = input_pars["nu"];
        tie(mass_, pos_) = gen_vertices_for_4bodies(rbase, kappa, mass_3, nu);
      }
      break;
    case 5: // regular polygon configuration
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

  if (config == 4)
  {
    mv = get_alpha_for_4bodies(pos_, mass_);
  }

  const double lambda = input_pars["lambda"];
  for (int i = 0; i != nbody; ++i)
  {
    const int ind_x = i + nbody * 0;
    const int ind_y = i + nbody * 1;
    const double md = sqrt(pos_[ind_x] * pos_[ind_x] + pos_[ind_y] * pos_[ind_y]);
    const double to_c_alpha = sqrt(mv[i] / md) / md * sin(lambda * Pi /180.0);
    // Note.- c_alpha[i] (i = 0, nbody - 1) should be equal
    c_alpha_[i] = to_c_alpha;
  }
  const double dtheta0 = sqrt(2.0) * c_alpha_[0];
  const double thetap = input_pars["angular_velocity"];
  const double thetap_min = -dtheta0;
  const double thetap_max =  dtheta0;
  if (!(thetap >= thetap_min && thetap <= thetap_max))
  {
    std::ostringstream os;
    os << "initial angular velocity, thetap not in range: ["
       << setw(8) << setprecision(3) << fixed << thetap_min << ","
       << setw(8) << setprecision(3) << fixed << thetap_max << "]";
    throw runtime_error(os.str());
  }
  vel_ = get_velocity_i(lambda, thetap, pos_);
}

vector<double> nbody::Nbody::get_velocity_i(
  const double lambda,
  const double thetap,
  const vector<double>& pos)
{
  const double ang = lambda * Pi / 180.0;
  const int nbody = pos.size() / 3;
  vector<double> v(nbody * 3);
  for (int i = 0; i < nbody; ++i)
  {
    const int ind_x = i + nbody * 0;
    const int ind_y = i + nbody * 1;
    const double ri = sqrt(pos[ind_x] * pos[ind_x] + pos[ind_y] * pos[ind_y]);
    const array<double, 2> u{pos[ind_x] / ri, pos[ind_y] / ri};
    const double mod_vel = ri * thetap / sin(ang);
    v[ind_x] = mod_vel * (u[0] * cos(ang) - u[1] * sin(ang));
    v[ind_y] = mod_vel * (u[1] * cos(ang) + u[0] * sin(ang));
  }

  return v;
}

void nbody::Nbody::rungekutta_for_nbody(
  vector<double>& r,
  vector<double>& v,
  const vector<double>& mass,
  const double dt,
  const int rk_order)
{
  const int nbody = r.size() / 3;
  array<double,4> c{0.0, 0.0, 0.0, 0.0};
  array<double,4> d{0.0, 0.0, 0.0, 0.0};
  switch(rk_order)
  {
    case 1: // rk_order=1; Euler's method
      d[0] = 0.0;
      c[0] = 1.0;
      break;
    case 2: // rk_order=2; Improved Euler's method
      d[0] = 0.0;
      d[1] = 1.0 / 2.0;
      c[0] = 0.0;
      c[1] = 1.0;
      break;
    case 3: // rk_order = 3
      d[0] = 0.0;
      d[1] = 1.0 / 2.0;
      d[2] = 3.0 / 4.0;
      c[0] = 2.0 / 9.0;
      c[1] = 3.0 / 9.0;
      c[2] = 4.0 / 9.0;
      break;
    case 4: // rk_order = 4
      d[0] = 0.0;
      d[1] = 1.0 / 2.0;
      d[2] = 1.0 / 2.0;
      d[3] = 1.0;
      c[0] = 1.0 / 6.0;
      c[1] = 2.0 / 6.0;
      c[2] = 2.0 / 6.0;
      c[3] = 1.0 / 6.0;
      break;
  };

  vector<double> r0(3 * nbody);
  vector<double> v0(3 * nbody);
  copy_n(r.data(), r.size(), r0.data());
  copy_n(v.data(), v.size(), v0.data());
  vector<double> dr(3 * nbody);
  vector<double> dv(3 * nbody);
  vector<double> k1(3 * nbody);
  vector<double> k2(3 * nbody);
  for (int k = 0; k < rk_order; ++k)
  {
    for (int i = 0; i < nbody; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        const int ind = i + nbody * j;
        k1[ind] = r0[ind] + dr[ind] * d[k];
        k2[ind] = v0[ind] + dv[ind] * d[k];
      }
    }
    const vector<double> acel = get_forces(mass, k1);
    for (int i = 0; i < nbody; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        const int ind = i + nbody * j;
        dr[ind] = dt * k2[ind];
        dv[ind] = dt * acel[ind];
        r[ind] += dr[ind] * c[k];
        v[ind] += dv[ind] * c[k];
      }
    }
  }
}

vector<double> nbody::Nbody::get_forces(
  const vector<double>& mass,
  const vector<double>& pos)
{
  const int nbody = pos.size() / 3;
  vector<double> forces(3 * nbody);
  for (int i = 0; i != nbody; ++i)
  {
    array<double, 3> ac{0.0, 0.0, 0.0};
    for (int j = 0; j != nbody; ++j)
    {
      if (i != j)
      {
        const array<double, 3> rij{pos[0 * nbody + j] - pos[0 * nbody + i],
                                   pos[1 * nbody + j] - pos[1 * nbody + i],
                                   pos[2 * nbody + j] - pos[2 * nbody + i]};
        const double rij_mod0 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
        const double rij3 = sqrt(rij_mod0 * rij_mod0 * rij_mod0);
        for (int k = 0; k != 3; ++k)
        {
          ac[k] += mass[j] * rij[k] / rij3;
        }
      }
    }
    for (int k = 0; k != 3; ++k)
    {
      forces[i + k * nbody] = Gconst * ac[k];
    }
  }

  return forces;
}

double nbody::Nbody::get_excentricity(
  const double lambda,
  const double c_alpha,
  const double theta_p)
{
  const double c1 = theta_p / sin(Pi / 180.0 * lambda) / 2.0;
  const double c3 = theta_p / c_alpha;

  return sqrt(1.0 + 2.0 * (c1 - c_alpha) * c3 * c3);
}

tuple<vector<double>, vector<double>> nbody::Nbody::check_angular_momentum(
  const vector<double>& pos,
  const vector<double>& vel,
  const vector<double>& mass)
{
  const int nbody = pos.size() / 3;
  vector<double> ang(3 * nbody);
  vector<double> total_ang(3);
  for (int i = 0; i != nbody; ++i)
  {
    ang[i + nbody * 0] = mass[i] * (pos[i + nbody * 1] * vel[i +nbody * 2] -
                                    pos[i + nbody * 2] * vel[i +nbody * 1]);
    ang[i + nbody * 1] = mass[i] * (pos[i + nbody * 2] * vel[i +nbody * 0] -
                                    pos[i + nbody * 0] * vel[i +nbody * 2]);
    ang[i + nbody * 2] = mass[i] * (pos[i + nbody * 0] * vel[i +nbody * 1] -
                                    pos[i + nbody * 1] * vel[i +nbody * 0]);
    for (int k = 0; k != 3; ++k)
    {
      total_ang[k] += ang[i + nbody * k];
    }
  }

  return make_tuple(total_ang, ang);
}

double nbody::Nbody::check_total_energy(
  const vector<double>& pos,
  const vector<double>& vel,
  const vector<double>& mass)
{
  double energy = 0.0;
  const int nbody = pos.size() / 3;
  for (int i = 0; i != nbody; ++i)
  {
    const double vel_i2 = vel[i + nbody * 0] * vel[i + nbody * 0] +
                          vel[i + nbody * 1] * vel[i + nbody * 1] +
                          vel[i + nbody * 2] * vel[i + nbody * 2];
    const double kinetic_i = mass[i] / 2.0 * vel_i2;
    double potential_energy_i = 0.0;
    for (int j = 0; j != i; ++j)
    {
      const array<double, 3> rij{pos[0 * nbody + j] - pos[0 * nbody + i],
                                 pos[1 * nbody + j] - pos[1 * nbody + i],
                                 pos[2 * nbody + j] - pos[2 * nbody + i]};
      const double rij_mod0 = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
      potential_energy_i -= mass[j] / rij_mod0;
    }
    potential_energy_i *= Gconst * mass[i];
    energy += kinetic_i + potential_energy_i;
  }
  return energy;
}

void nbody::Nbody::prn_vector(
  std::ofstream& ou,
  const string tag,
  const vector<double>& vec)
{
  ou << "  \"" + tag + "\" : [";
  for (int j = 0; j != (int)vec.size(); ++j)
  {
    if (j != (int)(vec.size() - 1))
    {
      ou << setw(15) << setprecision(8) << scientific << vec[j] << ",";
    } else
    {
      ou << setw(15) << setprecision(8) << scientific << vec[j] << "],\n";
    }
  }
}

void nbody::Nbody::nbody_simulator()
{
  get_ini_posvel(input_pars_);
  const double dt = input_pars_["dt"];
  const double tfinal = input_pars_["tfinal"];
  const int nsteps = tfinal / dt;
  std::ofstream ou("nbody_output.json");
  if (!ou.is_open())
  {
    throw runtime_error("file nbody_output.json could not be opened");
  }
  ou << "{\n";
  const int ndata = input_pars_["ndata"];
  if (ndata > nsteps)
  {
    const string tag = "ndata > nsteps, (ndata,nsteps) = " + to_string(ndata) + "," + to_string(nsteps) + ")";
    throw runtime_error(tag);
  }
  double t = 0.0;
  for (int i = 0; i != nsteps; ++i)
  {
    const string tag = "\"step_" + to_string(i) + "\" : {\n";
    rungekutta_for_nbody(pos_, vel_, mass_, dt, 4);
    vector<double> total_ang;
    vector<double> ang;
    const double energy = check_total_energy(pos_, vel_, mass_);
    tie(total_ang, ang) = check_angular_momentum(pos_, vel_, mass_);
    if (i < ndata)
    {
      energy_[i] = energy;
      copy_n(total_ang.data(), total_ang.size(), angular_momentum_.data() + i * total_ang.size());
      copy_n(ang.data(), ang.size(), angular_momentum_i_.data() + i * ang.size());
    }
    // output to a json file:
    ou << tag;
    prn_vector(ou, "positions", pos_);
    prn_vector(ou, "velocities", vel_);
    prn_vector(ou, "angular_momentum", total_ang);
    prn_vector(ou, "angular_momentum_i", ang);
    ou << "  \"energy\" : "
       << setw(15) << setprecision(8) << scientific << energy << "\n";
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
  input_pars_(input_pars0), nbody_(input_pars0["nbody"]), ndata_(input_pars0["ndata"]),
  c_alpha_(nbody_), pos_(3 * nbody_), vel_(3 * nbody_), mass_(nbody_),
  energy_(ndata_), angular_momentum_(3 * ndata_), angular_momentum_i_(3 * ndata_ * nbody_)
{
  nbody_simulator();
}
