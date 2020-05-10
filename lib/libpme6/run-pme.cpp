#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cfenv>

#include <complex>
#include <iomanip>
#include <iostream>

#include <armadillo>

#include <omp.h>

#include "pme.h"
#include "abort_unless.h"

using namespace std;
using namespace arma;

const complex<double> I(0.0, 1.0);

const double sigma = pow(2.0, 1.0 / 6.0);

int main(int argc, char* argv[]) {
  if (argc < 7) {
    cout << "Usage: " << argv[0]
         << " dimension cell-file positions-file charges-file cut-off tolerance\n";
    return EXIT_FAILURE;
  }

  const long dim0 = atoi(argv[1]);
  const long dim1 = dim0;
  const long dim2 = dim0;

  mat::fixed<3, 3> L;
  L.load(argv[2]);
  const mat::fixed<3, 3> LinvT = trans(inv(L));
  // const double V = det(L);
  L.print("L = ");

  mat r;
  r.load(argv[3]);
  abort_unless(r.n_rows == 3);
  // r.print("r = ");
  const mat x = inv(L) * r;

  vec qq;
  qq.load(argv[4]);
  std::vector<double> q(r.n_cols);
  {
    for (size_t k = 0; k < q.size(); ++k)
      q[k] = qq(k);
  }

  const size_t num_particles = q.size();
  clog << "Read " << num_particles << " particles into memory.\n";

  std::vector<double> coor_x(num_particles);
  std::vector<double> coor_y(num_particles);
  std::vector<double> coor_z(num_particles);
  {
    for (size_t k = 0; k < num_particles; ++k) {
      coor_x[k] = r(0, k);
      coor_y[k] = r(1, k); 
      coor_z[k] = r(2, k); 
    }
  }

  const double cut_off = atof(argv[5]);
  const double tolerance = atof(argv[6]);
  clog << "Real space cut-off: " << cut_off << "\n"
       << "Real space tolerance: " << tolerance << "\n\n";

  mat forces = zeros<mat>(3, num_particles);
  mat forces_cutoff = zeros<mat>(3, num_particles);
  mat forces_direct = zeros<mat>(3, num_particles);
  mat forces_recip = zeros<mat>(3, num_particles);

  long double ErecipPME = 0.0L;
  std::vector<double> forces_x(num_particles);
  std::vector<double> forces_y(num_particles);
  std::vector<double> forces_z(num_particles);
  mat forcesPME = zeros<mat>(3, num_particles);
    wall_clock timer;
    timer.tic();

    double LL[3 * 3];
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        LL[3*i+j] = L(i, j);

    ewald::pme p(dim0, dim1, dim2, LL, 
                 num_particles, &q[0], 
                 cut_off, tolerance, 
                 omp_get_max_threads());

    double seconds = timer.toc();
    clog << "PME initialization took " << seconds << " seconds.\n\n";

    timer.tic();

    ErecipPME = p.energy(&coor_x[0], &coor_y[0], &coor_z[0],
                         &forces_x[0], &forces_y[0], &forces_z[0]);

    for (size_t k = 0; k < num_particles; ++k) {
      forcesPME(0, k) = forces_x[k];
      forcesPME(1, k) = forces_y[k];
      forcesPME(2, k) = forces_z[k];
    }

    seconds = timer.toc();
    clog << "PME took " << seconds << " seconds.\n\n";

  long double Ebrute = 0.0L, Edirect = 0.0L;
  long double Ecutoff = 0.0L;
  const long double Eextra = p.energy_extra();
  const long double Eself  = p.energy_self();

  mat forces_ewald, forces_pme;
  
  double elapsed = 0.0;
  for (int n = 0; n < 1; n++) {
    wall_clock timer;
    timer.tic();
    for (int nx = -n; nx <= n; nx++) {
      for (int ny = -n; ny <= n; ny++) {
        for (int nz = -n; nz <= n; nz++) {
          if (abs(nx) != n && abs(ny) != n && abs(nz) != n)
            continue;

          const bool origin = (nx == 0 && ny == 0 && nz == 0);
          
          const vec::fixed<3> m = { double(nx), double(ny), double(nz) };
          const vec::fixed<3> Lm = L * m;
          const vec::fixed<3> LinvTm = LinvT * m;
          
          for (size_t i = 0; i < num_particles; i++) {
            for (size_t j = 0; j < num_particles; j++) {
              if (origin && i == j)
                continue;

              const vec::fixed<3> v = r.col(i) - r.col(j) + Lm;
              const double r2 = dot(v, v);
              const double r6 = r2 * r2 * r2;
              const double r8 = r6 * r2;

              assert(r2 >= 1e-16);

              Ebrute  += 0.5 * q[i] * q[j] / r6;
              forces.col(i) -= -6.0 * q[i] * q[j] / r8 * v;

              if (abs(nx) < 2 && abs(ny) < 2 && abs(nz) < 2) {
                forces_cutoff.col(i) -= -6.0 * q[i] * q[j] / r8 * v;
                Ecutoff = Ebrute;
              }

              vec::fixed<3> dv;
              Edirect += 0.5 * q[i] * q[j] * p.direct_convergence_term(v.memptr(), dv.memptr());
              forces_direct.col(i) -= q[i] * q[j] * dv;
            }
          }
        }
      }
    }

    elapsed += timer.toc();

    const long double EtotalPME = Edirect + ErecipPME + Eextra - Eself;

    forces_pme   = forces_direct + forcesPME;

    cout << "Real space computation over one shell took " << elapsed << " seconds." << "\n"
         << fixed << setprecision(16)
         << "Ecutoff      = " << Ecutoff << "\n"
         << "Etotal (PME) = " << EtotalPME  << "\n"
         << "Edirect + Erecip + Eextra - Eself = " 
         << Edirect << " + " << ErecipPME << " + " << Eextra << " - " << Eself << "\n";
  }
  
  return 0;
}
