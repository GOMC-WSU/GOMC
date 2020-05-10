#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cfenv>

#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

#include <armadillo>

#ifdef HAVE_OMP
#include <omp.h>
#endif

#include "pme.h"
#include "abort_unless.h"
#include "energies.h"

using namespace std;
using namespace arma;

const complex<double> I(0.0, 1.0);

ewald::exclusion_vector exclusions;

static bool excluded(size_t i, size_t j);

int main(int argc, char* argv[]) {
  if (argc < 8) {
    cout << "Usage: " << argv[0]
         << " dimension cell-file positions-file charges-file cut-off tolerance num-shells\n";
    return EXIT_FAILURE;
  }

  const long dim0 = atoi(argv[1]);
  const long dim1 = dim0;
  const long dim2 = dim0;

  mat::fixed<3, 3> L;
  L.load(argv[2]);
  const mat::fixed<3, 3> LinvT = trans(inv(L));
  L.print();

  mat r;
  r.load(argv[3]);
  abort_unless(r.n_rows == 3);
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
  const double cut_off2 = cut_off * cut_off;
  const double tolerance = atof(argv[6]);
  clog << "Real space cut-off: " << cut_off << "\n"
       << "Real space tolerance: " << tolerance << "\n\n";


  exclusions = ewald::exclusion_vector(num_particles);
  {                                     // Populate exclusion vector.
/*
    for (size_t i = 0; i < num_particles; ++i)
      for (size_t j = 0; j < num_particles; ++j)
        // if (i != j && ((i % 2) != (j % 2))) {
        if (abs(int(i) - int(j)) <= 2) {
          exclusions[i].push_back(j);
        }
*/
    /*
    for (size_t i = 0; i < num_particles; ++i) {
      std::cout << i << ": ";
      ewald::index_vector const& is = exclusions[i];
      for (size_t j = 0; j < is.size(); ++j)
        std::cout << is[j] << " ";
      std::cout << std::endl;
    }
    */
  }

  const double max_shell = int(atof(argv[7]));

  mat forces = zeros<mat>(3, num_particles);
  mat force_cutoff = zeros<mat>(3, num_particles);
  mat force_direct = zeros<mat>(3, num_particles);
  mat force_recip = zeros<mat>(3, num_particles);
  mat force_excluded = zeros<mat>(3, num_particles);
  mat force_excluded_pme = zeros<mat>(3, num_particles);

  std::vector<double> force_x(num_particles, 0.0);
  std::vector<double> force_y(num_particles, 0.0);
  std::vector<double> force_z(num_particles, 0.0);
  mat force_recip_pme = zeros<mat>(3, num_particles);
  double LL[3 * 3];
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      LL[3*i+j] = L(i, j);

  energies energy_pme, energy_ewald;

  wall_clock timer;
  timer.tic();

#ifdef HAVE_OMP
  ewald::pme p(dim0, dim1, dim2, LL, num_particles, &q[0], cut_off, tolerance, omp_get_max_threads());
#else
  ewald::pme p(dim0, dim1, dim2, LL, num_particles, &q[0], cut_off, tolerance, 1);
#endif
  double seconds = timer.toc();

  clog << "PME initialization took " << seconds << " seconds.\n\n";

  const double t = p.get_threshold();
  if (std::isnan(t)) {
    cerr << "Unable to run PME for the prescribed values of cut off and tolerance.\n";
    return EXIT_FAILURE;
  }

  timer.tic();

  energy_pme.recip = p.energy(&coor_x[0], &coor_y[0], &coor_z[0],
                              &force_x[0], &force_y[0], &force_z[0]);

  for (size_t k = 0; k < num_particles; ++k) {
    force_recip_pme(0, k) = force_x[k];
    force_recip_pme(1, k) = force_y[k];
    force_recip_pme(2, k) = force_z[k];
  }

  seconds = timer.toc();
  clog << "PME took " << seconds << " seconds.\n\n";

  std::vector<double> force_excl_x(num_particles, 0.0);
  std::vector<double> force_excl_y(num_particles, 0.0);
  std::vector<double> force_excl_z(num_particles, 0.0);

  energy_pme.masked = p.excluded(exclusions,
                                 &coor_x[0],
                                 &coor_y[0],
                                 &coor_z[0],
                                 &force_excl_x[0],
                                 &force_excl_y[0],
                                 &force_excl_z[0]);

  for (size_t k = 0; k < num_particles; ++k) {
    force_excluded_pme(0, k) = force_excl_x[k];
    force_excluded_pme(1, k) = force_excl_y[k];
    force_excluded_pme(2, k) = force_excl_z[k];
  }

  energy_ewald.extra = energy_pme.extra = p.energy_extra();
  energy_ewald.self = energy_pme.self = p.energy_self();

  double Eexhaustive = 0.0;
  double Ecutoff = 0.0;

  mat force_ewald, force_pme;

  double elapsed = 0.0;

  for (int n = 0; n <= max_shell; n++) {
    wall_clock timer;
    timer.tic();
    for (int nx = -n; nx <= n; nx++) {
      for (int ny = -n; ny <= n; ny++) {
        for (int nz = -n; nz <= n; nz++) {
          if (abs(nx) != n && abs(ny) != n && abs(nz) != n)
            continue;

          const bool origin = (n == 0);

          const vec::fixed<3> m = { double(nx), double(ny), double(nz) };
          const vec::fixed<3> Lm = L * m;
          const vec::fixed<3> LinvTm = LinvT * m;

          // Compute direct space contribution.
          for (size_t i = 0; i < num_particles; i++) {
            for (size_t j = 0; j < num_particles; j++) {
              if (origin && i == j)
                continue;

              const vec::fixed<3> v = r.col(i) - r.col(j) + Lm;
              const double r2 = dot(v, v);
              const double r6 = r2 * r2 * r2;
              const double r8 = r6 * r2;

              const double qij = q[i] * q[j];

              const bool is_excluded = excluded(i, j);

              const double ener = 0.5 * qij / r6;
              const vec::fixed<3> f = 6.0 * qij / r8 * v;

              if (origin) {
                if (!is_excluded) {     // Unmasked pair.
                  Eexhaustive  += ener;
                  forces.col(i) += f;

                  if (r2 < cut_off2) {
                    Ecutoff += ener;
                    force_cutoff.col(i) += f;

                    vec::fixed<3> dv;
                    const double pot = 0.5 * qij
                        * p.direct_convergence_term(v.memptr(), dv.memptr());
                    energy_pme.direct += pot;
                    energy_ewald.direct += pot;
                    force_direct.col(i) += -qij * dv;
                  }
                } else {                // Masked pair.
                  // Note that we don't do cut offs here because we
                  // want to cancel the Fourier space contribution,
                  // which does not involve cut offs.
                  vec::fixed<3> dv;
                  energy_ewald.masked += (ener - 0.5 * qij
                                          * p.direct_convergence_term(v.memptr(),
                                                                      dv.memptr()));
                  force_excluded.col(i) += f + qij * dv;
                }
              } else {                  // Outer shells.
                Eexhaustive  += ener;
                forces.col(i) += f;

                if (r2 < cut_off2) {
                  Ecutoff += ener;
                  force_cutoff.col(i) += f;

                  vec::fixed<3> dv;
                  const double pot = 0.5 * qij
                      * p.direct_convergence_term(v.memptr(), dv.memptr());
                  energy_pme.direct += pot;
                  energy_ewald.direct += pot;
                  force_direct.col(i) += -qij * dv;
                }
              }
            }
          }

          // Compute contribution from Fourier space.
          if (!origin) {
            complex<double> rho(0.0, 0.0);

            for (size_t i = 0; i < num_particles; i++)
              rho += q[i] * exp(2.0 * M_PI * I * dot(m, x.col(i)));

            const double h2 = dot(LinvTm, LinvTm);
            energy_ewald.recip += p.recip_convergence_term(h2, rho);

            // force_recip.col(0) = -4.0 * M_PI * q[0] * q[1] * p.psi(h2) * sin(2.0 * M_PI * dot(LinvTm, r.col(0) - r.col(1))) * (-m);
            // force_recip.col(1) = -4.0 * M_PI * q[0] * q[1] * p.psi(h2) * sin(2.0 * M_PI * dot(LinvTm, r.col(0) - r.col(1))) * m;
          }
        }
      }
    }

    elapsed += timer.toc();

    force_ewald = force_direct + force_recip;
    force_pme   = force_direct + force_recip_pme - force_excluded_pme;

    {
      const double pme_total   = energy_pme.total();
      const double ewald_total = energy_ewald.total();

      const double rel_error_ewald_vs_pme    = fabsl(ewald_total - pme_total) / fabsl(ewald_total);
      const double rel_error_ewald_vs_cutoff = fabsl(ewald_total - Ecutoff) / fabsl(ewald_total);

      const double rel_error_exhaustive_vs_cutoff = fabsl(Eexhaustive - Ecutoff) / fabsl(Eexhaustive);
      const double rel_error_exhaustive_vs_ewald  = fabsl(Eexhaustive - ewald_total) / fabsl(Eexhaustive);
      const double rel_error_exhaustive_vs_pme    = fabsl(Eexhaustive - pme_total) / fabsl(Eexhaustive);

      cout << "==============================================\n"
           << fixed << setprecision(4)
           << "Shell " << n << " (" << elapsed << " seconds)" << "\n"
           << "==============================================\n"
           << scientific << setprecision(9)
           << "\n"
           << "Eexhaustive = " << Eexhaustive << "\n"
           << "Ecutoff     = " << Ecutoff << "\n"
           << "Ewald       = " << energy_ewald << "\n"
           << "PME         = " << energy_pme << "\n"
           << "\n"
           << "Abs. error in excluded forces = " << norm(force_excluded - force_excluded_pme, "inf") << "\n"
           << "\n"
           << "Abs. error in forces (Exhaustive vs. cutoff) = " << norm(forces - force_cutoff, "inf") << "\n"
           << "Abs. error in forces (Exhaustive vs. PME)    = " << norm(forces - force_pme, "inf") << "\n"
           << "\n"
           << "Rel. error in total energies (Ewald vs. cutoff) = " << rel_error_ewald_vs_cutoff << "\n"
           << "Rel. error in total energies (Ewald vs. PME)    = " << rel_error_ewald_vs_pme << "\n"
           << "\n"
           << "Rel. error in total energies (Exhaustive vs. cutoff) = " << rel_error_exhaustive_vs_cutoff << "\n"
           << "Rel. error in total energies (Exhaustive vs. Ewald)  = " << rel_error_exhaustive_vs_ewald << "\n"
           << "Rel. error in total energies (Exhaustive vs. PME)    = " << rel_error_exhaustive_vs_pme << "\n"
           << endl;
    }
  }

  return 0;
}

bool excluded(size_t i, size_t j) {
  // assert(i < num_particles);
  // assert(j < num_particles);

  ewald::index_vector const& is = exclusions[i];
  for (size_t k = 0; k < is.size(); ++k)
    if (j == is[k])
      return true;

  return false;
}
