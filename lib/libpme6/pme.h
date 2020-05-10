#ifndef PME_H
#define PME_H

#include <vector>
#include <complex>

#include <fftw3.h>

#include "cube.h"

namespace ewald {

typedef std::vector<size_t> index_vector;
typedef std::vector<index_vector> exclusion_vector;

class pme {
 public:
  pme(long dimension0, long dimension1, long dimension2,
      double const* unit_cell,
      size_t num_particles,
      double const* q,
      double cut_off, double tolerance,
      int num_threads);

  ~pme();

  double get_threshold() const { return threshold; }

  double energy(double const* x, double const* y, double const* z,
                double* force_x, double* force_y, double* force_z);

  double energy_extra() const;

  double energy_self() const;

  double direct_convergence_term(double const* v,
                                 double* dv) const;

  // XXX Make the following two methods private at some point.
  double recip_convergence_term(double h2,
                                std::complex<double> const& rho) const;

  double psi(double h2) const;

  double excluded(exclusion_vector const& exclusions,
                  double const* x, double const* y, double const* z,
                  double* force_x, double* force_y, double* force_z) const;

 private:
  double compute_threshold(double cut_off, double tol) const;

  void initialize_table();

  void interpolate(double const* x, double const* y, double const* z);

 private:
  const size_t N;                       // Number of particles.

  double L[3 * 3];
  double Linv[3 * 3];
  double volume;

  double const* q;                      // B params. in OPLS LJ model.

  const double cut_off;

  const double threshold;               // Real-space threshold.
  const double threshold2;

  double K[3 * 3];
  double H[3 * 3];                      // Coordinate transformation.

  const long dim0, dim1, dim2;

  complex_cube Q;
  complex_cube Qhat;
  complex_cube BC;

  real_cube s;
  real_cube ds;

  fftw_plan forward_plan, backward_plan;
};

}

#endif
