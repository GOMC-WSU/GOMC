#include <cassert>
#include <cstdlib>
#include <cstring>

#include <limits>
#include <numeric>
#include <functional>
#include <iostream>

#include <fftw3.h>

#include "cube.cpp"
#include "pme.h"
#include "spline.h"

#include "linear-algebra.h"

namespace ewald {

const spline4 spline;
const int order = spline.order;

const double epsilon = std::numeric_limits<double>::epsilon();

pme::pme(long d0, long d1, long d2,
         double const unit_cell[],
         size_t N_,
         double const* q_,
         double cut_off_,
         double tolerance,
         int num_threads)
    : N(N_),
      q(q_),
      cut_off(cut_off_),
      threshold(compute_threshold(cut_off, tolerance)),
      threshold2(threshold * threshold),
      dim0(d0), dim1(d1), dim2(d2),
      Q(d0, d1, d2),
      Qhat(d0, d1, d2),
      BC(d0, d1, d2),
      s(3, N, order),
      ds(3, N, order)
{
  assert(dim0 > 0 && dim1 > 0 && dim2 > 0);

  std::memcpy(L, unit_cell, 3 * 3 * sizeof(double));
  volume = determinant(L);
  assert(volume > 1e2 * epsilon);
  matrix_inverse(L, Linv);

  std::memset(K, 0, 3 * 3 * sizeof(double));
  K[3*0+0] = dim0;
  K[3*1+1] = dim1;
  K[3*2+2] = dim2;

  matrix_matrix_product(K, Linv, H);

  initialize_table();

  fftw_init_threads();
  fftw_plan_with_nthreads(num_threads);

  {
    // TODO: Benchmark r2c.
    fftw_complex* in  = reinterpret_cast<fftw_complex*>(Q.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(Qhat.memptr());
    forward_plan = fftw_plan_dft_3d(dim0, dim1, dim2,
                                    in, out,
                                    FFTW_FORWARD, FFTW_MEASURE);
  }

  {
    fftw_complex* in  = reinterpret_cast<fftw_complex*>(Qhat.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(Qhat.memptr());
    backward_plan = fftw_plan_dft_3d(dim0, dim1, dim2,
                                     in, out,
                                     FFTW_BACKWARD, FFTW_MEASURE);
  }
}

pme::~pme() {
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(backward_plan);

  fftw_cleanup_threads();
}

void pme::initialize_table() {
  long k0, k1, k2;
  for (k0 = 0; k0 < dim0; ++k0) {
    const double h0 = k0 < dim0 / 2 ? k0 : k0 - dim0;

    for (k1 = 0; k1 < dim1; ++k1) {
      const double h1 = k1 < dim1 / 2 ? k1 : k1 - dim1;

      for (k2 = 0; k2 < dim2; ++k2) {

        if (k0 == 0 && k1 == 0 && k2 == 0)
          continue;

        const double h2 = k2 < dim2 / 2 ? k2 : k2 - dim2;

        double LinvTm[] = { h0, h1, h2 };
        matrix_vector_product_transpose(Linv, LinvTm);

        const double hh = dot_product(LinvTm, LinvTm);

        const double nrm =  norm(spline.factor(k0, dim0)
                                 * spline.factor(k1, dim1)
                                 * spline.factor(k2, dim2));

        BC(k0, k1, k2) = psi(hh) * nrm;
      }
    }
  }
}

static void compute_spline(double u,
                           double* __restrict__ s,
                           double* __restrict__ ds,
                           double& floor_u) {
  floor_u = floor(u);
  double x = u - floor_u + double(order) - 1.0;

  // Compute splines and their derivatives using Horner's method.
  const double a[] = {
    32.0 / 3.0, -8.0, 2.0, -1.0 / 6.0,
    -22.0 / 3.0, 10.0, -4.0, 1.0 / 2.0,
    2.0 / 3.0, -2.0, 2.0, -1.0 / 2.0,
    0.0, 0.0, 0.0, 1.0 / 6.0
  };

  for (size_t k = 0; k < 4; ++k, --x) {
     s[k] = a[4*k] + (a[4*k+1] + (a[4*k+2] + a[4*k+3] * x) * x) * x;
    ds[k] = a[4*k+1] + (2.0 * a[4*k+2] + 3.0 * a[4*k+3] * x) * x;
  }
}

void pme::interpolate(double const* __restrict__ x,
                      double const* __restrict__ y,
                      double const* __restrict__ z) {
  for (size_t n = 0; n < N; ++n) {
    double u[] = { x[n], y[n], z[n] };
    matrix_vector_product(H, u);

    const double q_n = q[n];

    double floor_u[3];
    for (size_t d = 0; d < 3; ++d)
      compute_spline(u[d], &s(d, n, 0), &ds(d, n, 0), floor_u[d]);

    for (long c0 = 0; c0 < order; ++c0) {
      const long l0 = floor_u[0] + c0 - order + 1;
      const long k0 = l0 < 0 ? l0 + dim0 : l0;

      const double q_n_s0 = q_n * s(0, n, c0);

      for (long c1 = 0; c1 < order; ++c1) {
        const long l1 = floor_u[1] + c1 - order + 1;
        const long k1 = l1 < 0 ? l1 + dim1 : l1;

        const double q_n_s0_s1 = q_n_s0 * s(1, n, c1);

        for (long c2 = 0; c2 < order; ++c2) {
          const long l2 = floor_u[2] + c2 - order + 1;
          const long k2 = l2 < 0 ? l2 + dim2 : l2;

          Q(k0, k1, k2) += q_n_s0_s1 * s(2, n, c2);
          // printf("cube (%s): %d %d %d %d %g\n", __func__, n, k0, k1, k2, std::real(Qhat(k0, k1, k2)));
        }
      }
    }
  }
}

double pme::energy(double const* __restrict__ x,
                   double const* __restrict__ y,
                   double const* __restrict__ z,
                   double* __restrict__ force_x,
                   double* __restrict__ force_y,
                   double* __restrict__ force_z) {
  Q.zeros();
  Qhat.zeros();
  s.zeros();
  ds.zeros();

  interpolate(x, y, z);
  /*
  {
    std::complex<double> const* Qptr = Q.memptr();
    std::complex<double> total { 0.0, 0.0 };
    for (long r = 0; r < dim0 * dim1 * dim2; ++r)
      total += Qptr[r];
    std::printf("%s: Interpolated a total \"charge\" of %g + %g i\n", __func__, real(total), imag(total));
  }
  */


  fftw_execute(forward_plan);

  Qhat *= BC;

  fftw_execute(backward_plan);

  double Erecip = 0.0;
  for (long k0 = 0; k0 < dim0; ++k0)
    for (long k1 = 0; k1 < dim1; ++k1)
      for (long k2 = 0; k2 < dim2; ++k2)
        Erecip += real(Q(k0, k1, k2)) * real(Qhat(k0, k1, k2));

  double Erecip1 = 0.0;
  for (size_t n = 0; n < N; ++n) {
    double u[] = { x[n], y[n], z[n] };
    matrix_vector_product(H, u);

    for (long c0 = 0; c0 < order; ++c0) {
      const long l0 = floor(u[0]) + c0 - order + 1;
      const long k0 = l0 < 0 ? l0 + dim0 : l0;

      for (long c1 = 0; c1 < order; ++c1) {
        const long l1 = floor(u[1]) + c1 - order + 1;
        const long k1 = l1 < 0 ? l1 + dim1 : l1;

        const double ds0_s1 = ds(0, n, c0) *  s(1, n, c1);
        const double s0_ds1 =  s(0, n, c0) * ds(1, n, c1);
        const double s0_s1  =  s(0, n, c0) *  s(1, n, c1);

        for (long c2 = 0; c2 < order; ++c2) {
          const long l2 = floor(u[2]) + c2 - order + 1;
          const long k2 = l2 < 0 ? l2 + dim2 : l2;

          const double Q_k = std::real(Qhat(k0, k1, k2));

          Erecip1 += Q_k * std::real(Q(k0, k1, k2));

          double w[] = {
            ds0_s1 *  s(2, n, c2),
            s0_ds1 *  s(2, n, c2),
             s0_s1 * ds(2, n, c2)
          };
          matrix_vector_product_transpose(H, w);

          // if (n == 1)
          //   printf("cube (%s): %3d\t%3d\t%3d\t%3d\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.9f\n",
          //          __func__, n, k0, k1, k2, 
          //          std::real(Q(k0, k1, k2)), std::imag(Q(k0, k1, k2)),
          //          std::real(Qhat(k0, k1, k2)), std::imag(Qhat(k0, k1, k2)),
          //          Q_k * std::real(Q(k0, k1, k2)));

          const double factor = -2.0 * Q_k * q[n];
          force_x[n] += factor * w[0];
          force_y[n] += factor * w[1];
          force_z[n] += factor * w[2];
        }
      }
    }
  }

  //printf("%s: Erecip = %3.16f vs. %3.16f\n", __func__, Erecip, Erecip1);
  return double(Erecip1);
}

double pme::energy_extra() const {
  const double sum_q = std::accumulate(q, q + N, 0.0, std::plus<double>());

  return 1.0 / 6.0 * pow(M_PI * threshold, 3.0) / volume * sum_q * sum_q;
}

double pme::energy_self() const {
  const double q2 = std::inner_product(q, q + N, q, 0.0);

  return 1.0 / 12.0 * pow(M_PI, 3.0) * pow(threshold, 6.0) * q2;
}

double pme::psi(double h2) const {
  const double b2 = M_PI * h2 / threshold2;
  const double b = sqrt(b2);
  const double b3 = b2 * b;

  const double h = sqrt(h2);
  const double h3 = h2 * h;

  return pow(M_PI, 9.0 / 2.0) / (3.0 * volume) * h3
      * (sqrt(M_PI) * erfc(b) + (1.0 / (2.0 * b3) - 1.0 / b) * exp(-b2));
}

double pme::recip_convergence_term(double h2,
                                   std::complex<double> const& rho) const {
  const double norm_rho2 = norm(rho);
  return norm_rho2 * psi(h2);
}

static inline double direct_term(double r2, double threshold2, double& fac) {
  const double r6 = r2 * r2 * r2;
  const double r8 = r6 * r2;
  const double a2 = M_PI * threshold2 * r2;
  const double a4 = a2 * a2;
  const double a6 = a4 * a2;
  const double exp_m_a2 = exp(-a2);

  fac = -(a6 + 3.0 * a4 + 6.0 * a2 + 6.0) * exp_m_a2 / r8;

  return (1.0 + a2 + 0.5 * a4) * exp_m_a2 / r6;
}

double pme::direct_convergence_term(double const* v, double* dv) const {
  const double r2 = dot_product(v, v);
  double fac;

  const double energy = direct_term(r2, threshold2, fac);

  for (size_t k = 0; k < 3; ++k)
    dv[k] = fac * v[k];

  return energy;
}

static double diff_direct_term(double threshold, double cut_off) {
  // Computes the derivative of direct_term() with respect to the
  // variable threshold. Note that the input is threshold, not
  // threshold squared (threshold2).

  const double threshold2 = threshold * threshold;
  const double threshold4 = threshold2 * threshold2;
  const double threshold5 = threshold4 * threshold;
  const double cut_off2 = cut_off * cut_off;

  return -threshold5 * exp(-threshold2 * cut_off2);
}

/** Obtains a suitable value of the threshold parameter. The threshold
 * parameter (beta) is obtained by solving the non-linear equation
 * direct_term(beta) = tol for a given tolerance tol.
 *
 * \param cut_off Cut off distance.
 *
 * \param tol Maximum value of the real space contribution for the
 * distance given by cut_off.
 *
 * \return The solution beta of the non-linear equation.
 */
double pme::compute_threshold(double cut_off, double tol) const {
  // XXX Maybe throw exception rather than return an integer status flag.
  const double r2 = cut_off * cut_off;
  double aux;

  // Our initial guess comes from imposing that the value of exp(-a2)
  // in direct_term should be greater than or equal to machine
  // epsilon.
  double beta_init = sqrt(-log(epsilon) / M_PI) / cut_off;
  double beta = beta_init;
  double beta2 = beta * beta;
  double val = direct_term(r2, beta2, aux);

  //std::clog << "Starting non-linear solver with beta = " << beta << std::endl;

  size_t iter;
  const size_t max_iters = size_t(1e6);
  for (iter = 0; iter < max_iters; ++iter) {
    beta -= (direct_term(r2, beta2, aux) - tol) / diff_direct_term(beta, cut_off);
    beta2 = beta * beta;
    if (std::isinf(beta) || std::isnan(beta)) {
      beta_init /= 2.0;
      beta = beta_init;

      if (beta_init == 0.0) {
        std::clog << "Non-linear solver failed to find a suitable value of beta."
                  << std::endl;
        // Red alert. XXX This should fail more gracefully.
        beta = NAN;
        break;
      }
    }

    val = direct_term(r2, beta2, aux);
    if (fabs(val - tol) < epsilon)
      break;
  }

  //std::clog << "Solver stopped after " << iter << " iterations.\n"
  //          << "Using beta = " << beta << " with value = " << val << "\n";

  /*if (iter == max_iters)
    std::clog << "Warning: Non-linear solver did not converge."
              << std::endl;*/

  return beta;
}

double pme::excluded(exclusion_vector const& exclusions,
                     double const* __restrict__ x,
                     double const* __restrict__ y,
                     double const* __restrict__ z,
                     double* __restrict__ force_x,
                     double* __restrict__ force_y,
                     double* __restrict__ force_z) const {
  double excluded_energy = 0.0;

  for (size_t i = 0; i < N; ++i) {
    const double q_i = q[i];

    index_vector const& excl = exclusions[i];
    for (size_t r = 0; r < excl.size(); ++r) {
      const size_t j = excl[r];
      if (j <= i)
        continue;

      const double v[] = { x[i] - x[j], y[i] - y[j], z[i] - z[j] };
      const double r2 = dot_product(v, v);
      // Note that we must NOT enforce a cut off here.  We are
      // compensating a masked contribution to the reciprocal part,
      // and that contribution does not involve cut offs.
      const double r6 = r2 * r2 * r2;
      const double invr6 = 1.0 / (r2 * r2 * r2);
      const double invr8 = 1.0 / (r6 * r2);

      double term;
      double f[3];
      const double qi_qj = q_i * q[j];
      excluded_energy += qi_qj * (invr6 - direct_term(r2, threshold2, term));
      const double factor = qi_qj * (6.0 * invr8 + term);
      for (size_t k = 0; k < 3; ++k)
        f[k] = factor * v[k];
      force_x[i] += f[0];
      force_y[i] += f[1];
      force_z[i] += f[2];
      force_x[j] -= f[0];
      force_y[j] -= f[1];
      force_z[j] -= f[2];
    }
  }

  return excluded_energy;
}

}
