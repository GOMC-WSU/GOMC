#include <cstdio>
#include <cstdlib>

#include <iostream>

#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " num-particles rand-seed\n";
    return EXIT_FAILURE;
  }

  const size_t num_particles = atoi(argv[1]);
  
  const ulong rand_seed = strtoul(argv[2], nullptr, 10);
  clog << "Using random seed: " << rand_seed << "\n";
  srand(rand_seed);

  const double sigma = pow(2.0, 1.0 / 6.0);
  mat::fixed<3, 3> L = 40.0 * sigma * eye(3, 3);
  L.save("cell.mat", arma::csv_ascii);
  
  mat x = randu<mat>(3, num_particles);
  mat r = L * x;
  r.save("positions.mat", arma::csv_ascii);
  
  vec q = 2e1 * randu<vec>(num_particles);
  // vec q = ones<vec>(num_particles);
  for (size_t n = 0; n < num_particles; ++n) if (n % 3 != 0) q[n] = 0.0;
  q.save("charges.mat", arma::csv_ascii);

  return EXIT_SUCCESS;
}
