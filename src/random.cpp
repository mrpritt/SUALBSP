#include "random.hpp"

#include <fstream>
using namespace std;

namespace sualbsp {
mt19937 rng;

unsigned setupRandom(unsigned seed) {
  if (seed == 0) {
    seed = time(0);
    ifstream f("/dev/urandom");
    if (f.good()) {
      f.read((char *)(&seed), sizeof(unsigned int));
    }
  }
  rng.seed(seed);
  srand48(seed);
  return seed;
}

int getRandom(int min, int max) {
  uniform_int_distribution<> U(min, max);
  return U(rng);
}

double getRandom() {
  uniform_real_distribution<> U;
  return U(rng);
}
double getRandom(double min, double max) {
  uniform_real_distribution<> U(min, max);
  return U(rng);
}
} // namespace sualbsp
