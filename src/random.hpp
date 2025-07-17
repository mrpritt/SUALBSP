#pragma once

#include <numeric>
#include <random>
#include <vector>

namespace sualbsp {
extern std::mt19937 rng;

unsigned setupRandom(unsigned seed = 0);
int getRandom(int min, int max);
double getRandom();
double getRandom(double min, double max);
} // namespace sualbsp
