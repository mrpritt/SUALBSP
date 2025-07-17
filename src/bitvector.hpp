
#pragma once

#include <vector>

#include "dimensions.hpp"

struct BitVector : public std::vector<char> {
  BitVector() {}
  BitVector(size_t n) : std::vector<char>(n) {}
  BitVector(size_t n, bool v) : std::vector<char>(n, v) {}
  void set() { fill(begin(), end(), true); }
  void reset() { fill(begin(), end(), false); }
  size_t count() {
    size_t result = 0;
    for (int i = size() - 1; i >= 0; i--)
      if (operator[](i))
        result++;
    return result;
  }
};

typedef BitVector BitRow;
