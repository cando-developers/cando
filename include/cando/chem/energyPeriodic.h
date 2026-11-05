#ifndef CANDO_CHEM_ENERGY_PERIODIC_H
#define CANDO_CHEM_ENERGY_PERIODIC_H

#include <cmath>

namespace chem {
// Generated energy expressions use Mathematica Mod: the remainder has the
// divisor's sign. fmod alone gives the dividend's sign and breaks periodic
// dihedral differences when angle - target + pi is negative.
inline double energy_mod(double value, double period) {
  double result = std::fmod(value, period);
  if (result != 0.0 && std::signbit(result) != std::signbit(period)) {
    result += period;
    // Rounding can land exactly on the excluded upper endpoint.
    if (result == period) result = 0.0;
  }
  return result;
}
}

#endif
