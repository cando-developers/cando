// Standalone regression (from the Cando extension directory):
// c++ -std=c++17 -Iinclude src/lisp/regression-tests/energy-periodic.cc -o /tmp/energy-periodic
#include <cando/chem/energyPeriodic.h>
#include <cassert>
#include <iostream>

int main() {
  const double pi = std::acos(-1.0);
  const double period = 2.0 * pi;
  // Match the generated dihedral energy's current Pi constant and formula.
  auto delta = [period](double angle, double target) {
    return -3.14159265 + chem::energy_mod(3.14159265 + angle - target, period);
  };
  auto energy = [&delta](double angle, double target) {
    const double d = delta(angle, target);
    return 100.0 * d * d;
  };
  assert(chem::energy_mod(-1.0, 4.0) == 3.0);
  assert(chem::energy_mod(1.0, -4.0) == -3.0);
  assert(chem::energy_mod(-8.0, 4.0) == 0.0);
  // REAL is float: storing +pi rounds it above double-precision pi.
  const double stored_pi = static_cast<float>(pi);
  const double rewritten_pi = static_cast<float>(stored_pi - period);
  assert(stored_pi > pi);
  for (double angle : {-3.0, -2.0, -0.5, 0.5, 2.0, 3.0}) {
    assert(std::abs(delta(angle, stored_pi)) <= pi + 1e-7);
    assert(std::abs(energy(angle, stored_pi) - energy(angle, rewritten_pi)) < 1e-3);
    for (int turns = -3; turns <= 3; ++turns) {
      assert(std::abs(energy(angle, stored_pi) - energy(angle, stored_pi + turns * period)) < 1e-9);
    }
    const double h = 1e-5;
    auto slope = [&](double target) {
      return (energy(angle + h, target) - energy(angle - h, target)) / (2.0 * h);
    };
    assert(std::abs(slope(stored_pi) - slope(rewritten_pi)) < 1e-3);
  }
  // The old fmod formula has the erroneous extra turn for this case.
  const double old_delta = -3.14159265 + std::fmod(3.14159265 - 0.5 - stored_pi, period);
  assert(std::abs(old_delta) > pi);
  std::cout << "ENERGY-PERIODIC-REGRESSION-PASSED\n";
}
