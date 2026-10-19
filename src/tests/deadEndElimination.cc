// Standalone numerical regression: g++ -std=c++17 -O2 this-file.cc -o dee-test
#define CANDO_DEE_STANDALONE
#include "../chem/deadEndElimination.cc"
#include <cassert>
#include <iostream>
#include <random>

int main() {
  std::mt19937 rng(71051);
  std::uniform_int_distribution<int> value(-12,12);
  for (size_t trial = 0; trial < 200; ++trial) {
    size_t n = 1 + trial%4, width = 2 + trial%2, slots = n*(width+1);
    std::vector<double> unary(slots), pairs(slots*(slots-1));
    std::vector<uint8_t> inter(slots*(slots-1)/2);
    for (auto& x : unary) x = value(rng)*0.25;
    for (auto& x : pairs) x = value(rng)*0.25;
    for (auto& x : inter) x = rng()%2;
    chem::dee::Domains domains(n);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < width; ++j) domains[i].push_back(i*(width+1)+j);
    chem::dee::Problem p{domains,unary.data(),pairs.data(),inter.data(),0.5,1.5};
    double exact = chem::dee::infinity;
    chem::dee::State state(n), optimum;
    std::function<void(size_t)> enumerate = [&](size_t i) {
      if (i == n) {
        double e = p.energy(state);
        if (e < exact) { exact = e; optimum = state; }
      } else for (auto r : domains[i]) { state[i] = r; enumerate(i+1); }
    };
    enumerate(0);
    while (p.desmet()) {}
    p.goldstein(100);
    for (size_t i = 0; i < n; ++i)
      assert(std::find(p.domains[i].begin(),p.domains[i].end(),optimum[i]) != p.domains[i].end());
    auto incumbent = p.icm({optimum});
    assert(incumbent.second == exact);
    for (size_t entries : {size_t(0),size_t(100000)}) {
      auto result = chem::dee::search(p,100000,entries,{},chem::dee::infinity);
      assert(result.optimal && result.energy == exact);
      assert(result.lowerBound <= exact);
      assert(p.energy(result.state) == exact);
      auto stopped = chem::dee::search(p,0,entries,incumbent.first,incumbent.second);
      assert(!stopped.optimal && stopped.nodes == 0 && stopped.energy == exact);
      assert(stopped.lowerBound <= exact);
    }
  }
  // Ties retain every candidate, and seeded ICM follows first-domain tie breaks.
  double unary[4] = {}, pairs[12] = {};
  uint8_t inter[6] = {};
  chem::dee::Problem tied{{{0,1},{2,3}},unary,pairs,inter,1,1};
  assert(tied.desmet() == 0 && tied.goldstein(100) == 0);
  auto icm = tied.icm({{1,3}});
  assert((icm.first == chem::dee::State{0,2}));
  std::cout << "DEE numerical tests passed: 200 exhaustive cases, ties, sparse domains, caps and both bound modes\n";
}
