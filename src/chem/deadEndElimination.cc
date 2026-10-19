// Dead-end elimination and bounded GMEC search over cached pairwise energies.
// Define CANDO_DEE_STANDALONE only for the dependency-free numerical tests.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
#ifndef CANDO_DEE_STANDALONE
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/array_double.h>
#include <clasp/core/array_int8.h>
#include <clasp/core/array_int32.h>
#include <clasp/core/cons.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/chemPackage.h>
#endif

namespace chem {
namespace dee {
using State = std::vector<uint32_t>;
using Domains = std::vector<State>;
constexpr double infinity = std::numeric_limits<double>::max();

// Borrow the Lisp arrays for the duration of one synchronous call. No energy
// table copies and no Lisp callbacks in the numerical loops.
struct Problem {
  Domains domains;
  const double *unary, *pairs;
  const uint8_t *inter;
  double lambda, rep;
  double pair(size_t a, size_t b) const {
    if (a > b) std::swap(a, b);
    size_t k = b * (b - 1) / 2 + a;
    double e = pairs[2*k] + rep * pairs[2*k+1];
    return inter[k] == 0 ? e : lambda * e;
  }
  size_t desmet() {
    size_t eliminated = 0;
    for (size_t i = 0; i < domains.size(); ++i) {
      auto alive = domains[i];
      if (alive.size() < 2) continue;
      std::vector<double> lo(alive.size()), hi(alive.size());
      for (size_t a = 0; a < alive.size(); ++a) {
        lo[a] = hi[a] = unary[alive[a]];
        for (size_t j = 0; j < domains.size(); ++j) if (i != j) {
          double mn = infinity, mx = -infinity;
          for (auto s : domains[j]) {
            double p = pair(alive[a], s);
            mn = std::min(mn, p); mx = std::max(mx, p);
          }
          lo[a] += mn; hi[a] += mx;
        }
      }
      State survivors;
      for (size_t a = 0; a < alive.size(); ++a) {
        double competitor = infinity;
        for (size_t b = 0; b < alive.size(); ++b)
          if (a != b) competitor = std::min(competitor, hi[b]);
        if (!(lo[a] > competitor)) survivors.push_back(alive[a]);
      }
      eliminated += alive.size() - survivors.size();
      domains[i] = std::move(survivors);
    }
    return eliminated;
  }
  size_t goldstein(size_t passes) {
    size_t total = 0;
    for (size_t pass = 0; pass < passes; ++pass) {
      size_t eliminated = 0;
      for (size_t i = 0; i < domains.size(); ++i) {
        const auto alive = domains[i];
        State survivors;
        for (auto r : alive) {
          bool dead = false;
          for (auto t : alive) if (t != r) {
            double margin = unary[r] - unary[t];
            for (size_t j = 0; j < domains.size(); ++j) if (j != i) {
              double best = infinity;
              for (auto s : domains[j]) best = std::min(best, pair(r,s) - pair(t,s));
              margin += best;
            }
            if (margin > 0.0) { dead = true; break; }
          }
          if (!dead) survivors.push_back(r);
        }
        eliminated += alive.size() - survivors.size();
        domains[i] = std::move(survivors);
      }
      total += eliminated;
      if (!eliminated) break;
    }
    return total;
  }
  double energy(const State& state) const {
    double e = 0;
    for (size_t i = 0; i < state.size(); ++i) {
      e += unary[state[i]];
      for (size_t j = i+1; j < state.size(); ++j) e += pair(state[i],state[j]);
    }
    return e;
  }
  std::pair<State,double> icm(const Domains& starts) const {
    State best;
    double bestEnergy = infinity;
    for (auto state : starts) {
      bool changed;
      do {
        changed = false;
        for (size_t i = 0; i < domains.size(); ++i) {
          auto bestR = state[i];
          double bestCost = infinity;
          for (auto r : domains[i]) {
            double cost = unary[r];
            for (size_t j = 0; j < domains.size(); ++j) if (j != i) cost += pair(r,state[j]);
            // Match Lisp's first-in-domain tie break exactly.
            if (cost < bestCost) {
              bestCost = cost; bestR = r;
            }
          }
          if (bestR != state[i]) { state[i] = bestR; changed = true; }
        }
      } while (changed);
      double e = energy(state);
      if (e < bestEnergy) { bestEnergy = e; best = state; }
    }
    return {best,bestEnergy};
  }
};

struct SearchResult {
  State state;
  double energy, lowerBound;
  size_t nodes = 0;
  bool optimal = true;
};

SearchResult search(const Problem& p, size_t maxNodes, size_t maxEntries,
                    State initial, double initialEnergy,
                    const std::function<void(size_t,double)>& progress = {}, size_t progressEvery = 0) {
  const size_t n = p.domains.size();
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](size_t a,size_t b) {
    return p.domains[a].size() < p.domains[b].size();
  });
  // Compact rows indexed by surviving rotamers, not the largest global index.
  size_t total = 0, largest = 0;
  for (const auto& d : p.domains) for (auto r : d) { ++total; largest = std::max(largest,size_t(r)); }
  std::vector<size_t> row;
  std::vector<double> minima;
  if (n && total <= maxEntries/n) {
    row.resize(largest+1);
    minima.resize(total*n);
    size_t a = 0;
    for (const auto& d : p.domains) for (auto r : d) {
      row[r] = a++;
      for (size_t j = 0; j < n; ++j) {
        double mn = infinity;
        // Same-locus minima are never used (nor are diagonal pairs defined).
        if (std::find(p.domains[j].begin(),p.domains[j].end(),r) == p.domains[j].end())
          for (auto s : p.domains[j]) mn = std::min(mn,p.pair(r,s));
        minima[row[r]*n+j] = mn;
      }
    }
  }
  auto minpair = [&](uint32_t r,size_t j) {
    if (!minima.empty()) return minima[row[r]*n+j];
    double mn = infinity;
    for (auto s : p.domains[j]) mn = std::min(mn,p.pair(r,s));
    return mn;
  };
  State assign(n);
  SearchResult result{std::move(initial),initialEnergy,0};
  auto bound = [&](size_t depth) {
    double val = 0;
    for (size_t i = 0; i < depth; ++i) {
      auto a = assign[order[i]];
      val += p.unary[a];
      for (size_t j = i+1; j < depth; ++j) val += p.pair(a,assign[order[j]]);
    }
    for (size_t i = depth; i < n; ++i) {
      double mn = infinity;
      for (auto r : p.domains[order[i]]) {
        double c = p.unary[r];
        for (size_t j = 0; j < depth; ++j) c += p.pair(assign[order[j]],r);
        for (size_t j = i+1; j < n; ++j) c += minpair(r,order[j]);
        mn = std::min(mn,c);
      }
      val += mn;
    }
    return val;
  };
  std::function<void(size_t)> visit = [&](size_t depth) {
    if (!result.optimal) return;
    if (result.nodes >= maxNodes) { result.optimal = false; return; }
    ++result.nodes;
    if (progress && progressEvery && result.nodes % progressEvery == 0) progress(result.nodes,result.energy);
    double b = bound(depth);
    if (depth == n) {
      if (b < result.energy) { result.energy = b; result.state = assign; }
    } else if (b < result.energy) {
      for (auto r : p.domains[order[depth]]) {
        assign[order[depth]] = r;
        visit(depth+1);
        if (!result.optimal) break;
      }
    }
  };
  result.lowerBound = bound(0);
  visit(0);
  return result;
}
} // namespace dee
} // namespace chem

// The binding scraper remembers the last named namespace, not C++ brace scopes.
// Reopen chem explicitly so the exported entry is registered in chem, not dee.
namespace chem {
#ifndef CANDO_DEE_STANDALONE
namespace {
dee::State deeState(core::SimpleVector_byte32_t_sp v) {
  dee::State state(v->length());
  for (size_t i = 0; i < state.size(); ++i) state[i] = (*v)[i];
  return state;
}
dee::Domains deeDomains(core::SimpleVector_sp input) {
  dee::Domains domains;
  for (size_t i = 0; i < input->length(); ++i)
    domains.push_back(deeState(gc::As<core::SimpleVector_byte32_t_sp>((*input)[i])));
  return domains;
}
core::T_sp deeLispState(const dee::State& s) {
  if (s.empty()) return nil<core::T_O>();
  auto v = core::SimpleVector_byte32_t_O::make(s.size());
  for (size_t i = 0; i < s.size(); ++i) (*v)[i] = s[i];
  return v;
}
core::SimpleVector_sp deeLispDomains(const dee::Domains& domains) {
  auto v = core::SimpleVector_O::make(domains.size());
  for (size_t i = 0; i < domains.size(); ++i) {
    core::List_sp list = nil<core::T_O>();
    for (auto r = domains[i].rbegin(); r != domains[i].rend(); ++r)
      list = core::Cons_O::create(core::make_fixnum(*r),list);
    (*v)[i] = list;
  }
  return v;
}
}

// Operations: 0 Desmet sweep, 1 Goldstein fixpoint, 2 ICM, 3 bounded search.
// Starts are generated in Lisp so resetting *random-state* gives both backends
// identical ICM starts. Unary terms are prepared once, not inside hot loops.
CL_DEFUN core::T_mv chem__dead_end_elimination(int operation, core::SimpleVector_sp domains,
    core::SimpleVector_double_sp unary, core::SimpleVector_double_sp pairs,
    core::SimpleVector_byte8_t_sp inter, double lambda, double repWeight,
    size_t limit, size_t maxEntries, core::SimpleVector_sp starts,
    core::T_sp initialState, double initialEnergy, size_t progressEvery) {
  auto ds = deeDomains(domains);
  if (ds.empty()) SIMPLE_ERROR("DEE requires at least one locus");
  if (!std::isfinite(lambda) || !std::isfinite(repWeight)) SIMPLE_ERROR("DEE weights must be finite");
  size_t largest = 0;
  std::vector<bool> seen(unary->length(),false);
  for (const auto& d : ds) {
    if (d.empty()) SIMPLE_ERROR("DEE requires nonempty domains");
    for (auto r : d) {
      if (r >= unary->length() || seen[r]) SIMPLE_ERROR("Invalid or duplicate DEE rotamer {}",r);
      seen[r] = true;
      if (!std::isfinite((*unary)[r])) SIMPLE_ERROR("Nonfinite DEE unary energy at {}",r);
      largest = std::max(largest,size_t(r));
    }
  }
  size_t required = largest*(largest-1)/2 + largest;
  if (inter->length() < required || pairs->length()/2 < required)
    SIMPLE_ERROR("DEE pair-energy tables are too short");
  dee::Problem p{std::move(ds),&(*unary)[0],pairs->length() ? &(*pairs)[0] : nullptr,
                 inter->length() ? &(*inter)[0] : nullptr,lambda,repWeight};
  // Validate only cross-locus terms actually used; padding may be uninitialized.
  for (size_t i = 0; i < p.domains.size(); ++i)
    for (size_t j = i+1; j < p.domains.size(); ++j)
      for (auto r : p.domains[i]) for (auto s : p.domains[j])
        if (!std::isfinite(p.pair(r,s))) SIMPLE_ERROR("Nonfinite DEE pair energy ({},{})",r,s);
  auto checkState = [&](const dee::State& state) {
    if (state.size() != p.domains.size()) SIMPLE_ERROR("DEE state has incorrect locus count");
    for (size_t i = 0; i < state.size(); ++i)
      if (std::find(p.domains[i].begin(),p.domains[i].end(),state[i]) == p.domains[i].end())
        SIMPLE_ERROR("DEE state selects inactive rotamer {} at locus {}",state[i],i);
  };
  if (operation == 0 || operation == 1) {
    size_t eliminated = operation == 0 ? p.desmet() : p.goldstein(limit);
    return Values(deeLispDomains(p.domains),core::make_fixnum(eliminated));
  }
  if (operation == 2) {
    auto ss = deeDomains(starts);
    if (ss.empty()) SIMPLE_ERROR("DEE ICM requires at least one restart");
    for (const auto& s : ss) checkState(s);
    auto result = p.icm(ss);
    return Values(deeLispState(result.first),core::DoubleFloat_O::create(result.second));
  }
  if (operation != 3) SIMPLE_ERROR("Unknown DEE operation {}",operation);
  dee::State initial;
  if (initialState.notnilp()) {
    initial = deeState(gc::As<core::SimpleVector_byte32_t_sp>(initialState));
    checkState(initial);
    // Never trust an independently supplied incumbent energy for pruning.
    initialEnergy = p.energy(initial);
  } else initialEnergy = dee::infinity;
  auto progress = [](size_t nodes,double best) {
    core::clasp_write_string(fmt::format("    C++ B&B {} nodes, incumbent E {:.3f}\n",nodes,best));
  };
  auto result = dee::search(p,limit,maxEntries,std::move(initial),initialEnergy,progress,progressEvery);
  return Values(deeLispState(result.state),core::DoubleFloat_O::create(result.energy),
                core::make_fixnum(result.nodes),result.optimal ? _lisp->_true() : nil<core::T_O>(),
                core::DoubleFloat_O::create(result.lowerBound));
}
#endif
} // namespace chem
