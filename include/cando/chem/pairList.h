/*
    File: pairList.h

    Template-based pair list building for nonbond energy components.

    Provides rebuildPairListImpl<Component> that factors out the heavy lifting
    shared by EnergyNonbond_O, EnergyRosettaNonbond_O, and EnergyRosettaElec_O:
      - Cell-list spatial hashing for O(N) neighbor finding
      - Distance screening
      - Bond/angle/1-4 exclusion checks
      - keepInteraction filter callbacks
      - Drift detection (rebuild only when the coordinates have actually moved)
      - Between-matters variant

    Component requirements (duck-typed):
      double          rpairlist() const
      double          rcut() const
      AtomTable_sp    atomTable() const
      core::T_sp      keepInteractionFactory() const
      core::T_sp      matter1() const
      core::T_sp      matter2() const
      void            clearTerms()
      bool            tryAddTerm(Atom_sp a1, Atom_sp a2,
                                 size_t i3x1, size_t i3x2,
                                 core::T_sp keepInteraction)
      bool            tryAddTermCached(Atom_sp a1, Atom_sp a2,
                                 size_t i1, size_t i2,
                                 size_t i3x1, size_t i3x2,
                                 core::T_sp keepInteraction)
      core::T_mv      rebuildPairList(core::T_sp coords)
      size_t          pairListEpoch() const
      double          pairListDrift() const
      void            notePairListBuilt(size_t epoch, double drift)
      void            invalidatePairList()
      static core::T_sp staticClass()      // already exists on all Lisp classes
*/

#ifndef PairList_H
#define PairList_H

#include <vector>
#include <unordered_map>
#include <cmath>
#include <clasp/core/array.h>            // SimpleBitVector_sp, for the query set / target mask
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>   // specializeKeepInteractionFactory
#include <cando/chem/nVector.h>
#include <cando/chem/loop.h>
#include <clasp/core/evaluator.h>

// Leave in for the new cached atom types
#define PAIRLIST_CACHED 1

namespace chem {

// ============================================================
//  CellGrid — spatial hash for O(N) neighbor pair enumeration
// ============================================================

class CellGrid {
public:
  struct CellKey {
    int32_t x, y, z;
    bool operator==(const CellKey& o) const {
      return x == o.x && y == o.y && z == o.z;
    }
  };

  struct CellKeyHash {
    size_t operator()(const CellKey& k) const {
      // Spatial hash using large primes
      return (size_t(uint32_t(k.x)) * 73856093u) ^
             (size_t(uint32_t(k.y)) * 19349663u) ^
             (size_t(uint32_t(k.z)) * 83492791u);
    }
  };

private:
  double _invCellSize;
  std::unordered_map<CellKey, std::vector<size_t>, CellKeyHash> _cells;

public:
  /*! Build the grid from coordinate data.
      atomIndices[i] is the local index (0..natoms-1).
      i3values[i]    is the coordinateIndexTimes3 for local index i.
      coords         is the flat coordinate array. */
  void build(const vecreal* coords,
             size_t natoms,
             const std::vector<size_t>& i3values,
             double cellSize)
  {
    _invCellSize = 1.0 / cellSize;
    _cells.clear();
    for (size_t i = 0; i < natoms; i++) {
      size_t i3 = i3values[i];
      CellKey key{
        (int32_t)std::floor(coords[i3]     * _invCellSize),
        (int32_t)std::floor(coords[i3 + 1] * _invCellSize),
        (int32_t)std::floor(coords[i3 + 2] * _invCellSize)
      };
      _cells[key].push_back(i);
    }
  }

  /*! Visit every unique pair of atom local-indices whose cells are
      the same or adjacent (26-neighborhood).  Uses the 13 forward
      neighbor offsets plus self to guarantee each pair appears once. */
  template <typename Visitor>
  void forEachPair(Visitor&& visit) const {
    static constexpr int fwd[14][3] = {
      {0,0,0},                                          // self
      {1,0,0},  {-1,1,0},  {0,1,0},  {1,1,0},          // 4 forward in y
      {-1,-1,1},{0,-1,1},  {1,-1,1},                    // 3 forward in z
      {-1,0,1}, {0,0,1},   {1,0,1},                     // 3 forward in z
      {-1,1,1}, {0,1,1},   {1,1,1}                      // 3 forward in z
    };
    for (const auto& [key, atoms] : _cells) {
      // Self pairs within this cell
      for (size_t i = 0; i < atoms.size(); i++)
        for (size_t j = i + 1; j < atoms.size(); j++)
          visit(atoms[i], atoms[j]);
      // Forward neighbor cells (avoids double-counting)
      for (int n = 1; n < 14; n++) {
        CellKey nkey{key.x + fwd[n][0], key.y + fwd[n][1], key.z + fwd[n][2]};
        auto it = _cells.find(nkey);
        if (it != _cells.end()) {
          const auto& natoms2 = it->second;
          for (size_t ai : atoms)
            for (size_t aj : natoms2)
              visit(ai, aj);
        }
      }
    }
  }

  /*! Visit every local-index in the 27 cells around the point at COORDS[I3].

      A POINT query, where forEachPair above is a full cross product.  The neighbour-list build
      needs this shape: it asks "what is near THIS atom" once per atom of the whole table, against
      a grid holding only the target atoms.

      All 27 offsets, not the 14 forward ones.  forEachPair can use forward-only because it emits
      each unordered pair once and the querying atoms ARE the gridded atoms; here the querying atom
      is generally not in the grid at all, so every surrounding cell has to be looked at.

      No distance test - the caller applies its own cutoff.  A cell is CELLSIZE on a side, so a
      27-cell neighbourhood covers everything within CELLSIZE and some things beyond it. */
  template <typename Visitor>
  void forEachNear(const vecreal* coords, size_t i3, Visitor&& visit) const {
    CellKey base{
      (int32_t)std::floor(coords[i3]     * _invCellSize),
      (int32_t)std::floor(coords[i3 + 1] * _invCellSize),
      (int32_t)std::floor(coords[i3 + 2] * _invCellSize)
    };
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          auto it = _cells.find(CellKey{base.x + dx, base.y + dy, base.z + dz});
          if (it == _cells.end()) continue;
          for (size_t local : it->second) visit(local);
        }
  }
};


// ============================================================
//  rebuildPairListBetweenMattersImpl
// ============================================================

// Helper to collect coordinate indexes from either a Matter_sp or a
// SimpleVector_sp of atoms.
//
// Only the coordinateIndexTimes3 is collected, not the Atom_sp: the atom is
// recoverable as atomTable->getVectorEnergyAtoms()[i3/3].atom(), because
// AtomTable_O::addAtomInfo assigns coordinateIndex = getNumberOfAtoms()*3, so
// i3/3 IS the atom-table index.  Keeping Atom_sp out of the cached vector means
// the cache holds no GC pointers to trace.
//
// Appends to I3VALUES, so calling it for matter1 then matter2 builds the
// concatenated index vector the CellGrid wants, with matter1 occupying
// [0, n1) and matter2 [n1, ntotal).
//
// Each getCoordinateIndexTimes3 is a Lisp hash-table lookup and this runs on
// every pair-list rebuild, so it is the hot spot in this file.  Caching it per
// component does NOT help: fill-piece-pair-energies sets new matters before
// every single evaluation, so any matter-keyed cache misses every time.
inline void collectAtomsFromMatter(core::T_sp tmatter, AtomTable_sp atomTable,
                                   std::vector<size_t>& i3values) {
  if (gc::IsA<core::SimpleVector_sp>(tmatter)) {
    core::SimpleVector_sp atomVec = gc::As<core::SimpleVector_sp>(tmatter);
    for (size_t i = 0; i < atomVec->length(); i++) {
      Atom_sp a = gc::As<Atom_sp>((*atomVec)[i]);
      i3values.push_back(atomTable->getCoordinateIndexTimes3(a));
    }
  } else if (gc::IsA<Matter_sp>(tmatter)) {
    Matter_sp mat = gc::As<Matter_sp>(tmatter);
    Loop l(mat, ATOMS);
    while (l.advanceLoopAndProcess()) {
      Atom_sp a = l.getAtom();
      i3values.push_back(atomTable->getCoordinateIndexTimes3(a));
    }
  } else {
    TYPE_ERROR(tmatter,core::Cons_O::createList(cl::_sym_or,cl::_sym_SimpleVector_O,_sym_Matter_O));
  };
}

template <typename Component>
core::T_mv rebuildPairListBetweenMattersImpl(Component* comp, core::T_sp tcoordinates) {
  core::T_sp keepInteractionFactory = comp->keepInteractionFactory();

  NVector_sp coords = gc::As<NVector_sp>(tcoordinates);
  vecreal* raw = &(*coords)[0];
  core::T_sp keepInteraction;
  bool hasKeepInteractionFunction = false;
  if (keepInteractionFactory.notnilp()) {
    keepInteraction = specializeKeepInteractionFactory(
        keepInteractionFactory, Component::staticClass());
    hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
  }

  core::T_sp tmat1 = comp->matter1();
  core::T_sp tmat2 = comp->matter2();
  auto atomTable = comp->atomTable();
  double rpairlist  = comp->rpairlist();
  double rpairlist2 = rpairlist * rpairlist;

  size_t interactionsKept = 0;
  size_t interactionsDiscarded = 0;
  comp->clearTerms();

  // Resolve matter1 then matter2 into one concatenated index vector.
  // Each matter can be a Matter_sp (Residue/Molecule/Aggregate) or a
  // SimpleVector_sp of atoms.  matter1 occupies [0, n1), matter2 [n1, ntotal).
  gctools::Vec0<EnergyAtom>& energyAtoms = atomTable->getVectorEnergyAtoms();
  std::vector<size_t> allI3;
  collectAtomsFromMatter(tmat1, atomTable, allI3);
  size_t n1 = allI3.size();
  collectAtomsFromMatter(tmat2, atomTable, allI3);
  size_t ntotal = allI3.size();
  size_t n2 = ntotal - n1;

  // At this point I know n1 and n2 the size of each pair
  // I can calculate the maximum number of pairs that we will generate and
  // decide if we want to use a grid or not.

#if 0
  // ---- TEMPORARY: measure how often a bounding-sphere test would reject ----
  if (n1 > 0 && n2 > 0) {
    static std::atomic<size_t> s_calls{0};
    static std::atomic<size_t> s_reject{0};

    // Calculate the centers c1 and c2 of the n1 and n2 atom sets
    double c1[3]={0,0,0}, c2[3]={0,0,0};
    for (size_t k = 0;  k < n1;     ++k) { size_t i3 = allI3[k];
      c1[0]+=raw[i3]; c1[1]+=raw[i3+1]; c1[2]+=raw[i3+2]; }
    for (size_t k = n1; k < ntotal; ++k) { size_t i3 = allI3[k];
      c2[0]+=raw[i3]; c2[1]+=raw[i3+1]; c2[2]+=raw[i3+2]; }
    for (int d = 0; d < 3; ++d) { c1[d] /= (double)n1; c2[d] /= (double)n2; }

    // Calculate the radius^2 for each of the spheres containing n1 and n2
    double r1sq = 0.0, r2sq = 0.0;
    for (size_t k = 0;  k < n1;     ++k) { size_t i3 = allI3[k];
      double dx=raw[i3]-c1[0], dy=raw[i3+1]-c1[1], dz=raw[i3+2]-c1[2];
      double d2=dx*dx+dy*dy+dz*dz; if (d2>r1sq) r1sq=d2;
    }
    for (size_t k = n1; k < ntotal; ++k) { size_t i3 = allI3[k];
      double dx=raw[i3]-c2[0], dy=raw[i3+1]-c2[1], dz=raw[i3+2]-c2[2];
      double d2=dx*dx+dy*dy+dz*dz; if (d2>r2sq) r2sq=d2;
    }

    // Reject if the spheres are too far apart so that no atoms
    // between them will ever be within a distance to form a pairlist
    double dx=c1[0]-c2[0], dy=c1[1]-c2[1], dz=c1[2]-c2[2];
    double reach = std::sqrt(r1sq) + std::sqrt(r2sq) + rpairlist;
    bool rejected = (dx*dx+dy*dy+dz*dz) > reach*reach;

    // fetch_add returns the PREVIOUS value; +1 gives this call's ordinal, so
    // exactly one thread observes each multiple of the report interval.
    size_t nth = s_calls.fetch_add(1, std::memory_order_relaxed) + 1;
    if (rejected) s_reject.fetch_add(1, std::memory_order_relaxed);
    if ((nth % 500000) == 0) {
      size_t rej = s_reject.load(std::memory_order_relaxed);
      fmt::print(stderr,
                 "sphere-test[{}]: {} / {} rejected ({:.1f}%)  n1={} n2={} rpairlist={:.2f}\n",
                 Component::staticClass()->_classNameAsString(),
                 rej, nth, 100.0*rej/nth, n1, n2, rpairlist);
    }
  }
#endif

  // ---- grid or brute force? -------------------------------------------
  // CellGrid cells are rpairlist on a side, so the grid can only prune when the
  // atoms span more than about one cell.  Two sidechains (~8 A across, 9 A cell)
  // land in one cell: forEachPair degenerates to all-pairs and the hashing plus
  // per-cell vectors are pure overhead.  A sidechain against a whole receptor
  // backbone spans tens of cells and the grid earns its keep.
  bool useGrid = false;
  if (ntotal > 1) {
    double lo[3] = { raw[allI3[0]], raw[allI3[0]+1], raw[allI3[0]+2] };
    double hi[3] = { lo[0], lo[1], lo[2] };
    for (size_t k = 1; k < ntotal; ++k) {
      size_t i3 = allI3[k];
      for (int d = 0; d < 3; ++d) {
        double v = raw[i3+d];
        if (v < lo[d]) lo[d] = v;
        if (v > hi[d]) hi[d] = v;
      }
    }
    double cells = std::max(1.0, std::ceil((hi[0]-lo[0])/rpairlist))
                       * std::max(1.0, std::ceil((hi[1]-lo[1])/rpairlist))
                       * std::max(1.0, std::ceil((hi[2]-lo[2])/rpairlist));
    useGrid = (cells > 8.0) && ((double)n1 * (double)n2 > 2048.0);
  }

#ifdef PAIRLIST_CACHED
  comp->ensureParameterCache();
#endif
  // Hoisted out of the per-pair lambda: set once by setMatters, constant for the build.
  const bool exclusionsPossible = comp->exclusionsPossible();
  auto consider = [&](size_t li, size_t lj) {
      // Only consider inter-set pairs (one from mat1, one from mat2)
      bool i_in_1 = (li < n1);
      bool j_in_1 = (lj < n1);
      if (i_in_1 == j_in_1) return;  // both in same set, skip

      // Ensure a1 is from mat1, a2 from mat2.  li/lj already index allI3.
      size_t i3x1 = allI3[i_in_1 ? li : lj];
      size_t i3x2 = allI3[i_in_1 ? lj : li];

      // Distance check
      double dx = raw[i3x1]     - raw[i3x2];
      double dy = raw[i3x1 + 1] - raw[i3x2 + 1];
      double dz = raw[i3x1 + 2] - raw[i3x2 + 2];
      double dist2 = dx*dx + dy*dy + dz*dz;
      if (dist2 >= rpairlist2) return;
      // Atoms are only needed past the cutoff - the distance test above reads raw[]
      // directly, so fetching them earlier costs two Vec0 lookups, two divisions and
      // two smart-pointer constructions on every rejected pair.  i3/3 is the
      // atom-table index, the same identity the exclusion check relies on.
      // NOTE: a1/a2 must stay in the lambda's scope - the keepInteraction filter
      // and tryAddTermCached below both use them.
      EnergyAtom& ea1 = energyAtoms[i3x1/3];
      Atom_sp a1 = ea1.atom();
      Atom_sp a2 = energyAtoms[i3x2/3].atom();

      // Exclusion check — skip 1-2, 1-3, 1-4 bonded pairs.  exclusionsPossible is
      // false when no bonded path can span the two matters (two sidechains at
      // different loci, two molecules, a sidechain vs an intermolecular partner);
      // the three container probes below would then always miss.
      if (exclusionsPossible && (ea1.inBondOrAngle(a2) || ea1.relatedBy14(a2))) return;

      // keepInteraction filter (only if factory was provided)
      if (hasKeepInteractionFunction) {
        core::T_sp result = core::eval::funcall(keepInteraction,
                                                a1, a2,
                                                core::make_fixnum(i3x1),
                                                core::make_fixnum(i3x2));
        if (result.nilp()) {
          ++interactionsDiscarded;
          return;
        }
      }
#ifdef PAIRLIST_CACHED
      if (comp->tryAddTermCached(a1, a2, i3x1/3, i3x2/3, i3x1, i3x2, keepInteraction)) {
        ++interactionsKept;
      }
#else
      if (comp->tryAddTerm(a1, a2, i3x1, i3x2, keepInteraction)) {
        ++interactionsKept;
      }
#endif
      };

  if (useGrid) {
    CellGrid grid;
    grid.build(raw, ntotal, allI3, rpairlist );
    grid.forEachPair(consider);
  } else {
    for (size_t li = 0; li < n1; ++li)
      for (size_t lj = n1; lj < ntotal; ++lj)
        consider(li, lj);
  }

  size_t totalInteractions = interactionsKept + interactionsDiscarded;
  return Values(core::clasp_make_fixnum(interactionsKept),
                core::clasp_make_fixnum(interactionsDiscarded),
                core::clasp_make_fixnum(totalInteractions));
}


// ============================================================
//  rebuildPairListImpl — main pair list builder
// ============================================================

// ============================================================
//  rebuildPairListFromNeighborsImpl - query-set driven
// ============================================================

/*! Build the pair list by walking the AtomTable's neighbour list from the component's QUERY SET,
    instead of gridding the whole table.

    The grid path examines every pair in the table and discards what the predicate rejects.  For a
    fan-out slot that is ~6.5 million candidates to keep ~4500 - the work scales with the protein
    while the answer scales with the slot.  Here the work scales with the slot.

    TWO LOOPS, and both are needed:

      CROSS   each query atom against its neighbour list.  The list holds only TARGET atoms, so a
              slot finds backbone partners and never another rotamer's.
      INTERNAL both atoms inside a RANGE query.  Those pairs are not in the lists - a slot's own
              atoms are not targets - and for a rotamer they are real energy: the 1-5-and-beyond
              contacts within the sidechain.  ~105 candidates over a contiguous range, brute force.

    DEDUPLICATION: when a query atom is ITSELF a target - the backbone querying the backbone - each
    pair is reachable from both ends, so it is emitted only when i3 < j3.  When the query atom is
    not a target, as for a slot, that cannot happen and every pair is emitted once. */
template <typename Component>
core::T_mv rebuildPairListFromNeighborsImpl(Component* comp, core::T_sp tcoordinates) {
  NVector_sp coordsObj = gc::As<NVector_sp>(tcoordinates);
  vecreal* raw = &(*coordsObj)[0];
  auto atomTable = comp->atomTable();
  double rpairlist = comp->rpairlist();
  double rpairlist2 = rpairlist * rpairlist;

  // One build serves every component and every slot at this fold: the cutoff test is >=, so a list
  // built at the widest rpairlist answers the narrower ones, and the coordinates do not change
  // across a fill.
  //
  // BUILD at queryCutoff when the caller set one - the max across the group - so the list does not
  // depend on which component happens to evaluate first.  TEST against this component's own
  // rpairlist2 below regardless: a wider list only offers more candidates, it does not widen the
  // interaction.
  double buildCutoff = (comp->queryCutoff() > 0.0) ? comp->queryCutoff() : rpairlist;
  if (buildCutoff < rpairlist) buildCutoff = rpairlist;
  atomTable->ensureNeighborList(tcoordinates, buildCutoff, comp->queryTargetMask());

  // AFTER ensureNeighborList, which brings the shared drift up to date and may re-snapshot it.
  // Stamping before that would record an epoch this list does not actually belong to.
  comp->notePairListBuilt(atomTable->refEpoch(), atomTable->maxDrift());
  comp->clearTerms();

#ifdef PAIRLIST_CACHED
  // NOT optional, and omitting it is why this path first failed with "LKSolvation term cache miss
  // for slot pair 8,0 of 0".  tryAddTermCached reads _TypeSlot and indexes _TermCache by
  // _NTypeSlots; ensureParameterCache is what fills both.  Without it _NTypeSlots stays 0 while
  // _TypeSlot may still hold slots from an earlier build, so the lookup runs off a zero-sized
  // table.  Both grid paths call it in the same position - see the two sites below.
  comp->ensureParameterCache();
#endif

  core::T_sp keepInteraction;
  bool hasKeepInteractionFunction = false;
  if (comp->keepInteractionFactory().notnilp()) {
    keepInteraction = specializeKeepInteractionFactory(comp->keepInteractionFactory(),
                                                       Component::staticClass());
    hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
  }
  const bool exclusionsPossible = comp->exclusionsPossible();
  auto& energyAtoms = atomTable->getVectorEnergyAtoms();
  auto& start = atomTable->_NeighborStart;
  auto& neighbors = atomTable->_Neighbors;

  size_t interactionsKept = 0;
  size_t interactionsDiscarded = 0;

  core::SimpleBitVector_sp targets;
  bool hasTargets = comp->queryTargetMask().notnilp();
  if (hasTargets) targets = gc::As_unsafe<core::SimpleBitVector_sp>(comp->queryTargetMask());

  auto consider = [&](size_t i3x1, size_t i3x2) {
      double dx = raw[i3x1]     - raw[i3x2];
      double dy = raw[i3x1 + 1] - raw[i3x2 + 1];
      double dz = raw[i3x1 + 2] - raw[i3x2 + 2];
      if (dx*dx + dy*dy + dz*dz >= rpairlist2) return;
      EnergyAtom& ea1 = energyAtoms[i3x1/3];
      Atom_sp a1 = ea1.atom();
      Atom_sp a2 = energyAtoms[i3x2/3].atom();
      if (exclusionsPossible && (ea1.inBondOrAngle(a2) || ea1.relatedBy14(a2))) return;
      if (hasKeepInteractionFunction) {
        core::T_sp result = core::eval::funcall(keepInteraction, a1, a2,
                                                core::make_fixnum(i3x1),
                                                core::make_fixnum(i3x2));
        if (result.nilp()) { ++interactionsDiscarded; return; }
      }
#ifdef PAIRLIST_CACHED
      if (comp->tryAddTermCached(a1, a2, i3x1/3, i3x2/3, i3x1, i3x2, keepInteraction))
        ++interactionsKept;
#else
      if (comp->tryAddTerm(a1, a2, i3x1, i3x2, keepInteraction)) ++interactionsKept;
#endif
    };

  // ---- CROSS: query atoms against their neighbour lists ----
  auto queryAtom = [&](size_t i3) {
      size_t i = i3 / 3;
      bool queryIsTarget = hasTargets && targets->testBit(i3);
      for (uint32_t k = start[i]; k < start[i+1]; ++k) {
        size_t j3 = neighbors[k];
        if (queryIsTarget && j3 < i3) continue;   // reachable from both ends - emit once
        consider(i3, j3);
      }
    };

  if (comp->queryHi3() != 0) {
    for (size_t i3 = comp->queryLo3(); i3 < comp->queryHi3(); i3 += 3) queryAtom(i3);
    // ---- INTERNAL: both atoms inside the range ----
    for (size_t a3 = comp->queryLo3(); a3 < comp->queryHi3(); a3 += 3)
      for (size_t b3 = a3 + 3; b3 < comp->queryHi3(); b3 += 3)
        consider(a3, b3);
  } else {
    core::SimpleBitVector_sp qmask = gc::As<core::SimpleBitVector_sp>(comp->queryMask());
    size_t n3 = atomTable->getNumberOfAtoms() * 3;
    for (size_t i3 = 0; i3 < n3; i3 += 3)
      if (qmask->testBit(i3)) queryAtom(i3);
    // No INTERNAL loop for a mask query: a mask-driven caller's own atoms are in the target set
    // (the backbone querying the backbone), so those pairs already came through the lists above.
  }

  return Values(core::make_fixnum(interactionsKept),
                core::make_fixnum(interactionsDiscarded));
}


template <typename Component>
core::T_mv rebuildPairListImpl(Component* comp, core::T_sp tcoordinates) {
  // Count every rebuild, whichever enumeration it ends up using.  All four components reach a
  // rebuild through here, so this is the one place that sees them all.  Tests read it to prove the
  // CACHED path was exercised - see AtomTable_O::_PairListRebuilds.
  {
    auto at = comp->atomTable();
    if (at.notnilp()) at->notePairListRebuild();
  }
  // A query set replaces the full-table enumeration entirely - see
  // rebuildPairListFromNeighborsImpl.  Checked before the matters branch because it is the more
  // specific instruction: matters say which two SETS to pair, a query set says which atoms to
  // start from and hands the rest to the neighbour list.
  if (comp->hasQuerySet()) {
    if (comp->matter1().notnilp()) {
      SIMPLE_ERROR("{} has both a query set and matters set - they are alternative ways to bound "
                   "the same enumeration and only one can apply",
                   _rep_(comp->asSmartPtr()));
    }
    return rebuildPairListFromNeighborsImpl(comp, tcoordinates);
  }

  // Delegate to between-matters variant if configured
  if (comp->matter1().notnilp()) {
    if (comp->matter2().notnilp()) {
      // Between-matters callers drive rebuilds via set-matters, which invalidates the list, so
      // drift detection is never used here.  Leave the stamp at "never built" so a list produced
      // here can never be mistaken for one that drift alone may keep alive - the matters can change
      // without a single coordinate moving, and drift would report all is well.
      comp->invalidatePairList();
      return rebuildPairListBetweenMattersImpl(comp, tcoordinates);
    }
    SIMPLE_ERROR("For {} matter1 is {} and matter2 is NIL", _rep_(comp->asSmartPtr()), _rep_(comp->matter1()));
  }

  // Stamp the shared drift for maybeRebuildPairList - see EnergyPairlistComponent_O::_PairListEpoch.
  // No atom table means no snapshot to stamp against; leaving the epoch at 0 makes every later
  // evaluation rebuild, which is the safe reading of "I cannot tell whether anything moved".
  {
    auto at = comp->atomTable();
    if (at.notnilp()) {
      at->driftFor(tcoordinates);
      comp->notePairListBuilt(at->refEpoch(), at->maxDrift());
    }
  }

  size_t interactionsKept = 0;
  size_t interactionsDiscarded = 0;
  size_t totalInteractions = 0;
  comp->clearTerms();

  if (comp->keepInteractionFactory().notnilp()) {
    core::T_sp keepInteraction = specializeKeepInteractionFactory(
        comp->keepInteractionFactory(), Component::staticClass());
    bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);

    double rpairlist  = comp->rpairlist();
    double rpairlist2 = rpairlist * rpairlist;
    auto atomTable = comp->atomTable();
    size_t natoms = atomTable->getNumberOfAtoms();

    if (natoms >= 2) {
      vecreal* coords = NULL;
      if (gc::IsA<NVector_sp>(tcoordinates)) {
        coords = &(*gc::As_unsafe<NVector_sp>(tcoordinates))[0];
      }
      if (coords) {
        // --- Cell-list accelerated path ---
        auto& energyAtoms = atomTable->getVectorEnergyAtoms();

        // Extract coordinate indices for grid construction
        std::vector<size_t> i3values(natoms);
        for (size_t i = 0; i < natoms; i++) {
          i3values[i] = energyAtoms[i].coordinateIndexTimes3();
        }

        // ---- grid or brute force? -------------------------------------
        // CellGrid cells are rpairlist on a side, so the grid can only prune
        // when the atoms span more than about one cell.  A compact group
        // (one residue, ~8 A across against a 9 A cell) lands in a single
        // cell: forEachPair degenerates to all-pairs and the hashing plus
        // per-cell vectors are pure overhead.
        bool useGrid = false;
        {
          double lo[3] = { coords[i3values[0]], coords[i3values[0]+1], coords[i3values[0]+2] };
          double hi[3] = { lo[0], lo[1], lo[2] };
          for (size_t k = 1; k < natoms; ++k) {
            size_t i3 = i3values[k];
            for (int d = 0; d < 3; ++d) {
              double v = coords[i3+d];
              if (v < lo[d]) lo[d] = v;
              if (v > hi[d]) hi[d] = v;
            } 
          }   
          double cells = std::max(1.0, std::ceil((hi[0]-lo[0])/rpairlist))
                             * std::max(1.0, std::ceil((hi[1]-lo[1])/rpairlist))
                             * std::max(1.0, std::ceil((hi[2]-lo[2])/rpairlist));
          double npairs = (double)natoms * (double)(natoms - 1) * 0.5;
          useGrid = (cells > 8.0) && (npairs > 2048.0);
        }             
            
#ifdef PAIRLIST_CACHED
        comp->ensureParameterCache();
#endif


        // Enumerate candidate pairs from neighboring cells
        auto consider = [&](size_t li, size_t lj) {
          EnergyAtom& ea1 = energyAtoms[li];
          EnergyAtom& ea2 = energyAtoms[lj];

          // Distance check first (cheap) — skip far pairs early
          size_t i3x1 = ea1.coordinateIndexTimes3();
          size_t i3x2 = ea2.coordinateIndexTimes3();
          double dx = coords[i3x1]     - coords[i3x2];
          double dy = coords[i3x1 + 1] - coords[i3x2 + 1];
          double dz = coords[i3x1 + 2] - coords[i3x2 + 2];
          double dist2 = dx*dx + dy*dy + dz*dz;
          if (dist2 >= rpairlist2) return;

          // Exclusion check (set lookups) — only for pairs within cutoff
          if (ea1.inBondOrAngle(ea2.atom()) || ea1.relatedBy14(ea2.atom())) return;

          totalInteractions++;

          // keepInteraction filter
          bool keep = true;
          if (hasKeepInteractionFunction) {
            core::T_sp result = core::eval::funcall(keepInteraction,
                                                    ea1.atom(), ea2.atom(),
                                                    core::make_fixnum(i3x1),
                                                    core::make_fixnum(i3x2));
            keep = result.notnilp();
          }
          if (keep) {
#ifdef PAIRLIST_CACHED
            if (comp->tryAddTermCached(ea1.atom(), ea2.atom(), li, lj, i3x1, i3x2, keepInteraction)) {
              ++interactionsKept;
            }
#else
            if (comp->tryAddTerm(ea1.atom(), ea2.atom(), i3x1, i3x2, keepInteraction)) {
              ++interactionsKept;
            }
#endif
          } else {
            ++interactionsDiscarded;
          }
        };

        if (useGrid) {
          CellGrid grid;
          grid.build(coords, natoms, i3values, rpairlist);
          grid.forEachPair(consider);
        } else {
          for (size_t li = 0; li + 1 < natoms; ++li)
            for (size_t lj = li + 1; lj < natoms; ++lj)
              consider(li, lj);
        } 
      } else {
        // --- Fallback: read positions from atoms (no coordinate vector) ---
#ifdef PAIRLIST_CACHED
        comp->ensureParameterCache();
#endif
        gctools::Vec0<EnergyAtom>::iterator iea1;
        gctools::Vec0<EnergyAtom>::iterator iea2;
        for (iea1 = atomTable->begin(); iea1 != atomTable->end() - 1; iea1++) {
          Vector3 pos1 = iea1->atom()->getPosition();
          for (iea2 = iea1 + 1; iea2 != atomTable->end(); iea2++) {
            if (iea1->inBondOrAngle(iea2->atom()) || iea1->relatedBy14(iea2->atom()))
              continue;
            totalInteractions++;
            Vector3 pos2 = iea2->atom()->getPosition();
            Vector3 delta = pos1 - pos2;
            double dist2 = delta.dotProduct(delta);
            if (dist2 < rpairlist2) {
              bool keep = true;
              if (hasKeepInteractionFunction) {
                core::T_sp result = core::eval::funcall(keepInteraction,
                                                        iea1->atom(), iea2->atom(),
                                                        core::make_fixnum(iea1->coordinateIndexTimes3()),
                                                        core::make_fixnum(iea2->coordinateIndexTimes3()));
                keep = result.notnilp();
              }
              if (keep) {
#ifdef PAIRLIST_CACHED
                size_t li = (size_t)(iea1 - atomTable->begin());
                size_t lj = (size_t)(iea2 - atomTable->begin());
                if (comp->tryAddTermCached(iea1->atom(), iea2->atom(), li, lj,
                                           iea1->coordinateIndexTimes3(),
                                           iea2->coordinateIndexTimes3(),
                                           keepInteraction)) {
                  ++interactionsKept;
                }
#else   
                if (comp->tryAddTerm(iea1->atom(), iea2->atom(),
                                     iea1->coordinateIndexTimes3(),
                                     iea2->coordinateIndexTimes3(),
                                     keepInteraction)) {
                  ++interactionsKept;
                }
#endif
              } else {
                ++interactionsDiscarded;
              }
            }
          }
        }
      }
    }
  }
  return Values(core::clasp_make_fixnum(interactionsKept),
                core::clasp_make_fixnum(interactionsDiscarded),
                core::clasp_make_fixnum(totalInteractions));
}


// ============================================================
//  maybeRebuildPairListImpl — only rebuild if atoms drifted
// ============================================================

/*! Rebuild only if the coordinates have drifted far enough to invalidate the pair list.
 *
 * O(1).  This used to walk a private copy of every coordinate - see
 * EnergyPairlistComponent_O::_PairListEpoch for what that cost across a blueprint's ~2829
 * components.  The walk now happens ONCE per evaluation, in AtomTable_O::updateDrift, and what is
 * left here is two loads and a compare.
 *
 * The bound is the triangle inequality: the component needs the drift since ITS OWN rebuild, and
 * what the atom table knows is the drift since a shared reference, so
 *
 *     |x_now - x_rebuilt| <= |x_now - ref| + |ref - x_rebuilt| <= maxDrift() + pairListDrift()
 *
 * Conservative in the safe direction: it can rebuild early, it cannot skip a rebuild that is due.
 * NOTHING IS LOST relative to the old walk, which was itself global - it asked whether ANY atom in
 * the table had moved, not just the ones this component enumerates. */
template <typename Component>
core::T_mv maybeRebuildPairListImpl(Component* comp, core::T_sp tcoordinates) {
  auto atomTable = comp->atomTable();
  if (atomTable.nilp()) {
    // Nothing to measure drift against.  A component in this state cannot enumerate anything
    // either - rebuildPairListImpl only reaches the grid through the atom table - so rebuilding
    // unconditionally is both correct and cheap.
    return comp->rebuildPairList(tcoordinates);
  }
  double movedTrigger = 0.5 * (comp->rpairlist() - comp->rcut());
  // Declare the threshold so the re-snapshot policy knows how far it may let the drift run before
  // resetting the reference - see AtomTable_O::updateDrift.  STAYS FIRST: the policy needs the
  // trigger registered even on the paths that return early below.
  atomTable->noteTrigger(movedTrigger);

  // BEFORE the drift read, because this branch rebuilds unconditionally and would discard it.  The
  // ordering note below governs the two branches that COMPARE epoch against drift; this one
  // compares nothing, so nothing can go stale by answering it early.  It matters during a pair
  // scan, where set-matters invalidates every component on every iteration.
  if (comp->pairListEpoch() == 0) {                       // never built, or explicitly invalidated
    return comp->rebuildPairList(tcoordinates);
  }

  // FIRST relative to the epoch comparisons below, because it may re-snapshot and so change
  // refEpoch.  Reading the epoch before this and the drift after would compare against two
  // different references and silently keep a stale list.
  double drift = atomTable->driftFor(tcoordinates);

  if (comp->pairListEpoch() != atomTable->refEpoch()) {    // reference moved out from under us
    return comp->rebuildPairList(tcoordinates);
  }
  if (drift + comp->pairListDrift() > movedTrigger) {
    return comp->rebuildPairList(tcoordinates);
  }
  return Values0<core::T_O>();
}


} // namespace chem

#endif // PairList_H
