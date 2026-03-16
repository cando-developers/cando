/*
    File: pairList.h

    Template-based pair list building for nonbond energy components.

    Provides rebuildPairListImpl<Component> that factors out the heavy lifting
    shared by EnergyNonbond_O, EnergyRosettaNonbond_O, and EnergyRosettaElec_O:
      - Cell-list spatial hashing for O(N) neighbor finding
      - Distance screening
      - Bond/angle/1-4 exclusion checks
      - keepInteraction filter callbacks
      - Displacement buffer management
      - Between-matters variant

    Component requirements (duck-typed):
      double          rpairlist() const
      AtomTable_sp    atomTable() const
      core::T_sp      keepInteractionFactory() const
      core::T_sp      matter1() const
      core::T_sp      matter2() const
      void            clearTerms()
      bool            tryAddTerm(Atom_sp a1, Atom_sp a2,
                                 size_t i3x1, size_t i3x2,
                                 core::T_sp keepInteraction)
      void            setDisplacementBuffer(NVector_sp)
      core::T_sp      displacementBuffer() const
      static core::T_sp staticClass()      // already exists on all Lisp classes
*/

#ifndef PairList_H
#define PairList_H

#include <vector>
#include <unordered_map>
#include <cmath>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>   // specializeKeepInteractionFactory
#include <cando/chem/nVector.h>
#include <cando/chem/loop.h>
#include <clasp/core/evaluator.h>

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
};


// ============================================================
//  rebuildPairListBetweenMattersImpl
// ============================================================

// Helper to collect atoms from either a Matter_sp or a SimpleVector_sp of atoms
template <typename Component>
void collectAtomsFromMatter(core::T_sp tmatter, AtomTable_sp atomTable,
                            std::vector<std::pair<Atom_sp, size_t>>& atoms) {
  if (gc::IsA<core::SimpleVector_sp>(tmatter)) {
    core::SimpleVector_sp atomVec = gc::As<core::SimpleVector_sp>(tmatter);
    for (size_t i = 0; i < atomVec->length(); i++) {
      Atom_sp a = gc::As<Atom_sp>((*atomVec)[i]);
      size_t i3 = atomTable->getCoordinateIndexTimes3(a);
      atoms.push_back({a, i3});
    }
  } else if (gc::IsA<core::SimpleVector_sp>(tmatter)) {
    Matter_sp mat = gc::As<Matter_sp>(tmatter);
    Loop l(mat, ATOMS);
    while (l.advanceLoopAndProcess()) {
      Atom_sp a = l.getAtom();
      size_t i3 = atomTable->getCoordinateIndexTimes3(a);
      atoms.push_back({a, i3});
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

  // Collect atoms and coordinates from mat1 and mat2
  // Each can be a Matter_sp (Residue/Molecule/Aggregate) or a SimpleVector_sp of atoms
  std::vector<std::pair<Atom_sp, size_t>> atoms1, atoms2;
  collectAtomsFromMatter<Component>(tmat1, atomTable, atoms1);
  collectAtomsFromMatter<Component>(tmat2, atomTable, atoms2);

  // Build cell grid over all atoms, tagging which set they belong to
  size_t ntotal = atoms1.size() + atoms2.size();
  std::vector<size_t> allI3(ntotal);
  for (size_t i = 0; i < atoms1.size(); i++) allI3[i] = atoms1[i].second;
  for (size_t i = 0; i < atoms2.size(); i++) allI3[atoms1.size() + i] = atoms2[i].second;

  CellGrid grid;
  grid.build(raw, ntotal, allI3, rpairlist);

  size_t n1 = atoms1.size();

  grid.forEachPair([&](size_t li, size_t lj) {
    // Only consider inter-set pairs (one from mat1, one from mat2)
    bool i_in_1 = (li < n1);
    bool j_in_1 = (lj < n1);
    if (i_in_1 == j_in_1) return;  // both in same set, skip

    // Ensure a1 is from mat1, a2 from mat2
    size_t idx1 = i_in_1 ? li : lj;
    size_t idx2 = i_in_1 ? lj - n1 : li - n1;
    Atom_sp a1 = atoms1[idx1].first;
    Atom_sp a2 = atoms2[idx2].first;
    size_t i3x1 = atoms1[idx1].second;
    size_t i3x2 = atoms2[idx2].second;

    // Distance check
    double dx = raw[i3x1]     - raw[i3x2];
    double dy = raw[i3x1 + 1] - raw[i3x2 + 1];
    double dz = raw[i3x1 + 2] - raw[i3x2 + 2];
    double dist2 = dx*dx + dy*dy + dz*dz;
    if (dist2 >= rpairlist2) return;
    // Exclusion check — skip 1-2, 1-3, 1-4 bonded pairs
    {
      EnergyAtom* ea1 = atomTable->getEnergyAtomPointer(a1);
      if (ea1->inBondOrAngle(a2) || ea1->relatedBy14(a2)) return;
    }
    
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
    if (comp->tryAddTerm(a1, a2, i3x1, i3x2, keepInteraction)) {
      ++interactionsKept;
    }
  });

  size_t totalInteractions = interactionsKept + interactionsDiscarded;
  return Values(core::clasp_make_fixnum(interactionsKept),
                core::clasp_make_fixnum(interactionsDiscarded),
                core::clasp_make_fixnum(totalInteractions));
}


// ============================================================
//  rebuildPairListImpl — main pair list builder
// ============================================================

template <typename Component>
core::T_mv rebuildPairListImpl(Component* comp, core::T_sp tcoordinates) {
  // Save displacement buffer for maybeRebuildPairList drift detection
  comp->setDisplacementBuffer(copy_nvector(gc::As<NVector_sp>(tcoordinates)));

  // Delegate to between-matters variant if configured
  if (comp->matter1().notnilp()) {
    if (comp->matter2().notnilp()) {
      return rebuildPairListBetweenMattersImpl(comp, tcoordinates);
    }
    SIMPLE_ERROR("For {} matter1 is {} and matter2 is NIL", _rep_(comp->asSmartPtr()), _rep_(comp->matter1()));
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

        // Build spatial hash grid with cell size = rpairlist
        CellGrid grid;
        grid.build(coords, natoms, i3values, rpairlist);

        // Enumerate candidate pairs from neighboring cells
        grid.forEachPair([&](size_t li, size_t lj) {
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
            if (comp->tryAddTerm(ea1.atom(), ea2.atom(), i3x1, i3x2, keepInteraction)) {
              ++interactionsKept;
            }
          } else {
            ++interactionsDiscarded;
          }
        });

      } else {
        // --- Fallback: read positions from atoms (no coordinate vector) ---
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
                if (comp->tryAddTerm(iea1->atom(), iea2->atom(),
                                     iea1->coordinateIndexTimes3(),
                                     iea2->coordinateIndexTimes3(),
                                     keepInteraction)) {
                  ++interactionsKept;
                }
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

template <typename Component>
core::T_mv maybeRebuildPairListImpl(Component* comp, core::T_sp tcoordinates) {
  auto coords = gc::As<NVector_sp>(tcoordinates);

  core::T_sp dispBuf = comp->displacementBuffer();
  if (dispBuf.nilp()) {
    return rebuildPairListImpl(comp, tcoordinates);
  } else if (gc::IsA<NVector_sp>(dispBuf)) {
    NVector_sp nvDisplacementBuffer = gc::As_unsafe<NVector_sp>(dispBuf);
    if (nvDisplacementBuffer->size() != coords->size()) {
      SIMPLE_ERROR("The size of the _DisplacementBuffer({}) MUST match the size of the coordinates({})",
                   nvDisplacementBuffer->size(), coords->size());
    }
    double rpairlist = comp->rpairlist();
    double rcut      = comp->rcut();
    double skinThickness = rpairlist - rcut;
    double movedTrigger = 0.5 * skinThickness;
    double movedTrigger2 = movedTrigger * movedTrigger;
    vecreal* raw_db = &(*nvDisplacementBuffer)[0];
    vecreal* raw_coords = &(*coords)[0];
    for (size_t ci = 0; ci < nvDisplacementBuffer->size(); ci += 3) {
      const vecreal& dx = raw_db[ci];
      const vecreal& dy = raw_db[ci + 1];
      const vecreal& dz = raw_db[ci + 2];
      const vecreal& cx = raw_coords[ci];
      const vecreal& cy = raw_coords[ci + 1];
      const vecreal& cz = raw_coords[ci + 2];
      if (dx == cx && dy == cy && dz == cz) continue;
      vecreal dist2 = (dx - cx) * (dx - cx) + (dy - cy) * (dy - cy) + (dz - cz) * (dz - cz);
      if (dist2 > movedTrigger2) {
        return rebuildPairListImpl(comp, tcoordinates);
      }
    }
    return Values0<core::T_O>();
  }
  SIMPLE_ERROR("{}: We should never get here", __FUNCTION__);
}


} // namespace chem

#endif // PairList_H
