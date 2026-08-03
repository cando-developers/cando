# `monteCarlo.cc` migration — interleaved `fa_rest`/`fa_rep` + `rep_weight`

## Context

The design `energies` object was refactored so its energy vectors store two channels
per entry — `fa_rest` (everything except the weighted repulsive term) and `fa_rep` (the
raw, unweighted repulsive sum) — **interleaved**: entry `i` lives at `[2i]=fa_rest`,
`[2i+1]=fa_rep`. The total at a given weight is `fa_rest + rep_weight·fa_rep`. Both
channels are rep-weight-independent, so one `energies` object is reused across weights and
the weight is applied at read time (a cheap re-score instead of a scan rebuild).

The Lisp readers (`mcstate-energy`, the DEE solver, the pose/mover plumbing) are already
migrated. `monteCarlo.cc` is a **parallel reader**: its `struct Energies` caches raw
pointers to those vectors and scores off them, still **single-stride**. Until it is
migrated it reads the wrong slots (and the backbone load throws — see Change 1), so every
C++ MC scorer — HREMC, simulated annealing, metropolis, Wang-Landau, and both
`mcstate-energy` variants — is currently broken.

## Design choice

- `lambda` varies **within** a run (per window / per state) → stays a method argument.
- `rep_weight` is **fixed for the whole run** → store it once on the transient `Energies`
  wrapper as `_RepWeight`. The wrapper is constructed fresh per entry-point call, so this
  is exactly the right lifetime, and it mirrors the Lisp side (the `energies` object stays
  weight-independent; the weight is supplied at read).

The payoff: the ~12 scorer **call sites** across HREMC/WL/SA/metropolis do **not change** —
only the three scorer method **bodies** read `this->_RepWeight`, plus the constructor and a
one-line-per-entry-point change.

---

## Change 1 — `struct Energies` members + constructor  (⚠️ hard crash until fixed)

Backbone is no longer a scalar — the accessors `intramolecular-backbone-energy` /
`intermolecular-backbone-energy` now return a 2-element NVector `[fa_rest, fa_rep]`. Store
each as a `double[2]` so the `term()` helper (below) evaluates it uniformly.

Members (lines ~146–147):
```cpp
// OLD
double _IntramolecularBackboneEnergy;
double _IntermolecularBackboneEnergy;
// NEW
double _IntramolecularBackbone[2];   // [0]=fa_rest, [1]=fa_rep
double _IntermolecularBackbone[2];   // [0]=fa_rest, [1]=fa_rep
double _RepWeight = 1.0;             // set by each entry point
```

Constructor load (lines ~161–164) — the old `gc::As<DoubleFloat_sp>(val)->get()` will
**throw** now that `val` is a `SimpleVector`:
```cpp
// OLD
val = core::eval::funcall( _sym_intramolecular_backbone_energy, energies );
this->_IntramolecularBackboneEnergy = gc::As<core::DoubleFloat_sp>(val)->get();
val = core::eval::funcall( _sym_intermolecular_backbone_energy, energies );
this->_IntermolecularBackboneEnergy = gc::As<core::DoubleFloat_sp>(val)->get();
// NEW
val = core::eval::funcall( _sym_intramolecular_backbone_energy, energies );   // NVector [fa_rest, fa_rep]
{ auto v = gc::As<core::SimpleVector_double_sp>(val);
  this->_IntramolecularBackbone[0] = (*v)[0];
  this->_IntramolecularBackbone[1] = (*v)[1]; }
val = core::eval::funcall( _sym_intermolecular_backbone_energy, energies );
{ auto v = gc::As<core::SimpleVector_double_sp>(val);
  this->_IntermolecularBackbone[0] = (*v)[0];
  this->_IntermolecularBackbone[1] = (*v)[1]; }
```

Add the combine helper on the struct (used by both scorer bodies and the backbone):
```cpp
// fa_rest + rep_weight * fa_rep, for entry i of an interleaved vector (or a double[2] with i=0)
inline double term(const double* v, size_t i) const { return v[2*i] + this->_RepWeight*v[2*i+1]; }
```

---

## Change 2 — the three scorer method bodies  (signatures unchanged)

### `Energies::physicalEnergy` (lines ~473–503)
```cpp
// singles (was 482–483)
for (int ii = 0; ii < state._State.size(); ++ii) {
  size_t mrk = state._State[ii];
  singleSum.Add( term(intramolecularSingleTerms, mrk) );
  singleSum.Add( term(intermolecularSingleTerms, mrk) * lambda );
}
// pairs (was 490–493) — the P flags stay SINGLE-stride
size_t lti = this->lowerTriangularIndex(xx,yy);
double pairTerm = term(pairTerms, lti);
int8_t intermolecular_p = intermolecularPPairTerms[lti];
pairSum.Add(intermolecular_p ? pairTerm*lambda : pairTerm);
// backbone (was 498–499)
total.Add( term(this->_IntramolecularBackbone, 0) );
total.Add( term(this->_IntermolecularBackbone, 0) * lambda );
```

### `Energies::deltaReducedEnergy` (lines ~520–544)
```cpp
// single delta (was 529–530)
d.Add( term(intra, newMrk) - term(intra, oldMrk) );
d.Add( (term(inter, newMrk) - term(inter, oldMrk)) * lambda );
// pair delta (was 535–538) — P flags stay SINGLE-stride
size_t ltiNew = this->lowerTriangularIndex(newMrk, yy);
size_t ltiOld = this->lowerTriangularIndex(oldMrk, yy);
double pNew = term(pair, ltiNew), pOld = term(pair, ltiOld);
d.Add( (pflag[ltiNew] ? pNew*lambda : pNew) - (pflag[ltiOld] ? pOld*lambda : pOld) );
// monomer-correction delta (541–543) — UNCHANGED (scalar per lmkey)
```

### `Energies::reducedEnergy` (lines ~505–517)
No edit — it calls `physicalEnergy` (now weight-aware) then adds monomer corrections,
which are **unchanged** (one scalar per lmkey, not interleaved).

---

## Change 3 — entry points: add `(rep-weight 1.0)` key + set `_RepWeight`

Mechanical and identical for each: add `(rep-weight 1.0)` to the `CL_LAMBDA`, add
`double repWeight` to the `CL_DEFUN` signature in the matching position, and immediately
after the `Energies energies(...)` construction add:
```cpp
energies._RepWeight = repWeight;
```
No scorer call inside these functions changes — they pick up the weight through `_RepWeight`.

| Line | Function | notes |
|---|---|---|
| 250 | `chem__mcstate_energy(mcstate, energies, lambda)` | this is `chem:mcstate-energy` (used by `graph.lisp:160`) |
| 416 | `chem__mcstate_energy(energies, beta, mcstate)` | ⚠️ shares the C++ name with 250 (overload) **and** `beta` is unused — it calls `reducedEnergy(state)` with default lambda. Confirm the `CL_LISPIFY_NAME` / whether this variant is still wanted before touching it. |
| 265 | `chem__simulatedAnnealing` | constructs `Energies energies(tenergies,temperatureScale)` |
| 364 | `chem__metropolis…` | its loop calls `reducedEnergy` |
| 568 | `chem__voelz_optimize_monomer_corrections_single_temperature` (Wang-Landau) | see note below |
| 740 | `chem__constantTemperatureHamiltonianReplicaExchangeMonteCarlo` (HREMC) | the main scoring loop |

Scorer call sites that now "just work" via `_RepWeight` (no edits):
HREMC 808, 823, 909, 912 (+`DEBUG_DELTA` 827/851); WL 614, 624 (+`WL_DELTA_CHECK` 618);
SA 296, 306; both `mcstate-energy` 254/419.

---

## Invariants & gotchas

- **Flags** (`_IntermolecularPPairTerms`) and **monomer corrections** (`_MonomerCorrections`)
  stay **single-stride** — one entry per pair / per lmkey. Only the three `double` energy
  vectors (two single-term + `_PairTerms`) interleave.
- **Lengths doubled:** `_PairTerms` and the two single-term vectors are now length `2N`.
  The wrapper sizes histograms/loops off `_MonomerCorrectionsLength` and `_NumberOfSlots`,
  not these — but confirm nothing downstream treats `cl__length(pairTerms)` as a count.
- Running the WL function (568) at `repWeight` makes the resulting monomer corrections
  **weight-specific** (they are a nonlinear log-sum over the DOS, not linearly reweightable).
  That is intended — corrections must be recomputed per weight — just note they are only
  valid at the weight they were built at.

## Lisp follow-ups (once the C++ accepts `:rep-weight`)

- `mover.lisp:460` (HREMC call) → add `:rep-weight (rep-weight pose)`
- `graph.lisp:160` → `(chem:mcstate-energy mcstate energies 1.0 <rep-weight>)`
