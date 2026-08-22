/*
    File: energyAtomTable.h
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyAtomTable.h
 *
 *	Maintains a database of atoms
 */

#ifndef EnergyAtomTable_H
#define	EnergyAtomTable_H //[
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/nativeVector.fwd.h>
#include <cando/units/quantity.fwd.h>
#include <cando/chem/indirectAtomCoordinateReference.h>
#include <cando/chem/molecule.h>
#include <cando/chem/energyComponent.fwd.h>
#include <cando/geom/vector3.h>


namespace chem
{

/*! Verlet skin for the shared neighbour list, in Angstroms.
 *
 * The list is built at (requested cutoff + this), and stays valid while every atom has drifted less
 * than half the leftover slack - the same criterion maybeRebuildPairListImpl uses per component,
 * where the skin is rpairlist - rcut (pairList.h:735).  Bigger means a longer-lived list and more
 * candidate pairs per query; 2.0 A is the usual compromise and matches the rpairlist/rcut gaps
 * already in the kernels (nonbond 9.0/6.0, LK 8.0/6.0, elec 7.5/5.5).
 *
 * Zero would be CORRECT but useless: with no slack any motion at all forces a rebuild. */
#define NEIGHBOR_LIST_SKIN 2.0

SMART(Atom);
SMART(ForceField);
SMART(Matter);
SMART(FFNonbondDb);

/*! Store a pointer to an Atom and an index into the coordinate vector array
 * which stores coordinates in a 1D array (x1,y1,z1,x2,y2,z2,x3,...,xN,yN,zN)
 */
class EnergyAtom : public IndirectAtomCoordinateReference {
  typedef IndirectAtomCoordinateReference Base;
public:
  MatterName		_AtomName;		//!<Atom name
  REAL                  _Charge;		//!<Atom charge in electrons
  REAL                  _Mass;			//!<Atom mass in Daltons
  INT                   _TypeIndex;		//!<Type index
  INT                   _AtomicNumber;           //!<Atomic number
  size_t                _Flag;
public:
	// Temporary variables, not necessary to store
	/*! Sets of all atoms that are bonded to this one at remove 0(bonded),
	  1(ends of angle) and 2(ends of dihedral) (indexed at 0, 1, 2 respectively */
  static const int max_remove = 2;
  gctools::SmallOrderedSet<Atom_sp> _AtomsAtRemoveBondAngle14[max_remove+1]; // s e t<Atom_sp>	_AtomsAtRemoveBondAngle14[3];
public:

  bool inBondOrAngle(Atom_sp a) { return (this->_AtomsAtRemoveBondAngle14[0].contains(a)) || (this->_AtomsAtRemoveBondAngle14[1].count(a)>0);};
  bool relatedBy14(Atom_sp a) { return (this->_AtomsAtRemoveBondAngle14[2].contains(a)); };

  EnergyAtom();
  EnergyAtom(Atom_sp atom, uint coordinateIndex);
  EnergyAtom(core::T_sp forceField, Atom_sp atom, uint coordinateIndex, core::HashTable_sp atomTypes, core::T_sp keepInteractionFactory );

		// methods
  core::List_sp encode() const;
  void		defineForAtom(core::T_sp forceField, Atom_sp atom, uint coordinateIndex, core::HashTable_sp atomTypes, core::T_sp keepInteractionFactory );
  string		getResidueAndName(Residue_sp res);
//	adapt::QDomNode_sp	asXmlRelativeToContainer(chem::Matter_sp parent);
//	void		parseFromXmlRelativeToContainer( adapt::QDomNode_sp xml, chem::Matter_sp parent );
};
};

namespace translate {

template <>
struct	to_object<chem::EnergyAtom >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::EnergyAtom& ea)
  {
    return ea.encode();
  }
};

template <>
struct	from_object<chem::EnergyAtom>
{
  typedef	chem::EnergyAtom	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR("Implement me");
        }
};
};

namespace chem {

EnergyAtom*	_findEnergyAtom(gctools::Vec0<EnergyAtom>& atoms, uint	idx3 );

SMART(AtomTable);
class AtomTable_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,AtomTable_O,"AtomTable",core::CxxObject_O);

 public:
  bool fieldsp() const { return true; }
  void fields(core::Record_sp node);
 public:
  void initialize();
public:
  gctools::Vec0<EnergyAtom>	      _Atoms;
  core::HashTable_sp                _AtomTableIndexes; // m a p<Atom_sp,uint>	_AtomTableIndexes;
  core::ComplexVector_int32_t_sp      _ResiduePointers;
  core::ComplexVector_T_sp            _ResidueNames;
  //! Store count of atoms in each molecule. The length of this vector is the number of molecules
  core::ComplexVector_int32_t_sp      _AtomsPerMolecule;
  //! Store the index of the first "solvent" molecule
  core::T_sp                     _firstSolventMoleculeNSPSOL;
  core::T_sp                     _finalSoluteResidueIPTRES;
  core::T_sp                     _totalNumberOfMoleculesNSPM;
  //! Stores actual residues from aggregate
  core::Vector_sp                _Residues;
  core::Vector_sp                _Molecules;
  core::T_sp                     _Matter;
  core::T_sp                     _AggregateName;
  core::T_sp                     _BoundingBox;
  core::T_sp                     _NonbondForceFieldForAggregate;
  core::T_sp                     _LKSolvationForceFieldForAggregate;
  /*! Shared LK-solvation type slots.
   *
   * An atom's LK solvation type is a property of the ATOM, not of the component asking for it, so
   * the atom-index -> parameter-slot mapping is identical for every EnergyRosettaLKSolvation_O
   * built over this atom table.  It used to be rebuilt per component, one Lisp funcall into
   * find-lksolvation-type per atom: with a blueprint's 471 slot-group components over ~6700 atoms
   * that is ~3.2 million funcalls before a single pair energy is computed.  Computed once here it
   * is ~6700.
   *
   * Only the SLOTS are shared.  Each component still builds its own _TermCache, because that
   * depends on its own _Parameters, and it is O(distinct-types^2) - a handful.
   *
   * INVALIDATION: keyed on the force field and the atom count via lkTypeSlotsValidFor.  Mutating
   * an atom's :lk-solvation-atom-type after the first component has cached will NOT be noticed -
   * call invalidateLKTypeSlots() if that ever happens.
   */
  gctools::Vec0<int>                    _LKTypeSlot;   // atom index -> slot, -1 = no params
  //! FLAT, 4 doubles per slot: [4s+0]=dGfree [4s+1]=lambda [4s+2]=radius [4s+3]=volume.
  //! Not Vec0<std::array<double,4>>, which reads better but would be the only std:: instantiation
  //! of Vec0 in the tree - the GC analyzer has precedent for Vec0 of builtins and of plain structs
  //! (Vec0<int>, Vec0<uint8_t>, Vec0<rosetta_lk_solvation_term>) and none for a std:: class
  //! template.  Not worth finding out during a build.
  gctools::Vec0<double>                 _LKUniq;
  core::T_sp                            _LKCachedForceField;
  /*! Bumped on every rebuild AND on every invalidation.
   *
   * Components keep their own copy of the slots and their own _TermCache derived from _LKUniq, so
   * "is the SHARED table valid" is not the question they need answered - it is "is the shared
   * table the same one I copied from".  A boolean cannot express that: after an invalidation the
   * first component rebuilds the shared table, and every other component would then see a valid
   * shared table plus its own already-set _CachedForAtomTable and return early holding stale
   * slots.  Comparing generations catches exactly that.
   */
  size_t                                _LKGeneration = 0;

  /*! The same sharing, for the ROSETTA-NONBOND (radius, epsilon) slots.
   *
   * An atom's nonbond type is a property of the ATOM and the force field, never of the component
   * asking - so this mapping is identical for all of a blueprint's slot-group components.
   *
   * Built per component it cost TWO hash lookups per atom: Atom_O::getType is a gethash
   * (atom.cc:1712) and chem__FFNonbond_findType is a second one keyed on the type STRING, so the
   * name is re-hashed every time.  Every slot group owns its own copies of the three pairlist
   * components and a copy starts with an invalid cache (energyRosettaNonbond.cc:284), so a
   * blueprint pays ~2820 cold starts of ~7000 atoms each - roughly 40 million lookups, measured at
   * 49% of the single-scan fill with the string hashing alone at 11%.  Built once here it is ~7000.
   *
   * The cache was never THRASHING: it is called from every pair-list rebuild (pairList.h:316 and
   * friends), n^2/2 times during a pair scan, and the guard returns immediately every time after
   * the first.  All of the cost is in the cold starts, and sharing is what removes them.
   *
   * Only the SLOTS are shared, exactly as for LK.  Each component still builds its own _TermCache,
   * which depends on its own _Parameters and is O(distinct-types^2) - a handful.
   *
   * INVALIDATION: same contract as the LK table.  Retyping an atom after the first component has
   * cached will NOT be noticed - call invalidateNBTypeSlots() if that ever happens.
   */
  gctools::Vec0<int>                    _NBTypeSlot;   // atom index -> slot, -1 = no params
  //! FLAT, 2 doubles per slot: [2s+0]=radius(Angstroms) [2s+1]=epsilon(kcal).  Flat for the reason
  //! _LKUniq is flat - Vec0 of a builtin is what the GC analyzer already has precedent for.
  gctools::Vec0<double>                 _NBUniq;
  core::T_sp                            _NBCachedForceField;
  //! Bumped on rebuild AND on invalidation - see _LKGeneration for why a boolean cannot do this.
  size_t                                _NBGeneration = 0;

  /*! NEIGHBOUR LIST, in CSR form: atom i's neighbours are
   *    _Neighbors[_NeighborStart[i] .. _NeighborStart[i+1])
   *
   * Stored as I3 values (coordinateIndexTimes3), because every consumer wants i3 and storing atom
   * indices would put a multiply by 3 in the innermost loop.  _NeighborStart is indexed by ATOM
   * index instead, which is touched once per query atom rather than once per pair.
   *
   * Only atoms in the TARGET MASK appear in the payload; every atom gets a list.  So the mask says
   * what may be FOUND, not what may be asked.  For a blueprint the mask is the backbone - about a
   * sixth of the table - because the other five sixths are fan-out rotamer atoms, and a pair drawn
   * from two different rotamers cannot contribute to a single-scan energy.
   *
   * DIRECTED, not symmetric: slot->backbone entries exist, backbone->slot ones do not.  That is
   * fine because the consumer always iterates from the query side, and it halves the storage.
   *
   * No Atom_sp anywhere - the atom is recoverable as _Atoms[i3/3].atom() when actually needed, so
   * the GC has nothing here to trace.
   *
   * WHY: rebuildPairListImpl grids the whole table and enumerates every pair within rpairlist, per
   * component, per evaluation.  Scoring one rotamer that way examined ~6.5 million candidate pairs
   * to keep ~4500, three times over, 470 times over - with the coordinates unchanged throughout.
   */
  gctools::Vec0<uint32_t>               _NeighborStart;      // length natoms+1
  gctools::Vec0<uint32_t>               _Neighbors;          // i3 values, concatenated
  double                                _NeighborCutoff = -1.0;
  core::T_sp                            _NeighborCoords;     // which coordinate vector
  //! Drift detection for the list, expressed against the SHARED snapshot below rather than against
  //! a private copy of the coordinates.  _MaxDrift at the moment the list was built, plus the epoch
  //! that value is measured in - exactly the pair every pair-list component carries.
  double                                _NeighborDrift = 0.0;
  size_t                                _NeighborEpoch = 0;  // 0 = never built
  core::T_sp                            _NeighborTargetMask; // which atoms may appear as neighbours
  size_t                                _NeighborGeneration = 0;

  // ============================================================
  //  SHARED DRIFT SNAPSHOT
  // ============================================================
  /*! ONE reference copy of the coordinates, for every consumer that needs "has anything moved?".
   *
   * Each pair-list component used to carry its own _DisplacementBuffer - a full copy of the
   * coordinates taken when that component last rebuilt - and walk it on every evaluation.  A
   * blueprint holds ~2829 of them over 6762 atoms: ~460 MB of copies and ~57 million coordinate
   * comparisons per evaluation cycle, every one of them recomputing the same answer.
   *
   * WHY ONE SNAPSHOT IS ENOUGH.  The per-component test was already GLOBAL: it walked the entire
   * coordinate vector, not just the atoms the component enumerates, so it asked "did ANY atom move
   * more than my threshold".  The only per-component parts were WHEN the snapshot was taken and
   * WHAT the threshold is.  Both are scalars.  Collapsing the buffer to one shared copy loses no
   * spatial resolution the old code had.
   *
   * HOW A COMPONENT ASKS ABOUT ITS OWN REBUILD TIME.  It stores _MaxDrift as of its last rebuild.
   * The triangle inequality then bounds what it needs:
   *
   *     |x_now - x_rebuilt| <= |x_now - ref| + |ref - x_rebuilt| <= maxDrift() + pairListDrift()
   *
   * That bound is conservative in the SAFE direction - it can rebuild early, never late.
   *
   * WHY IT STAYS TIGHT.  updateDrift re-snapshots once _MaxDrift exceeds _MaxTriggerSeen, the
   * largest threshold any consumer has asked about.  At that instant every consumer's bound already
   * fires, so they are all rebuilding regardless and invalidating the epoch costs nothing - while
   * resetting the reference to zero puts the conservatism back at zero too.
   *
   * COORDINATES ARE REUSED cycle to cycle, so the identity of the NVector proves nothing: the same
   * object comes back with new numbers in it.  _DriftOwner is therefore a SAFETY NET for callers
   * that bypass the funnel, not the validity test - see updateDrift. */
  core::T_sp                            _RefCoords;          // the one snapshot
  core::T_sp                            _DriftOwner;         // vector updateDrift last ran against
  size_t                                _RefEpoch = 1;       // 0 reserved for "never built"
  double                                _MaxDrift = 0.0;     // max_i |coords_i - ref_i|
  double                                _MaxTriggerSeen = 0.0;
  /*! What the last updateDrift MEASURED, before any re-snapshot reset it - reporting only.
   *
   * _MaxDrift is the live algorithmic value and a re-snapshot zeroes it, so reading it after a big
   * move reports 0.0 and makes a test table say nothing moved on exactly the rows where the most
   * moved.  This keeps the measurement. */
  double                                _LastDrift = 0.0;
  /*! Pair-list rebuilds since the counter was last reset - for tests, not for the algorithm.
   *
   * A drift regression test that never exercises the cached path passes for the wrong reason: if
   * every component rebuilt anyway there was nothing to be stale.  This is what lets the test
   * ASSERT that the fast path was actually taken. */
  size_t                                _PairListRebuilds = 0;
  /*! Skip the per-evaluation drift walk because the caller guarantees COORDS will not change.
   *
   * THE CALLER'S OBLIGATION, and it is not checkable from here: while this is set, updateDrift does
   * not run, so _MaxDrift is frozen at whatever the last walk measured.  A component that is NOT
   * invalidated will therefore read a stale drift and may skip a rebuild it needed.  Set this only
   * around a loop that (a) never mutates the coordinate vector and (b) invalidates every component
   * it evaluates - a pair scan does both, which is why it costs a quarter of that scan's runtime to
   * measure a number that is identically zero and discarded.
   *
   * Named for the promise rather than for the mechanism ("disable drift") so the obligation travels
   * with the call site. */
  bool                                  _CoordinatesPinned = false;
 public:
  /*! Recompute _MaxDrift against COORDS.  O(natoms), and meant to be called ONCE per evaluation.
   *
   * EnergyFunction_O::evaluateAll is the single funnel every evaluation path goes through -
   * evaluateEnergy, evaluateEnergyForce and evaluateEnergyIntoFaRestFaRepVector all route through
   * it - so one call there covers every component, including a blueprint evaluating one group with
   * the rest disabled.  Returns the new _MaxDrift.
   *
   * May RE-SNAPSHOT, which bumps refEpoch and so forces every consumer to rebuild.  That only
   * happens once the drift is past every threshold in play, where they were all rebuilding anyway. */
  double updateDrift(core::T_sp coords);
  /*! Read the drift, recomputing it if COORDS is not what updateDrift last ran against.
   *
   * The owner check catches a caller that never called updateDrift at all.  It CANNOT catch a
   * caller that mutated the same NVector in place and skipped updateDrift - nothing can, short of
   * the walk itself.  Correctness rests on evaluateAll doing the update first. */
  double driftFor(core::T_sp coords);
  /*! Replace the reference with COORDS, zero the drift, bump the epoch. */
  void resnapshotDrift(core::T_sp coords);
  size_t refEpoch() const { return this->_RefEpoch; };
  CL_LISPIFY_NAME("atom-table-ref-epoch");
  CL_DEFMETHOD size_t refEpochForLisp() const { return this->_RefEpoch; };
  double maxDrift() const { return this->_MaxDrift; };
  CL_LISPIFY_NAME("atom-table-last-drift");
  CL_DOCSTRING(R"dx(How far the furthest atom had moved from the shared reference snapshot when the drift was last measured.  Reporting only - unlike the live value this survives a re-snapshot, which zeroes it.)dx");
  CL_DEFMETHOD double lastDrift() const { return this->_LastDrift; };
  /*! Register a threshold, so the re-snapshot policy knows how far it can let the drift run. */
  void noteTrigger(double trigger) {
    if (trigger > this->_MaxTriggerSeen) this->_MaxTriggerSeen = trigger;
  };
  CL_LISPIFY_NAME("set-coordinates-pinned");
  CL_DOCSTRING(R"dx(Promise that the coordinate vector will not change, so the per-evaluation drift walk can be skipped.  The caller must also ensure every component it evaluates is invalidated; see WITH-PINNED-COORDINATES.)dx");
  CL_DEFMETHOD void setCoordinatesPinned(bool on) { this->_CoordinatesPinned = on; };
  CL_LISPIFY_NAME("coordinates-pinned-p");
  CL_DEFMETHOD bool coordinatesPinned() const { return this->_CoordinatesPinned; };
  /*! Throw away the reference so EVERY pair-list consumer must rebuild on its next evaluation.
   *
   * This is the regression-test lever: evaluate once to let the drift logic decide, call this,
   * evaluate again, and any difference is the drift logic having wrongly skipped a rebuild.  It
   * bumps the epoch rather than only nil'ing the snapshot, because the epoch is what components
   * and the neighbour list actually compare against. */
  CL_LISPIFY_NAME("invalidate-drift-snapshot");
  CL_DOCSTRING(R"dx(Force every pair-list component and the neighbour list to rebuild on their next evaluation.)dx");
  CL_DEFMETHOD void invalidateDriftSnapshot() {
    this->_RefCoords = nil<core::T_O>();
    this->_DriftOwner = nil<core::T_O>();
    this->_MaxDrift = 0.0;
    this->_RefEpoch++;
    if (this->_RefEpoch == 0) this->_RefEpoch = 1;
    this->invalidateNeighborList();
  };
  CL_LISPIFY_NAME("atom-table-rebuild-count");
  CL_DOCSTRING(R"dx(Pair-list rebuilds since the count was last reset.  For tests: a drift check whose count never rises never exercised the cached path.)dx");
  CL_DEFMETHOD size_t pairListRebuilds() const { return this->_PairListRebuilds; };
  CL_LISPIFY_NAME("reset-atom-table-rebuild-count");
  CL_DEFMETHOD void resetPairListRebuilds() { this->_PairListRebuilds = 0; };
  void notePairListRebuild() { this->_PairListRebuilds++; };
 public:
  /*! True when the cached neighbour list answers this request.
   *
   * CUTOFF is >=, not ==: a list built at 9.0 is a valid SUPERSET for a component asking 7.5,
   * which applies its own exact distance test per pair anyway.  That is what lets one list serve
   * rosetta nonbond (9.0), LK solvation (8.0) and elec (7.5) from a single build. */
  /*! True when the cached list still answers this request.
   *
   * CUTOFF is >=, not ==: a list built wider is a valid superset, and each component applies its
   * own exact distance test per pair anyway.
   *
   * The SLACK between what was built and what is asked for IS the Verlet skin.  A pair further
   * apart than CUTOFF can only come within it after moving, and two atoms each moving less than
   * half the slack cannot close more than the slack - so the list stays valid while every atom has
   * drifted less than (built - requested)/2.  Build with NEIGHBOR_LIST_SKIN of headroom and this
   * survives many cycles; build at exactly CUTOFF and the slack is zero, so any motion at all
   * forces a rebuild, which is correct but pointless.
   *
   * Reads the SHARED drift - see _RefCoords.  It does no walking of its own, so the caller must
   * have brought the drift up to date for these coordinates first; ensureNeighborList does. */
  bool neighborListValidFor(core::T_sp coords, double cutoff, core::T_sp targetMask) const;
  CL_DEFMETHOD void invalidateNeighborList() {
    this->_NeighborCoords = nil<core::T_O>();
    this->_NeighborTargetMask = nil<core::T_O>();
    this->_NeighborCutoff = -1.0;
    this->_NeighborEpoch = 0;
    this->_NeighborDrift = 0.0;
    this->_NeighborStart.clear();
    this->_Neighbors.clear();
    this->_NeighborGeneration++;
  }
  size_t neighborGeneration() const { return this->_NeighborGeneration; };
  /*! Build the neighbour list unless a valid one is already cached.
   *
   * TARGETMASK is a SIMPLE-BIT-VECTOR indexed by I3 - the same convention and the same object a
   * blueprint already builds for its backbone scope tests - or NIL meaning every atom is a target.
   * COORDS is the NVector the positions come from. */
  void ensureNeighborList(core::T_sp coords, double cutoff, core::T_sp targetMask);
 public:
  bool lkTypeSlotsValidFor(core::T_sp ff) const {
    return this->_LKCachedForceField == ff
        && this->_LKTypeSlot.size() == this->_Atoms.size();
  }
  CL_DEFMETHOD void invalidateLKTypeSlots() {
    this->_LKCachedForceField = nil<core::T_O>();
    this->_LKTypeSlot.clear();
    this->_LKUniq.clear();
    this->_LKGeneration++;
  }
  bool nbTypeSlotsValidFor(core::T_sp ff) const {
    return this->_NBCachedForceField == ff
        && this->_NBTypeSlot.size() == this->_Atoms.size();
  }
  CL_DEFMETHOD void invalidateNBTypeSlots() {
    this->_NBCachedForceField = nil<core::T_O>();
    this->_NBTypeSlot.clear();
    this->_NBUniq.clear();
    this->_NBGeneration++;
  }
 public:
  typedef gctools::Vec0<EnergyAtom>::iterator iterator;
public:
  static AtomTable_sp make();
  
 public:
  gctools::Vec0<EnergyAtom>&	getVectorEnergyAtoms() { return this->_Atoms;};
  CL_DEFMETHOD size_t	getNumberOfAtoms()	{ return this->_Atoms.size();};
  CL_DEFMETHOD size_t   getNumberOfMolecules()  { return this->_AtomsPerMolecule->length();};
  core::T_mv getAtomIndexOrNil(Atom_sp atom);

  core::List_sp atomsAtRemoveBondAngle14(Atom_sp atom) const;
  core::T_sp   firstSolventMoleculeNSPSOL() const;
  bool firstSolventMoleculeNSPSOL_boundP() const;
  void set_firstSolventMoleculeNSPSOL(size_t num);
  void makUnbound_firstSolventMoleculeNSPSOL();

  core::T_sp   finalSoluteResidueIPTRES() const;
  bool finalSoluteResidueIPTRES_boundP() const;
  void set_finalSoluteResidueIPTRES(size_t num);
  void makUnbound_finalSoluteResidueIPTRES();

  core::T_sp   totalNumberOfMoleculesNSPM() const;
  bool totalNumberOfMoleculesNSPM_boundP() const;
  void set_totalNumberOfMoleculesNSPM(size_t num);
  void makUnbound_totalNumberOfMoleculesNSPM();

  core::T_sp   aggregateName() const;
  bool aggregateNameBoundP() const;
  void setAggregateName(core::T_sp name);
  void makUnboundAggregateName();

  core::T_sp   boundingBox() const;
  bool boundingBoxBoundP() const;
  void setBoundingBox(core::T_sp name);
  void makUnboundBoundingBox();

  core::T_sp   nonbondForceFieldForAggregate() const;
  bool nonbondForceFieldForAggregateBoundP() const;
  void setNonbondForceFieldForAggregate(core::T_sp forceField);
  void makUnboundNonbondForceFieldForAggregate();

  core::T_sp   lksolvationForceFieldForAggregate() const;
  bool lksolvationForceFieldForAggregateBoundP() const;
  void setLKSolvationForceFieldForAggregate(core::T_sp forceField);
  void makUnboundLKSolvationForceFieldForAggregate();

  uint	getNVectorSize() const { return this->_Atoms.size()*3;};
  // Read coordinates from atoms into an NVector of coordinates
  void readAtomCoordinates(NVector_sp coords) const;
  // Write coordinates from atoms into an NVector of coordinates
  void writeAtomCoordinates(NVector_sp coords);
  
  EnergyAtom*	getEnergyAtomPointer(Atom_sp a);
  EnergyAtom&	energyAtomEntry(uint i) { return this->_Atoms[i]; };
  EnergyAtom*	findEnergyAtomWithCoordinateIndex(uint i) { return &this->_Atoms[i/3]; };

  
  void	add(EnergyAtom& aa) { SIMPLE_ERROR("AtomTable->add(EnergyAtom may cause recursive allocation"); this->_Atoms.push_back(aa);};

	/*! Add the info for a single atom and return the coordinate index (index*3) of the atom */
  int addAtomInfo(Atom_sp atom, units::Quantity_sp charge, units::Quantity_sp mass, int typeIndex, uint atomicNumber );

  iterator	begin() { return this->_Atoms.begin(); };
  iterator	end() { return this->_Atoms.end(); };

  void setAtomFlag(size_t index, size_t flag);
  size_t getAtomFlag(size_t index);

  void	dumpTerms(core::HashTable_sp atomTypes);
  void constructFromMolecule(Molecule_sp mol, core::T_sp nonbondForceField, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes );

  size_t push_back_excluded_atom_indexes_and_sort(core::ComplexVector_int32_t_sp excludedAtomIndexes, size_t atomIndex);
  /*! Calculate the excluded atom list with atom indexes starting at 0.
      Atoms with no excluded atoms get a -1 in the excluded_atom_list.
      This is different from AMBER, where it starts counting at 1.
      To generate an AMBER excluded atom list, add 1 to each entry.
      Return Values(number_excluded_atoms,excluded_atom_list).
  */
  CL_DEFMETHOD core::T_mv calculate_excluded_atom_list(core::T_sp keepInteractionsFactory);

  // ------------------------------------------------------------
  //
  // Access fields in AtomTable elements directly
  //
  core::HashTable_sp getAtomTableIndexes();
  CL_DEFMETHOD MatterName elt_atom_name(int index) { return this->_Atoms[index]._AtomName; };
  CL_DEFMETHOD core::T_sp elt_atom_type(int index,core::HashTable_sp atomTypes);
  CL_DEFMETHOD Atom_sp elt_atom(int index) { return this->_Atoms[index].atom(); };
  CL_DEFMETHOD int elt_atom_coordinate_index_times3(int index) { return this->_Atoms[index].coordinateIndexTimes3(); };
  CL_DEFMETHOD double elt_charge(int index) { return this->_Atoms[index]._Charge; };
  CL_DEFMETHOD double elt_mass(int index) { return this->_Atoms[index]._Mass; };
  CL_DEFMETHOD int elt_type_index(int index) { return this->_Atoms[index]._TypeIndex; };
  CL_DEFMETHOD int elt_atomic_number(int index) { return this->_Atoms[index]._AtomicNumber; };
  core::T_sp atom_table_residue_pointers() const;
  core::T_sp atom_table_residue_names() const;
  core::T_sp atom_table_atoms_per_molecule() const;
  core::T_sp atom_table_residues() const;
  core::T_sp atom_table_molecules() const;
  CL_DEFMETHOD void setf_atom_table_residue_pointers(core::ComplexVector_int32_t_sp val) {this->_ResiduePointers = val; };
  CL_DEFMETHOD void setf_atom_table_residue_names(core::ComplexVector_T_sp val) {this->_ResidueNames = val; };
  CL_DEFMETHOD void setf_atom_table_atoms_per_molecule(core::ComplexVector_int32_t_sp val) {this->_AtomsPerMolecule = val; };
  CL_DEFMETHOD void setf_atom_table_residues(core::Vector_sp val) {this->_Residues = val; };
  CL_DEFMETHOD void setf_atom_table_molecules(core::Vector_sp val) {this->_Molecules = val; };
  size_t getCoordinateIndexTimes3(Atom_sp atom);
  size_t getCoordinateIndexTimes3ForAtomAtIndex(size_t index);
  
 AtomTable_O() : _ResiduePointers(unbound<core::ComplexVector_int32_t_O>()),
                 _ResidueNames(unbound<core::ComplexVector_T_O>()),
                 _AtomsPerMolecule(unbound<core::ComplexVector_int32_t_O>()),
                 _firstSolventMoleculeNSPSOL(unbound<core::T_O>()),
                 _finalSoluteResidueIPTRES(unbound<core::T_O>()),
                 _totalNumberOfMoleculesNSPM(unbound<core::T_O>()),
                 _Residues(unbound<core::Vector_O>()),
                 _Molecules(unbound<core::Vector_O>()),
                 _AggregateName(unbound<core::T_O>()),
                 _BoundingBox(unbound<core::T_O>()),
                 _NonbondForceFieldForAggregate(unbound<core::T_O>()),
                 _LKSolvationForceFieldForAggregate(unbound<core::T_O>()),
                 // NIL, not unbound: lkTypeSlotsValidFor compares it with == against a real force
                 // field, so it only has to be something no force field is equal to.
                 _LKCachedForceField(nil<core::T_O>()),
                 _NBCachedForceField(nil<core::T_O>()),   // NIL for the same reason
                 _NeighborCoords(nil<core::T_O>()),
                 _NeighborTargetMask(nil<core::T_O>()),
                 _RefCoords(nil<core::T_O>()),
                 _DriftOwner(nil<core::T_O>())
  {};

  virtual void fill_atom_table_from_vectors(core::List_sp values);
//  int residue_index(int atom_index);
};


};


#endif //]
