/*
    File: energyFunction.h
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
       
       
#define	DEBUG_LEVEL_NONE


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyFunction.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyFunction_H
#define	EnergyFunction_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/hashTableEq.h>
#include <cando/geom/vector3.h>
#include <cando/chem/scoringFunction.h>
#include <cando/chem/energyComponent.fwd.h>
// Full definition, not just the fwd - ensureComponent's group branch below dereferences the
// group's _Components and calls addComponent.  No cycle: energyComponent.h does not include
// energyFunction.h.
#include <cando/chem/energyComponentGroup.h>
#include <cando/chem/atom.h>




#include <cando/adapt/quickDom.fwd.h>// energyFunction.h wants QDomNode needs quickDom.fwd.h
//#include "geom/render.fwd.h"// energyFunction.h wants DisplayList needs render.fwd.h
#include <clasp/core/iterator.fwd.h>// energyFunction.h wants Iterator needs iterator.fwd.h


#include <cando/chem/chemPackage.h>

namespace       chem
{
  FORWARD(FFStretchDb);
  FORWARD(FFAngleDb);
  FORWARD(FFPtorDb);
  FORWARD(FFItorDb);
  FORWARD(FFNonbondDb);
  
  class EnergyAtom;

  FORWARD(FFParameter);
  FORWARD(EnergyRosettaNonbond);
  FORWARD(AbstractLargeSquareMatrix);
  FORWARD(FFNonbondCrossTermTable);
  FORWARD(QDomNode);
  FORWARD(Atom);
  FORWARD(Matter);
  FORWARD(ForceField);
  FORWARD(CombinedForceField);
  FORWARD(AtomTable);
  FORWARD(EnergyStretch);
  FORWARD(EnergyAngle);
  FORWARD(EnergyDihedral);
  FORWARD(EnergyNonbond);
  FORWARD(EnergyPeriodicBoundaryConditionNonbond);
  FORWARD(EnergyDihedralRestraint);
  FORWARD(EnergyChiralRestraint);
  FORWARD(EnergyAnchorRestraint);
  FORWARD(EnergyFixedNonbondRestraint);
  FORWARD(Iterator);
  FORWARD(ForceMatchReport);
FORWARD(BoundingBox);


core::T_sp specializeKeepInteractionFactory( core::T_sp keepInteractionFactory, core::T_sp aclass );

bool skipInteraction_EnergyNonbond( core::T_sp keepInteractionFunction,
                                   Atom_sp atomA,
                                   Atom_sp atomB,
                                   core::T_sp ti3A,
                                   core::T_sp ti3B );

bool skipInteraction( core::T_sp keepInteractionFunction,
                      Atom_sp atomA = unbound<Atom_O>(),
                      Atom_sp atomB = unbound<Atom_O>(),
                      Atom_sp atomC = unbound<Atom_O>(),
                      Atom_sp atomD = unbound<Atom_O>() );
 /*! Throw this exception from Cando if you want the minimizer to deal with a bad interaction */
  class	InteractionCondition 
  {
  public:
    core::Symbol_sp  condition_name;
    core::List_sp    condition_arguments;
  };

  FORWARD(EnergyFunction);
};


template <>
struct gctools::GCInfo<chem::EnergyFunction_O> {
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {
  FORWARD(EnergyFunction);
  class EnergyFunction_O : public ScoringFunction_O
  {
    LISP_CLASS(chem,ChemPkg,EnergyFunction_O,"EnergyFunction",ScoringFunction_O);
  public:
    static EnergyFunction_sp make(Aggregate_sp matter,
                                  core::T_sp disableComponents=nil<core::T_O>(),
                                  core::List_sp enableComponents=nil<core::T_O>(),
                                  bool useExcludedAtoms=false,
                                  core::T_sp keepInteractionFactory=nil<core::T_O>(),
                                  bool assign_types=false,
                                  core::T_sp forceFieldOverride=nil<core::T_O>(),
                                  core::T_sp setup = nil<core::T_O>());
  public:
    void initialize();
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    Matter_sp				_Matter;	// Aggregate or Molecule
    /*! Stores cross terms for evaluating nonbond interactions
     */
    FFNonbondCrossTermTable_sp		_NonbondCrossTermTable;
    AtomTable_sp			_AtomTable;
    core::List_sp                       _EnergyComponents; // list of EnergyComponent_sp objects
    BoundingBox_sp                      _BoundingBox;
    /*! If true then secondary amides are
     * automatically restrainted to be trans
     */
    bool				_RestrainSecondaryAmides;
    core::T_sp                          _Message;
    core::Symbol_sp                     _ForceFieldName; // :amber or :rosetta
    /*! The FULL setup list that defineForAggregate was called with, kept verbatim.
     * _ForceFieldName above holds only its head.  generateIntoGroup has to re-run
     * generateNonbondEnergyFunctionTables for one rotamer slot, and the rosetta branch builds
     * SetupAccumulator(classSymbol,setup) from the WHOLE list - reconstructing (list
     * _ForceFieldName) would drop the tail and silently reparameterize.  _UseExcludedAtoms is
     * kept for the same reason: same call, same arguments.
     */
    core::T_sp                          _Setup;
    bool                                _UseExcludedAtoms;
    core::List_sp			_MissingParameters;
  public:
    void	_eraseMissingParameters() { this->_MissingParameters = nil<core::T_O>();};
    void	_addMissingParameter(FFParameter_sp p) { this->_MissingParameters = core::Cons_O::create(p,this->_MissingParameters);};
    void __createSecondaryAmideRestraints(VectorAtom& nitrogens, core::T_sp keepInteractionFactory );

    void	flagDihedralRestraintsAboveThreshold(NVector_sp nvPosition);
    core::T_mv enabledDisabled() const;
  private:
    int _applyRestraints(core::T_sp forceField, core::Iterator_sp restraintIterator, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes );
    void _addDihedralRestraintDegrees(EnergyDihedralRestraint_sp rest, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double phi0_degrees, double kdh, core::T_sp keepInteractionFactory );

  public:

    CL_LISPIFY_NAME("atomTable");
    CL_DEFMETHOD     AtomTable_sp atomTable() const { return this->_AtomTable;};

    core::List_sp allComponents() const;

    string	energyTermsEnabled() ;
    void	loadCoordinatesIntoVector(NVector_sp pos);
    void	saveCoordinatesFromVector(NVector_sp pos);
    void	saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force);
    size_t	getNVectorSize() const;
    double	evaluateRaw( NVector_sp pos, NVector_sp force );

    BoundingBox_sp boundingBox() const;
    bool boundingBoxBoundP() const;
    void setBoundingBox(BoundingBox_sp bounding_box);
    void makUnboundBoundingBox();

    void resolveMoleculeForceFields(Aggregate_sp aggregate,
                                    core::T_sp forceFieldOverrides,
                                    core::HashTable_sp molecule_force_fields,
                                    core::HashTable_sp molecule_force_field_names);

    /*! BONDED NIL skips the per-molecule pass entirely and generates only the nonbond components.
     *
     * Refusing the bonded classes in the keep-interaction-factory is NOT the same thing.  That
     * discards the TERMS, but generateForMolecule has already run the full force-field
     * parameterization by then - so a caller wanting only nonbond re-parameterizes every molecule
     * to produce zero bonded terms.  For a blueprint pair group that is pure waste, and worse: the
     * parameterization runs before any predicate is consulted, so a fully attached spiro NG with
     * 16 bonds trips SMIRNOFF's bond-order check and kills the run. */
    void generateIntoGroup(core::T_sp keepInteractionFactory, core::T_sp group, core::T_sp forceFieldOverrides, bool bonded=true, core::T_sp scopeAggregate=nil<core::T_O>());

    ForceMatchReport_sp checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, core::T_sp energyScale, NVector_sp force, core::T_sp activeAtomMask );

    CL_LISPIFY_NAME("getMatter");
    CL_DEFMETHOD     Matter_sp	getMatter() { return this->_Matter;};

    void	useDefaultSettings();


    /*! Set a single options */
    void	setOption( core::Symbol_sp option, core::T_sp val);


    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
    void	setOptions( core::List_sp options );

    CL_DEFMETHOD core::Symbol_sp forceFieldName() { return this->_ForceFieldName; };

    CL_DEFMETHOD core::T_sp findComponentOrNil(core::Symbol_sp className);
    CL_DEFMETHOD EnergyComponent_sp findComponent(core::Symbol_sp className);

    void addComponentIfMissing(EnergyComponent_sp comp);
    void ensureBaseComponents();

    CL_LISPIFY_NAME("getStretchComponent");
    CL_DEFMETHOD     EnergyStretch_sp	getStretchComponent();
    CL_LISPIFY_NAME("getAngleComponent");
    CL_DEFMETHOD     EnergyAngle_sp	getAngleComponent();
    CL_LISPIFY_NAME("getDihedralComponent");
    CL_DEFMETHOD EnergyDihedral_sp getDihedralComponent();
    CL_LISPIFY_NAME("getNonbondComponent");
    CL_DEFMETHOD EnergyNonbond_sp getNonbondComponent();
    CL_LISPIFY_NAME("getChiralRestraintComponent");
    CL_DEFMETHOD EnergyChiralRestraint_sp getChiralRestraintComponent();
    CL_LISPIFY_NAME("getAnchorRestraintComponent");
    CL_DEFMETHOD EnergyAnchorRestraint_sp getAnchorRestraintComponent();
    CL_LISPIFY_NAME("getDihedralRestraintComponent");
    CL_DEFMETHOD EnergyDihedralRestraint_sp getDihedralRestraintComponent();
    CL_LISPIFY_NAME("getFixedNonbondRestraintComponent");
    CL_DEFMETHOD EnergyFixedNonbondRestraint_sp getFixedNonbondRestraintComponent();

    core::List_sp allEnergyComponents() const;
    void          pushEnergyComponent(EnergyComponent_sp component);
    bool          removeEnergyComponent(EnergyComponent_sp component);

    CL_DEFMETHOD bool hasMissingParameters();
    CL_DEFMETHOD core::List_sp getMissingParameters();

    void	setupHessianPreconditioner( NVector_sp pos, AbstractLargeSquareMatrix_sp hessian,
                                            core::T_sp activeAtomMask );

    /*! Enable debugging on all energy components
     */
    void	enableDebug();
    /*! Disable debugging on all energy components
     */
    void	disableDebug();

    void	summarizeTerms();
    void	dumpTerms();
    CL_DEFMETHOD     core::T_sp	getMessage() { return this->_Message;};

    int	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos);
    int	compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition();


    void		writeForceToAtoms(NVector_sp f);
    EnergyAtom*     getEnergyAtomPointer(Atom_sp a);

    void assignAtomTypes(Matter_sp matter, bool show_progress);
    void defineForAggregate(Aggregate_sp agg, bool useExcludedAtoms, core::T_sp keepInteractionFactory=nil<core::T_O>(), bool assign_types=true, core::T_sp forceFieldOverrides=nil<core::T_O>(), core::T_sp setup=nil<core::T_O>(), core::T_sp maybe_energy_component_group=nil<core::T_O>() );
    /*! PASS 1 - assign MOLECULE's atom types and append its atoms to the atom table. */
    void constructAtomTableForMolecule(Molecule_sp molecule, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes, core::T_sp nonbondForceField, core::HashTable_sp force_fields, core::HashTable_sp force_field_names );
    /*! PASS 2 - generate MOLECULE's bonded terms into GROUP (or into the energy function's
        own components when GROUP is NIL).  The atom table is already complete. */
    void generateForMolecule(Molecule_sp molecule, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes, core::T_sp nonbondForceField, core::HashTable_sp force_fields, core::HashTable_sp force_field_names, core::T_sp group=nil<core::T_O>() );
    void generateStandardEnergyFunctionTables(Matter_sp mol,
                                              FFStretchDb_sp stretchDb,
                                              FFAngleDb_sp angleDb,
                                              FFPtorDb_sp ptorDb,
                                              FFItorDb_sp itorDb,
                                              core::T_sp keepInteractionFactory,
                                              core::HashTable_sp atomTypes,
                                              core::T_sp group=nil<core::T_O>());
    /*! GROUP NIL puts the nonbond components on the energy function, which is what
        defineForAggregate wants.  GROUP non-NIL collects them there instead - see
        collectComponent - which is how one rotamer slot gets its own nonbond terms. */
    /*! BONDED14 false suppresses the 1-4 component.
     *
     * A 1-4 interaction is a nonbond ENERGY whose MEMBERSHIP comes from the bond graph:
     * construct14InteractionTerms loops PROPERS over the matter and takes each proper dihedral's
     * terminal atoms.  Same Loop, same matter, same source as the dihedral terms - so it is
     * topology-determined in exactly the way the bonded terms are, and a caller that skipped the
     * bonded pass has no business paying for it either.
     *
     * Defaults true so defineForAggregate is unaffected. */
    void generateNonbondEnergyFunctionTables(bool useExcludedAtoms, Matter_sp agg, core::T_sp forceField, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes, core::T_sp setup, core::T_sp group=nil<core::T_O>(), bool bonded14=true );
    void generateRestraintEnergyFunctionTables(Matter_sp agg, core::T_sp nonbonds, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes );

    /*! Add the restraints to the energy function.
     * This allows restraints to be applied to the system
     * without having to add them to the molecule/aggregate.
     */
    void	addTermsForListOfRestraints( ForceField_sp forceField,  core::List_sp restraintList, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes );


    double	calculateNumericalDerivative(NVector_sp pos, core::T_sp energyScale, double delta, uint i, core::T_sp activeAtomMask );
    double	calculateNumericalSecondDerivative(NVector_sp pos, core::T_sp energyScale, double delta, uint i, uint j, core::T_sp activeAtomMask );
    double	evaluateAll(NVector_sp pos,
                            core::T_sp energyScale,
                            core::T_sp energyComponents,
                            bool calcForce,
                            gc::Nilable<NVector_sp> force,
                            bool calcDiagonalHessian,
                            bool calcOffDiagonalHessian,
                            gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                            gc::Nilable<NVector_sp> hdvec,
                            gc::Nilable<NVector_sp> dvec,
                            core::T_sp activeAtomMask,
                            core::T_sp debugInteractions,
                            bool disableRestraints );

    // faRepComponent is the ef's own EnergyRosettaNonbond, or NIL when this ef has no
    // rosetta-nonbond component (e.g. a ligand-only system) — then fa_rep is 0.
    void evaluateEnergyIntoFaRestFaRepVector(NVector_sp pos, NVector_sp faRestFaRepVector, size_t index, core::T_sp faRepComponent,
                                             core::T_sp energyScale, core::T_sp activeAtomMask, core::T_sp debugInteractions,
                                             bool disableRestraints ) ;

    /*! As above, but evaluate ONE NAMED component instead of scanning the component list.
     *
     * evaluateAll walks _EnergyComponents - a CONS LIST - testing isEnabled() and the virtual
     * restraintp() on every element.  A blueprint holds ~2829 components across ~471 groups and a
     * pair scan enables exactly one group, so every evaluation chases 2829 heap-scattered cells to
     * find the one that is on.  That search measured 84% of a scan against 12% for the physics it
     * was searching for.
     *
     * The component's own isEnabled() still gates it (energyComponentGroup.cc:181), so this changes
     * WHICH components are considered, not whether the named one runs.
     *
     * No DISABLE-RESTRAINTS: it only means anything while scanning a list, and here the caller has
     * named the component outright - including, if it likes, a restraint. */
    void evaluateComponentIntoFaRestFaRepVector(core::T_sp component, NVector_sp pos,
                                                NVector_sp faRestFaRepVector, size_t index,
                                                core::T_sp faRepComponent,
                                                core::T_sp energyScale, core::T_sp activeAtomMask,
                                                core::T_sp debugInteractions ) ;

    size_t runTestCalls(core::T_sp stream, NVector_sp pos);

    string	summarizeEnergyAsString();


    void	evaluateNumericalForce(NVector_sp pos, core::T_sp energyScale, NVector_sp numForce, double delta, core::T_sp activeAtomMask );
    void	evaluateNumericalHessian(NVector_sp pos, core::T_sp energyScale, AbstractLargeSquareMatrix_sp numHessian, bool calcOffDiagonalElements, double delta, core::T_sp activeAtomMask);

    string	debugLogAsString();

    void	dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments);

    EnergyFunction_sp copyFilter(core::T_sp keepInteractionFactory, core::List_sp setup, core::List_sp addComponents );

    EnergyFunction_O(BoundingBox_sp bounding_box) :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        ,_EnergyComponents(nil<core::T_O>())
        ,_BoundingBox(bounding_box)
        ,_ForceFieldName(nil<core::Symbol_O>())
        ,_Setup(nil<core::T_O>())
        ,_UseExcludedAtoms(false)
//      , _MissingParameters(unbound<core::List_O>())
    {};

    EnergyFunction_O() :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        ,_EnergyComponents(nil<core::T_O>())
        ,_BoundingBox(unbound<BoundingBox_O>())
        ,_ForceFieldName(nil<core::Symbol_O>())
        ,_Setup(nil<core::T_O>())
        ,_UseExcludedAtoms(false)
//      , _MissingParameters(unbound<core::List_O>())
    {};
    EnergyFunction_O( const EnergyFunction_O& ef ) :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        ,_EnergyComponents(nil<core::T_O>())
        ,_BoundingBox(unbound<BoundingBox_O>())
        ,_ForceFieldName(nil<core::Symbol_O>())
        ,_Setup(nil<core::T_O>())
        ,_UseExcludedAtoms(false)
//      , _MissingParameters(unbound<core::List_O>())
    {};
  };

    bool inAtomSet(core::T_sp atomSet, Atom_sp atom);

#define	FINITE_DIFFERENCE_TOLERANCE 0.10
  int	_areValuesClose(double numVal, double analVal, const char* funcName, const char* termName, int index );

  void energyFunction_initializeSmarts();

/*! Find-or-create the one component of class Type for this energy function.

    GROUP scopes the lookup.  NIL - the default, and every caller today - searches the energy
    function's own components, so there is exactly one EnergyStretch, one EnergyAngle and so on
    per energy function.  That is the right answer for the backbone.

    A blueprint fan-out needs the opposite: one set of components PER ROTAMER SLOT, so a slot's
    energy can be evaluated on its own.  Passing a GROUP searches only within that group, so a
    generation pass driven with a fresh group produces a fresh set of components - one stretch,
    one angle, one dihedral - holding only that pass's terms.  The group is both the lookup
    scope and the result: it is what gets stored on the ROTAMER-SCAN.

    GROUP is an EnergyComponentGroup_sp - itself an EnergyComponent_O, so it lives in the energy
    function's component list like any other and forwards evaluation to its children.

    Note the scoping is per-GROUP, not per-CALL.  ensureComponent is invoked once per TERM, so a
    fresh component per call would give one component per term.

    A grouped component is NOT also pushed onto the energy function.  The group itself is pushed,
    once, and forwards evaluation to its children; pushing both double-counts every term, and
    pushing neither makes enable/disable a no-op.  Once any group exists, findComponentOrNil is no
    longer meaningful for the grouped classes - it returns the FIRST, which is the backbone's - so
    consumers of slot components must hold direct references. */
template <typename Type>
gc::smart_ptr<Type> ensureComponent(EnergyFunction_sp mthis, core::T_sp group = nil<core::T_O>()) {
  if (group.notnilp()) {
    EnergyComponentGroup_sp egroup = gc::As<EnergyComponentGroup_sp>(group);
    for (auto& c : egroup->_Components) {
      if (gc::IsA<gc::smart_ptr<Type>>(c)) return gc::As<gc::smart_ptr<Type>>(c);
    }
    auto new_comp = Type::create();
    egroup->addComponent(new_comp);
    return new_comp;
  }
  auto comp = mthis->findComponentOrNil(Type::static_classSymbol());
  if (comp.nilp()) {
    auto new_comp = Type::create();
    mthis->pushEnergyComponent(new_comp);
    return new_comp;
  }
  return gc::As<gc::smart_ptr<Type>>(comp);
}

/*! Route an ALREADY-CONSTRUCTED component to a group, or to the energy function when GROUP is NIL.

    The nonbond counterpart of ensureComponent, and deliberately NOT a find-or-create: nonbond
    components are built fresh at each generation pass (constructNonbondTermsFromAtomTable fills a
    new one), so there is nothing to look up.  A slot wants its own EnergyNonbond regardless.

    Same group contract as ensureComponent: the child goes in the group and is NOT also pushed onto
    the energy function.  The group is pushed once, by the caller. */
inline void collectComponent(EnergyFunction_sp mthis, EnergyComponent_sp component, core::T_sp group) {
  if (group.notnilp()) {
    gc::As<EnergyComponentGroup_sp>(group)->addComponent(component);
    return;
  }
  mthis->pushEnergyComponent(component);
}

};

#endif
