/*
    File: energyFunction.cc
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
#define	DEBUG_LEVEL_FULL
//
// (C) 2004 Christian E. Schafmeister
//

/*
__BEGIN_DOC(classes.EnergyFunction,section,EnergyFunction)
EnergyFunction objects are used to calculate Amber or Gaff molecular
mechanics energies and first and second derivatives.
__END_DOC
*/

/*
 *	energyFunction.cc
 *
 */

#include <clasp/core/common.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/energyComponent.fwd.h>
#include <cando/chem/loop.h>
#include <cando/adapt/indexedObjectBag.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/minimizerLog.h>
#include <cando/chem/restraint.h>
#include <cando/chem/iterateRestraints.h>
#include <cando/adapt/iterateCons.h>
#include <cando/chem/ringFinder.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/atom.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyStretch.h>
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyLinearAngle.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyPeriodicBoundaryConditionsNonbond.h>
#include <cando/chem/energyDihedralRestraint.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyRosettaNonbond.h>
#include <cando/chem/energyRosettaElec.h>
#include <cando/chem/energyRosettaLKSolvation.h>
#include <cando/chem/energyNonbond14.h>
#include <cando/chem/energyFixedNonbond.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/ffItorDb.h>
#include <cando/chem/ffPtorDb.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/ffVdwDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/chemInfo.h>
#include <clasp/core/wrappers.h>


#if 0
#define	TURN_DEBUG_OFF
#define	MATHEMATICA_DEBUG_ON	1
#define	PARMDEBUG	1
#endif


namespace chem
{

SYMBOL_EXPORT_SC_(ChemPkg,nonbondTerm);
SYMBOL_EXPORT_SC_(ChemPkg,restraintAnchor);
SYMBOL_EXPORT_SC_(ChemPkg,parameterizationError);
SYMBOL_EXPORT_SC_(ChemPkg,missingType);
SYMBOL_EXPORT_SC_(ChemPkg,missingStretch);
SYMBOL_EXPORT_SC_(ChemPkg,missingAngle);
SYMBOL_EXPORT_SC_(ChemPkg,missingDihedral);
SYMBOL_EXPORT_SC_(ChemPkg,missingTypeError);
SYMBOL_EXPORT_SC_(ChemPkg,missingStretchError);
SYMBOL_EXPORT_SC_(ChemPkg,missingAngleError);
SYMBOL_EXPORT_SC_(ChemPkg,missingDihedralError);
SYMBOL_EXPORT_SC_(ChemPkg,skipMissingParameters);
SYMBOL_EXPORT_SC_(ChemPkg,define_for_molecule_using_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,construct_atom_table_for_molecule);
SYMBOL_EXPORT_SC_(ChemPkg,generate_for_molecule_using_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,with_perception);
SYMBOL_EXPORT_SC_(ChemPkg,amber_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,gaff_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,wrapped_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,construct_from_molecule);


//Maybe use this in place of InteractionError?
SYMBOL_EXPORT_SC_(ChemPkg,InteractionError);
/*! LinearAngleError has a argument :atoms (atom1 atom2 atom3) */
SYMBOL_EXPORT_SC_(ChemPkg,LinearAngleError);

SYMBOL_EXPORT_SC_(ChemPkg,LinearDihedralError);
SYMBOL_EXPORT_SC_(ChemPkg,LinearDihedralRestraintError);
SYMBOL_EXPORT_SC_(ChemPkg,OverlappingNonbondError);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);

bool energyFunctionInitialized = false;



void EnergyFunction_O::pushEnergyComponent(EnergyComponent_sp component) {
  // EnergyComponentGroups are deliberately MANY - one per blueprint rotamer slot, ~470 on the myc
  // system - so the one-per-class rule cannot apply to them.  It is narrowed rather than dropped
  // because it still catches a real bug for the singleton components: two EnergyStretch in one
  // energy function means every bonded term is counted twice, which shows up as a plausible
  // number rather than a failure.
  //
  // Skipping the scan for groups also keeps pushing N groups linear.  The check is O(n) per push,
  // so 470 groups through it would be ~110,000 comparisons.
  if (!gc::IsA<EnergyComponentGroup_sp>(component)) {
    auto componentClass = component->_instanceClass();
    for ( auto cur : this->_EnergyComponents ) {
      EnergyComponent_sp existing = gc::As<EnergyComponent_sp>(oCar(cur));
      if (existing->_instanceClass() == componentClass) {
        SIMPLE_ERROR("Attempted to add a duplicate EnergyComponent of class {}", _rep_(componentClass->_className()));
      }
    }
  }
  this->_EnergyComponents = core::Cons_O::create(component,this->_EnergyComponents);
}

CL_LISPIFY_NAME("remove-energy-component");
CL_DOCSTRING(R"dx(Drop COMPONENT from this energy function's component list.  Returns T if it was
  there.

  For scratch components - a group built only to be compared against another and then thrown away.
  Left in place such a group is evaluated by every later whole-function walk, which for a blueprint
  slot group means one rotamer's energy silently added twice.)dx")
CL_DEFMETHOD bool EnergyFunction_O::removeEnergyComponent(EnergyComponent_sp component) {
  core::List_sp kept = nil<core::T_O>();
  bool found = false;
  for ( auto cur : this->_EnergyComponents ) {
    core::T_sp one = oCar(cur);
    if (one == component) { found = true; continue; }
    kept = core::Cons_O::create(one, kept);
  }
  // The list is built by pushing, so it is already reverse-of-insertion; rebuilding by pushing
  // reverses it again.  Nothing depends on the order - evaluateAll sums - but preserving it keeps
  // a before/after inspection readable.
  core::List_sp restored = nil<core::T_O>();
  for ( auto cur : kept ) restored = core::Cons_O::create(oCar(cur), restored);
  this->_EnergyComponents = restored;
  return found;
}

int	_areValuesClose( double numVal, double analVal, const char* funcName, const char* termName, int index )
{
  double	rel = 0.0;
  if ( fabs(numVal)<FINITE_DIFFERENCE_TOLERANCE && fabs(analVal)<FINITE_DIFFERENCE_TOLERANCE ) {
    return true; // both close to zero
  }
  if ( numVal*analVal > 0.0 ) {
    double avg = (fabs(numVal)+fabs(analVal));
    rel = (fabs(numVal-analVal)/avg);
    if ( rel<FINITE_DIFFERENCE_TOLERANCE ) {
      return true;	// they are within a few percent of each other
    }
  }
  return false;
}

SYMBOL_EXPORT_SC_(ChemPkg,STARsecondaryAmideSmartsSTAR);
void energyFunction_initializeSmarts()
{
  if ( energyFunctionInitialized ) {
    printf("%s:%d:%s  energyFunction_initializeSmarts was already called\n", __FILE__, __LINE__, __FUNCTION__ );
    // return;
  }

  _sym_STARparserNodeHolderSTAR->defparameter(adapt::IndexedObjectBag_O::create());
  energyFunctionInitialized = true;
  SmartsRoot_sp secondaryAmide = chem__compile_smarts("[$([N:1]([C:2])(~[#1:3])~[C:4](=[O:5])[C:6])]",nil<core::T_O>());
  _sym_STARsecondaryAmideSmartsSTAR->defparameter(secondaryAmide);
}


core::List_sp EnergyFunction_O::allComponents() const {
  core::List_sp result = this->_EnergyComponents;
  return result;
}

CL_DOCSTRING(R"doc(Create an energy-scale object for an energy-function.)doc");
CL_LISPIFY_NAME(make_energy_scale);
CL_DEF_CLASS_METHOD EnergyScale_sp EnergyScale_O::make()
{
  auto energyScale = gctools::GC<EnergyScale_O>::allocate();
  return energyScale;
}



CL_DOCSTRING(R"doc(Create an energy function for the matter.
: disable-components - NIL (default), T or a list of component class names.  NIL means don't disable any components; T means disable all components;
and a list means disable those components. Use enable-components to in a second step enable components.
: enable-components - NIL (default) or a list of component class names. NIL means enable no components and a list of component class names enables
those components.
: use-excluded-atoms - T (default) uses excluded atoms for nonbond terms and NIL uses pairwise terms. Only small molecules can use pairwise terms.
: keep-interaction-factory - T (default) means keep all interactions; NIL means no interactions (an empty energy function that
evaluates nothing) and (lambda (component-class)) returns a (lambda (&rest atoms) ...) or specialized versions that take a fixed number of atom arguments
that will be called for each interaction and the function returns T or NIL if
each interaction should be added to the energy function.
: assign-types - T (default) assign atom types as part of generating the energy function.  [I don't know what will happen if assing-types is NIL.)doc");
CL_LAMBDA(&key matter disable-components enable-components (use-excluded-atoms nil) (keep-interaction-factory t) (assign-types t) force-field-overrides setup);
CL_LISPIFY_NAME("chem:%make-energy-function");
CL_DEF_CLASS_METHOD EnergyFunction_sp EnergyFunction_O::make(Aggregate_sp matter, core::T_sp disableComponents, core::List_sp enableComponents,
                                                             bool useExcludedAtoms, core::T_sp keepInteractionFactory, bool assign_types,
                                                             core::T_sp forceFieldOverrides, core::T_sp setup )
{
  auto  me  = gctools::GC<EnergyFunction_O>::allocate();
  if (disableComponents.notnilp()||enableComponents.notnilp()) {
    SIMPLE_ERROR("I think I stopped using disable-components and enable-components.  If this message appears then I'm wrong and we have to keep using the code below");
  }
  if ( matter.notnilp() ) me->defineForAggregate( matter, useExcludedAtoms, keepInteractionFactory, assign_types, forceFieldOverrides, setup );
#if 1
  //
  // Disable and then enable components
  //
  core::List_sp components = me->allComponents();
  if (disableComponents == _lisp->_true()) {
    for ( auto curcomp : components ) {
      auto onecomp = gc::As<EnergyComponent_sp>(CONS_CAR(curcomp));
      onecomp->disable();
    }
  } else if (disableComponents.consp()) {
    for ( auto cur : gc::As_unsafe<core::List_sp>(disableComponents)) {
      auto disname = gc::As<core::Symbol_sp>(CONS_CAR(cur));
      bool found = false;
      for ( auto curcomp : components ) {
        auto onecomp = gc::As<EnergyComponent_sp>(CONS_CAR(curcomp));
        core::Symbol_sp compname = onecomp->_instanceClass()->_className();
        if (compname == disname) {
          found = true;
          onecomp->disable();
        }
      }
      if (!found) {
        ql::list ll;
        for ( auto curcomp : components ) {
          auto onecomp = gc::As<EnergyComponent_sp>(CONS_CAR(curcomp));
          core::Symbol_sp compname = onecomp->_instanceClass()->_className();
          ll << compname;
        }
        SIMPLE_ERROR("The name {} is not a valid energy component class name - valid names: {}", _rep_(disname), _rep_(ll.cons()));
      }
    }
  }
  if (enableComponents.consp()) {
    for ( auto cur : enableComponents) {
      auto enname = gc::As<core::Symbol_sp>(CONS_CAR(cur));
      bool found = false;
      for ( auto curcomp : components ) {
        auto onecomp = gc::As<EnergyComponent_sp>(CONS_CAR(curcomp));
        core::Symbol_sp compname = onecomp->_instanceClass()->_className();
        if (compname == enname) {
          found = true;
          onecomp->enable();
        }
      }
      if (!found) {
        ql::list ll;
        for ( auto curcomp : components ) {
          auto onecomp = gc::As<EnergyComponent_sp>(CONS_CAR(curcomp));
          core::Symbol_sp compname = onecomp->_instanceClass()->_className();
          ll << compname;
        }
        SIMPLE_ERROR("The name {} is not a valid energy component class name - valid names: {}", _rep_(enname), _rep_(ll.cons()));
      }
    }
  }
#endif
  return me;
};


CL_DOCSTRING("Return (values enabled-component-class-names disabled-component-class-names)");
CL_DEFMETHOD core::T_mv EnergyFunction_O::enabledDisabled() const
{
  auto components = this->allComponents();
  ql::list lenabled;
  ql::list ldisabled;
  for ( auto curcomp : components ) {
    auto onecomp = gc::As<EnergyComponent_sp>(CONS_CAR(curcomp));
    core::Symbol_sp compname = onecomp->_instanceClass()->_className();
    if (onecomp->isEnabled()) {
      lenabled << compname;
    } else {
      ldisabled << compname;
    }
  }
  return Values(lenabled.cons(),ldisabled.cons());
}

bool inAtomSet(core::T_sp activeSet, Atom_sp a)
{
  if ( activeSet.nilp() ) return true;
  if ( core::Cons_sp cset = activeSet.asOrNull<core::Cons_O>() ) {
    for ( auto it : core::List_sp(cset) ) {
      if ( oCar(it) == a ) return true;
    }
    return false;
  } else if ( core::HashTable_sp htset = activeSet.asOrNull<core::HashTable_O>() ) {
    return htset->gethash(a).notnilp();
  } else if ( core::Vector_sp vset = activeSet.asOrNull<core::Vector_O>() ) {
    for ( size_t i(0), iEnd(core::cl__length(vset)); i<iEnd; ++i ) {
      if ( vset->rowMajorAref(i) == a ) return true;
    }
    return false;
  }
  SIMPLE_ERROR("Add support to search {} as a set\n" , _rep_(activeSet));
}



void	EnergyFunction_O::initialize()
{
  this->Base::initialize();
  this->_AtomTable = AtomTable_O::create();
  this->_EnergyComponents = nil<core::T_O>();
  this->setScoringFunctionName(nil<core::T_O>());
  this->_Message = nil<core::T_O>();
  this->useDefaultSettings();
}

void	EnergyFunction_O::useDefaultSettings()
{
  this->_RestrainSecondaryAmides = true;
}



void EnergyFunction_O::fields(core::Record_sp node)
{
  node->field_if_not_unbound(INTERN_(kw,AtomTable),this->_AtomTable);
  node->field_if_not_unbound(INTERN_(kw,BoundingBox),this->_BoundingBox);
  node->field(INTERN_(kw,OtherEnergyComponents),this->_EnergyComponents);
  // Needed to replay generateNonbondEnergyFunctionTables in generateIntoGroup - a deserialized
  // energy function that lost this could not build slot groups.
  //
  // _UseExcludedAtoms is NOT serialized: Record has no bool field overload (only
  // field_if_defined, whose bool is a presence flag, not the value), which is why the existing
  // bool _RestrainSecondaryAmides is absent here too.  It round-trips as false, and false is the
  // only value any live path uses - the true branch of generateNonbondEnergyFunctionTables
  // SIMPLE_ERRORs immediately.  Wrap it if that ever stops being true.
  node->field_if_not_nil(INTERN_(kw,setup),this->_Setup);
  this->Base::fields(node);
}

CL_DOCSTRING(R"dx(Return the bounding-box for the atom-table.)dx");
CL_LISPIFY_NAME(energyFunction-bounding-box);
CL_DEFMETHOD BoundingBox_sp EnergyFunction_O::boundingBox() const
{
  if (this->_BoundingBox.boundp()) {
    return this->_BoundingBox;
  }
  SIMPLE_ERROR("The bounding-box slot is unbound");
}

CL_DOCSTRING(R"dx(Return T if the bounding-box is bound)dx");
CL_LISPIFY_NAME(energyFunction-bounding-box-bound-p);
CL_DEFMETHOD bool EnergyFunction_O::boundingBoxBoundP() const
{
  return this->_BoundingBox.boundp();
}

CL_DOCSTRING(R"dx(Set the bounding-box)dx");
CL_LISPIFY_NAME(energyFunction-set-bounding-box);
CL_DEFMETHOD void EnergyFunction_O::setBoundingBox(BoundingBox_sp boundingBox) {
  this->_BoundingBox = boundingBox;
}

CL_DOCSTRING(R"dx(Make the bounding-box unbound)dx");
CL_LISPIFY_NAME(energyFunction-mak-unbound-bounding-box);
CL_DEFMETHOD void EnergyFunction_O::makUnboundBoundingBox() {
  this->_BoundingBox = unbound<BoundingBox_O>();
}




size_t EnergyFunction_O::getNVectorSize() const
{ 
  return this->_AtomTable->getNVectorSize();
};

#if 0
double EnergyFunction_O::getDihedralComponentEnergy()
{
  return this->_Dihedral->getEnergy();
}

double EnergyFunction_O::getDihedralRestraintComponentEnergy()
{
  return this->_DihedralRestraint->getEnergy();
}
#endif








void	EnergyFunction_O::setOption( core::Symbol_sp option, core::T_sp val)
{
  if ( option == _sym_nonbondTerm)
    {
      bool onval = val.isTrue();
      // Try AMBER nonbond component first
      core::T_sp amberNb = this->findComponentOrNil(EnergyNonbond_O::static_classSymbol());
      if (amberNb.notnilp()) {
        if ( onval ) gc::As<EnergyNonbond_sp>(amberNb)->enable();
        else gc::As<EnergyNonbond_sp>(amberNb)->disable();
      } else {
        // Try Rosetta components
        core::T_sp rnb = this->findComponentOrNil(EnergyRosettaNonbond_O::static_classSymbol());
        core::T_sp relec = this->findComponentOrNil(EnergyRosettaElec_O::static_classSymbol());
        core::T_sp nb14 = this->findComponentOrNil(EnergyNonbond14_O::static_classSymbol());
        if (rnb.nilp() && relec.nilp()) {
          SIMPLE_ERROR("No nonbond component found (neither AMBER nor Rosetta)");
        }
        if (rnb.notnilp()) {
          if ( onval ) gc::As<EnergyComponent_sp>(rnb)->enable();
          else gc::As<EnergyComponent_sp>(rnb)->disable();
        }
        if (relec.notnilp()) {
          if ( onval ) gc::As<EnergyComponent_sp>(relec)->enable();
          else gc::As<EnergyComponent_sp>(relec)->disable();
        }
        if (nb14.notnilp()) {
          if ( onval ) gc::As<EnergyComponent_sp>(nb14)->enable();
          else gc::As<EnergyComponent_sp>(nb14)->disable();
        }
      }
    }
  else if ( option == _sym_restraintAnchor )
    {
      bool onval = val.isTrue();
      //#if USE_ALL_ENERGY_COMPONENTS
      if ( onval ) this->getAnchorRestraintComponent()->enable();
      else	this->getAnchorRestraintComponent()->disable();
      //#endif
    } else
    {
      SIMPLE_ERROR("Unknown EnergyFunction setOption keyword[{}]" , _rep_(option) );
    }
}



/*
  __BEGIN_DOC(classes.EnergyFunction.setOptions,subsection,setOptions)
  \scriptCmd{setOptions}{energyFunction options}

  Set the \scriptArg{options} for the energy function. The options are represented by a list of list pairs
  of the form [: 'symbol value ]
  They include...\par
  [ 'nonbond true/false ]\par
  __END_DOC
*/
void	EnergyFunction_O::setOptions( core::List_sp options )
{
  while ( options.notnilp() )
    {
      core::Symbol_sp option = core::oCar(options).as<core::Symbol_O>();
      options = options.asCons()->cdr();
      core::T_sp val = core::oCar(options);
      options = options.asCons()->cdr();
      this->setOption(option,val);
    }
}




// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
//
//	This is where we include the code that interfaces to Mathematica
//	generated code.
//
//

class	FiniteDifferenceMismatch {
public:
  string	functionName;
  string	termName;
  int	index;
  string	getError() { return "FiniteDifferenceMismatch"; } //  @"+this->functionName()+" term:("+this->termName+")";
};



#define	delta2	0.00001



#undef	TEST_FORCE
#define	TEST_FORCE( func, delta, argLow, argHigh, term, idx ) {         \
    double eLow = this->func argLow;                                    \
    double eHigh = this->func argHigh;                                  \
    double numForce = -(eHigh-eLow)/(delta);                            \
    if ( !_areValuesClose( numForce, term, #func, #term, idx ) ) { fails++; }; \
  }

#undef	TEST_DIAGONAL_HESSIAN
#define	TEST_DIAGONAL_HESSIAN( func, delta, argLow, argMiddle, argHigh, term, idx ) { \
    double eLow = this->func argLow;                                    \
    double eMiddle = this->func argMiddle;                              \
    double eHigh = this->func argHigh;                                  \
    double numHessian = (eHigh+eLow-2.0*eMiddle)/(((delta)/2.0)*((delta)/2.0)); \
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; }; \
  }

#undef	TEST_OFF_DIAGONAL_HESSIAN
#define	TEST_OFF_DIAGONAL_HESSIAN( func, delta, argMM, argPM, argMP, argPP, term, idx ) { \
    double eMM = this->func argMM;                                      \
    double eMP = this->func argMP;                                      \
    double ePM = this->func argPM;                                      \
    double ePP = this->func argPP;                                      \
    double numHessian = ((ePP-ePM)-(eMP-eMM))/(delta*delta);            \
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; }; \
  }





SYMBOL_EXPORT_SC_(ChemPkg,monomer_corrections);
SYMBOL_EXPORT_SC_(ChemPkg,parse_energy_spec);
SYMBOL_EXPORT_SC_(ChemPkg,copy_energy_function);
SYMBOL_EXPORT_SC_(ChemPkg,make_energy_function);



//
double	EnergyFunction_O::evaluateAll( NVector_sp 	pos,
                                       core::T_sp       energyScale,
                                       core::T_sp       energyComponents,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec,
                                       gc::Nilable<NVector_sp> dvec,
                                       core::T_sp activeAtomMask,
                                       core::T_sp debugInteractions,
                                       bool disableRestraints )
{
  bool	hasForce = force.notnilp();
  bool   hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());

  if ( !calcForce && ( calcDiagonalHessian || calcOffDiagonalHessian ) ) {
    SIMPLE_ERROR("Inconsistant arguments: if you want to calcDiagonalHessian or calcOffDiagonalHessian you must calcForce");
  }
  if ( !calcDiagonalHessian & calcOffDiagonalHessian ) {
    SIMPLE_ERROR("Inconsistant arguments: if you want to calcOffDiagonalHessian you must calcDiagonalHessian");
  }

  double totalEnergy = 0.0;

  if ( hasForce ) force->zero();
  if ( hasHessian ) hessian->zero();
  if ( hasHdAndD ) {
    LOG("Zeroing hdvec" );
    hdvec->zero();	// Zero the result
  }

  LOG("Starting evaluation of energy" );

  // ONE drift walk for every pair-list component below.
  //
  // Each of them needs to know whether the coordinates have moved far enough to invalidate its
  // pair list, and each of them used to answer that by walking a private copy of all 3N
  // coordinates.  Across a blueprint's ~2829 components that was tens of millions of comparisons
  // per evaluation, all recomputing one number.  The walk happens here instead, once, and
  // maybeRebuildPairListImpl reads the result - see AtomTable_O::_RefCoords.
  //
  // THIS IS THE FUNNEL.  evaluateEnergy, evaluateEnergyForce and
  // evaluateEnergyIntoFaRestFaRepVector all route through evaluateAll, so this covers every
  // evaluation path including a blueprint scoring one group with the rest disabled.  A caller that
  // reaches a component directly is caught by driftFor's owner check only if it also brings a
  // different coordinate vector - mutating the same vector in place and skipping this would leave
  // stale pair lists in use.
  // SKIPPED while the coordinates are pinned - see AtomTable_O::_CoordinatesPinned.  A pair scan
  // re-points every component and evaluates once per pair without ever touching the coordinates, so
  // this walk runs n^2/2 times to produce a number that is identically zero and that every
  // component then discards, having been invalidated by set-matters.  It measured 24.8% of one such
  // scan.
  if (this->_AtomTable.notnilp() && !this->_AtomTable->coordinatesPinned()) {
    this->_AtomTable->updateDrift(pos);
  }

  for ( auto cur : this->_EnergyComponents ) {
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(CONS_CAR(cur));
    if (!(disableRestraints && component->restraintp()) && component->isEnabled()) {
      totalEnergy+= component->evaluateAllComponent(this->asSmartPtr(),
                                                    pos,
                                                    energyScale,
                                                    energyComponents,
                                                    calcForce,
                                                    force,
                                                    calcDiagonalHessian,
                                                    calcOffDiagonalHessian,
                                                    hessian,
                                                    hdvec,
                                                    dvec,
                                                    activeAtomMask,
                                                    debugInteractions );
    }
  }
  return totalEnergy;
}

CL_DEFMETHOD
size_t	EnergyFunction_O::runTestCalls( core::T_sp stream, NVector_sp 	pos) {
  size_t badCount = 0;
  for ( auto cur : this->_EnergyComponents ) {
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(CONS_CAR(cur));
    if (component->isEnabled()) {
      size_t oneBadCount = component->runTestCalls( stream, pos );
      if (stream.notnilp() && oneBadCount!=0) {
        core::print(fmt::format("EnergyFunction_O::runTestCalls found {} mismatches between analytical and finite difference for {}\n", oneBadCount, component->_instanceClass()->_classNameAsString() ), stream );
      }
      badCount += oneBadCount;
    }
  }
  if (stream.notnilp()) {
    core::print(fmt::format("EnergyFunction_O::runTestCalls found {} mismatches between analytical and finite difference terms\n", badCount ), stream );
  }
  return badCount;
}

/*!
 * Compare the analytical force and hessian components term by term with
 * numerical ones.  Print a message for every mismatch
 * Return the number of mismatches
 */
int	EnergyFunction_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm( NVector_sp 	pos)
{
  int	fails = 0;
  IMPLEMENT_ME();
#if 0
  {
    for ( auto cur : this->_EnergyComponents ) {
      EnergyComponent_sp component = gc::As<EnergyComponent_sp>(CONS_CAR(cur));
      component->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    }
  }
#endif
  return fails;
}


/*!
 * Compare the analytical force and hessian components term by term with
 * numerical ones at the current position.  Print a message for every mismatch
 * Return the number of mismatches
 */
int	EnergyFunction_O::compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition( )
{
  NVector_sp pos = NVector_O::create(this->getNVectorSize());
  this->loadCoordinatesIntoVector(pos);
  return this->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
}




SYMBOL_EXPORT_SC_(KeywordPkg,stretch_deviations);
SYMBOL_EXPORT_SC_(KeywordPkg,angle_deviations);


//
//
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------





string	EnergyFunction_O::energyTermsEnabled()
{
  stringstream ss;
  ss.str("");
  for ( auto cur : this->_EnergyComponents ) {
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(CONS_CAR(cur));
    ss << component->enabledAsString();
  }
  return ss.str();
}




#define	mathematicaSinNPhi(n,s,c) 0.0
#define	mathematicaCosNPhi(n,s,c) 0.0













#define	DELTA	0.00000001

double	EnergyFunction_O::calculateNumericalDerivative(NVector_sp pos, core::T_sp energyScale, double delta, uint i, core::T_sp activeAtomMask )
{
  double x, ylow, yhigh, fval;
  double	deltaDiv2 = delta/2.0;
  x = pos->element(i);
  pos->setElement(i,x-deltaDiv2);
  ylow = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
  pos->setElement(i,x+deltaDiv2);
  yhigh = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
  pos->setElement(i,x);
  fval = (yhigh-ylow)/delta;
  return fval;
}


double	EnergyFunction_O::calculateNumericalSecondDerivative(NVector_sp pos, core::T_sp energyScale, double delta, uint i, uint j, core::T_sp activeAtomMask )
{
  double	x, fxmh, fx, fxph, f2;
  double	y, fpipj, fpimj, fmipj, fmimj, fp, fm;
  if ( i==j ) {
    x = pos->element(i);
    pos->setElement(i,x-delta);
    fxmh = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
    pos->setElement(i,x+delta);
    fxph = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
    pos->setElement(i,x);
    fx = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
    f2 = (fxph+fxmh-2.0*(fx))/(delta*delta);
  } else {
    double	deltaDiv2 = delta/2.0;
    x = pos->element(i);
    y = pos->element(j);
    pos->setElement(i,x+deltaDiv2);
    pos->setElement(j,y+deltaDiv2);
    fpipj = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
    pos->setElement(i,x+deltaDiv2);
    pos->setElement(j,y-deltaDiv2);
    fpimj = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
    pos->setElement(i,x-deltaDiv2);
    pos->setElement(j,y+deltaDiv2);
    fmipj = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
    pos->setElement(i,x-deltaDiv2);
    pos->setElement(j,y-deltaDiv2);
    fmimj = this->evaluateEnergy(pos,energyScale,nil<core::T_O>(),activeAtomMask);
    pos->setElement(i,x);
    pos->setElement(j,y);
    LOG("fpipj = {}" , fpipj  );
    LOG("fpimj = {}" , fpimj  );
    LOG("fmipj = {}" , fmipj  );
    LOG("fmimj = {}" , fmimj  );
    fp = (fpipj-fpimj)/delta;
    fm = (fmipj-fmimj)/delta;
    LOG("fp = {}" , fp  );
    LOG("fm = {}" , fm  );
    f2 = (fp-fm)/delta;
    LOG("f2 = {}" , f2  );
  }
  return f2;
}

CL_DOCSTRING(R"doc(Fill SPHERES with one bounding sphere per piece, packed (x y z radius).
COORDS is the vector the pieces were placed into.  I3-VALUES holds every piece's atom
coordinate indices concatenated; piece P occupies [OFFSETS[P], OFFSETS[P+1]), so OFFSETS has
one more element than there are pieces.  SPHERES must be 4*npieces doubles.)doc");
CL_LAMBDA(coords i3-values offsets spheres);
CL_DEFUN void chem__calculatePieceSpheres(NVector_sp coords,
                                          core::SimpleVector_byte32_t_sp i3values,
                                          core::SimpleVector_byte32_t_sp offsets,
                                          NVector_sp spheres) {
  const double*   xyz = &(*coords)[0];
  const uint32_t* idx = &(*i3values)[0];
  const uint32_t* off = &(*offsets)[0];
  double*         out = &(*spheres)[0];
  size_t npieces = cl__length(offsets) - 1;
  if (cl__length(spheres) < 4*npieces)
    SIMPLE_ERROR("spheres vector is {} long, need {} for {} pieces", cl__length(spheres), 4*npieces, npieces);
  for (size_t p = 0; p < npieces; ++p) {
    uint32_t b = off[p], e = off[p+1];
    size_t n = e - b;
    if (n == 0) {   // no atoms: degenerate sphere at the origin, radius 0
      out[4*p] = out[4*p+1] = out[4*p+2] = out[4*p+3] = 0.0;
      continue;
    }
    double cx=0.0, cy=0.0, cz=0.0;
    for (uint32_t k = b; k < e; ++k) {
      uint32_t i3 = idx[k];
      cx += xyz[i3]; cy += xyz[i3+1]; cz += xyz[i3+2];
    }
    cx /= (double)n; cy /= (double)n; cz /= (double)n;
    double r2max = 0.0;
    for (uint32_t k = b; k < e; ++k) {
      uint32_t i3 = idx[k];
      double dx = xyz[i3]-cx, dy = xyz[i3+1]-cy, dz = xyz[i3+2]-cz;
      double d2 = dx*dx + dy*dy + dz*dz;
      if (d2 > r2max) r2max = d2;
    }
    out[4*p]   = cx;
    out[4*p+1] = cy;
    out[4*p+2] = cz;
    out[4*p+3] = std::sqrt(r2max);
  }
}


CL_DOCSTRING(R"doc(True when pieces I and J can have no atom pair within MAX-CUTOFF.
For a in I and b in J, |a-b| >= |ci-cj| - ri - rj; if that lower bound is already >=
MAX-CUTOFF the whole pair evaluation can be skipped.  Compares squared distances, no sqrt.)doc");
CL_LAMBDA(spheres i j max-cutoff);
CL_DEFUN bool chem__spheresTooFarApartP(NVector_sp spheres, size_t i, size_t j, double maxCutoff) {
  const double* s = &(*spheres)[0];
  const double* a = s + 4*i;
  const double* b = s + 4*j;
  double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
  double reach = a[3] + b[3] + maxCutoff;
  return (dx*dx + dy*dy + dz*dz) > (reach*reach);
}

CL_LISPIFY_NAME("evaluate-energy-into-fa-rest-fa-rep-vector");
CL_LAMBDA((energy-function chem:energy-function) positions fa-rest-fa-rep-vector index fa-rep-component &key energy-scale active-atom-mask debug-interactions disable-restraints);
CL_DEFMETHOD void EnergyFunction_O::evaluateEnergyIntoFaRestFaRepVector(NVector_sp pos, NVector_sp faRestFaRepVector, size_t index,
                                                                        core::T_sp faRepComponent,
                                                                        core::T_sp energyScale, core::T_sp activeAtomMask,
                                                                        core::T_sp debugInteractions,
                                                                        bool disableRestraints ) {
  // 1. total across all components (nil energyComponents => no allocation)
  double total = this->evaluateEnergy(pos, energyScale, nil<core::T_O>(),
                                      activeAtomMask, debugInteractions, disableRestraints );
  // 2. pull fa_rep + weight off the caller-supplied nonbond component, whose _LastFaRep
  //    was just set by evaluateEnergy(this,...) above.  The caller resolves this ONCE per
  //    ef and passes it in, so the hot pair loop never searches the component list.  It
  //    MUST be THIS ef's own rosetta-nonbond component, or getLastFaRep() is stale.
  //
  //    NIL means this ef has NO rosetta-nonbond component, which is a legitimate
  //    configuration — a ligand-only system has no nonbonded components at all (the
  //    backbone scan's intermolecular ef is an empty copy there, and evaluates to 0).
  //    With no repulsive term to separate out, fa_rep = 0 and fa_rest = total, which
  //    then reweights correctly at every rep_weight.  Callers derive this argument from
  //    find-component-or-nil, so NIL already means absent — we do NOT re-search the
  //    component list here (that search is exactly what hoisting the component avoids).
  double faRep = 0.0;
  double w     = 0.0;
  if (faRepComponent.notnilp()) {
    EnergyRosettaNonbond_sp nb = gc::As<EnergyRosettaNonbond_sp>(faRepComponent);
    faRep = nb->getLastFaRep();
    w     = nb->getRepWeight();
  }
  double faRest = total - w * faRep;          // everything except the weighted rep
  // 3. write raw doubles straight into the destination — no boxing
  double* v = &(*faRestFaRepVector)[0];
  v[2*index]     = faRest;
  v[2*index + 1] = faRep;
}

// Evaluate ONE named component, skipping the component-list scan.  See the header for why: a pair
// scan enables one group out of ~2829 components and evaluateAll re-finds it by walking a cons list
// on every one of its n^2/2 evaluations.
//
// Everything evaluateAll does BESIDES that loop is either preserved here or provably a no-op on
// this path: the calcForce/calcHessian consistency checks cannot fire with both false, and the
// force/hessian/hdvec zeroing is skipped because all three are NIL.  The drift update is NOT a
// no-op and is kept - evaluateAll is the funnel that keeps every pair list honest, and a second
// entry point that skipped it would leave stale lists in use with no way to notice.

CL_LISPIFY_NAME("evaluate-component-into-fa-rest-fa-rep-vector");
CL_LAMBDA((energy-function chem:energy-function) component positions fa-rest-fa-rep-vector index fa-rep-component &key energy-scale active-atom-mask debug-interactions);
CL_DEFMETHOD void EnergyFunction_O::evaluateComponentIntoFaRestFaRepVector(core::T_sp tcomponent,
                                                                          NVector_sp pos,
                                                                          NVector_sp faRestFaRepVector,
                                                                          size_t index,
                                                                          core::T_sp faRepComponent,
                                                                          core::T_sp energyScale,
                                                                          core::T_sp activeAtomMask,
                                                                          core::T_sp debugInteractions ) {
  if (this->_AtomTable.notnilp() && !this->_AtomTable->coordinatesPinned()) {
    this->_AtomTable->updateDrift(pos);
  }
  EnergyComponent_sp component = gc::As<EnergyComponent_sp>(tcomponent);
  double total = component->evaluateAllComponent(this->asSmartPtr(),
                                                 pos,
                                                 energyScale,
                                                 nil<core::T_O>(),   // energyComponents: no breakdown
                                                 false,              // calcForce
                                                 nil<NVector_O>(),   // force
                                                 false,              // calcDiagonalHessian
                                                 false,              // calcOffDiagonalHessian
                                                 nil<AbstractLargeSquareMatrix_O>(),
                                                 nil<NVector_O>(),   // hdvec
                                                 nil<NVector_O>(),   // dvec
                                                 activeAtomMask,
                                                 debugInteractions );
  // Same split as evaluateEnergyIntoFaRestFaRepVector: fa_rep comes off the caller-supplied nonbond
  // component, whose _LastFaRep the evaluation above just set.  NIL is legitimate - a system with no
  // rosetta-nonbond has no repulsive term to separate out, so fa_rep is 0 and fa_rest is the total.
  double faRep = 0.0;
  double w     = 0.0;
  if (faRepComponent.notnilp()) {
    EnergyRosettaNonbond_sp nb = gc::As<EnergyRosettaNonbond_sp>(faRepComponent);
    faRep = nb->getLastFaRep();
    w     = nb->getRepWeight();
  }
  double* v = &(*faRestFaRepVector)[0];
  v[2*index]     = total - w * faRep;
  v[2*index + 1] = faRep;
}


CL_LISPIFY_NAME("energy-total-fa-rest-fa-rep-rep-weight");
CL_LAMBDA(fa-rest-fa-rep-vector index rep-weight);
CL_DEFUN double chem__energyTotalFaRestFaRepRepWeight(NVector_sp faRestFaRepVector, size_t index2, double rep_weight)
{
  ASSERT((index2&1)==0); // index2 should ALWAYS be even
  ASSERT(index2<cl__length(faRestFaRepVector));
  double faRest = (*faRestFaRepVector)[index2];
  double faRep = (*faRestFaRepVector)[index2+1];
  return faRest + (rep_weight*faRep);
}




/*! Calculate the force numerically
 */
void	EnergyFunction_O::evaluateNumericalForce(NVector_sp pos, core::T_sp energyScale, NVector_sp numForce, double delta, core::T_sp activeAtomMask )
{
  double		fval;
  uint		i;

  for (i=0; i<pos->size(); i++ ) {
    fval = -this->calculateNumericalDerivative(pos,energyScale,delta,i,activeAtomMask);
    numForce->setElement(i,fval);
  }
}


/*! Calculate the hessian numerically
 */
void	EnergyFunction_O::evaluateNumericalHessian(NVector_sp pos, core::T_sp energyScale, AbstractLargeSquareMatrix_sp hessian, bool calcOffDiagonal, double delta, core::T_sp activeAtomMask )
{
  double		fval;
  uint		c, r;

  if ( hessian->columns() != pos->size() || hessian->rows()!=pos->size() ) {
    SIMPLE_ERROR("evaluateNumericalHessian must have the right size");
  }
  hessian->zero();
  for ( c=0; c<pos->size(); c++ ) {
    fval = this->calculateNumericalSecondDerivative(pos,energyScale,delta,c,c,activeAtomMask);
    hessian->setElement(c,c,fval);
  }
  if ( !calcOffDiagonal ) return;
  for ( c=0; c<pos->size(); c++ ) {
    for ( r=0; r<pos->size(); r++ ) {
      if ( c!=r) {
        fval = this->calculateNumericalSecondDerivative(pos,energyScale,delta,c,r,activeAtomMask);
        hessian->setElement(c,r,fval);
      }
    }
  }
}


/*!
 * Check if the force matches the numerical force
 *
 * Return an XML node that indicates what is going on.
 * If there is a mis-match then dump the EnergyFunction into the result.
 *
 */
ForceMatchReport_sp EnergyFunction_O::checkIfAnalyticalForceMatchesNumericalForce(NVector_sp pos, core::T_sp energyScale, NVector_sp analyticalForce, core::T_sp activeAtomMask )
{
  ForceMatchReport_sp report;
  NVector_sp	numForce, tempForce;
  stringstream	result;
  double	analyticalMag, numericalMag, dot;
  double	avg;

  report = ForceMatchReport_O::create();

  numForce = NVector_O::create(pos->size());
  this->evaluateNumericalForce(pos,energyScale,numForce,DELTA,activeAtomMask);
  dot = dotProductWithActiveAtomMask(numForce,analyticalForce,nil<core::T_O>());
  numericalMag = magnitudeWithActiveAtomMask(numForce,nil<core::T_O>());
  analyticalMag = magnitudeWithActiveAtomMask(analyticalForce,nil<core::T_O>());
  tempForce = NVector_O::create(pos->size());
  // Evaluate the force at pos again
  this->evaluateEnergyForce(pos,energyScale,true,tempForce,activeAtomMask);
  avg = (analyticalMag+numericalMag)/2.0;
  if ( analyticalMag < VERYSMALL && numericalMag < VERYSMALL ) {
    result.str("");
    result << "analyticalForce and numericalForces are both < VERYSMALL"<<std::endl;
    report->_Message = result.str();
    goto DONE;
  }
  if ( analyticalMag < VERYSMALL ) {
    result.str("");
    result << "analyticalForce < VERYSMALL and numericalForces is not"<<std::endl;
    report->_Message = result.str();
    goto DONE;
  }
  if ( numericalMag < VERYSMALL ) {
    result.str("");
    result << "numericalForce < VERYSMALL and analyticalForce is not"<<std::endl;
    report->_Message = result.str();
    goto DONE;
  }
  if ( avg < VERYSMALL ) {
    report->_Message = "average of Analytical & Numerical Forces is VERY small";
    goto DONE;
  }
  dot /= ( analyticalMag*numericalMag );

  // Dump the analytical and the numerical force vectors to the log
  //
  report->_AnalyticalForce = analyticalForce;
  report->_NumericalForce = numForce;
  if ( fabs(analyticalMag-numericalMag)/avg >0.1 ) {
    result.str("");
    result << "Lengths of Analytical and Numerical forces differ by more than 10%" << std::endl;
    result << "|analyticalForce| == " << analyticalMag << "  |numericalForce| == "<< numericalMag << std::endl;
    result << "(analyticalForce/|analyticalForce|).(numericalForce/|numericalForce|) = "<< dot << std::endl;
    report->_Message = result.str();
    this->saveCoordinatesAndForcesFromVectors(pos,analyticalForce);
    goto DONE;
  }
  if ( dot < 0.98 ) {
    result.str("");
    result << "The Lengths of Analytical and Numerical forces are very similar but they are not parallel%" << std::endl;
    result << "|analyticalForce| == " << analyticalMag << "  |numericalForce| == "<< numericalMag << std::endl;
    result << "(analyticalForce/|analyticalForce|).(numericalForce/|numericalForce|) = "<< dot << std::endl;
    report->_Message = result.str();
    this->saveCoordinatesAndForcesFromVectors(pos,analyticalForce);
    goto DONE;
  }
  report->_Message = "Analytical and Numerical forces are virtually identical";
 DONE:
  return report;
}



void	EnergyFunction_O::summarizeTerms()
{
  core::clasp_write_string(fmt::format("Number of atom terms: {}\n" , this->_AtomTable->getNumberOfAtoms() ));
};



void	EnergyFunction_O::dumpTerms()
{
  core::HashTable_sp atomTypes = this->_AtomTypes;
  this->_AtomTable->dumpTerms(atomTypes);
  for ( auto cur : this->_EnergyComponents ) {
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(CONS_CAR(cur));
    component->dumpTerms(atomTypes);
  }
}

/*
 * keepInteractionFactory can be one of ...
 * nil = discard all interactions
 * T   = keep all interactions
 * Fn  = pass up to 4 atoms to the function and if the function returns true then the interaction is kept
 */
core::T_sp specializeKeepInteractionFactory( core::T_sp keepInteractionFactory, core::T_sp aclass ) {
  if (keepInteractionFactory.nilp()) return keepInteractionFactory;
  if (keepInteractionFactory == _lisp->_true()) return keepInteractionFactory;
  core::T_sp keepInteractionFunction = core::eval::funcall(keepInteractionFactory,aclass);
  return keepInteractionFunction;
}

bool skipInteraction_EnergyNonbond( core::T_sp keepInteractionFunction,
                                    Atom_sp atomA,
                                    Atom_sp atomB,
                                    core::T_sp ti3A,
                                    core::T_sp ti3B )
{
  if (keepInteractionFunction.nilp()) return true;
  if (keepInteractionFunction == _lisp->_true()) return false;
  unlikely_if (!gc::IsA<core::Function_sp>(keepInteractionFunction)) {
    SIMPLE_ERROR("The object {} must be NIL, T or a function object", keepInteractionFunction );
  }
  core::T_sp result = core::eval::funcall( keepInteractionFunction, atomA, atomB, ti3A, ti3B );
  return result.nilp();
}

bool skipInteraction( core::T_sp keepInteractionFunction,
                      Atom_sp atomA,
                      Atom_sp atomB,
                      Atom_sp atomC,
                      Atom_sp atomD )
{
  if (keepInteractionFunction.nilp()) return true;
  if (keepInteractionFunction == _lisp->_true()) return false;
  unlikely_if (!gc::IsA<core::Function_sp>(keepInteractionFunction)) {
    SIMPLE_ERROR("The object {} must be NIL, T or a function object", keepInteractionFunction );
  }
  if (atomD.boundp()) {
    core::T_sp result = core::eval::funcall( keepInteractionFunction, atomA, atomB, atomC, atomD );
    return result.nilp();
  }
  if (atomC.boundp()) {
    core::T_sp result = core::eval::funcall( keepInteractionFunction, atomA, atomB, atomC );
    return result.nilp();
  }
  if (atomB.boundp()) {
    core::T_sp result = core::eval::funcall( keepInteractionFunction, atomA, atomB );
    return result.nilp();
  }
  if (atomA.boundp()) {
    core::T_sp result = core::eval::funcall( keepInteractionFunction, atomA );
    return result.nilp();
  }
  return false;
}

int EnergyFunction_O::_applyRestraints(core::T_sp nonbondDb, core::Iterator_sp restraintIterator, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes )
{
  int terms = 0;
  restraintIterator->first();
  if (keepInteractionFactory.nilp()) return 0;
  while ( restraintIterator->notDone() ) {
    Restraint_sp restraint = restraintIterator->current<Restraint_O>();
    if ( !restraint->isActive() ) goto CONT;
    if ( restraint.isA<RestraintDihedral_O>() )
      {
        RestraintDihedral_sp dih = (restraint).as<RestraintDihedral_O>();
        EnergyDihedralRestraint energyTerm;
        core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory,RestraintDihedral_O::staticClass());
        if (skipInteraction(keepInteraction, dih->getAtomA(), dih->getAtomB(), dih->getAtomC(), dih->getAtomD() )) goto CONT;
        energyTerm._Atom1 = dih->getAtomA();
        energyTerm._Atom2 = dih->getAtomB();
        energyTerm._Atom3 = dih->getAtomC();
        energyTerm._Atom4 = dih->getAtomD();
        EnergyAtom* ea1 = this->getEnergyAtomPointer(dih->getAtomA());
        EnergyAtom* ea2 = this->getEnergyAtomPointer(dih->getAtomB());
        EnergyAtom* ea3 = this->getEnergyAtomPointer(dih->getAtomC());
        EnergyAtom* ea4 = this->getEnergyAtomPointer(dih->getAtomD());
        energyTerm.term.I1 = ea1->coordinateIndexTimes3();
        energyTerm.term.I2 = ea2->coordinateIndexTimes3();
        energyTerm.term.I3 = ea3->coordinateIndexTimes3();
        energyTerm.term.I4 = ea4->coordinateIndexTimes3();
        energyTerm.term.phi0 = dih->getDegrees()*0.0174533;
        energyTerm.term.kdh = dih->getWeight();
        auto dihedralRestraint = ensureComponent<EnergyDihedralRestraint_O>(this->asSmartPtr());
        dihedralRestraint->addTerm(energyTerm);
        ++terms;
      } else if ( restraint.isA<RestraintAnchor_O>() ) {
      RestraintAnchor_sp anchor = (restraint).as<RestraintAnchor_O>();
      EnergyAnchorRestraint	iterm;
      Vector3		anchorPos;
      Atom_sp a1 = anchor->getAtom();
      core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory,RestraintAnchor_O::staticClass());
      if ( skipInteraction(keepInteraction,a1) ) goto CONT;
      EnergyAtom* ea1 = this->getEnergyAtomPointer(a1);
      LOG("Create an anchor restraint for {}" , a1->description()  );
      anchorPos = anchor->getAnchorPos();
      iterm.term.xa = anchorPos.getX();
      iterm.term.ya = anchorPos.getY();
      iterm.term.za = anchorPos.getZ();
      iterm.term.ka = DefaultAnchorRestraintWeight;
      iterm.term.i3x1 = ea1->coordinateIndexTimes3();
      auto anchorRestraint = ensureComponent<EnergyAnchorRestraint_O>(this->asSmartPtr());
      anchorRestraint->addTerm(iterm);
      ++terms;
    } else if ( restraint.isA<RestraintFixedNonbond_O>() ) {
      auto fixedNonbondRestraint = ensureComponent<EnergyFixedNonbondRestraint_O>(this->asSmartPtr());
      fixedNonbondRestraint->setupForEvaluation(this->_AtomTable,this->_NonbondCrossTermTable);
      RestraintFixedNonbond_sp fixedNonbond = restraint.as<RestraintFixedNonbond_O>();
      Matter_sp matter = fixedNonbond->getMatter();
      //	    EnergyAtom	energyAtom(_lisp);
      Loop loop;
      {
        loop.loopTopGoal(matter,ATOMS);
        while ( loop.advanceLoopAndProcess() ) 
          {
            Atom_sp a1 = loop.getAtom();
            core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory,RestraintFixedNonbond_O::staticClass());
            if ( skipInteraction( keepInteraction, a1 ) ) goto CONT;
            if ( a1.isA<VirtualAtom_O>() ) continue; // skip virtuals
            fixedNonbondRestraint->addFixedAtom(nonbondDb,a1,atomTypes);
            ++terms;
          }
      }
    } else if ( restraint.isA<RestraintDistance_O>() ) {
      RestraintDistance_sp rd = gc::As_unsafe<RestraintDistance_sp>(restraint);
      EnergyAtom* ea1 = this->getEnergyAtomPointer(rd->_A);
      EnergyAtom* ea2 = this->getEnergyAtomPointer(rd->_B);
      EnergyStretch   energyStretch(rd->_A,rd->_B,ea1->coordinateIndexTimes3(),ea2->coordinateIndexTimes3(),rd->_K,rd->_R0);
      auto stretchComponent = ensureComponent<EnergyStretch_O>(this->asSmartPtr());
      stretchComponent->addTerm(energyStretch);
    } else {
      SIMPLE_ERROR("Handle restraint: {}" , _rep_(restraint));
    }
  CONT:
    restraintIterator->next();
  }
  return terms;
}

void EnergyFunction_O::_addDihedralRestraintDegrees(EnergyDihedralRestraint_sp dihedralRestraintComponent, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double phi0_degrees, double kdh, core::T_sp keepInteraction)
{
  if ( skipInteraction( keepInteraction, a1, a2, a3, a4 ) ) return;
  dihedralRestraintComponent->addDihedralRestraint(this->asSmartPtr(),
                                                   kdh,
                                                   phi0_degrees*0.0174533,
                                                   a1,a2,a3,a4
                                                   );
}

void EnergyFunction_O::__createSecondaryAmideRestraints(gctools::Vec0<Atom_sp>& nitrogens, core::T_sp keepInteractionFactory )
{
  if (keepInteractionFactory.nilp()) return;
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, RestraintDihedral_O::staticClass());
  gctools::Vec0<Atom_sp>::iterator ni;
#if 0
  // right
  core::clasp_write_string("%s:%d:%s Generating correct dihedral restraint for amide.\n", __FILE__, __LINE__, __FUNCTION__);
  double trans  = -180.0;
  double cis    = 0.0;
#else
  // wrong/testing
  core::clasp_write_string(fmt::format("%s:%d:%s Generating WRONG dihedral restraint for amide.\n", __FILE__, __LINE__, __FUNCTION__));
  double trans  = 90.0;
  double cis    = -90.0;
#endif

  auto dihedralRestraintComponent = ensureComponent<EnergyDihedralRestraint_O>(this->asSmartPtr());
  double weight = 1.0;
  for ( ni=nitrogens.begin(); ni!=nitrogens.end(); ni++ ) {
    SmartsRoot_sp secondaryAmide = gctools::As<SmartsRoot_sp>(chem::_sym_STARsecondaryAmideSmartsSTAR->symbolValue());
    core::T_mv match_mv = chem__chem_info_match(secondaryAmide,*ni);
    core::MultipleValues &values = core::lisp_multipleValues();
    if ( match_mv.notnilp() ) {
      ChemInfoMatch_sp match = gc::As<ChemInfoMatch_sp>(values.second(match_mv.number_of_values()));
      Atom_sp ax = match->tag(core::make_fixnum(1));
      Atom_sp ax1 = match->tag(core::make_fixnum(2));
      Atom_sp ax2 = match->tag(core::make_fixnum(3));
      Atom_sp ay = match->tag(core::make_fixnum(4));	// Carbonyl carbon
      Atom_sp ay1 = match->tag(core::make_fixnum(5));
      Atom_sp ay2 = match->tag(core::make_fixnum(6));
      //
      // If amide carbonyl carbon is in ring
      // then we want a cis amide bond
      // otherwise trans
      //
      // Achieve this by swapping ay1 and ay2 if its
      // supposed to be a cis amide bond
      //
      if ( ay->getRingMembershipCount() > 0 )
        {
          Atom_sp azz = ay1;
          ay1 = ay2;
          ay2 = azz;
        }
      LOG("Applying a secondary amide restraint between {} and {}" , ax->description() , ay->description()  );
      //
      // H3(ax2) and O5(ay1) should be trans
      if (!skipInteraction(keepInteraction, ax1, ax, ay, ay1 ))
        this->_addDihedralRestraintDegrees(dihedralRestraintComponent,ax1,ax,ay,ay1,cis,weight,keepInteraction);
      LOG("Restrain cis {} - {} - {} -{}" , ax1->description() , ax->description() , ay->description() , ay1->description()  );
      if (!skipInteraction(keepInteraction, ax1, ax, ay, ay2 )) 
        this->_addDihedralRestraintDegrees(dihedralRestraintComponent,ax1,ax,ay,ay2,trans,weight,keepInteraction);
      LOG("Restrain trans {} - {} - {} -{}" , ax1->description() , ax->description() , ay->description() , ay2->description()  );
      if (!skipInteraction(keepInteraction, ax2, ax, ay, ay1 )) 
        this->_addDihedralRestraintDegrees(dihedralRestraintComponent,ax2,ax,ay,ay1,trans,weight,keepInteraction);
      LOG("Restrain trans {} - {} - {} -{}" , ax2->description() , ax->description() , ay->description() , ay1->description()  );
      if (!skipInteraction(keepInteraction, ax2, ax, ay, ay2 )) 
        this->_addDihedralRestraintDegrees(dihedralRestraintComponent,ax2,ax,ay,ay2,cis,weight,keepInteraction);
      LOG("Restrain cis {} - {} - {} -{}" , ax2->description() , ax->description() , ay->description() , ay2->description()  );
    }
  }
}





SYMBOL_EXPORT_SC_(ChemPkg,find_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,force_fields_as_list);
SYMBOL_EXPORT_SC_(ChemPkg, nonbond_force_field_name);
SYMBOL_EXPORT_SC_(ChemPkg,add_shadowing_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,assign_force_field_types);
SYMBOL_EXPORT_SC_(ChemPkg,compute_merged_nonbond_force_field_for_aggregate);
SYMBOL_EXPORT_SC_(ChemPkg,compute_merged_lksolvation_force_field_for_aggregate);
SYMBOL_EXPORT_SC_(ChemPkg,force_field_component_merge);
SYMBOL_EXPORT_SC_(ChemPkg,clear_combined_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,generate_molecule_energy_function_tables);
SYMBOL_EXPORT_SC_(ChemPkg,STARparameter_warningsSTAR);
SYMBOL_EXPORT_SC_(ChemPkg,report_parameter_warnings);

SYMBOL_EXPORT_SC_(ChemPkg,identify_aromatic_rings);
SYMBOL_EXPORT_SC_(ChemPkg,STARcurrent_aromaticity_informationSTAR);

SYMBOL_EXPORT_SC_(KeywordPkg,amber);
SYMBOL_EXPORT_SC_(KeywordPkg,rosetta);

core::T_sp EnergyFunction_O::findComponentOrNil(core::Symbol_sp classSymbol) {
  for ( auto cur : this->_EnergyComponents ) {
    EnergyComponent_sp comp = gc::As<EnergyComponent_sp>(oCar(cur));
    if (classSymbol == comp->_instanceClass()->_className()) {
      return comp;
    }
  }
  return nil<core::T_O>();
}

EnergyComponent_sp EnergyFunction_O::findComponent(core::Symbol_sp classSymbol) {
  for ( auto cur : this->_EnergyComponents ) {
    EnergyComponent_sp comp = gc::As<EnergyComponent_sp>(oCar(cur));
    if (classSymbol == comp->_instanceClass()->_className()) {
      return comp;
    }
  }
  SIMPLE_ERROR("Could not find component with name {}", _rep_(classSymbol));
}

void EnergyFunction_O::addComponentIfMissing(EnergyComponent_sp comp) {
  if (this->findComponentOrNil(comp->_instanceClass()->_className()).nilp()) {
    this->pushEnergyComponent(comp);
  }
}

EnergyStretch_sp	EnergyFunction_O::getStretchComponent() {
  auto comp = this->findComponentOrNil(EnergyStretch_O::static_classSymbol());
  return gc::As<EnergyStretch_sp>(comp);
};

EnergyAngle_sp	EnergyFunction_O::getAngleComponent() {
  auto comp = this->findComponentOrNil(EnergyAngle_O::static_classSymbol());
  return gc::As<EnergyAngle_sp>(comp);
};
EnergyDihedral_sp EnergyFunction_O::getDihedralComponent() {
  auto comp = this->findComponentOrNil(EnergyDihedral_O::static_classSymbol());
  return gc::As<EnergyDihedral_sp>(comp);
};

EnergyNonbond_sp EnergyFunction_O::getNonbondComponent() {
  auto comp = this->findComponentOrNil(EnergyNonbond_O::static_classSymbol());
  return gc::As<EnergyNonbond_sp>(comp);
};

EnergyChiralRestraint_sp EnergyFunction_O::getChiralRestraintComponent() {
  auto comp = this->findComponentOrNil(EnergyChiralRestraint_O::static_classSymbol());
  return gc::As<EnergyChiralRestraint_sp>(comp);
};

EnergyAnchorRestraint_sp EnergyFunction_O::getAnchorRestraintComponent() {
  auto comp = this->findComponentOrNil(EnergyAnchorRestraint_O::static_classSymbol());
  return gc::As<EnergyAnchorRestraint_sp>(comp);
};

EnergyDihedralRestraint_sp EnergyFunction_O::getDihedralRestraintComponent() {
  auto comp = this->findComponentOrNil(EnergyDihedralRestraint_O::static_classSymbol());
  return gc::As<EnergyDihedralRestraint_sp>(comp);
};

EnergyFixedNonbondRestraint_sp EnergyFunction_O::getFixedNonbondRestraintComponent() {
  auto comp = this->findComponentOrNil(EnergyFixedNonbondRestraint_O::static_classSymbol());
  return gc::As<EnergyFixedNonbondRestraint_sp>(comp);
};

void EnergyFunction_O::resolveMoleculeForceFields(Aggregate_sp aggregate,
                                                    core::T_sp forceFieldOverrides,
                                                    core::HashTable_sp molecule_force_fields,
                                                    core::HashTable_sp molecule_force_field_names)
{
  Loop moleculeLoop;
  moleculeLoop.loopTopGoal(aggregate,MOLECULES);
  while (moleculeLoop.advanceLoopAndProcess() ) {
    Molecule_sp molecule = moleculeLoop.getMolecule();
    core::T_sp force_field_name = molecule->force_field_name();
    if (forceFieldOverrides.notnilp()) {
      core::T_sp maybe_override = gc::As<core::HashTable_sp>(forceFieldOverrides)->gethash(molecule);
      if (maybe_override.notnilp()) force_field_name = maybe_override;
    }
    if (molecule_force_fields->gethash(molecule).nilp()) {
      molecule_force_fields->setf_gethash(molecule, core::eval::funcall(_sym_find_force_field,force_field_name));
      molecule_force_field_names->setf_gethash(molecule,force_field_name);
    } else {
      SIMPLE_ERROR("A molecule already has a force-field");
    }
  }
} 


CL_DOCSTRING(R"dx(Run pass 2 over the aggregate this energy function was already defined for,
  collecting every component into GROUP.

  The atom table is untouched: atom numbering, charges and nonbond parameters all stay as
  defineForAggregate left them, so the terms generated here index the same coordinate vector as
  the backbone's.  Call this once per rotamer slot with that slot's KEEP-INTERACTION-FACTORY and a
  fresh GROUP.

  BONDED NIL skips the per-molecule pass and generates only the nonbond components.  Refusing the
  bonded classes in the factory is NOT equivalent - that discards the TERMS, but generateForMolecule
  has already re-run the whole force-field parameterization by then.

  SCOPE-AGGREGATE runs this pass over a DIFFERENT aggregate than the one the energy function was
  defined for - a temporary one holding only what this pass needs to see.

  The atom table is NOT rebuilt and is not consulted for scoping: it still describes the real
  aggregate, and every atom reachable through SCOPE-AGGREGATE is the same object with the same
  coordinate index, so the terms generated index the same coordinate vector as always.  What changes
  is only WHAT GETS WALKED.

  This is what makes a blueprint's rotamer pass affordable.  Its aggregate holds every rotamer of
  every slot simultaneously, and a per-slot pass walks all of them to keep one slot's terms: the
  Loop hierarchy in map-bonds/angles/dihedrals/impropers, ring perception and SMIRKS matching in the
  real-SMIRNOFF path, and construct14InteractionTerms below - all of them scale with the fan-out
  rather than with the answer.  A scope aggregate holding backbone-plus-one-rotamer bounds every one
  of them at once.

  Its molecules must carry a FORCE-FIELD-NAME, like any other: resolveMoleculeForceFields runs over
  whatever aggregate this is given.

  BONDED NIL skips the per-molecule pass and generates only the nonbond components.)dx");
CL_LISPIFY_NAME(energy-function/generate-into-group);
CL_LAMBDA((energy-function chem:energy-function) keep-interaction-factory group &key force-field-overrides (bonded t) scope-aggregate);
CL_DEFMETHOD void EnergyFunction_O::generateIntoGroup(core::T_sp keepInteractionFactory,
                                                      core::T_sp group,
                                                      core::T_sp forceFieldOverrides,
                                                      bool bonded,
                                                      core::T_sp scopeAggregate)
{
  // ONE substitution, here, and every use below reads AGGREGATE - the force-field resolution, the
  // molecule loop and the nonbond tables.  Scoping applied at any one of those alone would leave
  // the others walking the full structure, which is exactly what a per-molecule restriction did.
  core::T_sp matterForPass = scopeAggregate;
  if (scopeAggregate.nilp()) {
    // not provided, so check _Matter and use it.
    matterForPass = this->_Matter;
    if (!matterForPass.isA<Aggregate_O>())
      SIMPLE_ERROR("generateIntoGroup requires an energy function already defined for an aggregate");
  } else if (!matterForPass.isA<Aggregate_O>())
    // scopeAggregate was provided but is not actually an aggregate.
    SIMPLE_ERROR("scope-aggregate must be an aggregate, got {}", _rep_(scopeAggregate));
  Aggregate_sp aggregate = matterForPass.as_unsafe<Aggregate_O>();
  if (this->_AtomTable->getNumberOfAtoms() == 0) {
    SIMPLE_ERROR("The atom table is empty - call defineForAggregate before generateIntoGroup");
  }   
  if (!this->_AtomTable->nonbondForceFieldForAggregateBoundP()) {
    SIMPLE_ERROR("The atom table has no aggregate nonbond force-field - call defineForAggregate before generateIntoGroup");
  }
  core::T_sp nonbondForceField = this->_AtomTable->nonbondForceFieldForAggregate();
  core::HashTable_sp ffs = core::HashTable_O::createEq();
  core::HashTable_sp ffnames = core::HashTable_O::createEq();
  this->resolveMoleculeForceFields(aggregate, forceFieldOverrides, ffs, ffnames);
  // The group joins the ef's component list; its children do not (see ensureComponent).
  if (group.notnilp()) this->pushEnergyComponent(gc::As<EnergyComponent_sp>(group));
  // Solute/solvent ordering mattered only for atom-table numbering, which is already fixed.
  //
  // Skipped entirely when BONDED is false.  generateForMolecule re-parameterizes the molecule
  // against its force field BEFORE any keep-interaction predicate is consulted, so a caller that
  // wants only nonbond cannot express that through the factory - it would pay for the whole
  // parameterization and then discard every term it produced.  For a blueprint's pair-scan group
  // that is not merely wasteful: the parameterization sees the aggregate as it stands, and a spiro
  // NG carrying all 16 of its fan-out bonds fails SMIRNOFF's bond-order check outright.
  //
  // Over AGGREGATE, which is the scope aggregate when one was given - so narrowing the pass needs no
  // special case here at all.  A scope aggregate holding one molecule visits one molecule.
  if (bonded) {
    Loop moleculeLoop;
    moleculeLoop.loopTopGoal(aggregate,MOLECULES);
    while (moleculeLoop.advanceLoopAndProcess()) {
      this->generateForMolecule( moleculeLoop.getMolecule(), keepInteractionFactory,
                                 this->_AtomTypes, nonbondForceField, ffs, ffnames, group );
    }
  }
  // Nonbond terms are AGGREGATE-wide, not per-molecule - the interactions that decide a rotamer
  // cross from the slot to the receptor - so this is one call outside the molecule loop, matching
  // defineForAggregate.  Without it a slot group holds only stretch/angle/dihedral and its energy
  // cannot rank rotamers at all.
  //
  // AGGREGATE here is the scope aggregate when one was given, and that matters: the components
  // themselves are built from _AtomTable, but construct14InteractionTerms walks the MATTER to find
  // 1-4 pairs.  Handed the full aggregate it walks every fanned-out rotamer to produce 1-4 terms the
  // factory then discards - 26% of the profiled scan.  Scoping loses nothing: a 1-4 pair is three
  // bonds apart, so both ends are in the same molecule, and for a slot both are in backbone-plus-
  // that-slot.
  //
  // _Setup and _UseExcludedAtoms are replayed verbatim from defineForAggregate.  If this energy
  // function was never defined for an aggregate they are NIL/false, which the setup dispatch below
  // rejects with "Provide a valid setup" - the guards at the top of this function should have
  // caught that first.
  if (this->_Setup.notnilp()) {
    // BONDED gates the 1-4 component too - it is derived from proper dihedrals, so a caller that
    // asked for no bonded terms was getting the bond-graph walk anyway.  A blueprint's remapped
    // rotamers rely on COPY-REMAPPED-BONDED-TERMS to supply 1-4 along with stretch/angle/dihedral;
    // if that ever stops handling EnergyNonbond14, those slots lose the terms silently.
    this->generateNonbondEnergyFunctionTables( this->_UseExcludedAtoms, aggregate, nonbondForceField,
                                               keepInteractionFactory, this->_AtomTypes,
                                               this->_Setup, group, bonded );
  }
}


// ============================================================
//  copy-remapped bonded terms - one parameterization per monomer, not per rotamer
// ============================================================
//
// Rotamers of the same monomer are THE SAME MOLECULE in different conformations.  Bonded parameters
// depend on atom types and connectivity; both are identical across rotamers.  Coordinates differ,
// but coordinates are not in the terms - only INDEXES INTO the coordinate vector are.  So one
// rotamer's stretch/angle/dihedral terms become another's by renaming atoms and shifting indexes,
// with every parameter carried across byte-for-byte.
//
// That turns the detach scan's bonded cost from one force-field parameterization per SLOT into one
// per (locus, monomer) - the pass that was 95% of the profiled run.
//
// NOT A UNIFORM SHIFT.  A slot's terms cross into the backbone: the CA-CB bond, the CA-CB-CG angle,
// backbone-CA-CB-CG dihedrals.  Backbone atoms are shared between rotamers and do not move.  Only
// indexes inside the source slot's range shift; everything else is copied unchanged.  This is what
// the contiguity check in BLUEPRINT-SLOT-I3-RANGE exists to guarantee.
//
// NONBOND IS DELIBERATELY NOT HANDLED HERE.  Which backbone atoms lie within a slot's cutoff is
// GEOMETRY, and geometry is exactly what differs between rotamers - copying those would be silently
// wrong.  The caller generates them per slot with :BONDED NIL, which is also what pushes the
// destination group onto the energy function.

namespace {

/*! Shift an index that lies inside the source slot; leave everything else alone. */
inline int remapSlotI3(int i3, size_t lo3, size_t hi3, ptrdiff_t delta) {
  size_t u = (size_t)i3;
  if (u >= lo3 && u < hi3) return (int)((ptrdiff_t)i3 + delta);
  return i3;
}

/*! The atom at an i3, straight from the atom table.
 *
 * Not from a per-slot atom vector: ATOM-SET-TO-VECTOR walks a hash set and its order is not i3
 * order, so indexing one by position would pair terms with the wrong atoms - and since the
 * INDEXES would still be right, the energy would be correct while every atom label was wrong. */
inline Atom_sp atomAtI3(AtomTable_sp at, int i3) {
  return at->getVectorEnergyAtoms()[i3/3].atom();
}

} // anonymous namespace

CL_LAMBDA(energy-function dest-group src-group src-lo3 src-hi3 delta);
CL_DOCSTRING(R"dx(Copy SRC-GROUP's bonded terms into DEST-GROUP, renaming the atoms of one rotamer
  slot to those of another.

  Indexes in [SRC-LO3, SRC-HI3) are shifted by DELTA - those are the source slot's own atoms.
  Everything outside that range is backbone, shared between rotamers, and is copied unchanged.
  Every force-field parameter is carried across untouched: this is a rename, not a
  reparameterization, and it is exact because bonded parameters depend on topology alone.

  Handles ENERGY-STRETCH, ENERGY-ANGLE and ENERGY-DIHEDRAL.  Any other component in SRC-GROUP -
  notably the nonbond components, whose terms depend on geometry - is SKIPPED, and the second return
  value counts them so a caller can assert it saw only what it expected.

  Returns (values TERMS-COPIED COMPONENTS-SKIPPED).)dx")
CL_DEFUN core::T_mv chem__copy_remapped_bonded_terms(EnergyFunction_sp ef,
                                                     core::T_sp tdestGroup,
                                                     core::T_sp tsrcGroup,
                                                     size_t srcLo3,
                                                     size_t srcHi3,
                                                     gc::Fixnum delta)
{
  EnergyComponentGroup_sp src = gc::As<EnergyComponentGroup_sp>(tsrcGroup);
  if (!gc::IsA<EnergyComponentGroup_sp>(tdestGroup)) {
    SIMPLE_ERROR("dest-group must be an energy-component-group, got {}", _rep_(tdestGroup));
  }
  if (srcHi3 <= srcLo3) {
    SIMPLE_ERROR("Empty source slot range lo3={} hi3={}", srcLo3, srcHi3);
  }
  AtomTable_sp at = ef->_AtomTable;
  size_t copied = 0;
  size_t skipped = 0;
  for (auto& c : src->_Components) {
    if (gc::IsA<EnergyStretch_sp>(c)) {
      EnergyStretch_sp s = gc::As_unsafe<EnergyStretch_sp>(c);
      EnergyStretch_sp d = ensureComponent<EnergyStretch_O>(ef, tdestGroup);
      // Carried, not defaulted.  Both happen to match a fresh component today - scale 1.0, enabled -
      // so leaving them out would be right by luck and wrong the moment either is set.
      d->setScale(s->getScale());
      if (s->isEnabled()) d->enable(); else d->disable();
      for (auto& t : s->_Terms) {
        EnergyStretch nt = t;                          // parameters copied verbatim
        nt.term.i3x1 = remapSlotI3(t.term.i3x1, srcLo3, srcHi3, delta);
        nt.term.i3x2 = remapSlotI3(t.term.i3x2, srcLo3, srcHi3, delta);
        nt._Atom1 = atomAtI3(at, nt.term.i3x1);
        nt._Atom2 = atomAtI3(at, nt.term.i3x2);
        d->_Terms.push_back(nt);
        copied++;
      }
    } else if (gc::IsA<EnergyAngle_sp>(c)) {
      EnergyAngle_sp s = gc::As_unsafe<EnergyAngle_sp>(c);
      EnergyAngle_sp d = ensureComponent<EnergyAngle_O>(ef, tdestGroup);
      d->setScale(s->getScale());
      if (s->isEnabled()) d->enable(); else d->disable();
      for (auto& t : s->_Terms) {
        EnergyAngle nt = t;
        nt.term.i3x1 = remapSlotI3(t.term.i3x1, srcLo3, srcHi3, delta);
        nt.term.i3x2 = remapSlotI3(t.term.i3x2, srcLo3, srcHi3, delta);
        nt.term.i3x3 = remapSlotI3(t.term.i3x3, srcLo3, srcHi3, delta);
        nt._Atom1 = atomAtI3(at, nt.term.i3x1);
        nt._Atom2 = atomAtI3(at, nt.term.i3x2);
        nt._Atom3 = atomAtI3(at, nt.term.i3x3);
        d->_Terms.push_back(nt);
        copied++;
      }
    } else if (gc::IsA<EnergyLinearAngle_sp>(c)) {
      // Same shape as EnergyAngle and equally topology-determined: a linear angle exists because
      // three atoms are colinear in the BOND GRAPH, not because of where they sit.
      EnergyLinearAngle_sp s = gc::As_unsafe<EnergyLinearAngle_sp>(c);
      EnergyLinearAngle_sp d = ensureComponent<EnergyLinearAngle_O>(ef, tdestGroup);
      d->setScale(s->getScale());
      if (s->isEnabled()) d->enable(); else d->disable();
      for (auto& t : s->_Terms) {
        EnergyLinearAngle nt = t;
        nt.term.i3x1 = remapSlotI3(t.term.i3x1, srcLo3, srcHi3, delta);
        nt.term.i3x2 = remapSlotI3(t.term.i3x2, srcLo3, srcHi3, delta);
        nt.term.i3x3 = remapSlotI3(t.term.i3x3, srcLo3, srcHi3, delta);
        nt._Atom1 = atomAtI3(at, nt.term.i3x1);
        nt._Atom2 = atomAtI3(at, nt.term.i3x2);
        nt._Atom3 = atomAtI3(at, nt.term.i3x3);
        d->_Terms.push_back(nt);
        copied++;
      }
    } else if (gc::IsA<EnergyDihedral_sp>(c)) {
      EnergyDihedral_sp s = gc::As_unsafe<EnergyDihedral_sp>(c);
      EnergyDihedral_sp d = ensureComponent<EnergyDihedral_O>(ef, tdestGroup);
      d->setScale(s->getScale());
      if (s->isEnabled()) d->enable(); else d->disable();
      for (auto& t : s->_Terms) {
        EnergyDihedral nt = t;
        nt.term.i3x1 = remapSlotI3(t.term.i3x1, srcLo3, srcHi3, delta);
        nt.term.i3x2 = remapSlotI3(t.term.i3x2, srcLo3, srcHi3, delta);
        nt.term.i3x3 = remapSlotI3(t.term.i3x3, srcLo3, srcHi3, delta);
        nt.term.i3x4 = remapSlotI3(t.term.i3x4, srcLo3, srcHi3, delta);
        nt._Atom1 = atomAtI3(at, nt.term.i3x1);
        nt._Atom2 = atomAtI3(at, nt.term.i3x2);
        nt._Atom3 = atomAtI3(at, nt.term.i3x3);
        nt._Atom4 = atomAtI3(at, nt.term.i3x4);
        d->_Terms.push_back(nt);
        copied++;
      }
    } else if (gc::IsA<EnergyNonbond14_sp>(c)) {
      // 1-4 is here rather than with the nonbond components because it is TOPOLOGY-determined:
      // construct14InteractionTerms loops PROPERS over the matter and takes each proper dihedral's
      // terminal atoms.  Same enumeration, same source as the dihedral terms - so it renames across
      // rotamers of a monomer for exactly the same reason they do.
      //
      // The reference's terms are already deduplicated on their (i3,i3) pairs, and the shift is
      // injective - a constant offset inside the source slot range, identity outside - mapping that
      // range onto the destination slot range, which is disjoint from the backbone.  So no two
      // reference pairs can collide after shifting and none can collide with a backbone pair; the
      // deduplication survives without being redone.
      EnergyNonbond14_sp s = gc::As_unsafe<EnergyNonbond14_sp>(c);
      EnergyNonbond14_sp d = ensureComponent<EnergyNonbond14_O>(ef, tdestGroup);
      d->setScale(s->getScale());
      if (s->isEnabled()) d->enable(); else d->disable();
      for (auto& t : s->_Terms14) {
        EnergyNonbond14 nt = t;
        nt.term.i3x1 = remapSlotI3(t.term.i3x1, srcLo3, srcHi3, delta);
        nt.term.i3x2 = remapSlotI3(t.term.i3x2, srcLo3, srcHi3, delta);
        nt._Atom1_enb = atomAtI3(at, nt.term.i3x1);
        nt._Atom2_enb = atomAtI3(at, nt.term.i3x2);
        d->_Terms14.push_back(nt);
        copied++;
      }
    } else {
      skipped++;
    }
  }
  return Values(core::make_fixnum(copied), core::make_fixnum(skipped));
}


void EnergyFunction_O::defineForAggregate(Aggregate_sp aggregate, bool useExcludedAtoms, core::T_sp keepInteractionFactory, bool assign_types, core::T_sp forceFieldOverrides, core::T_sp setup, core::T_sp maybe_energy_component_group )
{
  if (setup.nilp()) {
    setup = core::Cons_O::create(kw::_sym_amber,nil<core::T_O>());
  }
  if (!setup.consp()) {
    SIMPLE_ERROR("setup must be NIL or a list that starts with :amber or :rosetta");
  }
  if (oCar(setup) == kw::_sym_rosetta && useExcludedAtoms) {
    SIMPLE_ERROR("You cannot have a rosetta energy-function with use-excluded-atoms t");
  }
  this->_ForceFieldName = gc::As<core::Symbol_sp>(oCar(setup)); // we should call this something else like nonbond-function-type because it stores :amber or :rosetta
  // Keep the FULL setup and the useExcludedAtoms flag so generateIntoGroup can re-run
  // generateNonbondEnergyFunctionTables with exactly these arguments.  _ForceFieldName holds only
  // oCar(setup); the rosetta branch builds SetupAccumulator(classSymbol,setup) from the whole
  // list, so rebuilding it as (list _ForceFieldName) would drop the tail and reparameterize.
  this->_Setup = setup;
  this->_UseExcludedAtoms = useExcludedAtoms;
  core::DynamicScopeManager scope(_sym_STARparameter_warningsSTAR,nil<core::T_O>());

  this->_Matter = aggregate;
  //
  // Loop over the molecules
  //
  core::HashTable_sp molecule_force_fields = core::HashTable_O::createEq();
  core::HashTable_sp molecule_force_field_names = core::HashTable_O::createEq();

  this->resolveMoleculeForceFields( aggregate, forceFieldOverrides, molecule_force_fields, molecule_force_field_names );

  core::HashTable_sp atomTypes = core::HashTable_O::createEq();
  this->_AtomTypes = atomTypes;

  //
  // Setup the extra information needed by the energy function to generate a topology file
  //
  // Get the name to pass that to the atomTable
  core::T_sp matterName = aggregate->getName();
  this->_AtomTable->setAggregateName(matterName);
  core::T_sp boundingBox = unbound<core::T_O>();
  if (aggregate->boundingBoxBoundP()) {
    boundingBox = aggregate->boundingBox();
  }
  if (boundingBox.unboundp()) {
    this->_AtomTable->makUnboundBoundingBox();
  } else {
    this->_AtomTable->setBoundingBox(boundingBox);
  }

  if (chem__verbose(0)) core::clasp_write_string("Assembling aggregate nonbond force-field.\n");
  core::T_sp nonbondForceField = core::eval::funcall(chem::_sym_compute_merged_nonbond_force_field_for_aggregate,aggregate,atomTypes,molecule_force_field_names);
  this->_AtomTable->setNonbondForceFieldForAggregate(nonbondForceField);

  if (chem__verbose(0)) core::clasp_write_string("Assembling aggregate lksolvation force-field.\n");
  core::T_sp lksolvationForceField = core::eval::funcall(chem::_sym_compute_merged_lksolvation_force_field_for_aggregate,aggregate,molecule_force_field_names);
  this->_AtomTable->setLKSolvationForceFieldForAggregate(lksolvationForceField);

  // Separate the molecules for solute from the solvent and handle them solute first then solvent
  size_t final_solute_residue_iptres = 0;
  size_t number_of_molecules_nspm = 0;
  size_t first_solvent_molecule_nspsol = 0;
  bool solvent_exists = false;
  
  if (chem__verbose(0)) core::clasp_write_string("Classifying solute/solvent.\n");
  ql::list solute_molecules;
  ql::list solvent_molecules;
  {
    Loop moleculeLoop;
    moleculeLoop.loopTopGoal(aggregate,MOLECULES);
    while (moleculeLoop.advanceLoopAndProcess() ) {
      Molecule_sp molecule = moleculeLoop.getMolecule();
      if (molecule->molecule_type() == kw::_sym_solvent) {
        solvent_molecules << molecule;
        solvent_exists = true;
      } else {
        solute_molecules << molecule;
      }
    }
  }


  core::List_sp solute = solute_molecules.cons();
  core::List_sp solvent = solvent_molecules.cons();

  // ================================================================================
  // PASS 1 - atom types and the atom table, solute first then solvent.
  //
  // This runs for EVERY molecule before ANY molecule's bonded terms are generated.  Atom
  // numbering is fixed here and pass 2 only reads it, so pass 2 can be re-run later against a
  // different keep-interaction-factory and a fresh energy-component-group without disturbing
  // the coordinate vector - that is what generateIntoGroup does for a blueprint rotamer slot.
  //
  for ( auto cur_solute : solute ) {
    Molecule_sp onemol = gc::As_unsafe<Molecule_sp>(CONS_CAR(cur_solute));
    this->constructAtomTableForMolecule( onemol, keepInteractionFactory, atomTypes, nonbondForceField, molecule_force_fields, molecule_force_field_names );
    final_solute_residue_iptres += onemol->contentSize();
    ++number_of_molecules_nspm;
  }
  first_solvent_molecule_nspsol = number_of_molecules_nspm+1;
  if (chem__verbose(1)) {
    if (core::cl__length(solvent)>0) {
      Molecule_sp onemol = gc::As<Molecule_sp>(oCar(solvent));
      core::T_sp force_field_name = onemol->force_field_name();
      core::clasp_write_string(fmt::format("Generating parameters for solvent {} using {} force-field.\n" , _rep_(onemol->getName()) , _rep_(force_field_name) ));
    }
  }
  for ( auto cur_solvent : solvent ) {
    Molecule_sp onemol = gc::As_unsafe<Molecule_sp>(CONS_CAR(cur_solvent));
    this->constructAtomTableForMolecule( onemol, keepInteractionFactory, atomTypes, nonbondForceField, molecule_force_fields, molecule_force_field_names );
    ++number_of_molecules_nspm;
  }

  // ================================================================================
  // PASS 2 - stretch/angle/dihedral/improper terms, same order.
  //
  for ( auto cur_solute : solute ) {
    Molecule_sp onemol = gc::As_unsafe<Molecule_sp>(CONS_CAR(cur_solute));
    this->generateForMolecule( onemol, keepInteractionFactory, atomTypes, nonbondForceField, molecule_force_fields, molecule_force_field_names, maybe_energy_component_group );
  }
  for ( auto cur_solvent : solvent ) {
    Molecule_sp onemol = gc::As_unsafe<Molecule_sp>(CONS_CAR(cur_solvent));
    this->generateForMolecule( onemol, keepInteractionFactory, atomTypes, nonbondForceField, molecule_force_fields, molecule_force_field_names, maybe_energy_component_group );
  }
  if (solvent_exists) {
    printf("%s:%d  solvent_exists NSPSOL %lu  IPTRES %lu   NSM %lu \n", __FILE__, __LINE__,
           first_solvent_molecule_nspsol, final_solute_residue_iptres, number_of_molecules_nspm);
    this->_AtomTable->set_firstSolventMoleculeNSPSOL(first_solvent_molecule_nspsol);
    this->_AtomTable->set_finalSoluteResidueIPTRES(final_solute_residue_iptres);
    this->_AtomTable->set_totalNumberOfMoleculesNSPM(number_of_molecules_nspm);
  }



  if (keepInteractionFactory.notnilp()) {
    if (chem__verbose(1)) core::clasp_write_string("About to calculate nonbond and restraint terms");
    this->generateNonbondEnergyFunctionTables(useExcludedAtoms,aggregate,nonbondForceField,keepInteractionFactory,atomTypes,setup);
    this->generateRestraintEnergyFunctionTables(aggregate,nonbondForceField,keepInteractionFactory,atomTypes);
  }

  // Check if setup energy component names all match energy components in this force field.
  core::List_sp validNames = nil<core::T_O>();
  core::List_sp setupRest = oCdr(setup);
  // build a list of valid names
  for ( auto cur : this->_EnergyComponents ) {
    auto component = oCar(cur);
    validNames = core::Cons_O::create(component->static_classSymbol(),validNames);
  }
  if (validNames.consp()) {
    core::Cons_sp cvalidNames = gc::As_unsafe<core::Cons_sp>(validNames);
    for ( auto cur : setupRest ) {
      core::List_sp componentSetup = gc::As<core::Cons_sp>(oCar(cur));
      core::Symbol_sp sym = gc::As<core::Symbol_sp>(oCar(componentSetup));
      if (!cvalidNames->memberEq(sym)) {
        SIMPLE_ERROR("Illegal energy-component name {} - expected one of {}", _rep_(sym), _rep_(validNames));
      }
    }
  }
  core::eval::funcall(_sym_report_parameter_warnings);
}




// PASS 1 for one molecule - assign its atom types and append its atoms to the atom table.
void EnergyFunction_O::constructAtomTableForMolecule( Molecule_sp onemol, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes, core::T_sp nonbondForceField, core::HashTable_sp molecule_force_fields, core::HashTable_sp molecule_force_field_names )
{
  core::eval::funcall( _sym_construct_atom_table_for_molecule,
                       this->asSmartPtr(), onemol,
                       molecule_force_fields->gethash(onemol),
                       molecule_force_field_names->gethash(onemol),
                       atomTypes, nonbondForceField, keepInteractionFactory );
}

// PASS 2 for one molecule - generate its bonded terms.  GROUP is NIL or an
// EnergyComponentGroup that collects every component this pass creates.
void EnergyFunction_O::generateForMolecule( Molecule_sp onemol, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes, core::T_sp nonbondForceField, core::HashTable_sp molecule_force_fields, core::HashTable_sp molecule_force_field_names, core::T_sp group )
{
  core::eval::funcall( _sym_generate_for_molecule_using_force_field,
                       this->asSmartPtr(), onemol,
                       molecule_force_fields->gethash(onemol),
                       molecule_force_field_names->gethash(onemol),
                       atomTypes, nonbondForceField, keepInteractionFactory, group );
}



core::HashTable_sp createAtomToResidueHashTable(Matter_sp molecule)
{
  Loop lRes;
  Loop lAtoms;
  core::HashTable_sp ht = core::HashTable_O::createEq();
  lRes.loopTopGoal(molecule,RESIDUES);
  while (lRes.advanceLoopAndProcess()) {
    Residue_sp res = lRes.getResidue();
    lAtoms.loopTopGoal(res,ATOMS);
    while (lAtoms.advanceLoopAndProcess()) {
      Atom_sp atom = lAtoms.getAtom();
      ht->setf_gethash(atom,res);
    }
  }
  return ht;
}

CL_LAMBDA((energy-function chem:energy-function) molecule ffstretches ffangles ffptors ffitors &key (keep-interaction-factory t) atom-types group)
CL_DOCSTRING(R"dx(Generate the standard energy function tables. The atom types, and CIP priorities need to be precalculated.

GROUP is NIL or an EnergyComponentGroup.  NIL puts the stretch/angle/dihedral terms in the energy
function's own single component of each class.  A group scopes find-or-create to that group, so a
pass driven with a fresh group yields a fresh set of components holding only that pass's terms -
which is how one blueprint rotamer slot gets an energy of its own.)dx")
CL_DEFMETHOD void EnergyFunction_O::generateStandardEnergyFunctionTables(Matter_sp molecule,
                                                                         FFStretchDb_sp ffstretches,
                                                                         FFAngleDb_sp ffangles,
                                                                         FFPtorDb_sp ffptors,
                                                                         FFItorDb_sp ffitors,
                                                                         core::T_sp keepInteractionFactory,
                                                                         core::HashTable_sp atomTypes,
                                                                         core::T_sp group )
{
  if (keepInteractionFactory.nilp()) return;
  Loop loop;
  Atom_sp          a1, a2, a3, a4, aImproperCenter;
  core::Symbol_sp  t1, t2, t3, t4, t141, t144;
  EnergyAtom       *ea1, *ea2, *ea3, *ea4;
  FFPtor_sp        ffPtor;
  FFItor_sp        ffItor;
  FFNonbond_sp	   ffNonbond1, ffNonbond2;
  int              coordinateIndex;

  core::HashTable_sp atomToRes = createAtomToResidueHashTable(molecule);
  //
  // Initialize the energy components
  //
  this->_eraseMissingParameters();
  coordinateIndex = 0;
  ASSERTNOTNULL(forceField);
  // Search the stretch terms
  {
    size_t terms = 0;
    size_t missing_terms = 0;
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyStretch_O::staticClass() );
    if (keepInteraction.notnilp()) {
      auto stretchComponent = ensureComponent<EnergyStretch_O>(this->asSmartPtr(),group);
      loop.loopTopGoal(molecule,BONDS);
      while ( loop.advanceLoopAndProcess() ) {
        a1 = loop.getBondA1();
        a2 = loop.getBondA2();
        if ( skipInteraction( keepInteraction, a1, a2 ) ) continue;
        //      printf("%s:%d Looking at STRETCH term between %s - %s\n", __FILE__, __LINE__, _rep_(a1).c_str(), _rep_(a2).c_str());
        t1 = a1->getType(atomTypes).as<core::Symbol_O>();
        t2 = a2->getType(atomTypes).as<core::Symbol_O>();
        ea1 = this->getEnergyAtomPointer(a1);
        ea2 = this->getEnergyAtomPointer(a2);
        FFStretch_sp ffStretch = gc::As<FFStretch_sp>(ffstretches->findTermForTypes(t1,t2));
        if ( ffStretch->level() != parameterized ) {
          this->_addMissingParameter(ffStretch);
          ++missing_terms;
        }
        if ( ffStretch->level() != unknown ) {
          EnergyStretch   energyStretch;
          energyStretch.defineFrom(ffStretch,ea1,ea2,stretchComponent->getScale());
          stretchComponent->addTerm(energyStretch);
          ++terms;
        } else {
          Residue_sp res1 = gc::As<Residue_sp>(atomToRes->gethash(a1));
          Residue_sp res2 = gc::As<Residue_sp>(atomToRes->gethash(a2));
          SIMPLE_WARN("Could not find stretch parameter in molecule {} between {}/{} (atom type {}) and {}/{} (atom type {})" , _rep_(molecule) , _rep_(res1) , _rep_(a1) , _rep_(a1->getType(atomTypes)) , _rep_(res2) , _rep_(a2) , _rep_(a2->getType(atomTypes)));
        }
      }
      if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built stretch table with {} terms added and {} missing terms\n" , terms , missing_terms));
    }
  }
  // Search the angle terms
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyAngle_O::staticClass() );
    if (keepInteraction.notnilp()) {
      auto angleComponent = ensureComponent<EnergyAngle_O>(this->asSmartPtr(),group);
      size_t terms = 0;
      size_t missing_terms = 0;
      loop.loopTopGoal(molecule,ANGLES);
      while ( loop.advanceLoopAndProcess() )
        {
          a1 = loop.getAtom1();
          a2 = loop.getAtom2();
          a3 = loop.getAtom3();
          if ( skipInteraction( keepInteraction, a1, a2, a3 ) ) continue;
          ea1 = this->getEnergyAtomPointer(a1);
          ea2 = this->getEnergyAtomPointer(a2);
          ea3 = this->getEnergyAtomPointer(a3);
          FFAngle_sp ffAngle = ffangles->findTerm(ffstretches,a1,a2,a3,atomTypes);
          if ( ffAngle->level() != parameterized ) {
            LOG("Missing angle parameter between types: {}-{}-{}" , _rep_(a1->getType(atomTypes)) , _rep_(a2->getType(atomTypes)) , _rep_(a3->getType(atomTypes)) );
            this->_addMissingParameter(ffAngle);
            LOG("Added to missing parameters" );
            ++missing_terms;
          }
          if ( ffAngle->level() != unknown ) {
            EnergyAngle energyAngle;
            energyAngle.defineFrom(ffAngle,ea1,ea2,ea3,angleComponent->getScale());
            angleComponent->addTerm(energyAngle);
            ++terms;
          }
        }
      if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built angle table with {} terms and {} missing terms\n" , terms , missing_terms));
    }
  }
  // Search the ptor terms
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyDihedral_O::staticClass() );
    if (keepInteraction.notnilp()) {
      auto dihedralComponent = ensureComponent<EnergyDihedral_O>(this->asSmartPtr(),group);
      size_t terms = 0;
      size_t missing_terms = 0;
      loop.loopTopGoal(molecule,PROPERS);
      while ( loop.advanceLoopAndProcess() ) {
        a1 = loop.getAtom1();
        a2 = loop.getAtom2();
        a3 = loop.getAtom3();
        a4 = loop.getAtom4();
        if ( skipInteraction( keepInteraction, a1, a2, a3, a4 ) ) continue;
        t1 = a1->getType(atomTypes).as<core::Symbol_O>();
        t2 = a2->getType(atomTypes).as<core::Symbol_O>();
        t3 = a3->getType(atomTypes).as<core::Symbol_O>();
        t4 = a4->getType(atomTypes).as<core::Symbol_O>();
        //      core::clasp_write_string(fmt::format("atoms types: {}-{}-{}-{} \n" , t1 , t2 , t3 , t4));
        ea1 = this->getEnergyAtomPointer(a1);
        ea2 = this->getEnergyAtomPointer(a2);
        ea3 = this->getEnergyAtomPointer(a3);
        ea4 = this->getEnergyAtomPointer(a4);
        core::T_sp tffPtor = ffptors->findBestTerm(t1,t2,t3,t4);
        if (tffPtor.notnilp()) {
          FFPtor_sp ffPtor = gc::As<FFPtor_sp>(tffPtor);
          int numPtors = 0;
          for ( int n=1;n<=FFPtor_O::MaxPeriodicity; n++ ) {
            if ( ffPtor->hasPeriodicity(n) ) {
              ++numPtors;
              if (chem__verbose(1)) {
                Residue_sp res1 = gc::As<Residue_sp>(atomToRes->gethash(a1));
                Residue_sp res2 = gc::As<Residue_sp>(atomToRes->gethash(a2));
                Residue_sp res3 = gc::As<Residue_sp>(atomToRes->gethash(a3));
                Residue_sp res4 = gc::As<Residue_sp>(atomToRes->gethash(a4));
                core::clasp_write_string(fmt::format( "Adding proper term for atoms {}-{}-{}-{} types: {}-{}-{}-{} -> {}/n"
                                                      , ea1->getResidueAndName(res1)
                                                      , ea2->getResidueAndName(res2)
                                                      , ea3->getResidueAndName(res3)
                                                      , ea4->getResidueAndName(res4)
                                                      , core::_rep_(t1) , core::_rep_(t2) , core::_rep_(t3) , core::_rep_(t4)
                                                      , _rep_(ffPtor)
                                                      ));
              }
              EnergyDihedral energyDihedral;
              energyDihedral.defineFrom(n,ffPtor,ea1,ea2,ea3,ea4,dihedralComponent->getScale());
              dihedralComponent->addTerm(energyDihedral);
              ++terms;
            }
          }
          if ( numPtors == 0 ) {
            FFPtor_sp ptor = FFPtor_O::create_missing(t1,t2,t3,t4);
            this->_addMissingParameter(ptor);
            ++missing_terms;
          }
          if (ea1->inBondOrAngle(ea4->atom()) )
            {
#if 0
#ifdef	DEBUG_ON
              if ( _rep_(t1) < _rep_(t4) ) {
                t141 = t1;
                t144 = t4;
              }else{
                t141 = t4;
                t144 = t1;
              }
              //          LOG("Defining 1-4 interaction %-9s- %-9s   " , t1 , t4 );
              //          core::clasp_write_string(fmt::format("Defining 1-4 interaction {:<9}- {:<9}   \n" , t1 , t4 ));
#endif
#endif
            } else {
#if 0          
#ifdef	DEBUG_ON
            if ( _rep_(t1) < _rep_(t4) ) {
              t141 = t1;
              t144 = t4;
            }else{
              t141 = t4;
              t144 = t1;
            }
            //          LOG("Ignoring 1-4 interaction %-9s- %-9s    " , t1 , t4 );
            //          core::clasp_write_string(fmt::format("Ignoring 1-4 interaction {:<9}- {:<9}    \n" , t1 , t4 ));
#endif
#endif
          }
          //		ea1->_CloserThan15.insert(ea4->_Atom);
          //		ea4->_CloserThan15.insert(ea1->_Atom);
        }
      }
      if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built dihedral table with {} terms and {} missing terms\n" , terms , missing_terms));
      //
      // Search the itor terms
      //
      EnergyDihedral energyDihedral;
      terms = 0;
      loop.loopTopGoal(molecule,IMPROPERS);
      while ( loop.advanceLoopAndProcess() ) {
        a1 = loop.getAtom1();
        a2 = loop.getAtom2();
        a3 = loop.getAtom3();
        a4 = loop.getAtom4();
        FFItorDb_O::improperAtomSort(a1,a2,a4,atomTypes);
        if ( skipInteraction( keepInteraction, a1, a2, a3, a4 ) ) continue;
        t1 = a1->getType(atomTypes).as<core::Symbol_O>();
        t2 = a2->getType(atomTypes).as<core::Symbol_O>();
        t3 = a3->getType(atomTypes).as<core::Symbol_O>();
        t4 = a4->getType(atomTypes).as<core::Symbol_O>();
        ea1 = this->getEnergyAtomPointer(a1);
        ea2 = this->getEnergyAtomPointer(a2);
        ea3 = this->getEnergyAtomPointer(a3);
        ea4 = this->getEnergyAtomPointer(a4);
        core::T_sp itor = ffitors->findBestTerm(t1,t2,t3,t4);
        if ( itor.notnilp() ) {
          //
          // Only one improper per central atom.
          // We may not get exactly the same improper as AMBER does
          //
          if ( a3 != aImproperCenter ) {
            ffItor = gc::As<FFItor_sp>(ffitors->findBestTerm(t1,t2,t3,t4));
            for ( int n=1;n<=6; n++ ) {
              if ( ffItor->hasPeriodicity(n) ) {
                energyDihedral.defineFrom(n,ffItor,ea1,ea2,ea3,ea4,dihedralComponent->getScale());
                dihedralComponent->addTerm(energyDihedral);
                aImproperCenter = a3;
                ++terms;
              }
            }
          }
        }
      }
      if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built improper table for {} terms\n" , terms));
    }
  }
  //  this->summarizeTerms();
}

SYMBOL_EXPORT_SC_(ChemPkg,prepare_amber_energy_nonbond);

CL_DOCSTRING(R"dx(Generate the nonbond energy function tables. The atom types, and CIP priorities need to be precalculated.)dx");
CL_DEFMETHOD void EnergyFunction_O::generateNonbondEnergyFunctionTables(bool useExcludedAtoms, Matter_sp matter, core::T_sp nonbondForceField, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes, core::T_sp setup, core::T_sp group, bool bonded14 )
{
  if (keepInteractionFactory.nilp()) return;
  if (chem__verbose(0))
    core::clasp_write_string(fmt::format("Built atom table for {} atoms\n" , this->_AtomTable->getNumberOfAtoms()));
#ifdef	DEBUG_DEFINE_ENERGY
  core::clasp_write_string(fmt::format("{}:{} There were {} atoms\n" , __FILE__ , __LINE__ , this->_AtomTable.size() ));
#endif
  // Nonbonds here!!!!!!!!!!!!!!
  if (oCar(setup)==kw::_sym_amber) {
    EnergyNonbond_sp nonbond;
    if (this->_BoundingBox.boundp()) {
      nonbond = EnergyPeriodicBoundaryConditionsNonbond_O::create();
      nonbond->initialize();
    } else {
      nonbond = EnergyNonbond_O::create();
      nonbond->initialize();
    }
    collectComponent(this->asSmartPtr(),nonbond,group);
    if (useExcludedAtoms) {
      SIMPLE_ERROR(":use-excluded-atoms is T - Don't use excluded atoms for the time being");
      // The nonbond parameters are calculated in Common Lisp

      core::List_sp parts = core::eval::funcall(_sym_prepare_amber_energy_nonbond,this->asSmartPtr(),nonbondForceField);
      nonbond->constructNonbondTermsFromAList(parts);
      nonbond->constructExcludedAtomListFromAtomTable(this->_AtomTable, nonbondForceField,keepInteractionFactory);
      nonbond->construct14InteractionTerms(this->_AtomTable,matter,nonbondForceField,keepInteractionFactory,atomTypes);
    } else {
      nonbond->constructNonbondTermsFromAtomTable(this->_AtomTable, nonbondForceField,atomTypes, keepInteractionFactory );
      // BONDED14 - see the header.  The 1-4 set comes off the bond graph, not the atom table.
      if (bonded14) {
        nonbond->construct14InteractionTerms(this->_AtomTable,matter,nonbondForceField,keepInteractionFactory,atomTypes);
      }
    }
  } else if (oCar(setup)==kw::_sym_rosetta) {
    {
      core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyRosettaNonbond_O::staticClass() );
      if (keepInteraction.notnilp()) {
        SetupAccumulator setupAccNonbond(EnergyRosettaNonbond_O::static_classSymbol(),setup);
        // keepInteractionFactory, NOT the specialized keepInteraction above.  ::make runs
        // specializeKeepInteractionFactory itself (energyRosettaNonbond.cc:66) and stores the
        // result in _KeepInteractionFactory, so handing it the already-specialized predicate
        // specializes TWICE - the predicate gets funcalled with the component class as its only
        // argument.  A &rest predicate accepts that and returns NIL, which trips ::make's own
        // "Mismatch between keepInteractionFactory ... and ::make" guard.  The other three
        // components below already pass the factory; this line was the odd one out, and it only
        // ever worked because the design path's factory is T, for which specialize is idempotent.
        auto energyRosettaNonbond = EnergyRosettaNonbond_O::make(this->asSmartPtr(), keepInteractionFactory, setupAccNonbond );
        collectComponent(this->asSmartPtr(),gc::As<EnergyRosettaNonbond_sp>(energyRosettaNonbond),group);
      }
    }
    {
      core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyRosettaElec_O::staticClass() );
      if (keepInteraction.notnilp()) {
        SetupAccumulator setupAccElec(EnergyRosettaElec_O::static_classSymbol(),setup);
        auto energyRosettaElec = EnergyRosettaElec_O::make(this->asSmartPtr(), keepInteractionFactory, setupAccElec );
        collectComponent(this->asSmartPtr(),gc::As<EnergyRosettaElec_sp>(energyRosettaElec),group);
      }
    }
    {
      core::T_sp keepInteraction = specializeKeepInteractionFactory(
          keepInteractionFactory, EnergyRosettaLKSolvation_O::staticClass());
      if (keepInteraction.notnilp()) {
        SetupAccumulator setupAccLK(EnergyRosettaLKSolvation_O::static_classSymbol(), setup);
        auto energyLKSolvation = EnergyRosettaLKSolvation_O::make(
            this->asSmartPtr(), keepInteractionFactory, setupAccLK);
        collectComponent(this->asSmartPtr(),gc::As<EnergyRosettaLKSolvation_sp>(energyLKSolvation),group);
      }
    }
    // BONDED14 - see the header.  construct14InteractionTerms loops PROPERS over the MATTER and
    // takes each proper dihedral's terminal atoms, so this is the same enumeration, over the same
    // structure, that produces the dihedral terms.  It is topology-determined in exactly the way the
    // bonded terms are, which is why it is gated with them rather than with the nonbond components
    // above - those are built from _AtomTable and never touch the matter.
    if (bonded14) {
      core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyDihedral_O::staticClass() );
      if (keepInteraction.notnilp()) {
        SetupAccumulator setupAcc14(EnergyNonbond14_O::static_classSymbol(),setup);
        auto energyNonbond14 = EnergyNonbond14_O::make(this->asSmartPtr(),keepInteractionFactory,setupAcc14);
        energyNonbond14->construct14InteractionTerms(this->_AtomTable,matter,nonbondForceField,keepInteractionFactory,atomTypes);
        collectComponent(this->asSmartPtr(),energyNonbond14,group);
      }
    }
  } else {
    SIMPLE_ERROR("Provide a valid setup (:rosetta ...) or (:amber ...)");
  }
}


CL_LAMBDA((energy-function chem:energy-function) matter force-field &key (keep-interaction-factory t) atom-types);
CL_DOCSTRING(R"dx(Generate the restraint energy function tables. The atom types, and CIP priorities need to be precalculated.
This should be called after generateStandardEnergyFunctionTables.)dx")
CL_DEFMETHOD void EnergyFunction_O::generateRestraintEnergyFunctionTables(Matter_sp matter, core::T_sp ffNonbond, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes ) {
  Loop loop;
  Atom_sp          a1, a2, a3, a4, aImproperCenter;
  core::Symbol_sp  t1, t2, t3, t4, t141, t144;
  EnergyAtom      *eaCenter, *ea1, *ea2, *ea3, *ea4;
  FFAngle_sp       ffAngle;
  FFPtor_sp        ffPtor;
  FFItor_sp        ffItor;
  FFNonbond_sp	ffNonbond1, ffNonbond2;
  if (chem__verbose(1)) core::lisp_write(fmt::format("In generateRestraintEnergyFunctionTables keepInteractionFactory = {}\n", keepInteractionFactory) );
  if (keepInteractionFactory.nilp()) return;
  //
  // Setup the atom chiral restraints
  //
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyChiralRestraint_O::staticClass() );
    if (keepInteraction.notnilp()) {
      auto chiralRestraintComponent = ensureComponent<EnergyChiralRestraint_O>(this->asSmartPtr());
      auto cip_priorities = CipPrioritizer_O::assignPrioritiesHashTable(matter);
      if (!gc::IsA<core::HashTable_sp>(cip_priorities)) {
        SIMPLE_ERROR("You need to provide a hash-table of atoms to relative CIP priorities - see CipPrioritizer_O::assignPrioritiesHashTable(matter)");
      }
      core::HashTable_sp cip = gc::As_unsafe<core::HashTable_sp>(cip_priorities);
      EnergyChiralRestraint	ichiral;
      Atom_sp	n1,n2,n3,n4;
      string	s1,s2,s3,s4;
      double	side ;
      Loop residue_loop;
      residue_loop.loopTopGoal(matter,RESIDUES);
      while (residue_loop.advanceLoopAndProcess()) {
        Residue_sp res = residue_loop.getResidue();
        loop.loopTopGoal(res,ATOMS);
        while ( loop.advanceLoopAndProcess() ) {
          a1 = loop.getAtom();
          if (chem__verbose(1)) core::clasp_write_string(fmt::format("Looking to assign stereochemical restraint for {}\n" , _rep_(a1)));
          if ( a1->getStereochemistryType() != undefinedCenter ) {
            if (chem__verbose(1)) core::clasp_write_string(fmt::format("getStereochemistryType != undefinedCenter for {}\n" , _rep_(a1)));
            LOG("Create a chiral restraint for {}" , a1->description()  );
            //
            // Figure out what the desired configuration should be
            // If it has been set then use that
            // Otherwise if its pro-chiral center then set it to "R"
            // If its a chiral center and it hasn't been set then throw an exception
            //
            // Configuration is set in variable side
            // "R" -> side = 1.0
            // "S" -> side = -1.0
            //
            // Calculate the angle for the chiral restraint based on
            // the stereochemistry that we want to impose
            // default R-stereochemistry (1-center)x(2-center).(3-center) is POSITIVE
            //
            side = 1.0;
            core::List_sp priority;
            if ( a1->getConfiguration() != undefinedConfiguration )
              {
                if ( a1->getConfiguration() == R_Configuration ) {
                  side = 1.0;
                } else if ( a1->getConfiguration() == S_Configuration ) {
                  side = -1.0;
                } else if ( a1->getConfiguration() == RightHanded_Configuration ) {
                  side = -1.0;
                } else if ( a1->getConfiguration() == LeftHanded_Configuration ) {
                  side = 1.0;
                }
              } else {
              if ( a1->getStereochemistryType() == prochiralCenter ) {
                side = 1.0;
              } else {
                SIMPLE_WARN("Chiral center ({}:{}) with configuration settings[{}] doesn't have its configuration set"
                            , _rep_(res)
                            , _rep_(a1)
                            , a1->getConfigurationAsString() );
              }
            }
            if ( a1->getConfiguration() == R_Configuration
                 || a1->getConfiguration() == S_Configuration ) {
              priority = a1->getNeighborsByRelativePriority(cip);
            } else if (a1->getConfiguration() == RightHanded_Configuration
                       || a1->getConfiguration() == LeftHanded_Configuration) {
              priority = a1->getNeighborsForAbsoluteConfiguration();
            } else {
              priority = a1->getNeighborsForAbsoluteConfiguration();
            }
            if (core::cl__length(priority)!=4) {
              SIMPLE_ERROR("There must be 4 neighbors of {} - but there is only {}" , _rep_(a1) , _rep_(priority));
            }
            if (chem__verbose(1)) core::clasp_write_string(fmt::format("Assigning stereochemistry for central atom {} neighbors: {}\n" , _rep_(a1) , _rep_(priority)));
            n1 = gc::As<Atom_sp>(oFirst(priority));
            n2 = gc::As<Atom_sp>(oSecond(priority));
            n3 = gc::As<Atom_sp>(oThird(priority));
            n4 = gc::As<Atom_sp>(oFourth(priority));
#if 0
            s1 = a1->getConfigurationPriorityHighest();
            s2 = a1->getConfigurationPriorityHigh();
            s3 = a1->getConfigurationPriorityLow();
            s4 = a1->getConfigurationPriorityLowest();
            n1 = a1->bondedNeighborWithName(s1);
            n2 = a1->bondedNeighborWithName(s2);
            n3 = a1->bondedNeighborWithName(s3);
            n4 = a1->bondedNeighborWithName(s4);
            ASSERTNOTNULLP(n1, "Atom("+a1->getName()+") does not have neighbor1("+s1+")");
            ASSERTNOTNULLP(n2, "Atom("+a2->getName()+") does not have neighbor2("+s2+")");
            ASSERTNOTNULLP(n3, "Atom("+a3->getName()+") does not have neighbor3("+s3+")");
            ASSERTNOTNULLP(n4, "Atom("+a4->getName()+") does not have neighbor4("+s4+")");
#endif
            eaCenter = this->getEnergyAtomPointer(a1);
            ea1 = this->getEnergyAtomPointer(n1);
            ea2 = this->getEnergyAtomPointer(n2);
            ea3 = this->getEnergyAtomPointer(n3);
            ea4 = this->getEnergyAtomPointer(n4);

            //
            // Setup chiral restraints for 1->2->center->3
            //			and 1->2->center->4
            //
            if (chem__verbose(1)) core::clasp_write_string(fmt::format("Assigning stereochemistry for central atom {} neighbors: {}\n" , _rep_(a1) , _rep_(priority)));
            ichiral._Atom1 = ea1->atom();
            ichiral._Atom2 = ea2->atom();
            ichiral._Atom3 = eaCenter->atom();
            ichiral._Atom4 = ea3->atom();
            if ( !skipInteraction( keepInteraction, ichiral._Atom1, ichiral._Atom2, ichiral._Atom3, ichiral._Atom4 ) ) {
              ichiral.term.i3x1 =ea1->coordinateIndexTimes3();
              ichiral.term.i3x2 =ea2->coordinateIndexTimes3();
              ichiral.term.i3x3 =eaCenter->coordinateIndexTimes3();
              ichiral.term.i3x4 =ea3->coordinateIndexTimes3();
              ichiral.term.k =DefaultChiralRestraintWeight * side;
              ichiral.term.co =DefaultChiralRestraintOffset;
              chiralRestraintComponent->addTerm(ichiral);
              // Now apply it to the other atom
              // on the chiral center, just flip the sign
              // of K
            }
          
            ichiral._Atom4 = ea4->atom();
            if ( !skipInteraction( keepInteraction, ichiral._Atom1, ichiral._Atom2, ichiral._Atom3, ichiral._Atom4 ) ) {
              ichiral.term.i3x4 =ea4->coordinateIndexTimes3();
              // flip the sign of the chiral restraint
              ichiral.term.k =DefaultChiralRestraintWeight * side * -1.0;
              ichiral.term.co =DefaultChiralRestraintOffset;
              chiralRestraintComponent->addTerm(ichiral);
            }
            // To try and increase the number of molecules that
            // minimize into the correct configuration I'll add another
            // set of chiral restraints for this atom
            //
            // Setup chiral restraints for 2->4->center->3
            //			and 2->4->center->1
            //
            ichiral._Atom1 = ea2->atom();
            ichiral._Atom2 = ea4->atom();
            ichiral._Atom3 = eaCenter->atom();
            ichiral._Atom4 = ea3->atom();
            if ( !skipInteraction( keepInteraction, ichiral._Atom1, ichiral._Atom2, ichiral._Atom3, ichiral._Atom4 ) ) {
              ichiral.term.i3x1 =ea2->coordinateIndexTimes3();
              ichiral.term.i3x2 =ea4->coordinateIndexTimes3();
              ichiral.term.i3x3 =eaCenter->coordinateIndexTimes3();
              ichiral.term.i3x4 =ea3->coordinateIndexTimes3();
              ichiral.term.k =DefaultChiralRestraintWeight * side;
              ichiral.term.co =DefaultChiralRestraintOffset;
              chiralRestraintComponent->addTerm(ichiral);
            }
            // Now apply it to the other atom
            // on the chiral center, just flip the sign
            // of K
            ichiral._Atom4 = ea1->atom();
            if ( !skipInteraction( keepInteraction, ichiral._Atom1, ichiral._Atom2, ichiral._Atom3, ichiral._Atom4 ) ) {
              ichiral.term.i3x4 =ea1->coordinateIndexTimes3();
              // flip the sign of the chiral restraint
              ichiral.term.k =DefaultChiralRestraintWeight * side * -1.0;
              ichiral.term.co =DefaultChiralRestraintOffset;
              chiralRestraintComponent->addTerm(ichiral);
            }
          } else {
            LOG("There is no chiral restraint for: {}" , a1->description()  );
          }
        }
      }
      if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built chiral restraints table for {} terms\n" , chiralRestraintComponent->numberOfTerms()));
    }
  }

  //
  // Set up the anchor restraints
  //
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyAnchorRestraint_O::staticClass() );
    if (keepInteraction.notnilp()) {
      auto component = ensureComponent<EnergyAnchorRestraint_O>(this->asSmartPtr());
    }
  }

  //
  // Set up the dihedral restraints
  //
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyDihedralRestraint_O::staticClass() );
    if (keepInteraction.notnilp()) {
      auto component = ensureComponent<EnergyDihedralRestraint_O>(this->asSmartPtr());
    }
  }

  //
  // Set up the fixedNonbond restraints
  //
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyFixedNonbondRestraint_O::staticClass() );
    if (keepInteraction.notnilp()) {
      auto component = ensureComponent<EnergyFixedNonbondRestraint_O>(this->asSmartPtr());
    }
  }

  //
  // Define secondary amide restraints
  //
#if 0
  if ( this->_RestrainSecondaryAmides )
    {
      gctools::Vec0<Atom_sp>	nitrogens;
      for ( AtomTable_O::iterator it=this->_AtomTable->begin(); it!=this->_AtomTable->end(); it++ )
        {
          Atom_sp a = it->atom();
          if (a->getElement() == element_N ) {
            nitrogens.push_back(a);
          }
        }
      //int startTerms = this->getDihedralRestraintComponent()->numberOfTerms();
      // this->__createSecondaryAmideRestraints(nitrogens,keepInteraction);
      //if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built secondary amide restraints including {} terms\n" , (this->getDihedralRestraintComponent()->numberOfTerms() - startTerms)));
    } else
    {
      LOG("Skipping Secondary amide restraints because _RestrainSecondaryAmides = {}" , this->_RestrainSecondaryAmides );
    }
#endif
  //
  // Set up force-field restraints
  //
  {
    IterateRestraints_sp restraintIt = IterateRestraints_O::create(matter);
    int terms = this->_applyRestraints(ffNonbond,restraintIt,keepInteractionFactory,atomTypes);
    if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built restraints including {} terms\n" , terms ));
  }
  LOG("Done terms" );
}



CL_LISPIFY_NAME("addTermsForListOfRestraints");
CL_DEFMETHOD void	EnergyFunction_O::addTermsForListOfRestraints(ForceField_sp forceField, core::List_sp restraintList, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes)
{
  if (keepInteractionFactory.nilp()) return;
  adapt::IterateCons_sp	iterate;
  iterate = adapt::IterateCons_O::create(restraintList);
  this->_applyRestraints(forceField->getNonbondDb(),iterate,keepInteractionFactory,atomTypes);
}

CL_DOCSTRING(R"doc(Write the coordinates into the nvector from the atom positions of the structure that this energy-function is based on.)doc");
CL_LISPIFY_NAME(EnergyFunction/loadCoordinatesIntoVector);
CL_DEFMETHOD
void	EnergyFunction_O::loadCoordinatesIntoVector(NVector_sp pos)
{
  int                             ci;
  gctools::Vec0<EnergyAtom>::iterator    ai;
  if ( pos->size() != this->getNVectorSize()) {
    SIMPLE_ERROR("NVector is the incorrect length");
  }
  for ( ai=this->_AtomTable->begin(); ai!=this->_AtomTable->end(); ai++ ) {
    ci = ai->coordinateIndexTimes3();
    pos->setElement(ci,ai->atom()->getPosition().getX());
    pos->setElement(ci+1, ai->atom()->getPosition().getY());
    pos->setElement(ci+2, ai->atom()->getPosition().getZ());
  }
}


CL_DOCSTRING(R"dx(Return the coordinate vector of the energy function)dx");
DOCGROUP(cando);
CL_DEFUN NVector_sp chem__energy_function_coordinate_vector(EnergyFunction_sp energy_function)
{
  NVector_sp pos = NVector_O::create(energy_function->getNVectorSize());
  energy_function->loadCoordinatesIntoVector(pos);
  return pos;
}

CL_DOCSTRING(R"dx(Return an empty force vector of the energy function)dx");
DOCGROUP(cando);
CL_DEFUN NVector_sp chem__energy_function_empty_force_vector(EnergyFunction_sp energy_function)
{
  NVector_sp pos = NVector_O::make(energy_function->getNVectorSize(),0.0,true);
  return pos;
}

CL_DOCSTRING(R"doc(Write the coordinates in the nvector into the atom positions)doc");
CL_LISPIFY_NAME(EnergyFunction/saveCoordinatesFromVector);
CL_DEFMETHOD
void    EnergyFunction_O::saveCoordinatesFromVector(NVector_sp pos)
{
  int                             ci;
  double                          x,y,z;
  gctools::Vec0<EnergyAtom>::iterator    ai;
  Vector3                         v;
  for ( ai=this->_AtomTable->begin(); ai!=this->_AtomTable->end(); ai++ ) {
    ci = ai->coordinateIndexTimes3();
    x = pos->getElement(ci+0);
    y = pos->getElement(ci+1);
    z = pos->getElement(ci+2);
    LOG("Set atom({}) position ({},{},{})" , ci , x , y , z  );
    v.set(x,y,z);
    ai->atom()->setPosition(v);
  }
}

CL_LISPIFY_NAME("writeForceToAtoms");
CL_DEFMETHOD void    EnergyFunction_O::writeForceToAtoms(NVector_sp force)
{
  IMPLEMENT_ME();
  int                             ci;
  double                          x,y,z;
  gctools::Vec0<EnergyAtom>::iterator    ai;
  Vector3                         v;
  for ( ai=this->_AtomTable->begin(); ai!=this->_AtomTable->end(); ai++ ) {
    ci = ai->coordinateIndexTimes3();
    x = force->getElement(ci+0);
    y = force->getElement(ci+1);
    z = force->getElement(ci+2);
    v.set(x,y,z);
    //    ai->atom()->setForce(v);
  }
}


void    EnergyFunction_O::saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force)
{
  this->saveCoordinatesFromVector(pos);
  //  this->writeForceToAtoms(force);
}


EnergyAtom*	EnergyFunction_O::getEnergyAtomPointer(Atom_sp a)
{ 
  return this->_AtomTable->getEnergyAtomPointer(a);
};





void	EnergyFunction_O::dealWithProblem(core::Symbol_sp problem, core::T_sp error_args)
{
  core::List_sp atoms = nil<core::T_O>();
  if ( error_args.consp() && CONS_CAR(error_args) == kw::_sym_atoms ) {
    atoms = atoms;
  }
  for ( auto cur : atoms ) {
    Atom_sp a = gctools::As<Atom_sp>(oCar(cur));
    a->bumpPosition(0.1);
  }
}



CL_LISPIFY_NAME("hasMissingParameters");
CL_DEFMETHOD bool	EnergyFunction_O::hasMissingParameters()
{
  return this->_MissingParameters.notnilp();
}


CL_LISPIFY_NAME("getMissingParameters");
CL_DEFMETHOD core::List_sp EnergyFunction_O::getMissingParameters()
{
  return this->_MissingParameters;
}

CL_LISPIFY_NAME("debugLogAsString");
CL_DEFMETHOD string	EnergyFunction_O::debugLogAsString()
{
  stringstream ss;
  for ( auto cur : this->_EnergyComponents ) {
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(CONS_CAR(cur));
    ss << component->debugLogAsString();
    ss << std::endl;
  }
  return ss.str();
}

SYMBOL_EXPORT_SC_(KeywordPkg,atom_table);
SYMBOL_EXPORT_SC_(KeywordPkg,stretch);
SYMBOL_EXPORT_SC_(KeywordPkg,angle);
SYMBOL_EXPORT_SC_(KeywordPkg,dihedral);
SYMBOL_EXPORT_SC_(KeywordPkg,nonbond);

DOCGROUP(cando);
CL_DEFUN void chem__fill_energy_function_from_alist(EnergyFunction_sp energy, core::List_sp alist)
{
  energy->_AtomTable = (safe_alist_lookup<AtomTable_sp>(alist,kw::_sym_atom_table));
  energy->pushEnergyComponent(safe_alist_lookup<EnergyStretch_sp>(alist,kw::_sym_stretch));
  energy->pushEnergyComponent(safe_alist_lookup<EnergyAngle_sp>(alist,kw::_sym_angle));
  energy->pushEnergyComponent(safe_alist_lookup<EnergyDihedral_sp>(alist,kw::_sym_dihedral));
  energy->pushEnergyComponent(safe_alist_lookup<EnergyNonbond_sp>(alist,kw::_sym_nonbond));
}

CL_LAMBDA((energy-function chem:energy-function) &optional (keep-interaction-factory t) setup add-components);
CL_DEFMETHOD
EnergyFunction_sp EnergyFunction_O::copyFilter(core::T_sp keepInteractionFactory, core::List_sp setup, core::List_sp addComponents )
{
  auto  me  = gctools::GC<EnergyFunction_O>::allocate();
  me->_Matter = this->_Matter;
  me->_AtomTable = this->_AtomTable;
  me->_NonbondCrossTermTable = this->_NonbondCrossTermTable;
  me->_BoundingBox = this->_BoundingBox;
  // Carry the defining setup so a filtered copy can still generateIntoGroup.  This is THIS->_Setup
  // (what defineForAggregate was called with), deliberately not the SETUP argument above - that
  // one configures the per-component SetupAccumulator for this filtering pass and is a different
  // thing that merely tends to look the same.
  me->_Setup = this->_Setup;
  me->_UseExcludedAtoms = this->_UseExcludedAtoms;
  if (keepInteractionFactory.notnilp()) {
    ql::list ll;
    for ( auto cur : this->_EnergyComponents ) {
      EnergyComponent_sp component = gc::As<EnergyComponent_sp>(CONS_CAR(cur));
      SetupAccumulator setupAcc(component->static_classSymbol(),setup);
      ll << component->copyFilter(keepInteractionFactory,setupAcc);
    }
    me->_EnergyComponents = ll.cons();
  }
  // Add empty components for classes listed in addComponents
  // that are not already present in the result.
  for ( auto cur : addComponents ) {
    core::Symbol_sp className = gc::As<core::Symbol_sp>(CONS_CAR(cur));
    if (me->findComponentOrNil(className).nilp()) {
      // Use the class's make_instance() to create an empty C++ object.
      core::Instance_sp theClass = gc::As<core::Instance_sp>(cl__find_class(className, true, nil<core::T_O>()));
      core::T_sp newComp = theClass->make_instance();
      me->pushEnergyComponent(gc::As<EnergyComponent_sp>(newComp));
    }
  }
  return me;
}

void EnergyFunction_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                  AbstractLargeSquareMatrix_sp m,
                                                  core::T_sp activeAtomMask )
{
  m->fill(0.0);
  for ( auto cur : this->_EnergyComponents ) {
    EnergyComponent_sp comp = gc::As<EnergyComponent_sp>(oCar(cur));
    comp->setupHessianPreconditioner(nvPosition,m,activeAtomMask);
  }
}


};
