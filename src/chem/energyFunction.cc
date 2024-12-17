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
#if USE_ALL_ENERGY_COMPONENTS
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyPeriodicBoundaryConditionsNonbond.h>
#include <cando/chem/energyDihedralRestraint.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyFixedNonbond.h>
#endif
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

//Maybe use this in place of InteractionError?
SYMBOL_EXPORT_SC_(ChemPkg,InteractionError);
/*! LinearAngleError has a argument :atoms (atom1 atom2 atom3) */
SYMBOL_EXPORT_SC_(ChemPkg,LinearAngleError);

SYMBOL_EXPORT_SC_(ChemPkg,LinearDihedralError);
SYMBOL_EXPORT_SC_(ChemPkg,LinearDihedralRestraintError);
SYMBOL_EXPORT_SC_(ChemPkg,OverlappingNonbondError);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);

bool energyFunctionInitialized = false;

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
  core::List_sp result = this->_OtherEnergyComponents;
  result = core::Cons_O::create(this->_FixedNonbondRestraint,result);
  result = core::Cons_O::create(this->_AnchorRestraint,result);
  result = core::Cons_O::create(this->_ChiralRestraint,result);
  result = core::Cons_O::create(this->_DihedralRestraint,result);
  result = core::Cons_O::create(this->_Nonbond,result);
  result = core::Cons_O::create(this->_Dihedral,result);
  result = core::Cons_O::create(this->_Angle,result);
  result = core::Cons_O::create(this->_Stretch,result);
  return result;
}

CL_DOCSTRING(R"doc(Create an energy-scale object for an energy-function.)doc");
CL_LISPIFY_NAME(make_energy_scale);
CL_DEF_CLASS_METHOD EnergyScale_sp EnergyScale_O::make()
{
  auto energyScale = gctools::GC<EnergyScale_O>::allocate_with_default_constructor();
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
CL_LAMBDA(&key matter disable-components enable-components (use-excluded-atoms t) (keep-interaction-factory t) (assign-types t));
CL_LISPIFY_NAME(make_energy_function);
CL_DEF_CLASS_METHOD EnergyFunction_sp EnergyFunction_O::make(core::T_sp matter, core::T_sp disableComponents, core::List_sp enableComponents,
                                                             bool useExcludedAtoms, core::T_sp keepInteractionFactory, bool assign_types )
{
  auto  me  = gctools::GC<EnergyFunction_O>::allocate_with_default_constructor();
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
  if ( matter.notnilp() ) me->defineForMatter(gc::As<Matter_sp>(matter),useExcludedAtoms,keepInteractionFactory,assign_types);
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
  this->_Stretch = EnergyStretch_O::create();
#if USE_ALL_ENERGY_COMPONENTS
  this->_Angle = EnergyAngle_O::create();
  this->_Dihedral = EnergyDihedral_O::create();
  if (this->_BoundingBox.boundp()) {
    this->_Nonbond = EnergyPeriodicBoundaryConditionsNonbond_O::create();
  } else {
    this->_Nonbond = EnergyNonbond_O::create();
  }
  this->_ChiralRestraint = EnergyChiralRestraint_O::create();
  this->_AnchorRestraint = EnergyAnchorRestraint_O::create();
  this->_DihedralRestraint = EnergyDihedralRestraint_O::create();
  this->_FixedNonbondRestraint = EnergyFixedNonbondRestraint_O::create();
#endif
  this->setScoringFunctionName(nil<core::T_O>());
  this->_Message = nil<core::T_O>();
  this->useDefaultSettings();
}

void	EnergyFunction_O::useDefaultSettings()
{
  ASSERTNOTNULL(this->_Stretch);
  ASSERTNOTNULL(this->_Stretch);
  ASSERT(this->_Stretch.notnilp());
  this->_Stretch->initialize();
#if USE_ALL_ENERGY_COMPONENTS
  ASSERTNOTNULL(this->_Angle);
  ASSERT(this->_Angle.notnilp());
  ASSERTNOTNULL(this->_Dihedral);
  ASSERT(this->_Dihedral.notnilp());
  ASSERTNOTNULL(this->_Nonbond);
  ASSERT(this->_Nonbond.notnilp());
  ASSERTNOTNULL(this->_ChiralRestraint);
  ASSERT(this->_ChiralRestraint.notnilp());
  ASSERTNOTNULL(this->_AnchorRestraint);
  ASSERT(this->_AnchorRestraint.notnilp());
  ASSERTNOTNULL(this->_DihedralRestraint);
  ASSERT(this->_DihedralRestraint.notnilp());
  ASSERTNOTNULL(this->_FixedNonbondRestraint);
  ASSERT(this->_FixedNonbondRestraint.notnilp());
  this->_Angle->initialize();
  this->_Dihedral->initialize();
  this->_Nonbond->initialize();
  this->_ChiralRestraint->initialize();
  this->_AnchorRestraint->initialize();
  this->_DihedralRestraint->initialize();
  this->_FixedNonbondRestraint->initialize();
#endif
  this->_RestrainSecondaryAmides = true;
}



void EnergyFunction_O::fields(core::Record_sp node)
{
  node->field_if_not_unbound(INTERN_(kw,AtomTable),this->_AtomTable);
  node->field_if_not_unbound(INTERN_(kw,Stretch),this->_Stretch);
#if USE_ALL_ENERGY_COMPONENTS
  node->field_if_not_unbound(INTERN_(kw,Angle),this->_Angle);
  node->field_if_not_unbound(INTERN_(kw,Dihedral),this->_Dihedral);
  node->field_if_not_unbound(INTERN_(kw,Nonbond),this->_Nonbond);
  node->field_if_not_unbound(INTERN_(kw,DihedralRestraint),this->_DihedralRestraint);
  node->field_if_not_unbound(INTERN_(kw,ChiralRestraint),this->_ChiralRestraint);
  node->field_if_not_unbound(INTERN_(kw,AnchorRestraint),this->_AnchorRestraint);
  node->field_if_not_unbound(INTERN_(kw,FixedNonbondRestraint),this->_FixedNonbondRestraint);
#endif
  node->field_if_not_unbound(INTERN_(kw,BoundingBox),this->_BoundingBox);
  node->field(INTERN_(kw,OtherEnergyComponents),this->_OtherEnergyComponents);
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
#if USE_ALL_ENERGY_COMPONENTS
    if ( onval ) this->getNonbondComponent()->enable();
    else	this->getNonbondComponent()->disable();
#endif
  }
  else if ( option == _sym_restraintAnchor )
  {
    bool onval = val.isTrue();
#if USE_ALL_ENERGY_COMPONENTS
    if ( onval ) this->getAnchorRestraintComponent()->enable();
    else	this->getAnchorRestraintComponent()->disable();
#endif
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
#define	TEST_FORCE( func, delta, argLow, argHigh, term, idx ) {\
    double eLow = this->func argLow;\
    double eHigh = this->func argHigh;\
    double numForce = -(eHigh-eLow)/(delta);\
    if ( !_areValuesClose( numForce, term, #func, #term, idx ) ) { fails++; };\
  }

#undef	TEST_DIAGONAL_HESSIAN
#define	TEST_DIAGONAL_HESSIAN( func, delta, argLow, argMiddle, argHigh, term, idx ) {\
    double eLow = this->func argLow;\
    double eMiddle = this->func argMiddle;\
    double eHigh = this->func argHigh;\
    double numHessian = (eHigh+eLow-2.0*eMiddle)/(((delta)/2.0)*((delta)/2.0));\
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
  }

#undef	TEST_OFF_DIAGONAL_HESSIAN
#define	TEST_OFF_DIAGONAL_HESSIAN( func, delta, argMM, argPM, argMP, argPP, term, idx ) {\
    double eMM = this->func argMM;\
    double eMP = this->func argMP;\
    double ePM = this->func argPM;\
    double ePP = this->func argPP;\
    double numHessian = ((ePP-ePM)-(eMP-eMM))/(delta*delta);\
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
  }




//
// Interface to amberFunction.cc
//

void EnergyFunction_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                  AbstractLargeSquareMatrix_sp m,
                                                  core::T_sp activeAtomMask )
{
  m->fill(0.0);
  if (this->_Stretch->isEnabled())
    this->_Stretch->setupHessianPreconditioner(nvPosition, m, activeAtomMask );
  if (this->_Angle->isEnabled())
    this->_Angle->setupHessianPreconditioner(nvPosition, m, activeAtomMask );
  if (this->_Dihedral->isEnabled())
    this->_Dihedral->setupHessianPreconditioner(nvPosition, m, activeAtomMask );
		// Nonbond doesn't contribute to hessian preconditioner
//    this->_Nonbond->setupHessianPreconditioner(nvPosition, m );
  if (this->_ChiralRestraint->isEnabled())
    this->_ChiralRestraint->setupHessianPreconditioner(nvPosition, m, activeAtomMask );
  if (this->_AnchorRestraint->isEnabled())
    this->_AnchorRestraint->setupHessianPreconditioner(nvPosition, m, activeAtomMask );
  if (this->_DihedralRestraint->isEnabled())
    this->_DihedralRestraint->setupHessianPreconditioner(nvPosition, m, activeAtomMask );
//    this->_FixedNonbondRestraint->setupHessianPreconditioner(nvPosition, m );
}







uint	EnergyFunction_O::countTermsBeyondThreshold()
{
  int		terms;
  terms = 0;
  SIMPLE_ERROR("Should there be something here?");
//    node = rawAccumulateTermsBeyondThresholdAsXml(terms);
  return terms;
}






//
double	EnergyFunction_O::evaluateAll( NVector_sp 	pos,
                                       core::T_sp       energyScale,
                                       core::T_sp       componentEnergy,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec,
                                       gc::Nilable<NVector_sp> dvec,
                                       core::T_sp activeAtomMask,
                                       core::T_sp debugInteractions )
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
  if (this->_Stretch->isEnabled()) {
    totalEnergy += this->_Stretch->evaluateAllComponent( this->asSmartPtr(),
                                                         pos,
                                                         energyScale,
                                                         componentEnergy,
                                                         calcForce, force,
                                                         calcDiagonalHessian,
                                                         calcOffDiagonalHessian,
                                                         hessian, hdvec, dvec, activeAtomMask,
                                                         debugInteractions );
  }
  if (this->_Angle->isEnabled()) {
    totalEnergy += this->_Angle->evaluateAllComponent( this->asSmartPtr(),
                                                       pos,
                                                       energyScale,
                                                       componentEnergy,
                                                       calcForce, force,
                                                       calcDiagonalHessian,
                                                       calcOffDiagonalHessian,
                                                       hessian, hdvec, dvec, activeAtomMask,
                                                       debugInteractions);
  }
  if(this->_Dihedral->isEnabled()) {
    totalEnergy += this->_Dihedral->evaluateAllComponent( this->asSmartPtr(),
                                                          pos,
                                                          energyScale,
                                                          componentEnergy,
                                                          calcForce, force,
                                                          calcDiagonalHessian,
                                                          calcOffDiagonalHessian,
                                                          hessian, hdvec, dvec, activeAtomMask,
                                                          debugInteractions );
  }
  if(this->_Nonbond->isEnabled()) {
    totalEnergy += this->_Nonbond->evaluateAllComponent( this->asSmartPtr(),
                                                         pos,
                                                         energyScale,
                                                         componentEnergy,
                                                         calcForce, force,
                                                         calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec, activeAtomMask,
                                                         debugInteractions );
  }
  if(this->_DihedralRestraint.boundp() && this->_DihedralRestraint->isEnabled()) {
    totalEnergy += this->_DihedralRestraint->evaluateAllComponent( this->asSmartPtr(),
                                                                   pos,
                                                                   energyScale,
                                                                   componentEnergy,
                                                                   calcForce, force,
                                                                   calcDiagonalHessian,
                                                                   calcOffDiagonalHessian,
                                                                   hessian,
                                                                   hdvec,
                                                                   dvec,
                                                                   activeAtomMask,
                                                                   debugInteractions );
  }
  if(this->_ChiralRestraint->isEnabled()) {
    totalEnergy += this->_ChiralRestraint->evaluateAllComponent( this->asSmartPtr(),
                                                                 pos,
                                                                 energyScale,
                                                                 componentEnergy,
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
  if(this->_AnchorRestraint->isEnabled()) {
    totalEnergy += this->_AnchorRestraint->evaluateAllComponent( this->asSmartPtr(),
                                                                 pos,
                                                                 energyScale,
                                                                 componentEnergy,
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
  if(this->_FixedNonbondRestraint->isEnabled()) {
    totalEnergy += this->_FixedNonbondRestraint->evaluateAllComponent( this->asSmartPtr(),
                                                                       pos,
                                                                       energyScale,
                                                                       componentEnergy,
                                                                       calcForce, force,
                                                                       calcDiagonalHessian,
                                                                       calcOffDiagonalHessian,
                                                                       hessian,
                                                                       hdvec,
                                                                       dvec,
                                                                       activeAtomMask,
                                                                       debugInteractions );
  }
  for ( auto cur : this->_OtherEnergyComponents ) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(oCdr(pair));
    if (component->isEnabled()) {
      totalEnergy+= component->evaluateAllComponent(this->asSmartPtr(),
                                                    pos,
                                                    energyScale,
                                                    componentEnergy,
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
    this->_Stretch->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_Angle->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_Dihedral->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_Nonbond->compareAnalyticalAndNumericalForceAndHessianTermByTerm(this->asSmartPtr(),pos);
    this->_DihedralRestraint->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_ChiralRestraint->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_AnchorRestraint->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_FixedNonbondRestraint->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    for ( auto cur : this->_OtherEnergyComponents ) {
      core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
      EnergyComponent_sp component = gc::As<EnergyComponent_sp>(oCdr(pair));
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

/*!
 * Look for bad geometry
 *
 * Extract the coordinates of the atoms
 */
CL_DEFMETHOD
core::List_sp EnergyFunction_O::checkForBeyondThresholdInteractions(double threshold)
{
  NVector_sp	pos;
  stringstream	info;
  int	fails = 0;

  info.str("");
  pos = NVector_O::create(this->getNVectorSize());
  this->loadCoordinatesIntoVector(pos);
  fails = 0;
  ql::list result;
  result & this->_Stretch->checkForBeyondThresholdInteractionsWithPosition(pos,threshold);
  result & this->_Angle->checkForBeyondThresholdInteractionsWithPosition(pos,threshold);
  return result.cons();
}



//
//
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------





string	EnergyFunction_O::energyTermsEnabled()
{
  stringstream ss;
  ss.str("");
  ss << this->_Stretch->enabledAsString();
  ss << this->_Angle->enabledAsString();
  ss << this->_Dihedral->enabledAsString();
  ss << this->_Nonbond->enabledAsString();
  ss << this->_DihedralRestraint->enabledAsString();
  ss << this->_ChiralRestraint->enabledAsString();
  ss << this->_AnchorRestraint->enabledAsString();
  ss << this->_FixedNonbondRestraint->enabledAsString();
  for ( auto cur : this->_OtherEnergyComponents ) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(oCdr(pair));
    ss << _rep_(oCar(pair)) << "-" << component->enabledAsString();
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
  core::clasp_write_string(fmt::format("Number of Stretch terms: {}\n" , this->_Stretch->numberOfTerms() ));
  core::clasp_write_string(fmt::format("Number of Angle terms: {}\n" , this->_Angle->numberOfTerms() ));
  core::clasp_write_string(fmt::format("Number of Dihedral terms: {}\n" , this->_Dihedral->numberOfTerms() ));
};



void	EnergyFunction_O::dumpTerms()
{
  core::HashTable_sp atomTypes = this->_AtomTypes;
  this->_AtomTable->dumpTerms(atomTypes);
  this->_Stretch->dumpTerms(atomTypes);
  this->_Angle->dumpTerms(atomTypes);
  this->_Dihedral->dumpTerms(atomTypes);
  this->_Nonbond->dumpTerms(atomTypes);
  this->_DihedralRestraint->dumpTerms(atomTypes);
  this->_ChiralRestraint->dumpTerms(atomTypes);
  this->_AnchorRestraint->dumpTerms(atomTypes);
  this->_FixedNonbondRestraint->dumpTerms(atomTypes);
  for ( auto cur : this->_OtherEnergyComponents ) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(oCdr(pair));
    component->dumpTerms(atomTypes);
  }
}

core::T_sp specializeKeepInteractionFactory( core::T_sp keepInteractionFactory, core::T_sp aclass ) {
  if (keepInteractionFactory.nilp()) return keepInteractionFactory;
  if (keepInteractionFactory == _lisp->_true()) return keepInteractionFactory;
  core::T_sp keepInteractionFunction = core::eval::funcall(keepInteractionFactory,aclass);
  return keepInteractionFunction;
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
      energyTerm.term.phi0 = dih->getMaxDegrees()*0.0174533;
      energyTerm.term.kdh = dih->getWeight();
      this->_DihedralRestraint->addTerm(energyTerm);
      ++terms;
    } else if ( restraint.isA<RestraintAnchor_O>() )
    {
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
      iterm.term.I1 = ea1->coordinateIndexTimes3();
      this->_AnchorRestraint->addTerm(iterm);
      ++terms;
    } else if ( restraint.isA<RestraintFixedNonbond_O>() )
    {
      this->_FixedNonbondRestraint->setupForEvaluation(this->_AtomTable,this->_NonbondCrossTermTable);
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
          this->_FixedNonbondRestraint->addFixedAtom(nonbondDb,a1,atomTypes);
          ++terms;
        }
      }
    } else if ( restraint.isA<RestraintDistance_O>() ) {
      RestraintDistance_sp rd = gc::As_unsafe<RestraintDistance_sp>(restraint);
      EnergyAtom* ea1 = this->getEnergyAtomPointer(rd->_A);
      EnergyAtom* ea2 = this->getEnergyAtomPointer(rd->_B);
      EnergyStretch   energyStretch(rd->_A,rd->_B,ea1->coordinateIndexTimes3(),ea2->coordinateIndexTimes3(),rd->_K,rd->_R0);
      this->_Stretch->addTerm(energyStretch);
    } else {
      SIMPLE_ERROR("Handle restraint: {}" , _rep_(restraint));
    }
  CONT:
    restraintIterator->next();
  }
  return terms;
}

void EnergyFunction_O::_addDihedralRestraintDegrees(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double phi0_degrees, double kdh, core::T_sp keepInteraction)
{
  if ( skipInteraction( keepInteraction, a1, a2, a3, a4 ) ) return;
  this->_DihedralRestraint->addDihedralRestraint(this->asSmartPtr(),
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
  double trans  = 180.0;
  double cis    = 0.0;
  double weight = 1.0;
  for ( ni=nitrogens.begin(); ni!=nitrogens.end(); ni++ )
  {
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
        this->_addDihedralRestraintDegrees(ax1,ax,ay,ay1,cis,weight,keepInteraction);
      LOG("Restrain cis {} - {} - {} -{}" , ax1->description() , ax->description() , ay->description() , ay1->description()  );
      if (!skipInteraction(keepInteraction, ax1, ax, ay, ay2 )) 
        this->_addDihedralRestraintDegrees(ax1,ax,ay,ay2,trans,weight,keepInteraction);
      LOG("Restrain trans {} - {} - {} -{}" , ax1->description() , ax->description() , ay->description() , ay2->description()  );
      if (!skipInteraction(keepInteraction, ax2, ax, ay, ay1 )) 
        this->_addDihedralRestraintDegrees(ax2,ax,ay,ay1,trans,weight,keepInteraction);
      LOG("Restrain trans {} - {} - {} -{}" , ax2->description() , ax->description() , ay->description() , ay1->description()  );
      if (!skipInteraction(keepInteraction, ax2, ax, ay, ay2 )) 
        this->_addDihedralRestraintDegrees(ax2,ax,ay,ay2,cis,weight,keepInteraction);
      LOG("Restrain cis {} - {} - {} -{}" , ax2->description() , ax->description() , ay->description() , ay2->description()  );
    }
  }
}





SYMBOL_EXPORT_SC_(ChemPkg,find_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,force_fields_as_list);

SYMBOL_EXPORT_SC_(ChemPkg,add_shadowing_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,assign_force_field_types);
SYMBOL_EXPORT_SC_(ChemPkg,compute_merged_nonbond_force_field_for_aggregate);
SYMBOL_EXPORT_SC_(ChemPkg,force_field_component_merge);
SYMBOL_EXPORT_SC_(ChemPkg,clear_combined_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,generate_molecule_energy_function_tables);
SYMBOL_EXPORT_SC_(ChemPkg,STARparameter_warningsSTAR);
SYMBOL_EXPORT_SC_(ChemPkg,report_parameter_warnings);

SYMBOL_EXPORT_SC_(ChemPkg,identify_aromatic_rings);
SYMBOL_EXPORT_SC_(ChemPkg,STARcurrent_aromaticity_informationSTAR);

CL_LISPIFY_NAME("defineForMatter");
CL_LAMBDA((energy-function chem:energy-function) matter &key use-excluded-atoms (keep-interaction-factory t) (assign-types t));
CL_DEFMETHOD void EnergyFunction_O::defineForMatter(Matter_sp matter, bool useExcludedAtoms, core::T_sp keepInteractionFactory, bool assign_types )
{
  if ( !(matter.isA<Aggregate_O>() || matter.isA<Molecule_O>() ) )
  {
    SIMPLE_ERROR("You can only define energy functions for Aggregates or Molecules");
  }

  core::DynamicScopeManager scope(_sym_STARparameter_warningsSTAR,nil<core::T_O>());
	//
	// Identify rings
	//
  if (chem__verbose(0)) core::clasp_write_string("Searching for rings.\n");

  // If *current-rings* is bound then reuse it, otherwise, calculate rings and bind those
  core::T_sp rings = unbound<core::T_sp>();
  if (_sym_STARcurrent_ringsSTAR->boundP()) {
    rings = _sym_STARcurrent_ringsSTAR->symbolValue();
  } else {
    rings = RingFinder_O::identifyRings(matter);
  }
  core::DynamicScopeManager ring_scope(_sym_STARcurrent_ringsSTAR, rings );

  //
  // Assign relative Cahn-Ingold-Preylog priorities
  //
  if (chem__verbose(0)) core::clasp_write_string("Assigning CIP priorities.\n");
  core::HashTable_sp cip = CipPrioritizer_O::assignPrioritiesHashTable(matter);

    //
    // Assign atom types for each molecule
    //

  core::HashTable_sp force_fields = core::HashTable_O::createEq();
  Loop moleculeLoop;
  moleculeLoop.loopTopGoal(matter,MOLECULES);
  while (moleculeLoop.advanceLoopAndProcess() ) {
    Molecule_sp molecule = moleculeLoop.getMolecule();
    core::T_sp force_field_name = molecule->force_field_name();
    if (force_fields->gethash(force_field_name).nilp()) {
      core::T_sp combined_force_field = core::eval::funcall(_sym_find_force_field,force_field_name);
      force_fields->setf_gethash(force_field_name,combined_force_field);
    }
  }
  core::HashTable_sp atomTypes = core::HashTable_O::createEq();
  // Assign given atom-types
  Loop atom_loop;
  atom_loop.loopTopGoal(matter,ATOMS);
  while (atom_loop.advanceLoopAndProcess()) {
    Atom_sp atom = atom_loop.getAtom();
    core::T_sp type = atom->getPropertyOrDefault( kw::_sym_given_atom_type, nil<core::T_O>() );
    if (type.notnilp()) {
      atomTypes->setf_gethash(atom,type);
    }
  }
  this->_AtomTypes = atomTypes;
  if (assign_types) {
    if (chem__verbose(0)) core::clasp_write_string("Assigning atom types.\n");
    moleculeLoop.loopTopGoal(matter,MOLECULES);
    while (moleculeLoop.advanceLoopAndProcess() ) {
      Molecule_sp molecule = moleculeLoop.getMolecule();
      core::T_sp force_field_name = molecule->force_field_name();
      [[maybe_unused]]core::T_sp use_given_types = molecule->force_field_use_given_types();
//      if (use_given_types.nilp()) {
        core::T_sp combined_force_field = force_fields->gethash(force_field_name);
        if (chem__verbose(0)) core::clasp_write_string(fmt::format("Assigning atom types for molecule {} using {}.\n" , _rep_(molecule->getName()) , _rep_(force_field_name)));
  //
  // Calculate aromaticity using the rings we just calculated
  //
        core::T_sp aromaticity_info = core::eval::funcall(_sym_identify_aromatic_rings,matter,force_field_name);
        if (aromaticity_info.nilp()) SIMPLE_ERROR("The aromaticity-info was NIL when about to assign force field types - it should not be");
        core::DynamicScopeManager aromaticity_scope(_sym_STARcurrent_aromaticity_informationSTAR,aromaticity_info);
        core::eval::funcall(_sym_assign_force_field_types,combined_force_field,molecule,atomTypes);
#if 0
      } else {
        if (chem__verbose(0)) core::clasp_write_string(fmt::format("Assigning atom types for molecule {} using given-types for {}.\n" , _rep_(molecule->getName()) , _rep_(force_field_name)));
        Loop atom_loop;
        atom_loop.loopTopGoal(molecule,ATOMS);
        while (atom_loop.advanceLoopAndProcess()) {
          Atom_sp atom = atom_loop.getAtom();
          atomTypes->setf_gethash(atom,atom->atomType());
        }
      }
#endif
    }
  }
  this->defineForMatterWithAtomTypes(matter,useExcludedAtoms,keepInteractionFactory,cip,atomTypes);
}


CL_LAMBDA((energy-function chem:energy-function) matter &key use-excluded-atoms (keep-interaction-factory t) cip-priorities atom-types);
CL_DEFMETHOD void EnergyFunction_O::defineForMatterWithAtomTypes(Matter_sp matter, bool useExcludedAtoms, core::T_sp keepInteractionFactory, core::T_sp cip_priorities, core::HashTable_sp atomTypes )
{
  if (!gc::IsA<core::HashTable_sp>(cip_priorities)) {
    SIMPLE_ERROR("You need to provide a hash-table of atoms to relative CIP priorities - see CipPrioritizer_O::assignPrioritiesHashTable(matter)");
  }
  if (chem__verbose(0)) core::clasp_write_string("defineForMatterWithAtomTypes\n");
  this->_Matter= matter;
  if ( !(matter.isA<Aggregate_O>() || matter.isA<Molecule_O>() ) )
  {
    SIMPLE_ERROR("You can only define energy functions for Aggregates or Molecules");
  }

  //
  // Setup the extra information needed by the energy function to generate a topology file
  //
  // Get the name to pass that to the atomTable
  core::T_sp matterName = matter->getName();
  this->_AtomTable->setAggregateName(matterName);
  core::T_sp boundingBox = unbound<core::T_O>();
  if (gc::IsA<Aggregate_sp>(matter)) {
    Aggregate_sp agg = gc::As_unsafe<Aggregate_sp>(matter);
    if (agg->boundingBoxBoundP()) {
      boundingBox = agg->boundingBox();
    }
  }
  if (boundingBox.unboundp()) {
    this->_AtomTable->makUnboundBoundingBox();
  } else {
    this->_AtomTable->setBoundingBox(boundingBox);
  }
  if (chem__verbose(0)) core::clasp_write_string("Assembling aggregate nonbond force-field.\n");
  core::T_sp nonbondForceField = core::eval::funcall(chem::_sym_compute_merged_nonbond_force_field_for_aggregate,matter,atomTypes);
  this->_AtomTable->setNonbondForceFieldForAggregate(nonbondForceField);
  
  //

  // Separate the molecules for solute from the solvent and handle them solute first then solvent
  size_t final_solute_residue_iptres = 0;
  size_t number_of_molecules_nspm = 0;
  size_t first_solvent_molecule_nspsol = 0;
  bool solvent_exists = false;
  if (gc::IsA<Aggregate_sp>(matter)) {
    if (chem__verbose(0)) core::clasp_write_string("Classifying solute/solvent.\n");
    ql::list solute_molecules;
    ql::list solvent_molecules;
    {
      Loop moleculeLoop;
      moleculeLoop.loopTopGoal(matter,MOLECULES);
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
    for ( auto cur_solute : solute ) {
      Molecule_sp onemol = gc::As_unsafe<Molecule_sp>(CONS_CAR(cur_solute));
      core::T_sp force_field_name = onemol->force_field_name();
      core::T_sp forceField = core::eval::funcall(chem::_sym_find_force_field,force_field_name);
  //
  // Calculate aromaticity using the rings we just calculated
  //
      core::T_sp aromaticity_info = core::eval::funcall(_sym_identify_aromatic_rings,matter,force_field_name);
      if (aromaticity_info.nilp()) SIMPLE_ERROR("The aromaticity-info was NIL when about to call generate-molecule-energy-function-tables - it should not be");
      core::DynamicScopeManager aromaticity_scope(_sym_STARcurrent_aromaticity_informationSTAR,aromaticity_info);
      this->_AtomTable->constructFromMolecule(onemol,nonbondForceField,keepInteractionFactory,atomTypes);
      if (chem__verbose(0)) core::clasp_write_string(fmt::format("Generating parameters for {} using {} force-field.\n" , _rep_(onemol->getName()) , _rep_(force_field_name) ));
      core::eval::funcall(_sym_generate_molecule_energy_function_tables,this->asSmartPtr(),onemol,forceField,keepInteractionFactory);
      final_solute_residue_iptres += onemol->contentSize();
      ++number_of_molecules_nspm;
    }
    core::List_sp solvent = solvent_molecules.cons();
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
      core::T_sp force_field_name = onemol->force_field_name();
      core::T_sp forceField = core::eval::funcall(chem::_sym_find_force_field,force_field_name);
      this->_AtomTable->constructFromMolecule(onemol,nonbondForceField,keepInteractionFactory,atomTypes);
      core::eval::funcall(_sym_generate_molecule_energy_function_tables,this->asSmartPtr(),onemol,forceField,keepInteractionFactory);
      ++number_of_molecules_nspm;
    }
  } else {
    Molecule_sp molecule = gc::As<Molecule_sp>(matter);
    core::T_sp force_field_name = molecule->force_field_name();
    core::T_sp forceField = core::eval::funcall(chem::_sym_find_force_field,force_field_name);
  //
  // Calculate aromaticity using the rings we just calculated
  //
      core::T_sp aromaticity_info = core::eval::funcall(_sym_identify_aromatic_rings,matter,force_field_name);
      if (aromaticity_info.nilp()) SIMPLE_ERROR("The aromaticity-info was NIL when we were about to call generate-molecule-energy-function-tables for a single molecule - it should not be");
      core::DynamicScopeManager aromaticity_scope(_sym_STARcurrent_aromaticity_informationSTAR,aromaticity_info);
      this->_AtomTable->constructFromMolecule(molecule,nonbondForceField,keepInteractionFactory,atomTypes);
    if (chem__verbose(0)) core::clasp_write_string(fmt::format("Generating parameters for {} using {} force-field.\n" , _rep_(molecule->getName()) , _rep_(force_field_name) ));
    core::eval::funcall(_sym_generate_molecule_energy_function_tables,this->asSmartPtr(),molecule,forceField,keepInteractionFactory);
    final_solute_residue_iptres = molecule->contentSize();
    number_of_molecules_nspm = 1;
    first_solvent_molecule_nspsol = 2;
  }
  if (solvent_exists) {
    printf("%s:%d  solvent_exists NSPSOL %lu  IPTRES %lu   NSM %lu \n", __FILE__, __LINE__,
           first_solvent_molecule_nspsol, final_solute_residue_iptres, number_of_molecules_nspm);
    this->_AtomTable->set_firstSolventMoleculeNSPSOL(first_solvent_molecule_nspsol);
    this->_AtomTable->set_finalSoluteResidueIPTRES(final_solute_residue_iptres);
    this->_AtomTable->set_totalNumberOfMoleculesNSPM(number_of_molecules_nspm);
  }
  {
    if (chem__verbose(1)) core::clasp_write_string("About to calculate nonbond and restraint terms");
    core::T_sp nonbondForceField = this->_AtomTable->nonbondForceFieldForAggregate();
    this->generateNonbondEnergyFunctionTables(useExcludedAtoms,matter,nonbondForceField,keepInteractionFactory,atomTypes);
    this->generateRestraintEnergyFunctionTables(matter,nonbondForceField,keepInteractionFactory,cip_priorities,atomTypes);
  }
  core::eval::funcall(_sym_report_parameter_warnings);
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

CL_LAMBDA((energy-function chem:energy-function) molecule ffstretches ffangles ffptors ffitors &key (keep-interaction-factory t) atom-types)
CL_DOCSTRING(R"dx(Generate the standard energy function tables. The atom types, and CIP priorities need to be precalculated.)dx")
CL_DEFMETHOD void EnergyFunction_O::generateStandardEnergyFunctionTables(Matter_sp molecule,
                                                                         FFStretchDb_sp ffstretches,
                                                                         FFAngleDb_sp ffangles,
                                                                         FFPtorDb_sp ffptors,
                                                                         FFItorDb_sp ffitors,
                                                                         core::T_sp keepInteractionFactory,
                                                                         core::HashTable_sp atomTypes )
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

  //
  // Define a Nonbond cross term table
  //
  // FIXME: I think this should only be done once after all of the molecules have been added
  // FIXME: Why are we using the forcefield->getNonbondDb()?  We should be using the aggregate nonbond force-field
#if 0  
  if (chem__verbose(0)) core::clasp_write_string("Starting to build standard energy function tables\n");
  auto  temp  = gctools::GC<FFNonbondCrossTermTable_O>::allocate_with_default_constructor();
  this->_NonbondCrossTermTable = temp;
  this->_NonbondCrossTermTable->fillUsingFFNonbondDb(forceField->getNonbondDb());
#endif

  core::HashTable_sp atomToRes = createAtomToResidueHashTable(molecule);
    	//
	// Initialize the energy components
	//
  ALL_ENERGY_COMPONENTS(initialize());
  this->_eraseMissingParameters();
  coordinateIndex = 0;
  ASSERTNOTNULL(forceField);
  // Search the stretch terms
  {
    size_t terms = 0;
    size_t missing_terms = 0;
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyStretch_O::staticClass() );
    loop.loopTopGoal(molecule,BONDS);
    while ( loop.advanceLoopAndProcess() ) {
      a1 = loop.getBondA1();
      a2 = loop.getBondA2();
      if ( skipInteraction( keepInteraction, a1, a2 ) ) continue;
//      printf("%s:%d Looking at STRETCH term between %s - %s\n", __FILE__, __LINE__, _rep_(a1).c_str(), _rep_(a2).c_str());
      t1 = a1->getType(atomTypes);
      t2 = a2->getType(atomTypes);
      ea1 = this->getEnergyAtomPointer(a1);
      ea2 = this->getEnergyAtomPointer(a2);
      FFStretch_sp ffStretch = gc::As<FFStretch_sp>(ffstretches->findTermForTypes(t1,t2));
      if ( ffStretch->level() != parameterized ) {
        this->_addMissingParameter(ffStretch);
        ++missing_terms;
      }
      if ( ffStretch->level() != unknown ) {
        EnergyStretch   energyStretch;
        energyStretch.defineFrom(ffStretch,ea1,ea2,this->_Stretch->getScale());
        this->_Stretch->addTerm(energyStretch);
        ++terms;
      } else {
        Residue_sp res1 = gc::As<Residue_sp>(atomToRes->gethash(a1));
        Residue_sp res2 = gc::As<Residue_sp>(atomToRes->gethash(a2));
        SIMPLE_WARN("Could not find stretch parameter in molecule {} between {}/{} (atom type {}) and {}/{} (atom type {})" , _rep_(molecule) , _rep_(res1) , _rep_(a1) , _rep_(a1->getType(atomTypes)) , _rep_(res2) , _rep_(a2) , _rep_(a2->getType(atomTypes)));
      }
    }
    if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built stretch table with {} terms added and {} missing terms\n" , terms , missing_terms));
  }
#ifdef	DEBUG_DEFINE_ENERGY
  core::clasp_write_string(fmt::format("{}:{} There were {} stretch terms\n" , __FILE__ , __LINE__ , this->_Stretch.size() ));
#endif
  // Search the angle terms
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyAngle_O::staticClass() );
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
        energyAngle.defineFrom(ffAngle,ea1,ea2,ea3,this->_Angle->getScale());
        this->_Angle->addTerm(energyAngle);
        ++terms;
      }
    }
    if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built angle table with {} terms and {} missing terms\n" , terms , missing_terms));
  }
  // Search the ptor terms
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyDihedral_O::staticClass() );
    size_t terms = 0;
    size_t missing_terms = 0;
    loop.loopTopGoal(molecule,PROPERS);
    while ( loop.advanceLoopAndProcess() ) {
      a1 = loop.getAtom1();
      a2 = loop.getAtom2();
      a3 = loop.getAtom3();
      a4 = loop.getAtom4();
      if ( skipInteraction( keepInteraction, a1, a2, a3, a4 ) ) continue;
      t1 = a1->getType(atomTypes);
      t2 = a2->getType(atomTypes);
      t3 = a3->getType(atomTypes);
      t4 = a4->getType(atomTypes);
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
            energyDihedral.defineFrom(n,ffPtor,ea1,ea2,ea3,ea4,this->_Dihedral->getScale());
            this->_Dihedral->addTerm(energyDihedral);
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
  }
  // Search the itor terms
  {
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyDihedral_O::staticClass() );
    EnergyDihedral energyDihedral;
    size_t terms = 0;
    loop.loopTopGoal(molecule,IMPROPERS);
    while ( loop.advanceLoopAndProcess() ) {
      a1 = loop.getAtom1();
      a2 = loop.getAtom2();
      a3 = loop.getAtom3();
      a4 = loop.getAtom4();
      FFItorDb_O::improperAtomSort(a1,a2,a4,atomTypes);
      if ( skipInteraction( keepInteraction, a1, a2, a3, a4 ) ) continue;
      t1 = a1->getType(atomTypes);
      t2 = a2->getType(atomTypes);
      t3 = a3->getType(atomTypes);
      t4 = a4->getType(atomTypes);
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
              energyDihedral.defineFrom(n,ffItor,ea1,ea2,ea3,ea4,this->_Dihedral->getScale());
              this->_Dihedral->addTerm(energyDihedral);
              aImproperCenter = a3;
              ++terms;
            }
          }
        }
      }
    }
    if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built improper table for {} terms\n" , terms));
  }
//  this->summarizeTerms();
}

SYMBOL_EXPORT_SC_(ChemPkg,prepare_amber_energy_nonbond);

CL_DOCSTRING(R"dx(Generate the nonbond energy function tables. The atom types, and CIP priorities need to be precalculated.)dx");
CL_DEFMETHOD void EnergyFunction_O::generateNonbondEnergyFunctionTables(bool useExcludedAtoms, Matter_sp matter, core::T_sp nonbondForceField, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes )
{
  if (keepInteractionFactory.nilp()) return;
  if (chem__verbose(0))
    core::clasp_write_string(fmt::format("Built atom table for {} atoms\n" , this->_AtomTable->getNumberOfAtoms()));
#ifdef	DEBUG_DEFINE_ENERGY
  core::clasp_write_string(fmt::format("{}:{} There were {} atoms\n" , __FILE__ , __LINE__ , this->_AtomTable.size() ));
#endif
        // Nonbonds here!!!!!!!!!!!!!!
  if (useExcludedAtoms) {
// The nonbond parameters are calculated in Common Lisp
//    printf("atom number: %d\n", this->_AtomTable->getNumberOfAtoms());

//    printf("%s:%d in generateNonbondEnergyFunctionTables\n", __FILE__, __LINE__);
    core::List_sp parts = core::eval::funcall(_sym_prepare_amber_energy_nonbond,this->asSmartPtr(),nonbondForceField);
//      printf("%s:%d:%s    parts -> %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(parts).c_str());
//    printf("%s:%d in generateNonbondEnergyFunctionTables\n", __FILE__, __LINE__);
    this->_Nonbond->constructNonbondTermsFromAList(parts);
    // -----------
//    printf("%s:%d in generateNonbondEnergyFunctionTables\n", __FILE__, __LINE__);
      
    this->_Nonbond->constructExcludedAtomListFromAtomTable(this->_AtomTable, nonbondForceField,keepInteractionFactory);
    this->_Nonbond->construct14InteractionTerms(this->_AtomTable,matter,nonbondForceField,keepInteractionFactory,atomTypes);
    LOG("Done construct14InteractionTerms");
//    printf("%s:%d:%s    nonbond -> %d\n", __FILE__, __LINE__, __FUNCTION__, _Nonbond->numberOfTerms());

  } else {
    this->_Nonbond->constructNonbondTermsFromAtomTable(false,this->_AtomTable, nonbondForceField,atomTypes, keepInteractionFactory );
  }
  if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built nonbond table for {} terms\n" , this->_Nonbond->numberOfTerms()));
}


CL_LAMBDA((energy-function chem:energy-function) matter force-field &key (keep-interaction-factory t) cip-priorities atom-types);
CL_DOCSTRING(R"dx(Generate the restraint energy function tables. The atom types, and CIP priorities need to be precalculated.
This should be called after generateStandardEnergyFunctionTables.
You need to pass a hash-table of atoms to relative CIP priorities (calculated using CipPrioritizer_O::assignPrioritiesHashTable(matter) for stereochemical restraints.)dx")
CL_DEFMETHOD void EnergyFunction_O::generateRestraintEnergyFunctionTables(Matter_sp matter, core::T_sp ffNonbond, core::T_sp keepInteractionFactory, core::T_sp cip_priorities, core::HashTable_sp atomTypes )
{
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
    if (!gc::IsA<core::HashTable_sp>(cip_priorities)) {
      SIMPLE_ERROR("You need to provide a hash-table of atoms to relative CIP priorities - see CipPrioritizer_O::assignPrioritiesHashTable(matter)");
    }
    core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyChiralRestraint_O::staticClass() );
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
            ichiral.term.I1 = ea1->coordinateIndexTimes3();
            ichiral.term.I2 = ea2->coordinateIndexTimes3();
            ichiral.term.I3 = eaCenter->coordinateIndexTimes3();
            ichiral.term.I4 = ea3->coordinateIndexTimes3();
            ichiral.term.K = DefaultChiralRestraintWeight * side;
            ichiral.term.CO = DefaultChiralRestraintOffset;
            this->_ChiralRestraint->addTerm(ichiral);
				// Now apply it to the other atom
				// on the chiral center, just flip the sign
				// of K
          }
          
          ichiral._Atom4 = ea4->atom();
          if ( !skipInteraction( keepInteraction, ichiral._Atom1, ichiral._Atom2, ichiral._Atom3, ichiral._Atom4 ) ) {
            ichiral.term.I4 = ea4->coordinateIndexTimes3();
					// flip the sign of the chiral restraint
            ichiral.term.K = DefaultChiralRestraintWeight * side * -1.0;
            ichiral.term.CO = DefaultChiralRestraintOffset;
            this->_ChiralRestraint->addTerm(ichiral);
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
            ichiral.term.I1 = ea2->coordinateIndexTimes3();
            ichiral.term.I2 = ea4->coordinateIndexTimes3();
            ichiral.term.I3 = eaCenter->coordinateIndexTimes3();
            ichiral.term.I4 = ea3->coordinateIndexTimes3();
            ichiral.term.K = DefaultChiralRestraintWeight * side;
            ichiral.term.CO = DefaultChiralRestraintOffset;
            this->_ChiralRestraint->addTerm(ichiral);
          }
				// Now apply it to the other atom
				// on the chiral center, just flip the sign
				// of K
          ichiral._Atom4 = ea1->atom();
          if ( !skipInteraction( keepInteraction, ichiral._Atom1, ichiral._Atom2, ichiral._Atom3, ichiral._Atom4 ) ) {
            ichiral.term.I4 = ea1->coordinateIndexTimes3();
					// flip the sign of the chiral restraint
            ichiral.term.K = DefaultChiralRestraintWeight * side * -1.0;
            ichiral.term.CO = DefaultChiralRestraintOffset;
            this->_ChiralRestraint->addTerm(ichiral);
          }
        } else {
          LOG("There is no chiral restraint for: {}" , a1->description()  );
        }
      }
    }
    if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built chiral restraints table for {} terms\n" , this->_ChiralRestraint->numberOfTerms()));
  }

	//
	// Set up the anchor restraints
	//
	//
	// Define secondary amide restraints
	//
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
    int startTerms = this->_DihedralRestraint->numberOfTerms();
    // this->__createSecondaryAmideRestraints(nitrogens,keepInteraction);
    if (chem__verbose(0)) core::clasp_write_string(fmt::format("Built secondary amide restraints including {} terms\n" , (this->_DihedralRestraint->numberOfTerms() - startTerms)));
  } else
  {
    LOG("Skipping Secondary amide restraints because _RestrainSecondaryAmides = {}" , this->_RestrainSecondaryAmides );
  }

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



CL_LISPIFY_NAME("summarizeBeyondThresholdInteractionsAsString");
CL_DEFMETHOD string EnergyFunction_O::summarizeBeyondThresholdInteractionsAsString()
{
  stringstream	ss;
  ss.str("");
  ss << "Beyond threshold ";
  ss << this->_Stretch->beyondThresholdInteractionsAsString();
  ss << this->_Angle->beyondThresholdInteractionsAsString();
  ss << this->_Dihedral->beyondThresholdInteractionsAsString();
  ss << this->_Nonbond->beyondThresholdInteractionsAsString();
  ss << this->_DihedralRestraint->beyondThresholdInteractionsAsString();
  ss << this->_ChiralRestraint->beyondThresholdInteractionsAsString();
  ss << this->_AnchorRestraint->beyondThresholdInteractionsAsString();
  ss << this->_FixedNonbondRestraint->beyondThresholdInteractionsAsString();
#if 0
  for ( auto cur : this->_OtherEnergyComponents ) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(oCdr(pair));
    ss << component->beyondThreshholdInteractionsAsString();
  }
#endif
  return ss.str();
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
  ss << this->_Stretch->debugLogAsString() << std::endl;
  ss << this->_Angle->debugLogAsString() << std::endl;
  ss << this->_Dihedral->debugLogAsString() << std::endl;
  ss << this->_Nonbond->debugLogAsString() << std::endl;
  ss << this->_DihedralRestraint->debugLogAsString() << std::endl;
  ss << this->_ChiralRestraint->debugLogAsString() << std::endl;
  ss << this->_AnchorRestraint->debugLogAsString() << std::endl;
  ss << this->_FixedNonbondRestraint->debugLogAsString() << std::endl;
  for ( auto cur : this->_OtherEnergyComponents ) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    EnergyComponent_sp component = gc::As<EnergyComponent_sp>(oCdr(pair));
    ss << component->debugLogAsString();
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
  energy->_Stretch = (safe_alist_lookup<EnergyStretch_sp>(alist,kw::_sym_stretch));
  energy->_Angle = (safe_alist_lookup<EnergyAngle_sp>(alist,kw::_sym_angle));
  energy->_Dihedral = (safe_alist_lookup<EnergyDihedral_sp>(alist,kw::_sym_dihedral));
  energy->_Nonbond = (safe_alist_lookup<EnergyNonbond_sp>(alist,kw::_sym_nonbond));
}


void EnergyFunction_O::enableDebug()
{
  ALL_ENERGY_COMPONENTS(enableDebug());
}

void EnergyFunction_O::disableDebug()
{
  ALL_ENERGY_COMPONENTS(disableDebug());
}


CL_DEFMETHOD
EnergyFunction_sp EnergyFunction_O::copyFilter(core::T_sp keepInteractionFactory)
{
  auto  me  = gctools::GC<EnergyFunction_O>::allocate_with_default_constructor();
  me->_Matter = this->_Matter;
  me->_AtomTable = this->_AtomTable;
  me->_NonbondCrossTermTable = this->_NonbondCrossTermTable;
  me->_BoundingBox = this->_BoundingBox;
  if (keepInteractionFactory.notnilp()) {
    me->_Stretch = this->_Stretch->copyFilter(keepInteractionFactory);
    me->_Angle = this->_Angle->copyFilter(keepInteractionFactory);
    me->_Dihedral = this->_Dihedral->copyFilter(keepInteractionFactory);
    me->_Nonbond = this->_Nonbond->copyFilter(keepInteractionFactory);
    me->_DihedralRestraint = this->_DihedralRestraint->copyFilter(keepInteractionFactory);
    me->_ChiralRestraint = this->_ChiralRestraint->copyFilter(keepInteractionFactory);
    me->_AnchorRestraint = this->_AnchorRestraint->copyFilter(keepInteractionFactory);
    me->_FixedNonbondRestraint = this->_FixedNonbondRestraint->copyFilter(keepInteractionFactory);
    ql::list ll;
    for ( auto cur : this->_OtherEnergyComponents ) {
      core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
      core::T_sp name = oCar(pair);
      EnergyComponent_sp component = gc::As<EnergyComponent_sp>(oCdr(pair));
      ll << core::Cons_O::create(name, component->filterCopyComponent(keepInteractionFactory) );
    }
    me->_OtherEnergyComponents = ll.cons();
  }
  return me;
}



};
