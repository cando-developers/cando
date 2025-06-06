
/*
    File: sketchFunction.cc
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
__BEGIN_DOC(classes.SketchFunction,section,SketchFunction)
SketchFunction objects are used to calculate Amber or Gaff molecular 
mechanics energies and first and second derivatives.
__END_DOC
*/

/*
 *	sketchFunction.cc
 *
 */

#include <clasp/core/common.h>
#include <clasp/core/bformat.h>
#include <cando/chem/sketchFunction.h>
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
#include <cando/chem/energySketchStretch.h>
#include <cando/chem/energyPointToLineRestraint.h>
#include <cando/chem/energyOutOfZPlane.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/chemInfo.h>
#include <clasp/core/wrappers.h>


#if 0
#define	TURN_DEBUG_OFF
#define	MATHEMATICA_DEBUG_ON	1
#define	PARMDEBUG	1
#endif


#ifdef	darwin
#include <xmmintrin.h>
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
SYMBOL_EXPORT_SC_(ChemPkg,LinearImproperRestraintError);
SYMBOL_EXPORT_SC_(ChemPkg,OverlappingNonbondError);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);


SYMBOL_EXPORT_SC_(ChemPkg,make_node_table_from_graph);
SYMBOL_EXPORT_SC_(ChemPkg,node_table_size);
SYMBOL_EXPORT_SC_(ChemPkg,node_table_node_at_index);
SYMBOL_EXPORT_SC_(ChemPkg,node_table_coordinate_index_times3);
SYMBOL_EXPORT_SC_(ChemPkg,node_get_position);
SYMBOL_EXPORT_SC_(ChemPkg,node_set_position);

bool sketchFunctionInitialized = false;


CL_DOCSTRING(R"dx(Create a sketch-function from a molecule. Pass an instance of a dummy sketch-nonbond-force-field passed to chem:find-atom-type-position.)dx");
CL_LISPIFY_NAME(make-sketch-function);
CL_DEF_CLASS_METHOD
SketchFunction_sp SketchFunction_O::make(core::T_sp graph, core::T_sp sketchNonbondForceField,core::HashTable_sp atomTypes)
{
  auto sf = gctools::GC<SketchFunction_O>::allocate(graph);
  sf->_NodeTable = core::eval::funcall(chem::_sym_make_node_table_from_graph,graph,sketchNonbondForceField,atomTypes);
  return sf;
}

void	SketchFunction_O::initialize()
{
  this->Base::initialize();
  this->_Nonbond = EnergySketchNonbond_O::create();
  this->_Stretch = EnergySketchStretch_O::create();
  this->_PointToLineRestraint = EnergyPointToLineRestraint_O::create(this->_Stretch);
  this->_OutOfZPlane = EnergyOutOfZPlane_O::create();
  this->useDefaultSettings();
}

void	SketchFunction_O::useDefaultSettings()
{
  this->_Nonbond->initialize();
  this->_Stretch->initialize();
  this->_PointToLineRestraint->initialize();
  this->_OutOfZPlane->initialize();
  this->_PointToLineRestraintWeight = DefaultPointToLineRestraintWeight;
  this->_OutOfZPlaneWeight = DefaultOutOfZPlaneWeight;
  this->_VelocityScale.set(1.0, 1.0, 1.0);
}



void SketchFunction_O::fields(core::Record_sp node)
{
  node->field_if_not_unbound(INTERN_(kw,NodeTable),this->_NodeTable);
  node->field_if_not_unbound(INTERN_(kw,Nonbond),this->_Nonbond);
  node->field_if_not_unbound(INTERN_(kw,Stretch),this->_Stretch);
  node->field_if_not_unbound(INTERN_(kw,PointToLineRestraint),this->_PointToLineRestraint);
  node->field_if_not_unbound(INTERN_(kw,OutOfZPlane),this->_OutOfZPlane);
  this->Base::fields(node);
}

CL_DEFMETHOD AtomTable_sp SketchFunction_O::atomTable() const
{
  if (gc::IsA<AtomTable_sp>(this->_NodeTable)) {
    return gc::As_unsafe<AtomTable_sp>(this->_NodeTable);
  }
  SIMPLE_ERROR("An attempt was made to get the atom-table of a sketch-function but one hasn't been set - instead we have a node-table: {}" , _rep_(this->_NodeTable));
}

CL_DOCSTRING(R"dx(Return the instance graph slot if it is a Matter object. 
If it isn't then signal an error. 
This is used in situations where the graph slot contains a matter object.)dx")
CL_DEFMETHOD
Matter_sp SketchFunction_O::getMatter() {
  return gc::As<Matter_sp>(this->_Graph);
}

size_t SketchFunction_O::getNVectorSize() const
{
  core::T_sp tsize = core::eval::funcall(_sym_node_table_size,this->_NodeTable);
  if (tsize.fixnump()) {
    return tsize.unsafe_fixnum()*3;
  }
  SIMPLE_ERROR("chem:get-nvector-size must return a fixnum");
};


void	SketchFunction_O::setOption( core::Symbol_sp option, core::T_sp val)
{
  SIMPLE_ERROR("Unknown SketchFunction setOption keyword[{}]" , _rep_(option) );
}



/*
__BEGIN_DOC(classes.SketchFunction.setOptions,subsection,setOptions)
\scriptCmd{setOptions}{sketchFunction options}

Set the \scriptArg{options} for the energy function. The options are represented by a list of list pairs
of the form [: 'symbol value ]
  They include...\par
[ 'nonbond true/false ]\par
__END_DOC
*/
void	SketchFunction_O::setOptions( core::List_sp options )
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

void SketchFunction_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                  AbstractLargeSquareMatrix_sp m,
                                                  core::T_sp activeAtomMask )
{
  m->fill(0.0);
  if (this->_Stretch->isEnabled())
    this->_Stretch->setupHessianPreconditioner(nvPosition, m, activeAtomMask );
  if (this->_OutOfZPlane->isEnabled())
    this->_OutOfZPlane->setupHessianPreconditioner(nvPosition, m, activeAtomMask );
}







uint	SketchFunction_O::countTermsBeyondThreshold()
{
  int		terms;
  terms = 0;
  SIMPLE_ERROR("Should there be something here?");
//    node = rawAccumulateTermsBeyondThresholdAsXml(terms);
  return terms;
}



void maybe_dump_force(const string& msg, NVector_sp force)
{
  if (chem__verbose(2)) {
    core::clasp_write_string(fmt::format("Name: {}\n" , msg));
    for ( size_t idx=0; idx<force->length(); idx++ ) {
      core::clasp_write_string(fmt::format("force[{}] -> {}\n" , idx , (*force)[idx]));
    }
  }
}



//

double	SketchFunction_O::evaluateAll( NVector_sp 	pos,
                                       core::T_sp energyScale,
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
  double totalEnergy = 0.0;
#ifdef DEBUG_ENERGY_FUNCTION
    printf("%s:%d:%s Entered\n", __FILE__, __LINE__, __FUNCTION__ );
#endif
#if 0 //[
#define	ForceAcc(i,o,v) {\
    if ( hasForce ) {\
      force->setElement((i)+(o),(v)+force->getElement((i)+(o)));\
    }\
  }

//
// Accumulate an off diagonal Hessian element
//
#define	OffDiagHessAcc(i1,o1,i2,o2,v) {\
    if ( hasHessian ) {\
      hessian->addToElement((i1)+(o1),(i2)+(o2),v);\
    }\
    if ( hasHdAndD ) {\
      hdvec->addToElement((i1)+(o1),v*dvec->element((i2)+(o2)));\
      hdvec->addToElement((i2)+(o2),v*dvec->element((i1)+(o1)));\
    }\
  }

//
// Accumulate a diagonal Hessian element
//
#define	DiagHessAcc(i1,o1,i2,o2,v) {\
    if ( hasHessian ) {\
      hessian->addToElement((i1)+(o1),(i2)+(o2),v);\
    }\
    if ( hasHdAndD ) {\
      hdvec->addToElement((i1)+(o1),v*dvec->element((i1)+(o1)));\
    }\
  }
#endif //]

  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool   hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());

#ifdef	DEBUG_ON //[
	// Summarize entry state for debugging
  LOG("calcForce = {}" , calcForce  );
  LOG("calcDiagonalHessian = {}" , calcDiagonalHessian  );
  LOG("calcOffDiagonalHessian = {}" , calcOffDiagonalHessian  );
  LOG("hasForce = {}" , hasForce  );
  LOG("hasHdAndD = {}" , hasHdAndD  );
  if ( hasForce && force->size() < pos->size() ) {
    SIMPLE_ERROR("Force does not have the necessary dimensions");
  }
#endif //]

  if ( !calcForce && ( calcDiagonalHessian || calcOffDiagonalHessian ) ) {
    SIMPLE_ERROR("Inconsistant arguments: if you want to calcDiagonalHessian or calcOffDiagonalHessian you must calcForce");
  }
  if ( !calcDiagonalHessian & calcOffDiagonalHessian ) {
    SIMPLE_ERROR("Inconsistant arguments: if you want to calcOffDiagonalHessian you must calcDiagonalHessian");
  }

////    _lisp->profiler().pushTimerStates();

  if ( hasForce ) force->zero();
  if ( hasHessian ) hessian->zero();
  if ( hasHdAndD ) {
    LOG("Zeroing hdvec" );
    hdvec->zero();	// Zero the result
  }

  LOG("Starting evaluation of energy" );
	//
	// Evaluate the stretch term
	//
////	_lisp->profiler().timer(core::timerBond).start();
#define DUMP_FORCE 123
#ifdef DUMP_FORCE
  maybe_dump_force("start",force);
#endif
  if (this->_Stretch->isEnabled())
    totalEnergy += this->_Stretch->evaluateAllComponent( this->asSmartPtr(),
                                                         pos,
                                                         energyScale,
                                                         componentEnergy,
                                                         calcForce, force,
                                                         calcDiagonalHessian,
                                                         calcOffDiagonalHessian,
                                                         hessian, hdvec, dvec,
                                                         activeAtomMask,
                                                         debugInteractions );
////	_lisp->profiler().timer(core::timerBond).stop();
////	_lisp->profiler().timer(core::timerNonbond).start();
#ifdef DUMP_FORCE
  maybe_dump_force("stretch",force);
#endif
  if (this->_Nonbond->isEnabled())
    totalEnergy += this->_Nonbond->evaluateAllComponent( this->asSmartPtr(),
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
////	_lisp->profiler().timer(core::timerNonbond).stop();

//	_lisp->profiler().timer(core::timerAnchorRestraint).start();
#ifdef DUMP_FORCE
  maybe_dump_force("nonbond",force);
#endif
  if (this->_PointToLineRestraint->isEnabled())
    totalEnergy += this->_PointToLineRestraint->evaluateAllComponent( this->asSmartPtr(),
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
#ifdef DUMP_FORCE
  maybe_dump_force("pointToLineRestraint",force);
#endif
  if (this->_OutOfZPlane->isEnabled())
    totalEnergy += this->_OutOfZPlane->evaluateAllComponent( this->asSmartPtr(),
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
#ifdef DUMP_FORCE
  maybe_dump_force("oozp",force);
#endif

    	// More energy terms
  return totalEnergy;
}


SYMBOL_EXPORT_SC_(KeywordPkg,stretch_deviations);
SYMBOL_EXPORT_SC_(KeywordPkg,angle_deviations);


//
//
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------





string	SketchFunction_O::energyTermsEnabled()
{
  stringstream ss;
  ss.str("");
  ss << this->_Stretch->enabledAsString();
  ss << this->_PointToLineRestraint->enabledAsString();
  ss << this->_OutOfZPlane->enabledAsString();
  return ss.str();
}




#define	mathematicaSinNPhi(n,s,c) 0.0
#define	mathematicaCosNPhi(n,s,c) 0.0













#define	DELTA	0.00000001

double	SketchFunction_O::calculateNumericalDerivative(NVector_sp pos, core::T_sp energyScale, double delta, uint i, core::T_sp activeAtomMask )
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


double	SketchFunction_O::calculateNumericalSecondDerivative(NVector_sp pos, core::T_sp energyScale, double delta, uint i, uint j, core::T_sp activeAtomMask )
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
void	SketchFunction_O::evaluateNumericalForce(NVector_sp pos, core::T_sp energyScale, NVector_sp numForce, double delta, core::T_sp activeAtomMask )
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
void	SketchFunction_O::evaluateNumericalHessian(NVector_sp pos, core::T_sp energyScale, AbstractLargeSquareMatrix_sp hessian, bool calcOffDiagonal, double delta, core::T_sp activeAtomMask )
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
 * If there is a mis-match then dump the SketchFunction into the result.
 *
 */
ForceMatchReport_sp SketchFunction_O::checkIfAnalyticalForceMatchesNumericalForce(NVector_sp pos, core::T_sp energyScale, NVector_sp analyticalForce, core::T_sp activeAtomMask )
{
  ForceMatchReport_sp report;
  NVector_sp	numForce, tempForce;
  stringstream	result;
  double	analyticalMag, numericalMag, dot;
  double	avg;

  report = ForceMatchReport_O::create();

  numForce = NVector_O::create(pos->size());
  this->evaluateNumericalForce(pos,energyScale,numForce,DELTA,activeAtomMask);
  dot = dotProductWithActiveAtomMask(numForce,analyticalForce,activeAtomMask);
  numericalMag = magnitudeWithActiveAtomMask(numForce,activeAtomMask);
  analyticalMag = magnitudeWithActiveAtomMask(analyticalForce,activeAtomMask);
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


void	SketchFunction_O::dumpTerms()
{
  this->_Stretch->dumpTerms(this->_AtomTypes);
  this->_PointToLineRestraint->dumpTerms(this->_AtomTypes);
  this->_OutOfZPlane->dumpTerms(this->_AtomTypes);
}


void	SketchFunction_O::loadCoordinatesIntoVector(NVector_sp pos)
{
  if ( pos->size() != this->getNVectorSize()) {
    SIMPLE_ERROR("NVector is the incorrect length");
  }
  core::T_sp tsize = core::eval::funcall(_sym_node_table_size,this->_NodeTable);
  ASSERT(tsize.fixnump());
  for ( size_t ai=0; ai<tsize.unsafe_fixnum(); ai++ ) {
    core::T_sp node = core::eval::funcall(_sym_node_table_node_at_index,this->_NodeTable,core::make_fixnum(ai));
    core::T_sp tcoordinate_index_times3 = core::eval::funcall(_sym_node_table_coordinate_index_times3,this->_NodeTable,core::make_fixnum(ai));
    ASSERT(tcoordinate_index_times3.fixnump());
    size_t ci = tcoordinate_index_times3.unsafe_fixnum();
    core::T_sp tnode_pos = core::eval::funcall(_sym_node_get_position,node);
    geom::OVector3_sp node_pos = gc::As<geom::OVector3_sp>(tnode_pos);
    pos->setElement(ci,node_pos->getX());
    pos->setElement(ci+1, node_pos->getY());
    pos->setElement(ci+2, node_pos->getZ());
  }
}



void    SketchFunction_O::saveCoordinatesFromVector(NVector_sp pos)
{
  double                          x,y,z;
  Vector3                         v;
  core::T_sp tsize = core::eval::funcall(_sym_node_table_size,this->_NodeTable);
  ASSERT(tsize.fixnump());
  for ( size_t ai = 0; ai<tsize.unsafe_fixnum(); ai++ ) {
    core::T_sp node = core::eval::funcall(_sym_node_table_node_at_index,this->_NodeTable,core::make_fixnum(ai));
    core::T_sp tcoordinate_index_times3 = core::eval::funcall(_sym_node_table_coordinate_index_times3,this->_NodeTable,core::make_fixnum(ai));
    ASSERT(tcoordinate_index_times3.fixnump());
    size_t ci = tcoordinate_index_times3.unsafe_fixnum();
    x = pos->getElement(ci+0);
    y = pos->getElement(ci+1);
    z = pos->getElement(ci+2);
    LOG("Set atom({}) position ({},{},{})" , ci , x , y , z  );
    core::eval::funcall(_sym_node_set_position,node,geom::OVector3_O::make(x,y,z));
  }
}

#if 0
CL_LISPIFY_NAME("writeForceToAtoms");
CL_DEFMETHOD void    SketchFunction_O::writeForceToAtoms(NVector_sp force)
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
#endif


void    SketchFunction_O::saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force)
{
  this->saveCoordinatesFromVector(pos);
//  this->writeForceToAtoms(force);
}

#if 0
EnergyAtom*	SketchFunction_O::getEnergyAtomPointer(Atom_sp a)
{ 
  return this->_AtomTable->getEnergyAtomPointer(a);
};
#endif



void	SketchFunction_O::dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments)
{
  SIMPLE_ERROR("Handle dealWithProblem");
}


CL_LAMBDA(scoring-function position &key energy-scale velocity force force-dt delta-t-over-mass delta-t unfrozen)
CL_LISPIFY_NAME("sketch-function-velocity-verlet-step");
DOCGROUP(cando);
CL_DEFUN void chem__SketchFunction_velocity_verlet_step(SketchFunction_sp sketchFunc,
                                                        NVector_sp position,
                                                        core::T_sp energyScale,
                                                        NVector_sp velocity,
                                                        NVector_sp force,
                                                        NVector_sp force_dt,
                                                        NVector_sp delta_t_over_mass,
                                                        double delta_t,
                                                        core::T_sp tunfrozen)
{
//  SIMPLE_WARN("FIXactiveAtomMask - Is tunfrozen and activeAtomMask the same");
  core::SimpleBitVector_sp unfrozen;
  if (gc::IsA<core::SimpleBitVector_sp>(tunfrozen)) {
    unfrozen = gc::As_unsafe<core::SimpleBitVector_sp>(tunfrozen);
    if (unfrozen->length() != (position->length()/3)) {
      SIMPLE_ERROR("unfrozen must be a simple-bit-vector of length {} or NIL, it is a simple-bit-vector of length {}" , (position->length()), unfrozen->length());
    }
  } else if (tunfrozen.notnilp()) {
    SIMPLE_ERROR("unfrozen must be a simple-bit-vector or NIL");
  }
  size_t atom_idx = 0;
  for ( size_t idx = 0; idx<position->size(); idx += 3) {
    if (!unfrozen || unfrozen->testBit(atom_idx)==1) {
      double offsetx = delta_t*(*velocity)[idx+0] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+0];
      if (offsetx>1.5) offsetx = 1.5;
      if (offsetx<-1.5) offsetx = -1.5;
      (*position)[idx+0] = (*position)[idx+0] + offsetx;
      double offsety = delta_t*(*velocity)[idx+1] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+1];
      if (offsety>1.5) offsety = 1.5;
      if (offsety<-1.5) offsety = -1.5;
      (*position)[idx+1] = (*position)[idx+1] + offsety;
      double offsetz = delta_t*(*velocity)[idx+2] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+2];
      if (offsetz>1.5) offsetz = 1.5;
      if (offsetz<-1.5) offsetz = -1.5;
      (*position)[idx+2] = (*position)[idx+2] + offsetz;
    }
    atom_idx++;
  }
  sketchFunc->evaluateEnergyForce(position,energyScale,true,force_dt,tunfrozen);
  atom_idx = 0;
  for ( size_t idx = 0; idx<position->size(); idx+=3 ) {
    if (!unfrozen || unfrozen->testBit(atom_idx)==1) {
      (*velocity)[idx+0] = ((*velocity)[idx+0] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+0]+(*force_dt)[idx+0]))*sketchFunc->_VelocityScale.getX();
      (*force)[idx+0] = (*force_dt)[idx+0];
      (*velocity)[idx+1] = ((*velocity)[idx+1] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+1]+(*force_dt)[idx+1]))*sketchFunc->_VelocityScale.getY();
      (*force)[idx+1] = (*force_dt)[idx+1];
      (*velocity)[idx+2] = ((*velocity)[idx+2] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+2]+(*force_dt)[idx+2]))*sketchFunc->_VelocityScale.getZ();
      (*force)[idx+2] = (*force_dt)[idx+2];
    }
    atom_idx++;
  }
}


CL_DEFMETHOD void SketchFunction_O::resetSketchFunction()
{
  this->_Stretch->reset();
  this->_OutOfZPlane->reset();
  this->_PointToLineRestraint->reset();
  this->_Nonbond->reset();
}


core::List_sp SketchFunction_O::allComponents() const {
  core::List_sp result = nil<core::T_O>();
  result = core::Cons_O::create(this->_Nonbond,result);
  result = core::Cons_O::create(this->_PointToLineRestraint,result);
  result = core::Cons_O::create(this->_OutOfZPlane,result);
  result = core::Cons_O::create(this->_Stretch,result);
  return result;
}



};
