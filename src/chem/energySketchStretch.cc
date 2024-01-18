/*
    File: energySketchStretch.cc
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

#include <cando/chem/energySketchStretch.h>
#include <clasp/core/numerics.h>
#include <clasp/core/ql.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/array.h>
#include <cando/geom/ovector3.h>
#include <cando/chem/matter.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


namespace chem {

#define STRETCH_APPLY_ATOM_MASK(I1,I2) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         ) \
    ) goto SKIP_term;
#define STRETCH_DEBUG_INTERACTIONS(I1,I2) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,EnergySketchStretch_O::static_classSymbol(), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2)); \
    }

core::List_sp EnergySketchStretch::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kb),core::clasp_make_double_float(this->term.kb)),
                                  core::Cons_O::create(INTERN_(kw,r0),core::clasp_make_double_float(this->term.r0)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.I2)));
}

void EnergySketchStretch::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergySketchStretch");
}

num_real	_evaluateEnergyOnly_SketchStretch (
                                     num_real x1,
                                     num_real y1,
                                     num_real z1,
                                     num_real x2,
                                     num_real y2,
                                     num_real z2,
                                     num_real r0,
                                     num_real kb )
{
  IMPLEMENT_ME();
#if 0
#undef	STRETCH_SET_PARAMETER
#define	STRETCH_SET_PARAMETER(x)	{}
#undef	STRETCH_SET_POSITION
#define	STRETCH_SET_POSITION(x,ii,of)	{}
#undef	STRETCH_ENERGY_ACCUMULATE
#define	STRETCH_ENERGY_ACCUMULATE(e) {}
#undef	STRETCH_FORCE_ACCUMULATE
#define	STRETCH_FORCE_ACCUMULATE(i,o,v) {}
#undef	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	STRETCH_CALC_FORCE	// Don't calculate FORCE or HESSIAN


#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
#pragma clang diagnostic pop
#pragma clang diagnostic pop

#if !USE_EXPLICIT_DECLARES
  num_real fx1 = 0.0; 
  num_real fy1 = 0.0; 
  num_real fz1 = 0.0;
  num_real fx2 = 0.0;
  num_real fy2 = 0.0;
  num_real fz2 = 0.0;
#endif
#include <cando/chem/energy_functions/_Stretch_termCode.cc>

  return Energy;
#endif
}




void	EnergySketchStretch_O::setupHessianPreconditioner(
                                                    NVector_sp nvPosition,
                                                    AbstractLargeSquareMatrix_sp m,
                                                    core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  core::T_sp debugInteractions = nil<core::T_O>();
  MAYBE_SETUP_DEBUG_INTERACTIONS(false);
  bool		calcForce = true;
  bool		calcDiagonalHessian = true;
  bool		calcOffDiagonalHessian = true;
#undef	STRETCH_SET_PARAMETER
#define	STRETCH_SET_PARAMETER(x)	{x = si->term.x;}
#undef	STRETCH_SET_POSITION
#define	STRETCH_SET_POSITION(x,ii,of)	{x = nvPosition->getElement(ii+of);}
#undef	STRETCH_ENERGY_ACCUMULATE
#define	STRETCH_ENERGY_ACCUMULATE(e) {}
#undef	STRETCH_FORCE_ACCUMULATE
#define	STRETCH_FORCE_ACCUMULATE(i,o,v) {}
#undef	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {	\
    m->addToElement((i1)+(o1),(i2)+(o2),v);		\
  }
#undef	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {	\
    m->addToElement((i1)+(o1),(i2)+(o2),v);			\
  }
#define STRETCH_CALC_FORCE
#define STRETCH_CALC_DIAGONAL_HESSIAN
#define STRETCH_CALC_OFF_DIAGONAL_HESSIAN

  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
#pragma clang diagnostic pop
#if !USE_EXPLICIT_DECLARES
    num_real fx1 = 0.0;
    num_real fy1 = 0.0;
    num_real fz1 = 0.0;
    num_real fx2 = 0.0;
    num_real fy2 = 0.0;
    num_real fz2 = 0.0;
#endif
    num_real x1,y1,z1,x2,y2,z2,kb,r0;
    int I1, I2;
//	stretchScale = this->getScale();
    for ( gctools::Vec0<EnergySketchStretch>::iterator si=this->_Terms.begin();
          si!=this->_Terms.end(); si++ ) {
#include <cando/chem/energy_functions/_Stretch_termCode.cc>
    }
  }
}


num_real EnergySketchStretch_O::evaluateAllComponent( ScoringFunction_sp score,
                                                    NVector_sp 	pos,
                                                    core::T_sp componentEnergy,
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
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  this->_Evaluations++;
  if ( this->_DebugEnergy ) {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(("%s {\n") , this->className());
  }
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
  num_real totalEnergy = 0.0;
#define STRETCH_CALC_FORCE
#define STRETCH_CALC_DIAGONAL_HESSIAN
#define STRETCH_CALC_OFF_DIAGONAL_HESSIAN
#undef	STRETCH_SET_PARAMETER
#define	STRETCH_SET_PARAMETER(x)	{x = si->term.x;}
#undef	STRETCH_SET_POSITION
#define	STRETCH_SET_POSITION(x,ii,of)	{x = pos->getElement(ii+of);}
#undef	STRETCH_ENERGY_ACCUMULATE
#define	STRETCH_ENERGY_ACCUMULATE(e) totalEnergy += (e);
#undef	STRETCH_FORCE_ACCUMULATE
#undef	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
#undef	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_FORCE_ACCUMULATE 		ForceAcc
#define	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc

  { 
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    num_real x1,y1,z1,x2,y2,z2,kb,r0;
    int I1, I2,i;
    gctools::Vec0<EnergySketchStretch>::iterator si;
    for ( i=0,si=this->_Terms.begin();
          si!=this->_Terms.end(); si++,i++ ) {
#ifdef DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfTermsToCalculate ) {
          break;
        }
      }
#endif

		// ************* Evaluate the stretch energy/force/hessian
		// using code generated by Mathematica ***************
#include <cando/chem/energy_functions/_Stretch_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
      si->_calcForce = calcForce;
      si->_calcDiagonalHessian = calcDiagonalHessian;
      si->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define EVAL_SET(var,val) {si->eval.var=val;}
#include <cando/chem/energy_functions/_Stretch_debugEvalSet.cc>
#endif //]
      if ( this->_DebugEnergy ) {
        LOG_ENERGY(( "MEISTER stretch %d args cando\n") , (i+1) );
        LOG_ENERGY(( "MEISTER stretch %d address %x \n") , (i+1) , (&(*si)) );
        LOG_ENERGY(( "MEISTER stretch %d r0 %5.3lf\n") , (i+1) , r0 );
        LOG_ENERGY(( "MEISTER stretch %d kb %5.1lf\n") , (i+1) , kb );
        LOG_ENERGY(( "MEISTER stretch %d x1 %5.3lf %d\n") , (i+1) , x1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER stretch %d y1 %5.3lf %d\n") , (i+1) , y1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER stretch %d z1 %5.3lf %d\n") , (i+1) , z1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER stretch %d x2 %5.3lf %d\n") , (i+1) , x2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER stretch %d y2 %5.3lf %d\n") , (i+1) , y2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER stretch %d z2 %5.3lf %d\n") , (i+1) , z2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER stretch %d results\n") , (i+1) );
        LOG_ENERGY(( "MEISTER stretch %d Energy %lf\n") , (i+1) , Energy );
        if ( calcForce ) {
          LOG_ENERGY(( "MEISTER stretch %d fx1 %lf %d\n") , (i+1) , fx1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER stretch %d fy1 %lf %d\n") , (i+1) , fy1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER stretch %d fz1 %lf %d\n") , (i+1) , fz1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER stretch %d fx2 %lf %d\n") , (i+1) , fx2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER stretch %d fy2 %lf %d\n") , (i+1) , fy2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER stretch %d fz2 %lf %d\n") , (i+1) , fz2 , (I2/3+1) );
        }
        LOG_ENERGY(( "MEISTER stretch %d stop\n") , (i+1) );
      }
		/* Add the forces */

//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz2>10000.0);
    }
  }
  maybeSetEnergy( componentEnergy, EnergySketchStretch_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}


SYMBOL_EXPORT_SC_(KeywordPkg,stretch);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);
SYMBOL_EXPORT_SC_(KeywordPkg,r0);
SYMBOL_EXPORT_SC_(KeywordPkg,indices);
SYMBOL_EXPORT_SC_(KeywordPkg,stretch_deviation);
SYMBOL_EXPORT_SC_(KeywordPkg,force);




void EnergySketchStretch_O::initialize()
{
  this->Base::initialize();
  this->setErrorThreshold(0.05);
}

void EnergySketchStretch_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

void EnergySketchStretch_O::addTerm(const EnergySketchStretch& term)
{
  this->_Terms.push_back(term);
}

SYMBOL_EXPORT_SC_(KeywordPkg,kb);
SYMBOL_EXPORT_SC_(KeywordPkg,r0);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,atom1);
SYMBOL_EXPORT_SC_(KeywordPkg,atom2);

CL_DEFMETHOD core::List_sp EnergySketchStretch_O::extract_vectors_as_alist() const {
  size_t size = this->_Terms.size();
  printf("%s:%d The number of EnergySketchStretch terms -> %lu\n", __FILE__, __LINE__, size);
  NVector_sp kb_vec = NVector_O::make(size);
  NVector_sp r0_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  for (size_t i=0;i<size;++i) {
    const EnergySketchStretch& entry = this->_Terms[i];
    (*kb_vec)[i] = entry.term.kb;
    (*r0_vec)[i] = entry.term.r0;
    (*i1_vec)[i] = entry.term.I1;
    (*i2_vec)[i] = entry.term.I2;
  }
  return core::Cons_O::createList(core::Cons_O::create(kw::_sym_kb,kb_vec),
                                  core::Cons_O::create(kw::_sym_r0,r0_vec),
                                  core::Cons_O::create(kw::_sym_i1,i1_vec),
                                  core::Cons_O::create(kw::_sym_i2,i2_vec));
};

CL_DEFMETHOD void EnergySketchStretch_O::fill_from_vectors_in_alist(core::List_sp vectors)
{
  NVector_sp kb_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_kb));
  NVector_sp r0_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_r0));
  core::SimpleVector_int32_t_sp i1_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i1));
  core::SimpleVector_int32_t_sp i2_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i2));
  this->_Terms.resize(kb_vec->length());
  for (size_t i=0,iEnd(kb_vec->length());i<iEnd;++i) {
    EnergySketchStretch& entry = this->_Terms[i];
    entry.term.kb = (*kb_vec)[i];
    entry.term.r0 = (*r0_vec)[i];
    entry.term.I1 = (*i1_vec)[i];
    entry.term.I2 = (*i2_vec)[i];
  }
}

CL_DOCSTRING(R"dx(The indices I1 and I2 must be coordinat indices, in other words index*3.)dx");
CL_DEFMETHOD void EnergySketchStretch_O::addSketchStretchTerm(size_t i1, size_t i2, double kb, double r0) {
  EnergySketchStretch energySketchStretch(i1,i2,kb,r0);
  this->addTerm(energySketchStretch);
}

CL_DOCSTRING(R"dx(Invoke a callback for every term in the energy-stretch.)dx");
CL_DOCSTRING_LONG(R"dx(
The callback takes the lambda-list (index atom1 atom2 atom1-index*3 atom2-index*3 kb r0).
Both kb and r0 are single-floats to avoid consing.
Use the modify-stretch-term-kb and modify-stretch-term-r0 to modify the parameters.)dx")
CL_DEFMETHOD void EnergySketchStretch_O::walkSketchStretchTerms(core::T_sp callback)
{
  for (size_t i=0;i<this->_Terms.size();++i) {
    const EnergySketchStretch& entry = this->_Terms[i];
    core::eval::funcall(callback,core::make_fixnum(i),
                        core::make_fixnum(entry.term.I1),
                        core::make_fixnum(entry.term.I2),
                        core::make_single_float(entry.term.kb),
                        core::make_single_float(entry.term.r0));
  };
};

CL_DEFMETHOD void EnergySketchStretch_O::modifySketchStretchTermKb(size_t index, float kb) {
  if (index<this->_Terms.size()) {
    EnergySketchStretch& entry = this->_Terms[index];
    entry.term.kb = kb;
    return;
  }
  SIMPLE_ERROR("index {} is out of bounds as a energy-stretch term index (#entries {})" , index , this->_Terms.size() );
}

CL_DEFMETHOD void EnergySketchStretch_O::modifySketchStretchTermR0(size_t index, float r0) {
  if (index<this->_Terms.size()) {
    EnergySketchStretch& entry = this->_Terms[index];
    entry.term.r0 = r0;
    return;
  }
  SIMPLE_ERROR("index {} is out of bounds as a energy-stretch term index (#entries {})" , index , this->_Terms.size() );
}

void EnergySketchStretch_O::reset()
{
  this->_Terms.clear();
}

}
