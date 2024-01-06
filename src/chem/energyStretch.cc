/*
    File: energyStretch.cc
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

#include <cando/chem/energyStretch.h>
#include <clasp/core/numerics.h>
#include <clasp/core/ql.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/energyAtomTable.h>
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
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


namespace chem {

#define STRETCH_APPLY_ATOM_MASK(I1,I2) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3))) goto SKIP_term; \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,EnergyStretch_O::staticClass(), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2)); \
    }

core::List_sp EnergyStretch::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kb),core::clasp_make_double_float(this->term.kb)),
                                  core::Cons_O::create(INTERN_(kw,r0),core::clasp_make_double_float(this->term.r0)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.I2)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                  core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2));
}

void EnergyStretch::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergyStretch");
}

#if XML_ARCHIVE
void	EnergyStretch::archive(core::ArchiveP node)
{
  node->attribute("kb",this->term.kb);
  node->attribute("r0",this->term.r0);
  node->attribute("I1",this->term.I1);
  node->attribute("I2",this->term.I2);
  node->attribute("a1",this->_Atom1);
  node->attribute("a2",this->_Atom2);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
  node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
  node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
  node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include <cando/chem/energy_functions/_Stretch_debugEvalSerialize.cc>
#endif //]
}
#endif



void EnergyStretch::defineFrom( FFStretch_sp stretch, EnergyAtom *ea1, EnergyAtom *ea2, double scale)
{
  LOG("Defining EnergyStretch with kb={} r0={}"
      , this->term.kb
      , this->term.r0  );
//    this->_K3 = stretch->_K3;
//    this->_K4 = stretch->_K4;
//    this->_Bci = stretch->_Bci;
  this->_Atom1 = ea1->atom();
  this->_Atom2 = ea2->atom();
  LOG(" Adding stretch between {} - {}"
      , this->_Atom1->description() , this->_Atom2->description() );
  this->term.I1 = ea1->coordinateIndexTimes3();
  this->term.I2 = ea2->coordinateIndexTimes3();
  this->term.kb = stretch->getKb_kcalPerAngstromSquared();
  this->term.r0 = stretch->getR0_Angstrom();
}

#if 0
adapt::QDomNode_sp	EnergyStretch::asXml()
{
  adapt::QDomNode_sp	node;
  Vector3	vdiff;
  node = adapt::QDomNode_O::create(env,"EnergyStretch");
  node->addAttributeHex("address",(unsigned long)(this));
  node->addAttributeString("atom1Name",this->_Atom1->getName());
  node->addAttributeString("atom2Name",this->_Atom2->getName());
//    node->addAttributeInt("idx1",this->_Atom1->_CoordinateIndex);
//    node->addAttributeInt("idx2",this->_Atom2->_CoordinateIndex);
  node->addAttributeInt("I1",this->term.I1);
  node->addAttributeInt("I2",this->term.I2);
  node->addAttributeDoubleScientific("r0",this->term.r0);
  node->addAttributeDoubleScientific("kb",this->term.kb);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
  adapt::QDomNode_sp xml = adapt::QDomNode_O::create(env,"Evaluated");
  xml->addAttributeBool("calcForce",this->_calcForce );
  xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
  xml->addAttributeBool("calcOffDiagonalHessian",
                        this->_calcOffDiagonalHessian );


//#include <cando/chem/energy_functions/_Stretch_debugEvalXml.cc>


  node->addChild(xml);
#endif
  return node;
}
#endif

 
void	EnergyStretch::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml, AtomTable_sp atomTable )
{
  int	i1, i2;
  i1 = xml->getAttributeInt("idx1");
  i2 = xml->getAttributeInt("idx2");
  LOG("parse i1={} i2={}" , i1 , i2  );
  this->_Atom1 = atomTable->findEnergyAtomWithCoordinateIndex(i1)->atom();
  this->_Atom2 = atomTable->findEnergyAtomWithCoordinateIndex(i2)->atom();
  LOG("_Atom1 = {}  _Atom2 = {}" , i1 , i2  );
  this->term.r0 = xml->getAttributeDouble("r0");
  this->term.kb = xml->getAttributeDouble("kb");
  this->term.I1 = xml->getAttributeInt("I1");
  this->term.I2 = xml->getAttributeInt("I2");
}



double	EnergyStretch::getR()
{
  Vector3	pos1, pos2, diff;
  pos1 = this->getAtom1()->getPosition();
  pos2 = this->getAtom2()->getPosition();
  diff = pos1.sub(pos2);
  return diff.length();
}




num_real _evaluateEnergyOnly_Stretch (
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




void	EnergyStretch_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                    AbstractLargeSquareMatrix_sp m,
                                                    core::T_sp activeAtomMask)
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
    ASSERT_NOT_NAN(v); \
    m->addToElement((i1)+(o1),(i2)+(o2),v);		\
  }
#undef	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {	\
    ASSERT_NOT_NAN(v); \
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
    for ( gctools::Vec0<EnergyStretch>::iterator si=this->_Terms.begin();
          si!=this->_Terms.end(); si++ ) {
#include <cando/chem/energy_functions/_Stretch_termCode.cc>
    }
  }
}


CL_DEFMETHOD core::T_sp EnergyStretch_O::stretchTermBetweenAtoms(Atom_sp x, Atom_sp y)
{
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    if ( (si->_Atom1 == x && si->_Atom2 == y) ||
         (si->_Atom1 == y && si->_Atom2 == x ) ) {
      return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kb), core::clasp_make_double_float(si->term.kb)),
                                      core::Cons_O::create(INTERN_(kw,r0), core::clasp_make_double_float(si->term.r0)));
    }
  }
  return nil<core::T_O>();
}

num_real EnergyStretch_O::evaluateAllComponent( ScoringFunction_sp score,
                                                NVector_sp 	pos,
                                                core::T_sp componentEnergy,
                                                bool 		calcForce,
                                                gc::Nilable<NVector_sp> 	force,
                                                bool		calcDiagonalHessian,
                                                bool		calcOffDiagonalHessian,
                                                gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                gc::Nilable<NVector_sp>	hdvec,
                                                gc::Nilable<NVector_sp> dvec,
                                                core::T_sp activeAtomMask ,
                                                core::T_sp debugInteractions )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  num_real totalEnergy = 0.0;
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

#define STRETCH_CALC_FORCE
#define STRETCH_CALC_DIAGONAL_HESSIAN
#define STRETCH_CALC_OFF_DIAGONAL_HESSIAN
#undef	STRETCH_SET_PARAMETER
#define	STRETCH_SET_PARAMETER(x)	{x = si->term.x;}
#undef	STRETCH_SET_POSITION
#define	STRETCH_SET_POSITION(x,ii,of)	{x = pos->getElement(ii+of);}
#undef	STRETCH_ENERGY_ACCUMULATE
#define	STRETCH_ENERGY_ACCUMULATE(e) {ASSERT_NOT_NAN(e);  totalEnergy += (e); }
#undef	STRETCH_FORCE_ACCUMULATE
#undef	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
#undef	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_FORCE_ACCUMULATE 		ForceAcc
#define	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
#pragma clang diagnostic pop
  fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
  fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
  num_real x1,y1,z1,x2,y2,z2,kb,r0;
  int I1, I2,i;
  gctools::Vec0<EnergyStretch>::iterator si;
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
  }
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY(("%s }") , this->className());
  }
  maybeSetEnergy( componentEnergy, EnergyStretch_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}

void	EnergyStretch_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm( NVector_sp 	pos)
{
  IMPLEMENT_ME();
#if 0
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;
  {

#define STRETCH_CALC_FORCE
#define STRETCH_CALC_DIAGONAL_HESSIAN
#define STRETCH_CALC_OFF_DIAGONAL_HESSIAN
#undef	STRETCH_SET_PARAMETER
#define	STRETCH_SET_PARAMETER(x)	{x = si->term.x;}
#undef	STRETCH_SET_POSITION
#define	STRETCH_SET_POSITION(x,ii,of)	{x = pos->getElement(ii+of);}
#undef	STRETCH_ENERGY_ACCUMULATE
#define	STRETCH_ENERGY_ACCUMULATE(e) {}
#undef	STRETCH_FORCE_ACCUMULATE
#define	STRETCH_FORCE_ACCUMULATE(i,o,v) {}
#undef	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
    {
      
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
#pragma clang diagnostic pop
      fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
      fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
      num_real x1,y1,z1,x2,y2,z2,kb,r0;
      int I1, I2,i;
//	    stretchScale = this->getScale();
      gctools::Vec0<EnergyStretch>::iterator si;
      for ( i=0,si=this->_Terms.begin();
            si!=this->_Terms.end(); si++,i++ ) {
		    // ************* Evaluate the stretch energy/force/hessian
		    // using code generated by Mathematica ***************
#include <cando/chem/energy_functions/_Stretch_termCode.cc>
        LOG("Energy = {}" , Energy  );
        int index = i;
        LOG("x1 = {}" , x1  );
        LOG("y1 = {}" , y1  );
        LOG("z1 = {}" , z1  );
        LOG("x2 = {}" , x2  );
        LOG("y2 = {}" , y2  );
        LOG("z2 = {}" , z2  );
#include <cando/chem/energy_functions/_Stretch_debugFiniteDifference.cc>
      }
    }
  }
#endif
}


SYMBOL_EXPORT_SC_(KeywordPkg,stretch);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);
SYMBOL_EXPORT_SC_(KeywordPkg,r0);
SYMBOL_EXPORT_SC_(KeywordPkg,indices);
SYMBOL_EXPORT_SC_(KeywordPkg,stretch_deviation);
SYMBOL_EXPORT_SC_(KeywordPkg,force);


core::List_sp	EnergyStretch_O::checkForBeyondThresholdInteractionsWithPosition(NVector_sp pos, double threshold)
{
  IMPLEMENT_ME();
#if 0
  ql::list result;
  bool calcForce = true;
//	bool	calcForce = false;
//	bool	calcDiagonalHessian = false;
//	bool	calcOffDiagonalHessian = false;

#define STRETCH_CALC_FORCE
#undef STRETCH_CALC_DIAGONAL_HESSIAN
#undef STRETCH_CALC_OFF_DIAGONAL_HESSIAN
#undef	STRETCH_SET_PARAMETER
#define	STRETCH_SET_PARAMETER(x)	{x = si->term.x;}
#undef	STRETCH_SET_POSITION
#define	STRETCH_SET_POSITION(x,ii,of)	{x = pos->getElement(ii+of);}
#undef	STRETCH_ENERGY_ACCUMULATE
#define	STRETCH_ENERGY_ACCUMULATE(e) {}
#undef	STRETCH_FORCE_ACCUMULATE
#define	STRETCH_FORCE_ACCUMULATE(i,o,v) {}
#undef	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

  {
    
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    num_real x1,y1,z1,x2,y2,z2,kb,r0;
    int I1, I2,i;
//	    stretchScale = this->getScale();
    gctools::Vec0<EnergyStretch>::iterator si;
    for ( i=0,si=this->_Terms.begin();
          si!=this->_Terms.end(); si++,i++ ) {
		// ************* Evaluate the stretch energy/force/hessian
		// using code generated by Mathematica ***************
#include <cando/chem/energy_functions/_Stretch_termCode.cc>
      ASSERT(r0 != 0.0);
      if ( fabs(StretchDeviation)/r0 > threshold ) {
        Atom_sp a1, a2;
        a1 = (*si)._Atom1;
        a2 = (*si)._Atom2;
        ql::list one_deviation;
        one_deviation << kw::_sym_stretch
                      << kw::_sym_atoms << core::Cons_O::createList(a1,a2)
                      << kw::_sym_indices << core::Cons_O::createList(core::make_fixnum(I1/3),core::make_fixnum(I2/3))
                      << kw::_sym_r0 << core::clasp_make_double_float(r0)
                      << kw::_sym_stretch_deviation << core::clasp_make_double_float(StretchDeviation)
                      << kw::_sym_force << (ql::list()
                                            << geom::OVector3_O::make(fx1,fy1,fz1)
                                            << geom::OVector3_O::make(fx2,fy2,fz2)).result();
        result << one_deviation.result();
      }
    }
  }
  return result.result();
#endif
}




void	EnergyStretch_O::dumpTerms(core::HashTable_sp atomTypes)
{
  gctools::Vec0<EnergyStretch>::iterator	esi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4;
  LOG("Dumping EnergyStretch terms");
  uint idx;
  for ( idx = 0, esi=this->_Terms.begin(); esi!=this->_Terms.end(); esi++, idx++ )
  {
    as1 = atomLabel(esi->_Atom1);
    as2 = atomLabel(esi->_Atom2);
    if ( as1 < as2 ) {
      str1 = as1;
      str2 = as2;
    } else {
      str2 = as1;
      str1 = as2;
    }
    core::clasp_write_string(fmt::format("TERM 1BND {:<9} - {:<9} {:8.2f} {:8.2f}\n" , str1 , str2 , esi->term.kb , esi->term.r0 ));
    LOG("TERM[#{}] 1BND %-9s - %-9s %8.2lf %8.2lf ; I1={} I2={}"
        , idx
        , str1 
        , str2
        , esi->term.kb
        , esi->term.r0
        , esi->term.I1
        , esi->term.I2 );
  }
}










void EnergyStretch_O::initialize()
{
  this->Base::initialize();
  this->setErrorThreshold(0.05);
}

void EnergyStretch_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

void EnergyStretch_O::addTerm(const EnergyStretch& term)
{
  this->_Terms.push_back(term);
}

string EnergyStretch_O::beyondThresholdInteractionsAsString()
{
  return component_beyondThresholdInteractionsAsString<EnergyStretch_O,EnergyStretch>(*this);
}


SYMBOL_EXPORT_SC_(KeywordPkg,kb);
SYMBOL_EXPORT_SC_(KeywordPkg,r0);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,atom1);
SYMBOL_EXPORT_SC_(KeywordPkg,atom2);

CL_DEFMETHOD core::List_sp EnergyStretch_O::extract_vectors_as_alist() const {
  size_t size = this->_Terms.size();
  printf("%s:%d The number of EnergyStretch terms -> %lu\n", __FILE__, __LINE__, size);
  NVector_sp kb_vec = NVector_O::make(size);
  NVector_sp r0_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_sp atom1_vec    = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom2_vec    = core::SimpleVector_O::make(size);
  for (size_t i=0;i<size;++i) {
    const EnergyStretch& entry = this->_Terms[i];
    (*kb_vec)[i] = entry.term.kb;
    (*r0_vec)[i] = entry.term.r0;
    (*i1_vec)[i] = entry.term.I1;
    (*i2_vec)[i] = entry.term.I2;
    (*atom1_vec)[i] = entry._Atom1 ;
    (*atom2_vec)[i] = entry._Atom2;
  }
  return core::Cons_O::createList(core::Cons_O::create(kw::_sym_kb,kb_vec),
                                  core::Cons_O::create(kw::_sym_r0,r0_vec),
                                  core::Cons_O::create(kw::_sym_i1,i1_vec),
                                  core::Cons_O::create(kw::_sym_i2,i2_vec),
                                  core::Cons_O::create(kw::_sym_atom1,atom1_vec),
                                  core::Cons_O::create(kw::_sym_atom2,atom2_vec));
};

CL_DEFMETHOD void EnergyStretch_O::fill_from_vectors_in_alist(core::List_sp vectors)
{
  NVector_sp kb_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_kb));
  NVector_sp r0_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_r0));
  core::SimpleVector_int32_t_sp i1_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i1));
  core::SimpleVector_int32_t_sp i2_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i2));
  core::SimpleVector_sp atom1_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom1));
  core::SimpleVector_sp atom2_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom2));
  this->_Terms.resize(kb_vec->length());
  for (size_t i=0,iEnd(kb_vec->length());i<iEnd;++i) {
    EnergyStretch& entry = this->_Terms[i];
    entry.term.kb = (*kb_vec)[i];
    entry.term.r0 = (*r0_vec)[i];
    entry.term.I1 = (*i1_vec)[i];
    entry.term.I2 = (*i2_vec)[i];
    entry._Atom1 = gc::As_unsafe<Atom_sp>((*atom1_vec)[i]);
    entry._Atom2 = gc::As_unsafe<Atom_sp>((*atom2_vec)[i]);
  }
}

CL_DEFMETHOD void EnergyStretch_O::addStretchTerm(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, double kb, double r0) {
  EnergyAtom* ea1 = atomTable->getEnergyAtomPointer(a1);
  EnergyAtom* ea2 = atomTable->getEnergyAtomPointer(a2);
  EnergyStretch energyStretch(a1,a2,ea1->coordinateIndexTimes3(),ea2->coordinateIndexTimes3(),kb,r0);
  this->addTerm(energyStretch);
}

CL_DOCSTRING(R"dx(Invoke a callback for every term in the energy-stretch.
The callback takes the lambda-list (index atom1 atom2 atom1-index*3 atom2-index*3 kb r0).
Both kb and r0 are single-floats to avoid consing.
Use the modify-stretch-term-kb and modify-stretch-term-r0 to modify the parameters.)dx")
CL_DEFMETHOD void EnergyStretch_O::walkStretchTerms(core::T_sp callback)
{
  for (size_t i=0;i<this->_Terms.size();++i) {
    const EnergyStretch& entry = this->_Terms[i];
    core::eval::funcall(callback,core::make_fixnum(i),
                        entry._Atom1,
                        entry._Atom2,
                        core::make_fixnum(entry.term.I1),
                        core::make_fixnum(entry.term.I2),
                        core::make_single_float(entry.term.kb),
                        core::make_single_float(entry.term.r0));
  };
};

CL_DEFMETHOD void EnergyStretch_O::modifyStretchTermKb(size_t index, float kb) {
  if (index<this->_Terms.size()) {
    EnergyStretch& entry = this->_Terms[index];
    entry.term.kb = kb;
    return;
  }
  SIMPLE_ERROR("index {} is out of bounds as a energy-stretch term index (#entries {})" , index , this->_Terms.size() );
}

CL_DEFMETHOD void EnergyStretch_O::modifyStretchTermR0(size_t index, float r0) {
  if (index<this->_Terms.size()) {
    EnergyStretch& entry = this->_Terms[index];
    entry.term.r0 = r0;
    return;
  }
  SIMPLE_ERROR("index {} is out of bounds as a energy-stretch term index (#entries {})" , index , this->_Terms.size() );
}

void EnergyStretch_O::reset()
{
  this->_Terms.clear();
}

CL_DEFMETHOD
core::List_sp	EnergyStretch_O::lookupStretchTerms(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, core::HashTable_sp atomTypes)
{
  ql::list  result;
  core::T_sp tia1 = atomTable->_AtomTableIndices->gethash(a1);
  core::T_sp tia2 = atomTable->_AtomTableIndices->gethash(a2);
  if (!tia1.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a1));
  if (!tia2.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a2));
  for (auto edi=this->_Terms.begin();edi!=this->_Terms.end();edi++) {
    if ((edi->_Atom1==a1 &&
         edi->_Atom2==a2)
        || (edi->_Atom2==a1 &&
            edi->_Atom1==a2)) {
      ql::list oneResult;
      oneResult << INTERN_(kw,type1) << edi->_Atom1->getType(atomTypes)
                << INTERN_(kw,type2) << edi->_Atom2->getType(atomTypes)
                << INTERN_(kw,kb) << core::clasp_make_double_float(edi->term.kb)
                << INTERN_(kw,r0) << core::clasp_make_double_float(edi->term.r0);
      result << oneResult.result();
    }
  }
  return result.result();
}

EnergyStretch_sp EnergyStretch_O::copyFilter(core::T_sp keepInteraction) {
  EnergyStretch_sp copy = EnergyStretch_O::create();
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    if ( skipInteraction( keepInteraction, EnergyStretch_O::staticClass(), a1, a2 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}


};
