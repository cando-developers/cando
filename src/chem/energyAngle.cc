/*
    File: energyAngle.cc
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

#include <clasp/core/foundation.h>
#include <clasp/core/numerics.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/ql.h>
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


namespace chem {

#define ENERGY_FUNCTION I1, I2, I3, activeAtomMask

#define ANGLE_APPLY_ATOM_MASK(I1,I2,I3) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         ) \
    ) goto SKIP_term;


//
// Copy this from implementAmberFunction.cc
//
double	_evaluateEnergyOnly_Angle(
    int I1, int I2, int I3, core::T_sp activeAtomMask,
		num_real x1, num_real y1, num_real z1,
		num_real x2, num_real y2, num_real z2,
		num_real x3, num_real y3, num_real z3,
		num_real t0, num_real kt )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();

#undef	ANGLE_SET_PARAMETER
#define	ANGLE_SET_PARAMETER(x)	{}
#undef	ANGLE_SET_POSITION
#define	ANGLE_SET_POSITION(x,ii,of)	{}
#undef	ANGLE_ENERGY_ACCUMULATE
#define	ANGLE_ENERGY_ACCUMULATE(e) {}
#undef	ANGLE_FORCE_ACCUMULATE
#define	ANGLE_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANGLE_CALC_FORCE	// Don't calculate FORCE or HESSIAN


#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
#if !USE_EXPLICIT_DECLARES
	num_real fx1;
	num_real fy1;
	num_real fz1;
	num_real fx2;
	num_real fy2;
	num_real fz2;
	num_real fx3;
	num_real fy3;
	num_real fz3;
#endif
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
    bool IllegalAngle = false;
#define ANGLE_DEBUG_INTERACTIONS(I1,I2,I3)
#include <cando/chem/energy_functions/_Angle_termCode.cc>
#undef ANGLE_DEBUG_INTERACTIONS

    return Energy;
}


CL_LAMBDA((energy-angle chem:energy-angle) pos &optional activeAtomMask);
CL_DEFMETHOD void	EnergyAngle_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(chem::NVector_sp 	pos, core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//
#define ANGLE_CALC_FORCE
#define ANGLE_CALC_DIAGONAL_HESSIAN
#define ANGLE_CALC_OFF_DIAGONAL_HESSIAN
#undef ANGLE_SET_PARAMETER
#define ANGLE_SET_PARAMETER(x)	{x=ai->term.x;}
#undef ANGLE_SET_POSITION
#define ANGLE_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef ANGLE_ENERGY_ACCUMULATE
#define ANGLE_ENERGY_ACCUMULATE(e) {}
#undef ANGLE_FORCE_ACCUMULATE
#define ANGLE_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0;
    int I1, I2, I3,i;
    gctools::Vec0<EnergyAngle>::iterator ai;
    for ( i=0,ai=this->_Terms.begin();
          ai!=this->_Terms.end(); ai++,i++ ) {
      LOG("ai->term.kt = {}" , ai->term.kt  );
      LOG("ai->term.t0 = {}" , ai->term.t0  );
      LOG("angleScale = {}" , angleScale  );
      LOG("ai->x1 = {}" , pos->element(ai->term.I1 ) );
      LOG("ai->y1 = {}" , pos->element(ai->term.I1+1 ) );
      LOG("ai->z1 = {}" , pos->element(ai->term.I1+2 ) );
      LOG("ai->x2 = {}" , pos->element(ai->term.I2 ) );
      LOG("ai->y2 = {}" , pos->element(ai->term.I2+1 ) );
      LOG("ai->z2 = {}" , pos->element(ai->term.I2+2 ) );
      LOG("ai->x3 = {}" , pos->element(ai->term.I3 ) );
      LOG("ai->y3 = {}" , pos->element(ai->term.I3+1 ) );
      LOG("ai->z3 = {}" , pos->element(ai->term.I3+2 ) );
      bool IllegalAngle = false;
#define ANGLE_DEBUG_INTERACTIONS(I1,I2,I3)
#include	<cando/chem/energy_functions/_Angle_termCode.cc>
#undef ANGLE_DEBUG_INTERACTIONS
      LOG("Energy = {}" , Energy  );
      LOG("x1 = {}" , x1  );
      LOG("y1 = {}" , y1  );
      LOG("z1 = {}" , z1  );
      LOG("x2 = {}" , x2  );
      LOG("y2 = {}" , y2  );
      LOG("z2 = {}" , z2  );
      LOG("x3 = {}" , x3  );
      LOG("y3 = {}" , y3  );
      LOG("z3 = {}" , z3  );
      int index = i;
#include <cando/chem/energy_functions/_Angle_debugFiniteDifference.cc>
    }
    core::lisp_write(fmt::format("There are {} fails out of {} terms\n", fails, i ));
  }
}





#define ANGLE_DEBUG_INTERACTIONS(I1,I2,I3) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,EnergyAngle_O::static_classSymbol(), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3)); \
    }

core::List_sp EnergyAngle::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kt),core::clasp_make_double_float(this->term.kt)),
                                  core::Cons_O::create(INTERN_(kw,t0),core::clasp_make_double_float(this->term.t0)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.I2)),
                                  core::Cons_O::create(INTERN_(kw,i3), core::make_fixnum(this->term.I3)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                  core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2),
                                  core::Cons_O::create(INTERN_(kw,atom3), this->_Atom3));
}

void EnergyAngle::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergyAngle");
}

#ifdef XML_ARCHIVE
    void	EnergyAngle::archive(core::ArchiveP node)
{
//    node->attribute("_Term",this->_Term);
//    node->attribute("_Type1",this->_Type1);
//    node->attribute("_Type2",this->_Type2);
//    node->attribute("_Type3",this->_Type3);
//    node->attribute("_K3",this->_K3);
//    node->attribute("_K4",this->_K4);
//    node->attribute("_Ub_k",this->_Ub_k);
//    node->attribute("_Ub_len",this->_Ub_len);
    node->attribute("kt",this->term.kt);
    node->attribute("t0",this->term.t0);
    node->attribute("I1",this->term.I1);
    node->attribute("I2",this->term.I2);
    node->attribute("I3",this->term.I3);
    node->attribute("a1",this->_Atom1);
    node->attribute("a2",this->_Atom2);
    node->attribute("a3",this->_Atom3);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
    node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
    node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
    node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include <cando/chem/energy_functions/_Angle_debugEvalSerialize.cc>
#endif //]
}
#endif




void EnergyAngle::defineFrom( FFAngle_sp term , EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, double scale)
{

//    this->_Term = term;
//    this->_Type1 = term->_Type1;
//    this->_Type2 = term->_Type2;
//    this->_Type3 = term->_Type3;
    this->term.t0 = term->_AngRad;
    this->term.kt = term->getK2_kcalPerRadianSquared()*scale;
//    this->_K3 = term->_K3;
//    this->_K4 = term->_K4;
//    this->_Ub_k = term->_Ub_k;
//    this->_Ub_len = term->_Ub_len;
    this->_Atom1 = ea1->atom();
    this->_Atom2 = ea2->atom();
    this->_Atom3 = ea3->atom();
    this->term.I1 = ea1->coordinateIndexTimes3();
    this->term.I2 = ea2->coordinateIndexTimes3();
    this->term.I3 = ea3->coordinateIndexTimes3();
}



#if 0 //[
void	EnergyAngle::defineMissing( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3 )
{
    this->_Type1 = ea1->_Atom->getType(atomTypes);
    this->_Type2 = ea2->_Atom->getType(atomTypes);
    this->_Type3 = ea3->_Atom->getType(atomTypes);
}
#endif //]

double EnergyAngle::getT()
{
    Vector3	pos1, pos2, pos3;
num_real	angle;
    pos1 = this->_Atom1->getPosition();
    pos2 = this->_Atom2->getPosition();
    pos3 = this->_Atom3->getPosition();
    angle = geom::calculateAngle( pos1, pos2, pos3 );
    return angle;
}

#if 0
    adapt::QDomNode_sp	EnergyAngle::asXml()
{
    adapt::QDomNode_sp	node;
    Vector3	vdiff;

    node = adapt::QDomNode_O::create(env,"EnergyAngle");
    node->addAttributeString("atom1Name",this->_Atom1->getName());
    node->addAttributeString("atom2Name",this->_Atom2->getName());
    node->addAttributeString("atom3Name",this->_Atom3->getName());
//    node->addAttributeString("atom1Type",this->_Type1 );
//    node->addAttributeString("atom2Type",this->_Type2 );
//    node->addAttributeString("atom3Type",this->_Type3 );
    node->addAttributeInt("I1",this->term.I1);
    node->addAttributeInt("I2",this->term.I2);
    node->addAttributeInt("I3",this->term.I3);
    node->addAttributeDoubleScientific("kt",this->term.kt);
    node->addAttributeDoubleScientific("t0",this->term.t0);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
    adapt::QDomNode_sp xml = adapt::QDomNode_O::create(env,"Evaluated");
    xml->addAttributeBool("calcForce",this->_calcForce );
    xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
    xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include <_Angle_debugEvalXml.cc>
    node->addChild(xml);
#endif
    return node;
}

    void	EnergyAngle::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
					AtomTable_sp atomTable )
{
    int i1 = xml->getAttributeInt("idx1");
    int i2 = xml->getAttributeInt("idx2");
    int i3 = xml->getAttributeInt("idx3");
    this->_Atom1 = atomTable->findEnergyAtomWithCoordinateIndex(i1)->atom();
    this->_Atom2 = atomTable->findEnergyAtomWithCoordinateIndex(i2)->atom();
    this->_Atom3 = atomTable->findEnergyAtomWithCoordinateIndex(i3)->atom();
    this->term.kt = xml->getAttributeDouble("kt");
    this->term.t0 = xml->getAttributeDouble("t0");
    this->term.I1 = xml->getAttributeInt("I1");
    this->term.I2 = xml->getAttributeInt("I2");
    this->term.I3 = xml->getAttributeInt("I3");
}
#endif


void EnergyAngle_O::addTerm(const EnergyAngle& term)
{
    this->_Terms.push_back(term);
}


string EnergyAngle_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyAngle_O,EnergyAngle>(*this);
}



void	EnergyAngle_O::dumpTerms(core::HashTable_sp atomTypes)
{
gctools::Vec0<EnergyAngle>::iterator	eai;
string				as1,as2,as3,as4;
string				str1, str2, str3, str4;
    for (eai=this->_Terms.begin();eai!=this->_Terms.end();eai++)
    {
	as1 = atomLabel(eai->_Atom1);
	as2 = atomLabel(eai->_Atom2);
	as3 = atomLabel(eai->_Atom3);
	if ( as1 < as3 ) {
	    str1 = as1;
	    str2 = as2;
	    str3 = as3;
	} else {
	    str3 = as1;
	    str2 = as2;
	    str1 = as3;
	}
        core::clasp_write_string(fmt::format("TERM 2ANG {:<9} - {:<9} - {:<9} {:8.2f} {:8.2f}\n" , str1.c_str() , str2.c_str() , str3.c_str() , eai->term.kt , core::degrees(eai->term.t0) ));
    }
}





void	EnergyAngle_O::setupHessianPreconditioner(
					chem::NVector_sp nvPosition,
					chem::AbstractLargeSquareMatrix_sp m,
                                        core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  core::T_sp debugInteractions = nil<core::T_O>();
  MAYBE_SETUP_DEBUG_INTERACTIONS(false);
bool		calcForce = true;
bool		calcDiagonalHessian = true;
bool		calcOffDiagonalHessian = true;

//
// Copy from implementAmberFunction::setupHessianPreconditioner

#undef	ANGLE_SET_PARAMETER
#define	ANGLE_SET_PARAMETER(x)	{x=ai->term.x;}
#undef	ANGLE_SET_POSITION
#define	ANGLE_SET_POSITION(x,ii,of)	{x=nvPosition->element(ii+of);}
#undef	ANGLE_ENERGY_ACCUMULATE
#define	ANGLE_ENERGY_ACCUMULATE(e) {}
#undef	ANGLE_FORCE_ACCUMULATE
#define	ANGLE_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
    ENSURE_NOT_NAN(v); \
    m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
    ENSURE_NOT_NAN(v); \
    m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define ANGLE_CALC_FORCE
#define ANGLE_CALC_DIAGONAL_HESSIAN
#define ANGLE_CALC_OFF_DIAGONAL_HESSIAN

 {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include	<cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
	fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0;
	int I1, I2, I3;
	for ( gctools::Vec0<EnergyAngle>::iterator ai=this->_Terms.begin();
		    ai!=this->_Terms.end(); ai++ )
	{
          bool IllegalAngle = false;
#include	<cando/chem/energy_functions/_Angle_termCode.cc>
	}
    }


}





double EnergyAngle_O::evaluateAllComponent( ScoringFunction_sp score,
                                            chem::NVector_sp 	pos,
                                            core::T_sp componentEnergy,
                                            bool 		calcForce,
                                            gc::Nilable<chem::NVector_sp> 	force,
                                            bool		calcDiagonalHessian,
                                            bool		calcOffDiagonalHessian,
                                            gc::Nilable<chem::AbstractLargeSquareMatrix_sp>	hessian,
                                            gc::Nilable<chem::NVector_sp> hdvec,
                                              gc::Nilable<chem::NVector_sp> dvec,
                                              core::T_sp activeAtomMask,
                                              core::T_sp debugInteractions )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  num_real termEnergy = 0.0;
  this->_Evaluations++;
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(("%s {\n") , this->className());
  }

//bool	fail = false;
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());


//
// Copy from implementAmberFunction::evaluateAll
//
#define ANGLE_CALC_FORCE
#define ANGLE_CALC_DIAGONAL_HESSIAN
#define ANGLE_CALC_OFF_DIAGONAL_HESSIAN
#undef ANGLE_SET_PARAMETER
#define ANGLE_SET_PARAMETER(x)	{x=ai->term.x;}
#undef ANGLE_SET_POSITION
#define ANGLE_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef ANGLE_ENERGY_ACCUMULATE
#define ANGLE_ENERGY_ACCUMULATE(e) { ENSURE_NOT_NAN(e); termEnergy += (e); }
#undef	ANGLE_FORCE_ACCUMULATE
#undef	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
#undef	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_FORCE_ACCUMULATE 		ForceAcc
#define	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc



  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0; //,angleScale;
//	double DotAbCb;
    int I1, I2, I3,i;
    gctools::Vec0<EnergyAngle>::iterator ai;
    for ( i=0,ai=this->_Terms.begin();
          ai!=this->_Terms.end(); ai++,i++ ) {
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfTermsToCalculate ) {
          break;
        }
      }
#endif
      bool IllegalAngle = false;
#include	<cando/chem/energy_functions/_Angle_termCode.cc>
                
      if ( IllegalAngle ) {
        ERROR(chem::_sym_LinearAngleError,core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(ai->_Atom1,ai->_Atom2,ai->_Atom3),
                                                                   kw::_sym_coordinates,pos,
                                                                   kw::_sym_indices,core::Cons_O::createList(core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3))));
      }
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
      ai->_calcForce = calcForce;
      ai->_calcDiagonalHessian = calcDiagonalHessian;
      ai->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef	EVAL_SET
#define	EVAL_SET(var,val)	{ ai->eval.var=val;};
#include	<cando/chem/energy_functions/_Angle_debugEvalSet.cc>
#endif //]
      if ( this->_DebugEnergy ) 
      {
        LOG_ENERGY(( "MEISTER angle %d args cando\n") , (i+1) );
        LOG_ENERGY(( "MEISTER angle %d t0 %lf\n") , (i+1) , t0 );
        LOG_ENERGY(( "MEISTER angle %d kt %lf\n") , (i+1) , kt);
        LOG_ENERGY(( "MEISTER angle %d x1 %5.3lf %d\n") , (i+1) , x1 %(I1/3+1) );
        LOG_ENERGY(( "MEISTER angle %d y1 %5.3lf %d\n") , (i+1) , y1 %(I1/3+1) );
        LOG_ENERGY(( "MEISTER angle %d z1 %5.3lf %d\n") , (i+1) , z1 %(I1/3+1) );
        LOG_ENERGY(( "MEISTER angle %d x2 %5.3lf %d\n") , (i+1) , x2 %(I2/3+1) );
        LOG_ENERGY(( "MEISTER angle %d y2 %5.3lf %d\n") , (i+1) , y2 %(I2/3+1) );
        LOG_ENERGY(( "MEISTER angle %d z2 %5.3lf %d\n") , (i+1) , z2 %(I2/3+1) );
        LOG_ENERGY(( "MEISTER angle %d x3 %5.3lf %d\n") , (i+1) , x3 %(I3/3+1) );
        LOG_ENERGY(( "MEISTER angle %d y3 %5.3lf %d\n") , (i+1) , y3 %(I3/3+1) );
        LOG_ENERGY(( "MEISTER angle %d z3 %5.3lf %d\n") , (i+1) , z3 %(I3/3+1) );
        LOG_ENERGY(( "MEISTER angle %d results\n") , (i+1) );
        LOG_ENERGY(( "MEISTER angle %d Energy %lf\n") , (i+1) , Energy );
        if ( calcForce ) 
        {
          LOG_ENERGY(( "MEISTER angle %d fx1 %8.5lf %d\n") , (i+1) , fx1 %(I1/3+1) );
          LOG_ENERGY(( "MEISTER angle %d fy1 %8.5lf %d\n") , (i+1) , fy1 %(I1/3+1) );
          LOG_ENERGY(( "MEISTER angle %d fz1 %8.5lf %d\n") , (i+1) , fz1 %(I1/3+1) );
          LOG_ENERGY(( "MEISTER angle %d fx2 %8.5lf %d\n") , (i+1) , fx2 %(I2/3+1) );
          LOG_ENERGY(( "MEISTER angle %d fy2 %8.5lf %d\n") , (i+1) , fy2 %(I2/3+1) );
          LOG_ENERGY(( "MEISTER angle %d fz2 %8.5lf %d\n") , (i+1) , fz2 %(I2/3+1) );
          LOG_ENERGY(( "MEISTER angle %d fx3 %8.5lf %d\n") , (i+1) , fx3 %(I3/3+1) );
          LOG_ENERGY(( "MEISTER angle %d fy3 %8.5lf %d\n") , (i+1) , fy3 %(I3/3+1) );
          LOG_ENERGY(( "MEISTER angle %d fz3 %8.5lf %d\n") , (i+1) , fz3 %(I3/3+1) );
        }
        LOG_ENERGY(( "MEISTER angle %d stop\n") , (i+1) );
      }
		    /* Add the forces */
      if ( calcForce ) {
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx3>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy3>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz3>10000.0);
      }
    }
  }
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY(("%s }\n") , this->className());
  }
  maybeSetEnergy( componentEnergy, EnergyAngle_O::static_classSymbol(), termEnergy );
  return termEnergy;
}







SYMBOL_EXPORT_SC_(KeywordPkg,angle);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);
SYMBOL_EXPORT_SC_(KeywordPkg,t0);
SYMBOL_EXPORT_SC_(KeywordPkg,angle_deviation);


core::List_sp	EnergyAngle_O::checkForBeyondThresholdInteractionsWithPosition(chem::NVector_sp pos, double threshold)
{
  IMPLEMENT_ME();
#if 0
  ql::list result;
#if 0
  bool	calcForce = false;
  bool	calcDiagonalHessian = false;
  bool	calcOffDiagonalHessian = false;
#endif
  this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//

#undef ANGLE_CALC_FORCE
#undef ANGLE_CALC_DIAGONAL_HESSIAN
#undef ANGLE_CALC_OFF_DIAGONAL_HESSIAN
#undef ANGLE_SET_PARAMETER
#define ANGLE_SET_PARAMETER(x)	{x=ai->term.x;}
#undef ANGLE_SET_POSITION
#define ANGLE_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef ANGLE_ENERGY_ACCUMULATE
#define ANGLE_ENERGY_ACCUMULATE(e) {}
#undef ANGLE_FORCE_ACCUMULATE
#define ANGLE_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

  {
    
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0; //,angleScale;
//	    double DotAbCb;
    int I1, I2, I3,i;
    gctools::Vec0<EnergyAngle>::iterator ai;
    for ( i=0,ai=this->_Terms.begin();
          ai!=this->_Terms.end(); ai++,i++ ) {
      bool IllegalAngle = false;
#include	<cando/chem/energy_functions/_Angle_termCode.cc>
      if ( fabs(AngleDeviation)/t0 > threshold ) {
        chem::Atom_sp a1, a2, a3;
        a1 = (*ai)._Atom1;
        a2 = (*ai)._Atom2;
        a3 = (*ai)._Atom3;
        ql::list one_deviation;
        one_deviation << kw::_sym_angle
                      << kw::_sym_atoms << core::Cons_O::createList(a1,a2,a3)
                      << kw::_sym_r0 << core::clasp_make_double_float(t0)
                      << kw::_sym_angle_deviation << core::clasp_make_double_float(AngleDeviation);
        result << one_deviation.result();
      }
    }
  }
  return result.cons();
#endif
}











void EnergyAngle_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(0.2);
}

void EnergyAngle_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}


SYMBOL_EXPORT_SC_(KeywordPkg,kt);
SYMBOL_EXPORT_SC_(KeywordPkg,t0);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,i3);
SYMBOL_EXPORT_SC_(KeywordPkg,atom1);
SYMBOL_EXPORT_SC_(KeywordPkg,atom2);
SYMBOL_EXPORT_SC_(KeywordPkg,atom3);

CL_DEFMETHOD core::List_sp EnergyAngle_O::extract_vectors_as_alist() const {
  size_t size = this->_Terms.size();
  NVector_sp kt_vec = NVector_O::make(size);
  NVector_sp t0_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i3_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_sp atom1_vec    = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom2_vec    = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom3_vec    = core::SimpleVector_O::make(size);

  for (size_t i=0;i<size;++i) {
    const EnergyAngle& entry = this->_Terms[i];
    (*kt_vec)[i] = entry.term.kt;
    (*t0_vec)[i] = entry.term.t0;
    (*i1_vec)[i] = entry.term.I1;
    (*i2_vec)[i] = entry.term.I2;
    (*i3_vec)[i] = entry.term.I3;
    (*atom1_vec)[i] = entry._Atom1;
    (*atom2_vec)[i] = entry._Atom2;
    (*atom3_vec)[i] = entry._Atom3;
  }
  return core::Cons_O::createList(core::Cons_O::create(kw::_sym_kt,kt_vec),
                                  core::Cons_O::create(kw::_sym_t0,t0_vec),
                                  core::Cons_O::create(kw::_sym_i1,i1_vec),
                                  core::Cons_O::create(kw::_sym_i2,i2_vec),
                                  core::Cons_O::create(kw::_sym_i3,i3_vec),
                                  core::Cons_O::create(kw::_sym_atom1,atom1_vec),
                                  core::Cons_O::create(kw::_sym_atom2,atom2_vec),
                                  core::Cons_O::create(kw::_sym_atom3,atom3_vec));

}

CL_DEFMETHOD void EnergyAngle_O::fill_from_vectors_in_alist(core::List_sp vectors)
{
  NVector_sp kt_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_kt));
  NVector_sp t0_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_t0));
  core::SimpleVector_int32_t_sp i1_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i1));
  core::SimpleVector_int32_t_sp i2_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i2));
  core::SimpleVector_int32_t_sp i3_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i3));
  core::SimpleVector_sp atom1_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom1));
  core::SimpleVector_sp atom2_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom2));
  core::SimpleVector_sp atom3_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom3));
  this->_Terms.resize(kt_vec->length());
  for (size_t i=0, iEnd(kt_vec->length());i<iEnd;++i){
       EnergyAngle& entry = this->_Terms[i];
       entry.term.kt = (*kt_vec)[i];
       entry.term.t0 = (*t0_vec)[i];
       entry.term.I1 = (*i1_vec)[i];
       entry.term.I2 = (*i2_vec)[i];
       entry.term.I3 = (*i3_vec)[i];
       entry._Atom1 = gc::As_unsafe<Atom_sp>((*atom1_vec)[i]);
       entry._Atom2 = gc::As_unsafe<Atom_sp>((*atom2_vec)[i]);
       entry._Atom3 = gc::As_unsafe<Atom_sp>((*atom3_vec)[i]);
  }
}

CL_DEFMETHOD void EnergyAngle_O::addAngleTerm(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, double kt, double t0) {
  EnergyAtom* ea1 = atomTable->getEnergyAtomPointer(a1);
  EnergyAtom* ea2 = atomTable->getEnergyAtomPointer(a2);
  EnergyAtom* ea3 = atomTable->getEnergyAtomPointer(a3);
  EnergyAngle energyAngle(a1,a2,a3,
                          ea1->coordinateIndexTimes3(),
                          ea2->coordinateIndexTimes3(),
                          ea3->coordinateIndexTimes3(),
                          kt,t0);
  this->addTerm(energyAngle);
}




CL_DEFMETHOD
core::List_sp	EnergyAngle_O::lookupAngleTerms(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, core::HashTable_sp atomTypes )
{
  ql::list  result;
  core::T_sp tia1 = atomTable->_AtomTableIndices->gethash(a1);
  core::T_sp tia2 = atomTable->_AtomTableIndices->gethash(a2);
  core::T_sp tia3 = atomTable->_AtomTableIndices->gethash(a3);
  if (!tia1.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a1));
  if (!tia2.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a2));
  if (!tia3.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a3));
  for (auto edi=this->_Terms.begin();edi!=this->_Terms.end();edi++) {
    if ((edi->_Atom1==a1 &&
         edi->_Atom2==a2 &&
         edi->_Atom3==a3)
        || (edi->_Atom3==a1 &&
            edi->_Atom2==a2 &&
            edi->_Atom1==a3)) {
      ql::list oneResult;
      oneResult << INTERN_(kw,type1) << edi->_Atom1->getType(atomTypes)
                << INTERN_(kw,type2) << edi->_Atom2->getType(atomTypes)
                << INTERN_(kw,type3) << edi->_Atom3->getType(atomTypes)
                << INTERN_(kw,kt) << core::clasp_make_double_float(edi->term.kt)
                << INTERN_(kw,t0) << core::clasp_make_double_float(edi->term.t0/0.0174533);
      result << oneResult.result();
    }
  }
  return result.result();
}

EnergyAngle_sp EnergyAngle_O::copyFilter(core::T_sp keepInteraction) {
  EnergyAngle_sp copy = EnergyAngle_O::create();
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    Atom_sp a3 = edi->_Atom3;
    if ( skipInteraction( keepInteraction, EnergyAngle_O::staticClass(), a1, a2, a3 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}

};
