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

#include <cando/chem/energyAngle.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/profiler.h>
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
    this->term.kt = term->getK2_kCalPerRadianSquared()*scale;
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
    this->_Type1 = ea1->_Atom->getType();
    this->_Type2 = ea2->_Atom->getType();
    this->_Type3 = ea3->_Atom->getType();
}
#endif //]

double	EnergyAngle::getT()
{
    Vector3	pos1, pos2, pos3;
double	angle;
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

//
// Copy this from implementAmberFunction.cc
//
double	_evaluateEnergyOnly_Angle(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3,
		double t0, double kt )
{
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
	double fx1;
	double fy1;
	double fz1;
	double fx2;
	double fy2;
	double fz2;
	double fx3;
	double fy3;
	double fz3;
#endif
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
#include <cando/chem/energy_functions/_Angle_termCode.cc>

    return Energy;
}



void EnergyAngle_O::addTerm(const EnergyAngle& term)
{
    this->_Terms.push_back(term);
}


string EnergyAngle_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyAngle_O,EnergyAngle>(*this);
}



void	EnergyAngle_O::dumpTerms()
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
	_lisp->print(BF("TERM 2ANG %-9s - %-9s - %-9s %8.2lf %8.2lf")
			    % str1.c_str()
			    % str2.c_str()
			    % str3.c_str()	
			    % eai->term.kt
			    % core::degrees(eai->term.t0) );
    }
}





void	EnergyAngle_O::setupHessianPreconditioner(
					chem::NVector_sp nvPosition,
					chem::AbstractLargeSquareMatrix_sp m )
{
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
	m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define ANGLE_CALC_FORCE
#define ANGLE_CALC_DIAGONAL_HESSIAN
#define ANGLE_CALC_OFF_DIAGONAL_HESSIAN

    if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include	<cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
	fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0;
	int I1, I2, I3;
	for ( gctools::Vec0<EnergyAngle>::iterator ai=this->_Terms.begin();
		    ai!=this->_Terms.end(); ai++ )
	{
#include	<cando/chem/energy_functions/_Angle_termCode.cc>
	}
    }


}





double EnergyAngle_O::evaluateAll(chem::NVector_sp 	pos,
                                  bool 		calcForce,
                                  gc::Nilable<chem::NVector_sp> 	force,
                                  bool		calcDiagonalHessian,
                                  bool		calcOffDiagonalHessian,
                                  gc::Nilable<chem::AbstractLargeSquareMatrix_sp>	hessian,
                                  gc::Nilable<chem::NVector_sp> hdvec,
                                  gc::Nilable<chem::NVector_sp> dvec)
{
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(BF("%s {\n")% this->className());
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
#define ANGLE_ENERGY_ACCUMULATE(e) this->_TotalEnergy += (e);
#undef	ANGLE_FORCE_ACCUMULATE
#undef	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
#undef	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_FORCE_ACCUMULATE 		ForceAcc
#define	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc



  if ( this->isEnabled() ) 
  {_BLOCK_TRACE("Angle");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0; //,angleScale;
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
#include	<cando/chem/energy_functions/_Angle_termCode.cc>
                
      if ( IllegalAngle ) {
        ERROR(chem::_sym_LinearAngleError,core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(ai->_Atom1,ai->_Atom2,ai->_Atom3)));
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
        LOG_ENERGY(BF( "MEISTER angle %d args cando\n")% (i+1) );
        LOG_ENERGY(BF( "MEISTER angle %d t0 %lf\n")%(i+1) % t0 );
        LOG_ENERGY(BF( "MEISTER angle %d kt %lf\n")%(i+1) % kt);
        LOG_ENERGY(BF( "MEISTER angle %d x1 %5.3lf %d\n")%(i+1) % x1 %(I1/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d y1 %5.3lf %d\n")%(i+1) % y1 %(I1/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d z1 %5.3lf %d\n")%(i+1) % z1 %(I1/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d x2 %5.3lf %d\n")%(i+1) % x2 %(I2/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d y2 %5.3lf %d\n")%(i+1) % y2 %(I2/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d z2 %5.3lf %d\n")%(i+1) % z2 %(I2/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d x3 %5.3lf %d\n")%(i+1) % x3 %(I3/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d y3 %5.3lf %d\n")%(i+1) % y3 %(I3/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d z3 %5.3lf %d\n")%(i+1) % z3 %(I3/3+1) );
        LOG_ENERGY(BF( "MEISTER angle %d results\n")% (i+1) );
        LOG_ENERGY(BF( "MEISTER angle %d Energy %lf\n")%(i+1) % Energy );
        if ( calcForce ) 
        {
          LOG_ENERGY(BF( "MEISTER angle %d fx1 %8.5lf %d\n")%(i+1) % fx1 %(I1/3+1) );
          LOG_ENERGY(BF( "MEISTER angle %d fy1 %8.5lf %d\n")%(i+1) % fy1 %(I1/3+1) );
          LOG_ENERGY(BF( "MEISTER angle %d fz1 %8.5lf %d\n")%(i+1) % fz1 %(I1/3+1) );
          LOG_ENERGY(BF( "MEISTER angle %d fx2 %8.5lf %d\n")%(i+1) % fx2 %(I2/3+1) );
          LOG_ENERGY(BF( "MEISTER angle %d fy2 %8.5lf %d\n")%(i+1) % fy2 %(I2/3+1) );
          LOG_ENERGY(BF( "MEISTER angle %d fz2 %8.5lf %d\n")%(i+1) % fz2 %(I2/3+1) );
          LOG_ENERGY(BF( "MEISTER angle %d fx3 %8.5lf %d\n")%(i+1) % fx3 %(I3/3+1) );
          LOG_ENERGY(BF( "MEISTER angle %d fy3 %8.5lf %d\n")%(i+1) % fy3 %(I3/3+1) );
          LOG_ENERGY(BF( "MEISTER angle %d fz3 %8.5lf %d\n")%(i+1) % fz3 %(I3/3+1) );
        }
        LOG_ENERGY(BF( "MEISTER angle %d stop\n")% (i+1) );
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
    LOG_ENERGY(BF("%s }\n")% this->className());
  }
  return this->_TotalEnergy;
}







void	EnergyAngle_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
		chem::NVector_sp 	pos)
{
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

	if ( this->isEnabled() ) {
		_BLOCK_TRACE("AngleEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
	    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	    double x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0;
	    int I1, I2, I3,i;
	    gctools::Vec0<EnergyAngle>::iterator ai;
	    for ( i=0,ai=this->_Terms.begin();
			ai!=this->_Terms.end(); ai++,i++ ) {
		LOG(BF("ai->term.kt = %le") % ai->term.kt  );
		LOG(BF("ai->term.t0 = %le") % ai->term.t0  );
		LOG(BF("angleScale = %le") % angleScale  );
		LOG(BF("ai->x1 = %le") % pos->element(ai->term.I1 ) );
		LOG(BF("ai->y1 = %le") % pos->element(ai->term.I1+1 ) );
		LOG(BF("ai->z1 = %le") % pos->element(ai->term.I1+2 ) );
		LOG(BF("ai->x2 = %le") % pos->element(ai->term.I2 ) );
		LOG(BF("ai->y2 = %le") % pos->element(ai->term.I2+1 ) );
		LOG(BF("ai->z2 = %le") % pos->element(ai->term.I2+2 ) );
		LOG(BF("ai->x3 = %le") % pos->element(ai->term.I3 ) );
		LOG(BF("ai->y3 = %le") % pos->element(ai->term.I3+1 ) );
		LOG(BF("ai->z3 = %le") % pos->element(ai->term.I3+2 ) );
#include	<cando/chem/energy_functions/_Angle_termCode.cc>
		LOG(BF("Energy = %le") % Energy  );
		LOG(BF("x1 = %le") % x1  );
		LOG(BF("y1 = %le") % y1  );
		LOG(BF("z1 = %le") % z1  );
		LOG(BF("x2 = %le") % x2  );
		LOG(BF("y2 = %le") % y2  );
		LOG(BF("z2 = %le") % z2  );
		LOG(BF("x3 = %le") % x3  );
		LOG(BF("y3 = %le") % y3  );
		LOG(BF("z3 = %le") % z3  );
		int index = i;
#include <cando/chem/energy_functions/_Angle_debugFiniteDifference.cc>
	    }
	}


}



int	EnergyAngle_O::checkForBeyondThresholdInteractions(
			stringstream& info, chem::NVector_sp pos )
{
int	fails = 0;
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

	if ( this->isEnabled() ) {
		_BLOCK_TRACE("AngleEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Angle_termDeclares.cc>
#pragma clang diagnostic pop
	    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	    double x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0; //,angleScale;
//	    double DotAbCb;
	    int I1, I2, I3,i;
	    gctools::Vec0<EnergyAngle>::iterator ai;
	    for ( i=0,ai=this->_Terms.begin();
			ai!=this->_Terms.end(); ai++,i++ ) {
#include	<cando/chem/energy_functions/_Angle_termCode.cc>
		if ( fabs(AngleDeviation)/t0 > this->_ErrorThreshold ) {
		    chem::Atom_sp a1, a2, a3;
		    a1 = (*ai)._Atom1;
		    a2 = (*ai)._Atom2;
		    a3 = (*ai)._Atom3;
		    info<< "AngleDeviation ";
//		    info<< a1->getAbsoluteIdPath() << " ";
//		    info<< a2->getAbsoluteIdPath() << " ";
//		    info<< a3->getAbsoluteIdPath() << " ";
		    info<< "value " << fabs(AngleDeviation)/t0 << " ";
		    info<<"threshold " << this->_ErrorThreshold << " Atoms(";
		    info << a1->getName() << " ";
		    info << a2->getName() << " ";
		    info << a3->getName() << ") ";
		    info << std::endl;
		    this->_BeyondThresholdTerms.push_back(*ai);
		    fails++;
		}
	    }
	}

    return fails;
}











void EnergyAngle_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(0.2);
}

#ifdef XML_ARCHIVE
void EnergyAngle_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    archiveEnergyComponentTerms<EnergyAngle_O,EnergyAngle>(node,*this);
}
#endif

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
  core::SimpleVectorDouble_sp kt_vec = core::SimpleVectorDouble_O::make(size);
  core::SimpleVectorDouble_sp t0_vec = core::SimpleVectorDouble_O::make(size);
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


};
