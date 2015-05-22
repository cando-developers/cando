       
#define	DEBUG_LEVEL_NONE


#include "largeSquareMatrix.h"
#include "energyImproperRestraint.h"
#include "energyAtomTable.h"
#include "energyFunction.h"
#include "core/profiler.h"
#include "bond.h"
#include "addon/vector3.h"
#include "matter.h"
#include "atom.h"
#include "residue.h"
#include "aggregate.h"
#include "nVector.h"
#include "ffBaseDb.h"
#include "ffTypesDb.h"
#include "core/wrappers.h"


namespace chem {


    EnergyImproperRestraint::EnergyImproperRestraint()
{
    this->_Atom1.reset();
    this->_Atom2.reset();
    this->_Atom3.reset();
    this->_Atom4.reset();
}


EnergyImproperRestraint::~EnergyImproperRestraint()
{
}

#ifdef XML_ARCHIVE
    void	EnergyImproperRestraint::archive(core::ArchiveP node)
{
    node->attribute("_AboveThreshold",this->_AboveThreshold);
    node->attribute("_AboveThreshold_Phi",this->_AboveThreshold_Phi);
    node->attribute("U",this->term.U);
    node->attribute("L",this->term.L);
    node->attribute("K",this->term.K);
    node->attribute("I1",this->term.I1);
    node->attribute("I2",this->term.I2);
    node->attribute("I3",this->term.I3);
    node->attribute("I4",this->term.I4);
    node->attribute("a1",this->_Atom1);
    node->attribute("a2",this->_Atom2);
    node->attribute("a3",this->_Atom3);
    node->attribute("a4",this->_Atom4);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
    node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
    node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
    node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include "_ImproperRestraint_debugEvalSerialize.cc"
#endif //]
}
#endif



double	EnergyImproperRestraint::getAngle()
{
    Vector3	pos1, pos2, pos3, pos4;
    pos1 = this->_Atom1->getPosition();
    pos2 = this->_Atom2->getPosition();
    pos3 = this->_Atom3->getPosition();
    pos4 = this->_Atom4->getPosition();
    return calculateDihedral(pos1,pos2,pos3,pos4,_lisp);
}


#if 0
    addon::QDomNode_sp	EnergyImproperRestraint::asXml(core::Lisp_sp env)
{
    addon::QDomNode_sp	node,child;
    Vector3	vdiff;

    node = addon::QDomNode_O::create(env,"EnergyImproperRestraint");
    node->addAttributeString("atom1Name",this->_Atom1->getName());
    node->addAttributeString("atom2Name",this->_Atom2->getName());
    node->addAttributeString("atom3Name",this->_Atom3->getName());
    node->addAttributeString("atom4Name",this->_Atom4->getName());
    node->addAttributeInt("I1",this->term.I1);
    node->addAttributeInt("I2",this->term.I2);
    node->addAttributeInt("I3",this->term.I3);
    node->addAttributeInt("I4",this->term.I4);
    node->addAttributeDouble("UDeg",this->term.U/0.0174533,6,2);
    node->addAttributeDouble("LDeg",this->term.L/0.0174533,6,2);
    node->addAttributeDoubleScientific("U",this->term.U);
    node->addAttributeDoubleScientific("L",this->term.L);
    node->addAttributeDoubleScientific("K",this->term.K);
    if ( this->_AboveThreshold ) {
	child = addon::QDomNode_O::create(env,"AboveThreshold");
	child->addAttributeDoubleScientific("Phi",this->_AboveThreshold_Phi );
	child->addAttributeDouble("PhiDeg",this->_AboveThreshold_Phi/0.0174533,6,2 );
	stringstream ss;
	ss << std::endl << "This improper has the value("<<this->_AboveThreshold_Phi/0.0174533<<")"
	    << " outside of bounds ["<<this->term.U/0.0174533<<","
           <<this->term.L/0.0174533<<"]" << std::endl;
	child->setCharacters(ss.str());
	node->addChild(child);
    }
#if TURN_ENERGY_FUNCTION_DEBUG_ON
    addon::QDomNode_sp xml = addon::QDomNode_O::create(env,"Evaluated");
    xml->addAttributeBool("calcForce",this->_calcForce );
    xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
    xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include "_ImproperRestraint_debugEvalXml.cc"
    node->addChild(xml);
#endif
    return node;
}

    void	EnergyImproperRestraint::parseFromXmlUsingAtomTable(addon::QDomNode_sp	xml,
					AtomTable_sp at)
{
    this->term.U = xml->getAttributeDouble("U");
    this->term.L = xml->getAttributeDouble("L");
    this->term.K = xml->getAttributeDouble("K");
    this->term.I1 = xml->getAttributeInt("I1");
    this->term.I2 = xml->getAttributeInt("I2");
    this->term.I3 = xml->getAttributeInt("I3");
    this->term.I4 = xml->getAttributeInt("I4");
    this->_Atom1 = at->findEnergyAtomWithCoordinateIndex(this->term.I1)->atom();
    this->_Atom2 = at->findEnergyAtomWithCoordinateIndex(this->term.I2)->atom();
    this->_Atom3 = at->findEnergyAtomWithCoordinateIndex(this->term.I3)->atom();
    this->_Atom4 = at->findEnergyAtomWithCoordinateIndex(this->term.I4)->atom();
}

#endif




//
// Copy this from implementAmberFunction.cc
//
double	_evaluateEnergyOnly_ImproperRestraint(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3,
		double x4, double y4, double z4,
		double K, double L, double U )
{
double	EraseLinearDihedral;
double	UShift, PhiShift;
bool	RestraintActive;

#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of)	{}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x)	{}
#undef	IMPROPER_RESTRAINT_CALC_FORCE	// Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "_ImproperRestraint_termDeclares.cc"
#pragma clang diagnostic pop
#include "_ImproperRestraint_termCode.cc"

    return Energy;
}


void EnergyImproperRestraint_O::addTerm(const EnergyImproperRestraint& e)
{
    this->_Terms.push_back(e);
}



void	EnergyImproperRestraint_O::dumpTerms()
{
}


string EnergyImproperRestraint_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyImproperRestraint_O,EnergyImproperRestraint>(*this);
}




void	EnergyImproperRestraint_O::setupHessianPreconditioner(
					chem::NVector_sp nvPosition,
					chem::AbstractLargeSquareMatrix_sp m )
{
bool		calcForce = true;
bool		calcDiagonalHessian = true;
bool		calcOffDiagonalHessian = true;

//
// Copy from implementAmberFunction::setupHessianPreconditioner
//
// -----------------------

#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{x=iri->term.x;}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of) {x=nvPosition->element(ii+of);}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x) {}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define IMPROPER_RESTRAINT_CALC_FORCE
#define IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN

    if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	#include	"_ImproperRestraint_termDeclares.cc"
#pragma clang diagnostic pop
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
	double K,U,L;
	double EraseLinearDihedral, UShift, PhiShift;
	bool RestraintActive;
	int	I1, I2, I3, I4;
	for ( gctools::Vec0<EnergyImproperRestraint>::iterator iri=this->_Terms.begin();
		    iri!=this->_Terms.end(); iri++ ) {
	    #include	"_ImproperRestraint_termCode.cc"
	}
    }
}





void	EnergyImproperRestraint_O::evaluateAll(
		chem::NVector_sp 	pos,
		bool 		calcForce,
		chem::NVector_sp 	force,
       		bool		calcDiagonalHessian,
		bool		calcOffDiagonalHessian,
		chem::AbstractLargeSquareMatrix_sp	hessian,
		chem::NVector_sp	hdvec, chem::NVector_sp dvec)
{_G();
    if ( this->_DebugEnergy ) 
    {
	LOG_ENERGY_CLEAR();
	LOG_ENERGY(BF("%s {")% this->className());
    }

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
// -----------------------
#define IMPROPER_RESTRAINT_CALC_FORCE
#define IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{x = iri->term.x;}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x) {}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) this->_TotalEnergy += (e);
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE 		ForceAcc
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc

    if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	#include "_ImproperRestraint_termDeclares.cc"
#pragma clang diagnostic pop
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
	double K,U,L;
	double EraseLinearDihedral, UShift, PhiShift;
	bool RestraintActive;
	int	I1, I2, I3, I4,i;
        gctools::Vec0<EnergyImproperRestraint>::iterator iri;
	for ( i=0,iri =this->_Terms.begin();
		    iri!=this->_Terms.end(); iri++,i++ ) {
	    #ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
		if ( this->_Debug_NumberOfImproperRestraintTermsToCalculate > 0 ) {
		    if ( i>= this->_Debug_NumberOfImproperRestraintTermsToCalculate ) {
			break;
		    }
		}
	    #endif
	    #include "_ImproperRestraint_termCode.cc"
	    if ( EraseLinearDihedral == 0.0 ) {
	       InteractionProblem problem;
	       problem._Atom1 = iri->_Atom1;
	       problem._Atom2 = iri->_Atom2;
	       problem._Atom3 = iri->_Atom3;
	       problem._Atom4 = iri->_Atom4;
	       problem._Message = "Found linear improper restraint";
	       problem._Type = linearImproperRestraint;
	       LOG(BF("%s") % problem._Message.c_str()  );
	       throw(problem);
	    }
	    #if TURN_ENERGY_FUNCTION_DEBUG_ON //[
		iri->_calcForce = calcForce;
		iri->_calcDiagonalHessian = calcDiagonalHessian;
		iri->_calcOffDiagonalHessian = calcOffDiagonalHessian;
		#undef EVAL_SET
		#define EVAL_SET(var,val)	{ iri->eval.var=val;};
		#include "_ImproperRestraint_debugEvalSet.cc"
	    #endif //]

	    if ( this->_DebugEnergy ) 
	    {
		LOG_ENERGY(BF( "MEISTER improperRestraint %d args cando\n")% (i+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d K %lf\n")% (i+1) % K );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d U %lf\n")% (i+1) % U );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d L %lf\n")% (i+1) % L );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d x1 %5.3lf %d\n")%(i+1) % x1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d y1 %5.3lf %d\n")%(i+1) % y1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d z1 %5.3lf %d\n")%(i+1) % z1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d x2 %5.3lf %d\n")%(i+1) % x2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d y2 %5.3lf %d\n")%(i+1) % y2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d z2 %5.3lf %d\n")%(i+1) % z2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d x3 %5.3lf %d\n")%(i+1) % x3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d y3 %5.3lf %d\n")%(i+1) % y3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d z3 %5.3lf %d\n")%(i+1) % z3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d x4 %5.3lf %d\n")%(i+1) % x4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d y4 %5.3lf %d\n")%(i+1) % y4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d z4 %5.3lf %d\n")%(i+1) % z4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d results\n")% (i+1) );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d CosPhi %lf\n")% (i+1) % CosPhi );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d SinPhi %lf\n")% (i+1) % SinPhi );
		if ( CosPhi>0.1 ) {
		    Phi = asin(SinPhi);
		} else {
		    Phi = acos(CosPhi)*SIGN(SinPhi);
		}
		LOG_ENERGY(BF( "MEISTER improperRestraint %d Phi %lf\n")% (i+1) % Phi );
		LOG_ENERGY(BF( "MEISTER improperRestraint %d Energy %lf\n")% (i+1) % Energy);
		if ( calcForce ) 
		{
//			LOG_ENERGY(BF( "MEISTER improperRestraint %d DePhi %lf\n")% (i+1) % DePhi);
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fx1 %8.5lf %d\n")%(i+1) % fx1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fy1 %8.5lf %d\n")%(i+1) % fy1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fz1 %8.5lf %d\n")%(i+1) % fz1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fx2 %8.5lf %d\n")%(i+1) % fx2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fy2 %8.5lf %d\n")%(i+1) % fy2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fz2 %8.5lf %d\n")%(i+1) % fz2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fx3 %8.5lf %d\n")%(i+1) % fx3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fy3 %8.5lf %d\n")%(i+1) % fy3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fz3 %8.5lf %d\n")%(i+1) % fz3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fx4 %8.5lf %d\n")%(i+1) % fx4 % (I4/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fy4 %8.5lf %d\n")%(i+1) % fy4 % (I4/3+1) );
		    LOG_ENERGY(BF( "MEISTER improperRestraint %d fz4 %8.5lf %d\n")%(i+1) % fz4 % (I4/3+1) );
		}
		LOG_ENERGY(BF( "MEISTER improperRestraint %d stop\n")% (i+1) );
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
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx4>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy4>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz4>10000.0);
	    }
	}
    }
    LOG_ENERGY(BF( "ImproperRestraint energy = %lf\n")% (double)(this->_TotalEnergy) );
    if ( this->_DebugEnergy ) 
    {
	LOG_ENERGY(BF("%s }\n")% this->className());
    }
}







void	EnergyImproperRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
		chem::NVector_sp 	pos)
{
int	fails = 0;
bool	calcForce = true;
bool	calcDiagonalHessian = true;
bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//------------------
#define IMPROPER_RESTRAINT_CALC_FORCE
#define IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{x = iri->term.x;}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of) {x = pos->element(ii+of);}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x) {}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


	if ( this->isEnabled() ) {
		_BLOCK_TRACE("ImproperRestraintEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	    #include "_ImproperRestraint_termDeclares.cc"
#pragma clang diagnostic pop
	    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
	    double K,U,L;
	    double EraseLinearDihedral, UShift, PhiShift;
	    bool RestraintActive;
	    int	I1, I2, I3, I4,i;
            gctools::Vec0<EnergyImproperRestraint>::iterator iri;
	    for ( i=0,iri =this->_Terms.begin();
			iri!=this->_Terms.end(); iri++,i++ ) {
		#include "_ImproperRestraint_termCode.cc"
		int index = i;
		#include "_ImproperRestraint_debugFiniteDifference.cc"
	    }
	}


}



int	EnergyImproperRestraint_O::checkForBeyondThresholdInteractions(
			stringstream& info, chem::NVector_sp pos )
{
int	fails = 0;

    this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------
#undef IMPROPER_RESTRAINT_CALC_FORCE
#undef IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN
#undef IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{x = iri->term.x;}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of) {x = pos->element(ii+of);}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x) {}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


	if ( this->isEnabled() ) {
		_BLOCK_TRACE("ImproperRestraintEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	    #include "_ImproperRestraint_termDeclares.cc"
#pragma clang diagnostic pop
	    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
	    double K,U,L;
	    double EraseLinearDihedral, UShift, PhiShift;
	    bool RestraintActive;
	    int	I1, I2, I3, I4,i;
            gctools::Vec0<EnergyImproperRestraint>::iterator iri;
	    for ( i=0,iri =this->_Terms.begin();
			iri!=this->_Terms.end(); iri++,i++ ) {
		#include "_ImproperRestraint_termCode.cc"
		if ( RestraintActive ) {
		    chem::Atom_sp a1, a2, a3, a4;
		    a1 = (*iri)._Atom1;
		    a2 = (*iri)._Atom2;
		    a3 = (*iri)._Atom3;
		    a4 = (*iri)._Atom4;
		    info<< "ImproperRestraintDeviation ";
		    info<< "Phi[degrees](" << Phi/0.0174533 << ") ";
		    info<< "U[degrees](" << U/0.0174533 << ") ";
		    info<< "L[degrees](" << L/0.0174533 << ") ";
		    info << a1->description() << " ";
		    info << a2->description() << " ";
		    info << a3->description() << " ";
		    info << a4->description() << " ";
		    info << std::endl;
		    this->_BeyondThresholdTerms.push_back(*iri);
		    fails++;
		}
	    }
	}


    return fails;
}







void EnergyImproperRestraint_O::exposeCando(core::Lisp_sp e)
{
    core::class_<EnergyImproperRestraint_O>()
    ;
}

    void EnergyImproperRestraint_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,EnergyImproperRestraint,"","",_lisp)
	;
#endif
}


void EnergyImproperRestraint_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(0.001);
}

#ifdef XML_ARCHIVE
void EnergyImproperRestraint_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    IMPLEMENT_ME();
}
#endif


EXPOSE_CLASS(chem,EnergyImproperRestraint_O);
};










