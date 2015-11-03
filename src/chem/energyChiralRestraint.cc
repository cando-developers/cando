       
#define	DEBUG_LEVEL_NONE


#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/profiler.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>



//#define	CHIRAL_RESTRAINT_WEIGHT 1000000.0	// 10000.0 is used in MOE


namespace chem {

    EnergyChiralRestraint::EnergyChiralRestraint()
{
  this->_Atom1 = _Nil<core::T_O>();
  this->_Atom2 = _Nil<core::T_O>();
  this->_Atom3 = _Nil<core::T_O>();
  this->_Atom4 = _Nil<core::T_O>();
}

EnergyChiralRestraint::~EnergyChiralRestraint()
{
}


#ifdef XML_ARCHIVE
    void	EnergyChiralRestraint::archive(core::ArchiveP node)
{
    node->attribute("K",this->term.K);
    node->attribute("CO",this->term.CO);
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
#include <cando/chem/energy_functions/_ChiralRestraint_debugEvalSerialize.cc>
#endif //]
}
#endif



string	EnergyChiralRestraint::description()
{
    stringstream ss;
    ss << "EnergyChiralRestraint[";
    ss << "atoms: ";
    ss << this->_Atom1->description() << "-";
    ss << this->_Atom2->description() << "-";
    ss << this->_Atom3->description() << "-";
    ss << this->_Atom4->description() ;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
    ss << " eval.Energy=" << this->eval.Energy;
#endif
    return ss.str();
}


#if 0
    adapt::QDomNode_sp	EnergyChiralRestraint::asXml()
{
    adapt::QDomNode_sp	node,child;
    Vector3	vdiff;

    node = adapt::QDomNode_O::create(env,"EnergyChiralRestraint");
    node->addAttributeString("atom1Name",this->_Atom1->getName());
    node->addAttributeString("atom2Name",this->_Atom2->getName());
    node->addAttributeString("atom3Name",this->_Atom3->getName());
    node->addAttributeString("atom4Name",this->_Atom4->getName());
    node->addAttributeInt("I1",this->term.I1);
    node->addAttributeInt("I2",this->term.I2);
    node->addAttributeInt("I3",this->term.I3);
    node->addAttributeInt("I4",this->term.I4);
    node->addAttributeDoubleScientific("K",this->term.K);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
    adapt::QDomNode_sp xml = adapt::QDomNode_O::create(env,"Evaluated");
    xml->addAttributeBool("calcForce",this->_calcForce );
    xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
    xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include <_ChiralRestraint_debugEvalXml.cc>
    node->addChild(xml);
#endif
    return node;
}

    void	EnergyChiralRestraint::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
				AtomTable_sp at )
{
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
double	_evaluateEnergyOnly_ChiralRestraint(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3,
		double x4, double y4, double z4,
		double K, double CO )
{
#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of)	{}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	CHIRAL_RESTRAINT_CALC_FORCE	// Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_ChiralRestraint_termCode.cc>

    return Energy;
}



void EnergyChiralRestraint_O::addTerm(const EnergyChiralRestraint& e)
{
    this->_Terms.push_back(e);
}

void EnergyChiralRestraint_O::dumpTerms()
{
    gctools::Vec0<EnergyChiralRestraint>::iterator cri;
    int i;
    string as1, as2, as3, as4;
    for ( i=0,cri=this->_Terms.begin();
	  cri!=this->_Terms.end(); cri++,i++ ) {
	as1 = atomLabel(cri->_Atom1);
	as2 = atomLabel(cri->_Atom2);
	as3 = atomLabel(cri->_Atom3);
	as4 = atomLabel(cri->_Atom4);
	_lisp->print(BF("TERM 7CHIRAL %-9s %-9s %-9s %-9s %8.2lf %8.2lf ; I1=%d I2=%d I3=%d I4=%d")
		     % as1 % as2 % as3 % as4 % cri->term.K % cri->term.CO
		     % cri->term.I1 % cri->term.I2 % cri->term.I3 % cri->term.I4 );
    }
}


string EnergyChiralRestraint_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyChiralRestraint_O,EnergyChiralRestraint>(*this);
}




void	EnergyChiralRestraint_O::setupHessianPreconditioner(
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


#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{x=cri->term.x;}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of) {x=nvPosition->element(ii+of);}
#undef	CHIRAL_RESTRAINT_PHI_SET
#define	CHIRAL_RESTRAINT_PHI_SET(x) {}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define CHIRAL_RESTRAINT_CALC_FORCE
#define CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN

    if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include	<cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,K, CO;
	int	I1, I2, I3, I4, i;
	for ( gctools::Vec0<EnergyChiralRestraint>::iterator cri=this->_Terms.begin();
		    cri!=this->_Terms.end(); cri++ ) {
#include	<cando/chem/energy_functions/_ChiralRestraint_termCode.cc>
	}
    }



}





void	EnergyChiralRestraint_O::evaluateAll(
		chem::NVector_sp 	pos,
		bool 		calcForce,
		gc::Nilable<chem::NVector_sp> 	force,
       		bool		calcDiagonalHessian,
		bool		calcOffDiagonalHessian,
		gc::Nilable<chem::AbstractLargeSquareMatrix_sp>	hessian,
		chem::NVector_sp	hdvec, chem::NVector_sp dvec)
{_G();
    if ( this->_DebugEnergy ) 
    {
	LOG_ENERGY_CLEAR();
	LOG_ENERGY(BF("%s {\n")% this->className());
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


#define CHIRAL_RESTRAINT_CALC_FORCE
#define CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	CHIRAL_RESTRAINT_PHI_SET
#define	CHIRAL_RESTRAINT_PHI_SET(x) {}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) this->_TotalEnergy += (e);
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE 		ForceAcc
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc


    if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,K, CO;
	int	I1, I2, I3, I4, i;
        gctools::Vec0<EnergyChiralRestraint>::iterator cri;
	for ( i=0,cri=this->_Terms.begin();
		    cri!=this->_Terms.end(); cri++,i++ ) {
	    #ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
		if ( this->_Debug_NumberOfChiralRestraintTermsToCalculate > 0 ) {
		    if ( i>= this->_Debug_NumberOfChiralRestraintTermsToCalculate ) {
			break;
		    }
		}
	    #endif

			/* Obtain all the parameters necessary to calculate */
			/* the amber and forces */
#include <cando/chem/energy_functions/_ChiralRestraint_termCode.cc>

	    #if TURN_ENERGY_FUNCTION_DEBUG_ON //[
		cri->_calcForce = calcForce;
		cri->_calcDiagonalHessian = calcDiagonalHessian;
		cri ->_calcOffDiagonalHessian = calcOffDiagonalHessian;
		    // Now write all of the calculated values into the eval structure
		#undef EVAL_SET
		#define	EVAL_SET(var,val)	{ cri->eval.var=val;};
#include	<cando/chem/energy_functions/_ChiralRestraint_debugEvalSet.cc>
	    #endif //]

	    if ( this->_DebugEnergy ) {
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d args cando\n")% (i+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d K %lf\n")% (i+1) % K );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d x1 %5.3lf %d\n")%(i+1) % x1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d y1 %5.3lf %d\n")%(i+1) % y1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d z1 %5.3lf %d\n")%(i+1) % z1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d x2 %5.3lf %d\n")%(i+1) % x2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d y2 %5.3lf %d\n")%(i+1) % y2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d z2 %5.3lf %d\n")%(i+1) % z2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d x3 %5.3lf %d\n")%(i+1) % x3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d y3 %5.3lf %d\n")%(i+1) % y3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d z3 %5.3lf %d\n")%(i+1) % z3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d x4 %5.3lf %d\n")%(i+1) % x4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d y4 %5.3lf %d\n")%(i+1) % y4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d z4 %5.3lf %d\n")%(i+1) % z4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d results\n")% (i+1) );
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d Energy %lf\n")% (i+1) % Energy);
		if ( calcForce ) {
//			LOG_ENERGY(BF( "MEISTER chiralRestraint %d DePhi %lf\n")% (i+1) % DePhi);
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fx1 %8.5lf %d\n")%(i+1) % fx1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fy1 %8.5lf %d\n")%(i+1) % fy1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fz1 %8.5lf %d\n")%(i+1) % fz1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fx2 %8.5lf %d\n")%(i+1) % fx2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fy2 %8.5lf %d\n")%(i+1) % fy2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fz2 %8.5lf %d\n")%(i+1) % fz2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fx3 %8.5lf %d\n")%(i+1) % fx3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fy3 %8.5lf %d\n")%(i+1) % fy3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fz3 %8.5lf %d\n")%(i+1) % fz3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fx4 %8.5lf %d\n")%(i+1) % fx4 % (I4/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fy4 %8.5lf %d\n")%(i+1) % fy4 % (I4/3+1) );
		    LOG_ENERGY(BF( "MEISTER chiralRestraint %d fz4 %8.5lf %d\n")%(i+1) % fz4 % (I4/3+1) );
		}
		LOG_ENERGY(BF( "MEISTER chiralRestraint %d stop\n")% (i+1) );
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
    if ( this->_DebugEnergy ) 
    {
	LOG_ENERGY(BF("%s }\n")% this->className());
    }
}







void	EnergyChiralRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
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

#define CHIRAL_RESTRAINT_CALC_FORCE
#define CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	CHIRAL_RESTRAINT_PHI_SET
#define	CHIRAL_RESTRAINT_PHI_SET(x) {}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


	if ( this->isEnabled() ) {
		_BLOCK_TRACE("ChiralRestraintEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
	    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,K, CO;
	    int	I1, I2, I3, I4, i;
            gctools::Vec0<EnergyChiralRestraint>::iterator cri;
	    for ( i=0,cri=this->_Terms.begin();
			cri!=this->_Terms.end(); cri++,i++ ) {
			  /* Obtain all the parameters necessary to calculate */
			  /* the amber and forces */
#include <cando/chem/energy_functions/_ChiralRestraint_termCode.cc>
		LOG(BF("fx1 = %le") % fx1 );
		LOG(BF("fy1 = %le") % fy1 );
		LOG(BF("fz1 = %le") % fz1 );
		LOG(BF("fx2 = %le") % fx2 );
		LOG(BF("fy2 = %le") % fy2 );
		LOG(BF("fz2 = %le") % fz2 );
		LOG(BF("fx3 = %le") % fx3 );
		LOG(BF("fy3 = %le") % fy3 );
		LOG(BF("fz3 = %le") % fz3 );
		LOG(BF("fx4 = %le") % fx4 );
		LOG(BF("fy4 = %le") % fy4 );
		LOG(BF("fz4 = %le") % fz4 );
		int index = i;
#include <cando/chem/energy_functions/_ChiralRestraint_debugFiniteDifference.cc>

	    }
	}
}



int	EnergyChiralRestraint_O::checkForBeyondThresholdInteractions(
			stringstream& info, chem::NVector_sp pos )
{_G();
int	fails = 0;

    this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------


#undef CHIRAL_RESTRAINT_CALC_FORCE
#undef CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN
#undef CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	CHIRAL_RESTRAINT_PHI_SET
#define	CHIRAL_RESTRAINT_PHI_SET(x) {}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


    if ( this->isEnabled() ) {
	    _BLOCK_TRACE("ChiralRestraintEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,K, CO;
	int	I1, I2, I3, I4, i;
        gctools::Vec0<EnergyChiralRestraint>::iterator cri;
	LOG(BF("Entering checking loop, there are %d terms") % this->_Terms.end()-this->_Terms.begin() );
	for ( i=0,cri=this->_Terms.begin();
		    cri!=this->_Terms.end(); cri++,i++ ) {
		      /* Obtain all the parameters necessary to calculate */
		      /* the amber and forces */
	    LOG(BF("Checking term# %d") % i  );
#include <cando/chem/energy_functions/_ChiralRestraint_termCode.cc>
	    LOG(BF("Status") );
	    if ( ChiralTest>0.0 ) {
		chem::Atom_sp a1, a2, a3, a4;
		a1 = (*cri)._Atom1;
		a2 = (*cri)._Atom2;
		a3 = (*cri)._Atom3;
		a4 = (*cri)._Atom4;
	        LOG(BF("Status") );
		info<< "ChiralRestraintDeviation ";
		info << "value " << ChiralTest << " Atoms(";
		info << a1->getName() << " ";
		info << a2->getName() << " ";
		info << a3->getName() << " ";
		info << a4->getName() << ")";
		info << std::endl;
	        LOG(BF("Info: %s") % info.str().c_str()  );
		this->_BeyondThresholdTerms.push_back(*cri);
	        LOG(BF("Status") );
		fails++;
	    }
	    LOG(BF("Status") );
	}
	LOG(BF("Status") );
    }

    return fails;
}






void EnergyChiralRestraint_O::exposeCando(core::Lisp_sp e)
{
    core::class_<EnergyChiralRestraint_O>()
    ;
}

    void EnergyChiralRestraint_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,EnergyChiralRestraint,"","",_lisp)
    ;
#endif
}

void EnergyChiralRestraint_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(0.2);
}


#ifdef XML_ARCHIVE
void EnergyChiralRestraint_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    IMPLEMENT_ME();
}
#endif

    EXPOSE_CLASS(chem,EnergyChiralRestraint_O);

};




