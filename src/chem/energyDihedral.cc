#define	DEBUG_LEVEL_FULL


#include "energyDihedral.h"
#include "energyAtomTable.h"
#include "energyFunction.h"
#include "core/profiler.h"
#include "bond.h"
#include "matter.h"
#include "atom.h"
#include "residue.h"
#include "aggregate.h"
#include "nVector.h"
#include "ffBaseDb.h"
#include "ffTypesDb.h"
#include "ffPtorDb.h"
#include "ffItorDb.h"
#include "ffAngleDb.h"
#include "largeSquareMatrix.h"
#include "core/wrappers.h"

namespace chem {

    EnergyDihedral::EnergyDihedral()
{
    this->_Atom1.reset();
    this->_Atom2.reset();
    this->_Atom3.reset();
    this->_Atom4.reset();
}

EnergyDihedral::~EnergyDihedral()
{
}

#ifdef XML_ARCHIVE
    void	EnergyDihedral::archive(core::ArchiveP node)
{
//    node->attribute("_Type1",this->_Type1);
//    node->attribute("_Type2",this->_Type2);
//    node->attribute("_Type3",this->_Type3);
//    node->attribute("_Type4",this->_Type4);
    node->attribute("_Proper",this->_Proper);
    node->attribute("_PhaseRad",this->_PhaseRad);
    node->attribute("sinPhase",this->term.sinPhase);
    node->attribute("cosPhase",this->term.cosPhase);
    node->attribute("V",this->term.V);
    node->attribute("DN",this->term.DN);
    node->attribute("IN",this->term.IN);
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
#include "_Dihedral_debugEvalSerialize.cc"
#endif //]
}
#endif


void	sinNPhiCosNPhi(int n, double* sinNPhi, double* cosNPhi,
			 double sinPhi, double cosPhi )
{
    double	sinNm1Phi, cosNm1Phi;
    if ( n==1 ) {
	*sinNPhi = sinPhi;
	*cosNPhi = cosPhi;
	return;
    };
    sinNPhiCosNPhi(n-1,&sinNm1Phi,&cosNm1Phi,sinPhi,cosPhi);
    *sinNPhi = cosPhi*sinNm1Phi+sinPhi*cosNm1Phi;
    *cosNPhi = cosPhi*cosNm1Phi-sinPhi*sinNm1Phi;
    return;
}


void EnergyDihedral::defineFrom( int n, FFPtor_sp ffterm , EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale)
{
//    this->_Term = term;
    this->_Proper = true;
//    this->_Type1 = term->_T1;
//    this->_Type2 = term->_T2;
//    this->_Type3 = term->_T3;
//    this->_Type4 = term->_T4;
    this->_Atom1 = ea1->atom();
    this->_Atom2 = ea2->atom();
    this->_Atom3 = ea3->atom();
    this->_Atom4 = ea4->atom();
    this->_PhaseRad = ffterm->getPhaseRad(n);
    this->term.sinPhase = sin(this->_PhaseRad);
    this->term.cosPhase = cos(this->_PhaseRad);
    this->term.V = ffterm->getV_kCal(n)*scale;
    this->term.DN = n;
    this->term.IN = n;
    this->term.I1 = ea1->coordinateIndexTimes3();
    this->term.I2 = ea2->coordinateIndexTimes3();
    this->term.I3 = ea3->coordinateIndexTimes3();
    this->term.I4 = ea4->coordinateIndexTimes3();
}

double	EnergyDihedral::getDihedral()
{
    Vector3	pos1, pos2, pos3, pos4;
double	angle;
    pos1 = this->_Atom1->getPosition();
    pos2 = this->_Atom2->getPosition();
    pos3 = this->_Atom3->getPosition();
    pos4 = this->_Atom4->getPosition();
    angle = calculateDihedral( pos1, pos2, pos3, pos4, _lisp);
    return angle;
}

double	EnergyDihedral::getDihedralDeviation()
{
    Vector3	pos1, pos2, pos3, pos4;
double	phi, dev;
    pos1 = this->_Atom1->getPosition();
    pos2 = this->_Atom2->getPosition();
    pos3 = this->_Atom3->getPosition();
    pos4 = this->_Atom4->getPosition();
    phi = calculateDihedral( pos1, pos2, pos3, pos4, _lisp);
    dev = 1.0+cos(this->term.DN*phi-this->_PhaseRad);
    return dev;
}


void EnergyDihedral::defineFrom( int n, FFItor_sp term , EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale)
{
//    this->_Term = term;
    this->_Proper = false;
//    this->_Type1 = term->_T1;
//    this->_Type2 = term->_T2;
//    this->_Type3 = term->_T3;
//    this->_Type4 = term->_T4;
    this->_Atom1 = ea1->atom();
    this->_Atom2 = ea2->atom();
    this->_Atom3 = ea3->atom();
    this->_Atom4 = ea4->atom();
    this->term.DN = n;
    this->term.IN = n;
    this->term.V = term->getV_kCal(n)*scale;
    this->_PhaseRad = term->getPhaseRad(n);
    this->term.cosPhase = cos(this->_PhaseRad);
    this->term.sinPhase = sin(this->_PhaseRad);
    this->term.I1 = ea1->coordinateIndexTimes3();
    this->term.I2 = ea2->coordinateIndexTimes3();
    this->term.I3 = ea3->coordinateIndexTimes3();
    this->term.I4 = ea4->coordinateIndexTimes3();
}

#if 0 //[
void	EnergyDihedral::defineMissingProper( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4 )
{
    this->_Proper = true;
//    this->_Type1 = ea1->_Atom->getType();
//    this->_Type2 = ea2->_Atom->getType();
/    this->_Type3 = ea3->_Atom->getType();
    this->_Type4 = ea4->_Atom->getType();
}
#endif //]


#if 0
addon::QDomNode_sp	EnergyDihedral::asXml(core::Lisp_sp env)
{
addon::QDomNode_sp	node;
Vector3	vdiff;

    node = addon::QDomNode_O::create(env,"EnergyDihedral");
    node->addAttributeBool("proper",this->_Proper );
    node->addAttributeString("atom1Name",this->_Atom1->getName());
    node->addAttributeString("atom2Name",this->_Atom2->getName());
    node->addAttributeString("atom3Name",this->_Atom3->getName());
    node->addAttributeString("atom4Name",this->_Atom4->getName());
//    node->addAttributeString("atom1Type",this->_Type1 );
//    node->addAttributeString("atom2Type",this->_Type2 );
//    node->addAttributeString("atom3Type",this->_Type3 );
//    node->addAttributeString("atom4Type",this->_Type4 );
    node->addAttributeInt("I1",this->term.I1);
    node->addAttributeInt("I2",this->term.I2);
    node->addAttributeInt("I3",this->term.I3);
    node->addAttributeInt("I4",this->term.I4);
    node->addAttributeDoubleScientific("DN",this->term.DN);
    node->addAttributeInt("IN",this->term.IN);
    node->addAttributeDoubleScientific("V",this->term.V);
    node->addAttributeDoubleScientific("PhaseRad",this->_PhaseRad);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
    addon::QDomNode_sp xml = addon::QDomNode_O::create(env,"Evaluated");
    xml->addAttributeBool("calcForce",this->_calcForce );
    xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
    xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include "_Dihedral_debugEvalXml.cc"
    node->addChild(xml);
#endif
    return node;
}

void	EnergyDihedral::parseFromXmlUsingAtomTable(addon::QDomNode_sp	xml,
					AtomTable_sp at )
{
    this->term.I1 = xml->getAttributeInt("I1");
    this->term.I2 = xml->getAttributeInt("I2");
    this->term.I3 = xml->getAttributeInt("I3");
    this->term.I4 = xml->getAttributeInt("I4");
    this->_Atom1 = at->findEnergyAtomWithCoordinateIndex(this->term.I1)->atom();
    this->_Atom2 = at->findEnergyAtomWithCoordinateIndex(this->term.I2)->atom();
    this->_Atom3 = at->findEnergyAtomWithCoordinateIndex(this->term.I3)->atom();
    this->_Atom4 = at->findEnergyAtomWithCoordinateIndex(this->term.I4)->atom();
    this->term.DN = xml->getAttributeDouble("DN");
    this->term.IN = xml->getAttributeInt("IN");
    this->term.V = xml->getAttributeDouble("V");
    this->_PhaseRad = xml->getAttributeDouble("PhaseRad");
    this->_Proper = xml->getAttributeBool("proper");
}
#endif


//
// Copy this from implementAmberFunction.cc
//
double	_evaluateEnergyOnly_Dihedral(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3,
		double x4, double y4, double z4,
		double V, double DN, int IN,
		double cosPhase, double sinPhase )
{
double	EraseLinearDihedral;
double	SinNPhi, CosNPhi;

#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) {}
#undef	DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(i,o,v) {}
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	DIHEDRAL_CALC_FORCE	// Don't calculate FORCE or HESSIAN


#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "_Dihedral_termDeclares.cc"
#pragma clang diagnostic pop
	fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
#include "_Dihedral_termCode.cc"

    return Energy;
}


void EnergyDihedral_O::addTerm(const EnergyDihedral& term)
{
    this->_Terms.push_back(term);
}


string EnergyDihedral_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyDihedral_O,EnergyDihedral>(*this);
}



#if 0
void	EnergyDihedral_O::dumpTerms()
{
    gctools::Vec0<EnergyDihedral>::iterator	edi;
string				as1,as2,as3,as4;
string				str1, str2, str3, str4, type;
    for (edi=this->_Terms.begin();edi!=this->_Terms.end();edi++)
    {
	if ( edi->term.V == 0.0 ) continue;
	if ( edi->_Proper ) {
	    type = "1PROPER  ";
	} else {
	    type = "9IMPROPER";
	}
	as1 = atomLabel(edi->_Atom1);
	as2 = atomLabel(edi->_Atom2);
	as3 = atomLabel(edi->_Atom3);
	as4 = atomLabel(edi->_Atom4);
	if ( edi->_Proper ) {
	    if ( as1 < as4 ) {
		str1 = as1;
		str2 = as2;
		str3 = as3;
		str4 = as4;
	    } else {
		str4 = as1;
		str3 = as2;
		str2 = as3;
		str1 = as4;
	    }
	} else {
	    stringstream ss;
	    ss.str("");
	    ss << as3 << "#Cr";
	    str1 = ss.str();
	    ss.str("");
	    ss << as1 << "#A1";
	    str2 = ss.str();
	    ss.str("");
	    ss << as2 << "#A2";
	    str3 = ss.str();
	    ss.str("");
	    ss << as4 << "#A4";
	    str4 = ss.str();
	}
	_lisp->print(BF("TERM 3DIH %s %-9s - %-9s - %-9s - %-9s %8.2lf %8.2lf %2.0lf")
			% type
			% str1
			% str2
			% str3
			% str4
			% edi->term.V
			% edi->_PhaseRad
			% edi->term.DN);
    }
}
#endif




void	EnergyDihedral_O::setupHessianPreconditioner(
					NVector_sp nvPosition,
					AbstractLargeSquareMatrix_sp m )
{
bool		calcForce = true;
bool		calcDiagonalHessian = true;
bool		calcOffDiagonalHessian = true;

//
// Copy from implementAmberFunction::setupHessianPreconditioner
//
// -----------------------

#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x=di->term.x;}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x=nvPosition->element(ii+of);}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) {}
#undef	DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(i,o,v) {}
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
    if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	#include	"_Dihedral_termDeclares.cc"
#pragma clang diagnostic pop
	fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
	double EraseLinearDihedral;
	int	I1, I2, I3, I4, IN;
	double sinPhase, cosPhase, SinNPhi, CosNPhi;
	for ( gctools::Vec0<EnergyDihedral>::iterator di=this->_Terms.begin();
		    di!=this->_Terms.end(); di++ ) {
	    #include	"_Dihedral_termCode.cc"
	}
    }


}





void	EnergyDihedral_O::evaluateAll(
		NVector_sp 	pos,
		bool 		calcForce,
		NVector_sp 	force,
       		bool		calcDiagonalHessian,
		bool		calcOffDiagonalHessian,
		AbstractLargeSquareMatrix_sp	hessian,
		NVector_sp	hdvec, NVector_sp dvec)
{ _G();
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
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x=di->term.x;}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) this->_TotalEnergy += (e);
#undef	DIHEDRAL_FORCE_ACCUMULATE
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE 		ForceAcc
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc



    if ( this->isEnabled() ) 
    {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	#include "_Dihedral_termDeclares.cc"
#pragma clang diagnostic pop
	fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
	double EraseLinearDihedral;
	int	I1, I2, I3, I4, IN, i;
	double sinPhase, cosPhase, SinNPhi, CosNPhi;
        gctools::Vec0<EnergyDihedral>::iterator di;
	for ( i=0,di=this->_Terms.begin();
		    di!=this->_Terms.end(); di++,i++ ) {
	    #ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
		if ( this->_Debug_NumberOfDihedralTermsToCalculate > 0 ) {
		    if ( i>= this->_Debug_NumberOfDihedralTermsToCalculate ) {
			break;
		    }
		}
	    #endif
	    #include "_Dihedral_termCode.cc"
	    if ( EraseLinearDihedral == 0.0 ) {
	       LOG(BF("Found linear dihedral") );
	       InteractionProblem problem;
	       problem._Atom1 = di->_Atom1;
	       problem._Atom2 = di->_Atom2;
	       problem._Atom3 = di->_Atom3;
	       problem._Atom4 = di->_Atom4;
	       problem._Message = "Found linear dihedral";
	       problem._Type = linearDihedral;
	       LOG(BF("%s") % problem._Message.c_str()  );
	       throw(problem);
	    }
	    #if TURN_ENERGY_FUNCTION_DEBUG_ON //[
		di->_calcForce = calcForce;
		di->_calcDiagonalHessian = calcDiagonalHessian;
		di->_calcOffDiagonalHessian = calcOffDiagonalHessian;
		#undef EVAL_SET
		#define	EVAL_SET(var,val)	{ di->eval.var=val;};
		#include	"_Dihedral_debugEvalSet.cc"
	    #endif //]

	    if ( this->_DebugEnergy ) 
	    {
		LOG_ENERGY(BF( "MEISTER dihedral %d args cando\n")% (i+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d V %lf\n")% (i+1) % V );
		LOG_ENERGY(BF( "MEISTER dihedral %d DN %lf\n")% (i+1) % DN );
		LOG_ENERGY(BF( "MEISTER dihedral %d IN %d\n")% (i+1) % IN );
//		    LOG_ENERGY(BF( "MEISTER dihedral %d phase %lf\n")% (i+1) % phase );
		LOG_ENERGY(BF( "MEISTER dihedral %d sinPhase %lf\n")% (i+1) % sinPhase );
		LOG_ENERGY(BF( "MEISTER dihedral %d cosPhase %lf\n")% (i+1) % cosPhase );
		LOG_ENERGY(BF( "MEISTER dihedral %d x1 %5.3lf %d\n")%(i+1) % x1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d y1 %5.3lf %d\n")%(i+1) % y1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d z1 %5.3lf %d\n")%(i+1) % z1 % (I1/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d x2 %5.3lf %d\n")%(i+1) % x2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d y2 %5.3lf %d\n")%(i+1) % y2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d z2 %5.3lf %d\n")%(i+1) % z2 % (I2/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d x3 %5.3lf %d\n")%(i+1) % x3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d y3 %5.3lf %d\n")%(i+1) % y3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d z3 %5.3lf %d\n")%(i+1) % z3 % (I3/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d x4 %5.3lf %d\n")%(i+1) % x4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d y4 %5.3lf %d\n")%(i+1) % y4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d z4 %5.3lf %d\n")%(i+1) % z4 % (I4/3+1) );
		LOG_ENERGY(BF( "MEISTER dihedral %d results\n")% (i+1) );
	//	LOG_ENERGY(BF( "MEISTER dihedral %d Phi %lf (%lf degrees)\n")% (i+1) %
	//		    	Phi, Phi/0.0174533 );
		LOG_ENERGY(BF( "MEISTER dihedral %d Energy %lf\n")% (i+1) % Energy);
		if ( calcForce ) {
		    LOG_ENERGY(BF( "MEISTER dihedral %d fx1 %8.5lf %d\n")%(i+1) % fx1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fy1 %8.5lf %d\n")%(i+1) % fy1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fz1 %8.5lf %d\n")%(i+1) % fz1 % (I1/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fx2 %8.5lf %d\n")%(i+1) % fx2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fy2 %8.5lf %d\n")%(i+1) % fy2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fz2 %8.5lf %d\n")%(i+1) % fz2 % (I2/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fx3 %8.5lf %d\n")%(i+1) % fx3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fy3 %8.5lf %d\n")%(i+1) % fy3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fz3 %8.5lf %d\n")%(i+1) % fz3 % (I3/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fx4 %8.5lf %d\n")%(i+1) % fx4 % (I4/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fy4 %8.5lf %d\n")%(i+1) % fy4 % (I4/3+1) );
		    LOG_ENERGY(BF( "MEISTER dihedral %d fz4 %8.5lf %d\n")%(i+1) % fz4 % (I4/3+1) );
		}
		LOG_ENERGY(BF( "MEISTER dihedral %d stop\n")% (i+1) );
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







void	EnergyDihedral_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
		NVector_sp 	pos)
{_OF();
int	fails = 0;
bool	calcForce = true;
bool	calcDiagonalHessian = true;
bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm
//
//------------------
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x=di->term.x;}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) {}
#undef	DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(i,o,v) {}
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}



	if ( this->isEnabled() ) {
		_BLOCK_TRACE("DihedralEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	    #include "_Dihedral_termDeclares.cc"
#pragma clang diagnostic pop
	fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
#if !USE_EXPLICIT_DECLARES
	double dhx1x1 = 0.0;
	double ohx1y1 = 0.0;
	double ohx1z1 = 0.0;
	double ohx1x2 = 0.0;
	double ohx1y2 = 0.0;
	double ohx1z2 = 0.0;
	double ohx1x3 = 0.0;
	double ohx1y3 = 0.0;
	double ohx1z3 = 0.0;
	double ohx1x4 = 0.0;
	double ohx1y4 = 0.0;
	double ohx1z4 = 0.0;
	double dhy1y1 = 0.0;
	double ohy1z1 = 0.0;
	double ohy1x2 = 0.0;
	double ohy1y2 = 0.0;
	double ohy1z2 = 0.0;
	double ohy1x3 = 0.0;
	double ohy1y3 = 0.0;
	double ohy1z3 = 0.0;
	double ohy1x4 = 0.0;
	double ohy1y4 = 0.0;
	double ohy1z4 = 0.0;
	double dhz1z1 = 0.0;
	double ohz1x2 = 0.0;
	double ohz1y2 = 0.0;
	double ohz1z2 = 0.0;
	double ohz1x3 = 0.0;
	double ohz1y3 = 0.0;
	double ohz1z3 = 0.0;
	double ohz1x4 = 0.0;
	double ohz1y4 = 0.0;
	double ohz1z4 = 0.0;
	double dhx2x2 = 0.0;
	double ohx2y2 = 0.0;
	double ohx2z2 = 0.0;
	double ohx2x3 = 0.0;
	double ohx2y3 = 0.0;
	double ohx2z3 = 0.0;
	double ohx2x4 = 0.0;
	double ohx2y4 = 0.0;
	double ohx2z4 = 0.0;
	double dhy2y2 = 0.0;
	double ohy2z2 = 0.0;
	double ohy2x3 = 0.0;
	double ohy2y3 = 0.0;
	double ohy2z3 = 0.0;
	double ohy2x4 = 0.0;
	double ohy2y4 = 0.0;
	double ohy2z4 = 0.0;
	double dhz2z2 = 0.0;
	double ohz2x3 = 0.0;
	double ohz2y3 = 0.0;
	double ohz2z3 = 0.0;
	double ohz2x4 = 0.0;
	double ohz2y4 = 0.0;
	double ohz2z4 = 0.0;
	double dhx3x3 = 0.0;
	double ohx3y3 = 0.0;
	double ohx3z3 = 0.0;
	double ohx3x4 = 0.0;
	double ohx3y4 = 0.0;
	double ohx3z4 = 0.0;
	double dhy3y3 = 0.0;
	double ohy3z3 = 0.0;
	double ohy3x4 = 0.0;
	double ohy3y4 = 0.0;
	double ohy3z4 = 0.0;
	double dhz3z3 = 0.0;
	double ohz3x4 = 0.0;
	double ohz3y4 = 0.0;
	double ohz3z4 = 0.0;
	double dhx4x4 = 0.0;
	double ohx4y4 = 0.0;
	double ohx4z4 = 0.0;
	double dhy4y4 = 0.0;
	double ohy4z4 = 0.0;
	double dhz4z4 = 0.0;
#endif
	    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
	    double EraseLinearDihedral;
	    int	I1, I2, I3, I4, IN, i;
	    double sinPhase, cosPhase, SinNPhi, CosNPhi;
            gctools::Vec0<EnergyDihedral>::iterator di;
	    for ( i=0,di=this->_Terms.begin();
			di!=this->_Terms.end(); di++,i++ ) {
		#include "_Dihedral_termCode.cc"
		int index = i;
		#include "_Dihedral_debugFiniteDifference.cc"
	    }
	}

}



int	EnergyDihedral_O::checkForBeyondThresholdInteractions(
			stringstream& info, NVector_sp pos )
{_OF();
int	fails = 0;

    this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------

#undef DIHEDRAL_CALC_FORCE
#undef DIHEDRAL_CALC_DIAGONAL_HESSIAN
#undef DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x=di->term.x;}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) {}
#undef	DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(i,o,v) {}
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


	if ( this->isEnabled() ) {
		_BLOCK_TRACE("DihedralEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
	    #include "_Dihedral_termDeclares.cc"
#pragma clang diagnostic pop
	fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
	fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
	fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
	fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
	    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
	    double EraseLinearDihedral;
	    int	I1, I2, I3, I4, IN, i;
	    double sinPhase, cosPhase, SinNPhi, CosNPhi;
            gctools::Vec0<EnergyDihedral>::iterator di;
	    for ( i=0,di=this->_Terms.begin();
			di!=this->_Terms.end(); di++,i++ ) {
		#include "_Dihedral_termCode.cc"
		if ( fabs(DihedralDeviation)>this->_ErrorThreshold) {
		    Atom_sp a1, a2, a3, a4;
		    a1 = (*di)._Atom1;
		    a2 = (*di)._Atom2;
		    a3 = (*di)._Atom3;
		    a4 = (*di)._Atom4;
		    info<< "DihedralDeviation ";
//		    info<< a1->getAbsoluteIdPath() << " ";
//		    info<< a2->getAbsoluteIdPath() << " ";
//		    info<< a3->getAbsoluteIdPath() << " ";
//		    info<< a4->getAbsoluteIdPath() << " ";
		    info<< "value " << fabs(DihedralDeviation) << " ";
		    info<<"threshold " << this->_ErrorThreshold;
		    info << a1->getName() << " ";
		    info << a2->getName() << " ";
		    info << a3->getName() << " ";
		    info << a4->getName() << " ";
		    info << std::endl;
		    EnergyDihedral ed = *di;
		    ed._CalculatedDihedralDeviation = DihedralDeviation;
		    this->_BeyondThresholdTerms.push_back(ed);
		    fails++;
		}
	    }
	}


    return fails;
}





void EnergyDihedral_O::exposeCando(core::Lisp_sp e)
{
    core::class_<EnergyDihedral_O>()
    ;
}

void EnergyDihedral_O::exposePython(core::Lisp_sp lisp)
{_G();
#ifdef USEBOOSTPYTHON
    PYTHON_CLASS(ChemPkg,EnergyDihedral,"","",_lisp)
	;
#endif
}


void EnergyDihedral_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(3.0);
}

#ifdef XML_ARCHIVE
void EnergyDihedral_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    archiveEnergyComponentTerms<EnergyDihedral_O,EnergyDihedral>(node,*this);
}
#endif


EXPOSE_CLASS(chem,EnergyDihedral_O);
};







