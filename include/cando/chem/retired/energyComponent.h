#define	DEBUG_LEVEL_FULL

//
// (C) 2004 Christian E. Schafmeister
//

/*
 *	energyComponent.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyComponent_H 
#define	EnergyComponent_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
//#include "atom.h"
//#include "bond.h"
//#include "largeSquareMatrix.h"
//#include "ffBaseDb.h"
//#include "ffTypesDb.h"
//#include "ffStretchDb.h"
//#include "ffAngleDb.h"
#include <cando/geom/quickDom.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/nVector.h>
#include <cando/chem/residue.h>

		//
		// SCALE FACTOR FOR ELECTROSTIC INTERACTIONS???
		//
const double ELECTROSTATIC_MODIFIER = (93865.10/3.53*93865.10/93958.78);

#define	DECLARE_FLOAT(x) double x = 0.0

	// If we declare mathematica generated term variables
	// on the fly then set ON_THE_FLY_DECLARES
	// and variable declarations won't be included
	//
	// This should avoid a lot of variables being declared but
	// not used and eliminate a lot of compiler warnings
	//
#define	USE_EXPLICIT_DECLARES 1



#ifndef	delta2
#define		delta2	0.00001
#endif
#include "geom/quickDom.fwd.h"// energyComponent.h wants QDomNode needs quickDom.fwd.h

namespace       chem {


SMART(QDomNode);
SMART(AbstractLargeSquareMatrix);



#if TURN_ENERGY_FUNCTION_DEBUG_ON == 1
#define	LOG_ENERGY_CLEAR() {this->_DebugLog.str("");};
#define LOG_ENERGY(s) {this->_DebugLog << __FILE__ << ":" << __LINE__ << " " << s;}
#else
#define	LOG_ENERGY_CLEAR() {}
#define LOG_ENERGY(s) {}
#endif

/*! Store a pointer to an Atom and an index into the coordinate vector array
 * which stores coordinates in a 1D array (x1,y1,z1,x2,y2,z2,x3,...,xN,yN,zN)
 */
#if 0
template <class TermClass>
class	Dumb_EnergyComponent {

public:

public:
public:

	typedef	typename std::vector<TermClass>::iterator	viter;

virtual	void	initialize() 
	{
	    this->customInitialize();
	};
virtual	void	customInitialize();


	int		numberOfTerms() { return this->_Terms.size();};
	viter 	begin_Terms() {
			return this->_Terms.begin(); };
	viter	end_Terms() {
			return this->_Terms.end(); };

	viter	begin_BeyondThresholdTerms() {
			return this->_BeyondThresholdTerms.begin(); };
	viter	end_BeyondThresholdTerms() {
			return this->_BeyondThresholdTerms.end(); };

//	int	getNumberOfMissingTerms();

	void	addTerm(TermClass& term) 
	{ _G();
	    this->_Terms.push_back(term); 
	};

	void	asXmlFillEnergy(geom::QDomNode_sp node) {
		stringstream ss;
		ss << this->getName()<<"TotalEnergy";
		node->addAttributeDoubleScientific(ss.str(),this->_TotalEnergy);
	};

	void	parseXmlEnergy(geom::QDomNode_sp node) {
		stringstream ss;
		ss << this->getName()<<"TotalEnergy";
		this->_TotalEnergy = node->getAttributeDouble(ss.str());
	};


geom::QDomNode_sp	asXmlLimitedTerms(core::Lisp_sp lisp, int numTerms)
{
geom::QDomNode_sp tnode, child;
int		t;
    viter si;
    tnode = geom::QDomNode_O::create(this->getName());
    tnode->addAttributeInt("numberOfTerms",this->_Terms.size());
    for (t=0,si=this->_Terms.begin();
	    si!=this->_Terms.end();si++,t++){
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
	if ( this->_Debug_NumberOfTermsToCalculate> 0 ) {
	    if ( t>= this->_Debug_NumberOfTermsToCalculate ) break;
	}
#endif
	child = (*si).asXml(lisp);
	tnode->addChild(child);

    }
    return tnode;
}



void	parseFromXmlUsingAtomTable(geom::QDomNode_sp node, AtomTable_sp atomTable )
{_G();
VectorQDomNodes	children;
VectorQDomNodes::iterator	child;
    this->_Terms.clear();
    children = node->getChildren();
    for ( child=children.begin(); child!=children.end(); child++ ) {
	TermClass st;
	if ( (*child)->isNamed("Term") ) {
	    st.parseFromXmlUsingAtomTable(*child, atomTable );
	    this->_Terms.push_back(st);
	} else {
	    THROW(core::LispError_O::create(BF("Unknown xml file(%s) line(%d)") % (*child)->getFileName() % (*child)->getLineNumber() ,_lisp));
	}
    }
}

geom::QDomNode_sp	asXmlMissingTerms();



virtual void setupHessianPreconditioner(NVector_sp nvPosition,
					AbstractLargeSquareMatrix_sp m ) {_OF();SUBCLASS_MUST_IMPLEMENT();};

virtual	void evaluateAll ( NVector_sp 	pos,
		bool 		calcForce,
		NVector_sp 	force,
       		bool		calcDiagonalHessian,
		bool		calcOffDiagonalHessian,
		AbstractLargeSquareMatrix_sp	hessian,
		NVector_sp	hdvec,
		NVector_sp dvec);

	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
				NVector_sp pos );

	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

	string	summarizeEnergyAsString() {
			stringstream	ss;
			    ss.str("");
			    ss<< this->getName();
			    ss << " energy: " << this->_TotalEnergy << std::endl;
			    return ss.str();
	};
	string	enabledAsString() {
	    stringstream	ss;
		ss.str("");
		if ( this->isEnabled() ) {
		    ss<<this->getName();
	    #ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
		    if ( this->_Debug_NumberOfTermsToCalculate < 0 ) {
			ss<<"(debug_All)";
		    } else {
			ss<<"(debug_";
			ss<<this->_Debug_NumberOfTermsToCalculate<<")";
		    }
	    #endif
		    ss<<" ";
		}
		return ss.str();
	};

	string debugLogAsString()
	{
#if TURN_ENERGY_FUNCTION_DEBUG_ON == 1
	    return this->_DebugLog.str();
#else
	    return this->getName()+"-no_log";
#endif
	}


	Dumb_EnergyComponent() 
	{
	    REF_CREATE(__FUNCTION__);
	};
virtual	~Dumb_EnergyComponent() 
	{
	    REF_DELETE(__FUNCTION__);
	};


};
#endif

class EnergyTerm 
{
private:
public:
    virtual string	className() {THROW_HARD_ERROR(BF("Subclass must implement"));};

    EnergyTerm() {};
};


class BeyondThresholdEnergyTerm
{
    public:
};



inline string	atomLabel(Atom_sp a)
{
Residue_sp	res;
stringstream	sstr;
    res= a->getResidueContainedBy().as<Residue_O>();
    sstr.str("");
    sstr <<":"<<res->getId()<<"@"<<a->getName();
    return sstr.str();
}






#define	REAL	double
#define	INT	int



#ifndef VERYSMALL
#define VERYSMALL       0.000000000001
#endif
#define VERYSMALLSQUARED       (VERYSMALL*VERYSMALL)
#define	TENM3		0.001


#ifndef	SIGN
#define	SIGN(x)	((x)<0.0?-1.0:1.0)
#endif
#define	myPow(x,y)	(pow(x,y))
#define	ArcCos(x)	(acos(x))
#define	Cos(x)		(cos(x))
#define	Sin(x)		(sin(x))

#define	TWOPI		(3.14159265*2.0)
#define	reciprocal(x)	(1.0/(x))
#define	mysqrt(x)	sqrt(x)
#define	power1(x)	(x)
#define	power2(x)	((x)*(x))
#define	two		2.0
#define	six		6.0
#define	oneHalf		0.5
#define	twelve		12.0
#define	True		true
#define	False		false

	// ChiralOffset is a term added in the chiral restraint energy function
	// to try and prevent the four atoms in a chiral restraint from
	// settling into the same plane
#define	ChiralOffset	0.0



#define	TEST_FORCE( func, delta, argLow, argHigh, term, idx ) {	\
	double eLow = func argLow;\
	double eHigh = func argHigh;\
	double numForce = -(eHigh-eLow)/(delta);\
	LOG(BF("eHigh = %le eLow = %le delta = %le")% (eHigh) % (eLow) % (delta) ); \
	LOG(BF("numForce = %le")% numForce );			\
	if ( !_areValuesClose( numForce, term, #func, #term, idx ) ) { fails++; };\
}

#define	TEST_DIAGONAL_HESSIAN( func, delta, argLow, argMiddle, argHigh, term, idx ) {\
	double eLow = func argLow;\
	double eMiddle = func argMiddle;\
	double eHigh = func argHigh;\
	double numHessian = (eHigh+eLow-2.0*eMiddle)/(((delta)/2.0)*((delta)/2.0));\
	if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
}

#define	TEST_OFF_DIAGONAL_HESSIAN( func, delta, argMM, argPM, argMP, argPP, term, idx ) {\
	double eMM = func argMM;\
	double eMP = func argMP;\
	double ePM = func argPM;\
	double ePP = func argPP;\
	double numHessian = ((ePP-ePM)-(eMP-eMM))/(delta*delta);\
	if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
}


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



#define	ALL_ENERGY_COMPONENTS(msg) {		\
	this->_Stretch->msg; 			\
	this->_Angle->msg; 			\
	this->_Dihedral->msg; 			\
	this->_Nonbond->msg; 			\
	this->_ImproperRestraint->msg; 		\
	this->_ChiralRestraint->msg; 		\
	this->_AnchorRestraint->msg; 		\
	this->_FixedNonbondRestraint->msg; 	\
}

template <class ComponentType, class EntryType>
string	component_beyondThresholdInteractionsAsString(ComponentType& component) 
{
    int	bt;
    stringstream	ss;
	bt = component._BeyondThresholdTerms.end()-
		component._BeyondThresholdTerms.begin();
	ss << component.className() << "(#"<< bt <<") ";
	return ss.str();
};


#ifdef XML_ARCHIVE
template <class ComponentType, class EntryType>
    void	archiveEnergyComponentTerms(core::ArchiveP node,ComponentType& component)
{_G();
    if ( node->loading() )
    {
	core::VectorNodes::iterator	vi;
	EntryType term(node->lisp());
	component._Terms.clear();
	component._BeyondThresholdTerms.clear();
	for ( vi=node->begin_Children(); vi!=node->end_Children(); vi++ )
	{
	    if ( (*vi)->isNamed(component.className()+"Term") )
	    {
		term.archive(*vi);
		component._Terms.push_back(term);
	    } else if ( (*vi)->isNamed(component.className()+"BeyondThresholdTerm") )
	    {
		term.archive(*vi);
		component._BeyondThresholdTerms.push_back(term);
	    } else
	    {
		ARCHIVE_ERROR(BF("Illegal node"), *vi);
	    }
	}
    } else
    {
	typename vector<EntryType>::iterator ti;
	core::ArchiveP	child;
	for (ti=component._Terms.begin(); ti!=component._Terms.end(); ti++ )
	{
	    child = node->createChildNode(component.className()+"Term");
	    ti->archive(child);
	}
	for (ti=component._BeyondThresholdTerms.begin(); ti!=component._BeyondThresholdTerms.end(); ti++ )
	{
	    child = node->createChildNode(component.className()+"BeyondThresholdTerm");
	    ti->archive(child);
	}
    }
}
#endif


SMART(EnergyComponent );
class EnergyComponent_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,EnergyComponent_O,"EnergyComponent");
public:
public: // virtual functions inherited from Object
	void initialize();
#if XML_ARCHIVE
	void	archiveBase(core::ArchiveP node);
#endif

protected: // instance variables
	bool		_Enabled;
	double		_Scale;
	double		_ErrorThreshold;
	bool		_DebugEnergy;
	double		_TotalEnergy;
	int		_Debug_NumberOfTermsToCalculate;
	stringstream	_DebugLog;
//protected:		// Define these in subclasses
//	vector<TermClass>	_Terms;
//	vector<TermClass>	_BeyondThresholdTerms;
public:
	virtual int numberOfTerms() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
void setScale(double s) {this->_Scale = s; };
	double getScale() { return this->_Scale ; };
	void enable() {this->_Enabled = true; };
	void disable() {this->_Enabled = false; };
	void enableDebug() {this->_DebugEnergy = true; };
	void disableDebug() {this->_DebugEnergy = false; };
	void setDebug_NumberOfTermsToCalculate(int i) {this->_Debug_NumberOfTermsToCalculate=i;}
	bool isEnabled() { return this->_Enabled; };
	void	setErrorThreshold(double t) { this->_ErrorThreshold = t; };
	double	getErrorThreshold() { return this->_ErrorThreshold; };

	string summarizeEnergyAsString();
	string enabledAsString();
	string debugLogAsString();

public:	// Virtual methods
    virtual double getEnergy() { return this->_TotalEnergy; };

    virtual void zeroEnergy();
virtual void dumpTerms() {_OF();SUBCLASS_MUST_IMPLEMENT();};

virtual	void evaluateAll ( NVector_sp 	pos,
		bool 		calcForce,
		NVector_sp 	force,
       		bool		calcDiagonalHessian,
		bool		calcOffDiagonalHessian,
		AbstractLargeSquareMatrix_sp	hessian,
		NVector_sp	hdvec,
		NVector_sp dvec) 
{_OF();SUBCLASS_MUST_IMPLEMENT();};


virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos ) {_OF();SUBCLASS_MUST_IMPLEMENT();};

virtual string	beyondThresholdInteractionsAsString() {_OF();SUBCLASS_MUST_IMPLEMENT();};

virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
    NVector_sp pos ) {_OF();SUBCLASS_MUST_IMPLEMENT();};
public:
	EnergyComponent_O( const EnergyComponent_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(EnergyComponent_O);
};

};


TRANSLATE(chem::EnergyComponent_O);
#endif




