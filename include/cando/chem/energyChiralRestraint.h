       
       

//
// (C) 2004 Christian E. Schafmeister
//



/*
 *	energyChiralRestraint.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyChiralRestraint_H  //[
#define	EnergyChiralRestraint_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>


namespace chem
{


struct TermChiralRestraint
{
	double		K;
	double		CO;
	int		I1;
	int		I2;
	int		I3;
	int		I4;
};

inline	string	XmlTag_ChiralRestraint() { return "ChiralRestraint"; };



class	EnergyChiralRestraint : public EnergyTerm {
public:
	string	className()	{ return "EnergyChiralRestraint"; };
public:
		// Variables
        gc::Nilable<Atom_sp>      _Atom1;
        gc::Nilable<Atom_sp>      _Atom2;
        gc::Nilable<Atom_sp>      _Atom3;
        gc::Nilable<Atom_sp>      _Atom4;
		// Threshold
	TermChiralRestraint	term;

#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_ChiralRestraint_debugEvalDeclares.cc>
#endif
        gc::Nilable<Atom_sp>	getAtom1() { return this->_Atom1; };
        gc::Nilable<Atom_sp>	getAtom2() { return this->_Atom2; };
        gc::Nilable<Atom_sp>	getAtom3() { return this->_Atom3; };
        gc::Nilable<Atom_sp>	getAtom4() { return this->_Atom4; };
public:
//	void	archive(core::ArchiveP node);
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );
public:
	string	description();
    EnergyChiralRestraint();
	virtual ~EnergyChiralRestraint();
};




double	_evaluateEnergyOnly_ChiralRestraint(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3,
		double x4, double y4, double z4,
		double K, double CO );




class EnergyChiralRestraint_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyChiralRestraint_O,"EnergyChiralRestraint",EnergyComponent_O);
public:
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

public:
    typedef EnergyChiralRestraint	TermType;
public: // instance variables
    gctools::Vec0<TermType>	_Terms;
    gctools::Vec0<TermType>	_BeyondThresholdTerms;

public:	// Creation class functions
public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
public:
    virtual int numberOfTerms() { return this->_Terms.size();};
    

public:
    void addTerm(const TermType& term);
    virtual void dumpTerms();

    virtual void setupHessianPreconditioner(NVector_sp nvPosition,
					    AbstractLargeSquareMatrix_sp m );
    virtual	void evaluateAll ( NVector_sp 	pos,
				   bool 		calcForce,
				   gc::Nilable<NVector_sp> 	force,
				   bool		calcDiagonalHessian,
				   bool		calcOffDiagonalHessian,
				   gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
				   NVector_sp	hdvec,
				   NVector_sp dvec);

    virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
	NVector_sp pos );

    virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

    virtual string	beyondThresholdInteractionsAsString();

public:
    EnergyChiralRestraint_O( const EnergyChiralRestraint_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyChiralRestraint_O);
};

};

TRANSLATE(chem::EnergyChiralRestraint_O);
#endif //]


