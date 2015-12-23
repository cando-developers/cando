//
// (C) 2004 Christian E. Schafmeister
//

/*
 *	energyDihedral.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyDihedral_H  //[
#define	EnergyDihedral_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>



#include <cando/chem/chemPackage.h>


namespace chem {

class EnergyAtom;

SMART(FFPtor);
SMART(FFItor);

inline	string	XmlTag_Dihedral() { return "Dihedral"; };
inline	string	XmlTag_EnergyDihedral() { return "EnergyDihedral"; };

struct TermDihedral
{
	REAL		sinPhase;	//!< Sine of phase, this must match Mathematica code!
	REAL		cosPhase;	//!< Cosine of phase, this must match Mathematica code!
	REAL		V;	//!< Dihedral force constant, this must match Mathematica code!
        REAL          	DN;	//!< Dihedral multiplicity as double, this must match Mathematica code!
        INT		IN;	//!< Dihedral multiplicity as int, this must match Mathematica code!
	INT		I1;	//!< i*3 index into coordinate array, this must match Mathematica code!
	INT		I2;	//!< i*3 index into coordinate array, this must match Mathematica code!
	INT		I3;	//!< i*3 index into coordinate array, this must match Mathematica code!
	INT		I4;	//!< i*3 index into coordinate array, this must match Mathematica code!
};

/*! Dihedral energy term
 */
class EnergyDihedral : public EnergyTerm 
{
public:
	string	className()	{ return "EnergyDihedral"; };
public:
                // Parameters
//	string		_Type1;
//	string		_Type2;
//	string		_Type3;
//	string		_Type4;
	bool		_Proper;
                // Variables
        gc::Nilable<Atom_sp>      _Atom1;
        gc::Nilable<Atom_sp>      _Atom2;
        gc::Nilable<Atom_sp>      _Atom3;
        gc::Nilable<Atom_sp>      _Atom4;
	double		_PhaseRad;
	TermDihedral	term;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_Dihedral_debugEvalDeclares.cc>
#endif
	double		_CalculatedDihedralDeviation;

        gc::Nilable<Atom_sp>		getAtom1() { return this->_Atom1; };
        gc::Nilable<Atom_sp>		getAtom2() { return this->_Atom2; };
        gc::Nilable<Atom_sp>		getAtom3() { return this->_Atom3; };
        gc::Nilable<Atom_sp>		getAtom4() { return this->_Atom4; };
	double		getDihedral();
	double		getDihedralDeviation();
	double		getCalculatedDihedralDeviation() { return this->_CalculatedDihedralDeviation;};

public:
//	void	archive(core::ArchiveP node);
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );

        void defineFrom( int n, FFPtor_sp term, EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale);
	void defineMissingProper( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4);
        void defineFrom( int n, FFItor_sp term, EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale);
    EnergyDihedral();
	virtual ~EnergyDihedral();
};


double  _evaluateEnergyOnly_Dihedral(
                double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3,
                double x4, double y4, double z4,
                double V, double DN, int IN,
                double cosPhase, double sinPhase );




class EnergyDihedral_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyDihedral_O,"EnergyDihedral",EnergyComponent_O);
public:
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);

public:
    typedef EnergyDihedral	TermType;
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
    EnergyDihedral_O( const EnergyDihedral_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyDihedral_O);

};


};

TRANSLATE(chem::EnergyDihedral_O);
#endif //]


