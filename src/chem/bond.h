       
       
//
// (C) 2004 Christian E. Schafmeister
//




#ifndef	BOND_H
#define	BOND_H

#include "core/common.h"
#include "bond.fwd.h"
#include "matter.h"

#include "constitutionAtoms.fwd.h"

#include "chemPackage.h"


namespace chem {






SMART(Bond);
class Bond_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,Bond_O,"Bond");
//    DECLARE_SERIALIZE();
    public:
static void lisp_initGlobals(core::Lisp_sp lisp);

	friend	class Loop;
private:
    BondOrder	order;
    Atom_sp		_Atom1;
    Atom_sp		_Atom2;
    /*! Store Symbol keyed properties of bond in alist */
    core::Cons_sp	_Properties;
public:
	void initialize();
//	void	archiveBase(core::ArchiveP node);


public:
	/*! Copy the bond but don't redirect atoms yet */
	Bond_sp copyDontRedirectAtoms();

	/*! Redirect atoms to their copy atoms */
	void redirectAtoms();

public:
	static Bond_sp create(Atom_sp a1, Atom_sp a2, BondOrder o);
public:
	/*! Return true if atom is Atom1, false if its Atom2 and throw an exception if its neither */
	bool isAtom1(Atom_sp atom) const;



    Atom_sp getAtom1() const { return this->_Atom1; };
    Atom_sp getAtom2() const { return this->_Atom2; };

	BondOrder getOrder()	{ return this->order; };
	int getOrderAsInt()	{ return this->order; };
	string	getOrderAsString();
	void	setOrder(BondOrder o)	{ this->order=o; };
	void	setOrderFromInt(int o)	{ this->order=(BondOrder)(o); };
	Atom_sp	getOtherAtom(Atom_sp atom ) const;
	bool	isInterResidueBond();
	bool	isIntraResidueBond() {return !this->isInterResidueBond();};
	void	redirectToAtomCopies();

	/*! Dump all of the properties associated with this matter */
	string	propertiesAsString() const;
	    
		/*! Remove the property from this Matters property list
		 */
	void	clearProperty(core::Symbol_sp propertySymbol);
		/*! Define/set the value of the property.
		 * If it already exists it is overwritten.
		 */
	void	pushProperty(core::Symbol_sp propertySymbol, core::T_sp value);
		/*! Define/set the value of the property to TRUE.
		 * If it already exists it is overwritten.
		 */
	void	pushPropertyTrue(core::Symbol_sp propertySymbol);
		/*! Return the value of the property or the
		 * default if it isn't defined.
		 */
	core::T_sp getProperty(core::Symbol_sp propertySymbol, core::T_sp defVal );
		/*! Return true if the property exists.
		 */
	bool	hasProperty(core::Symbol_sp propertySymbol );


	bool isBondBetween(Atom_sp a, Atom_sp b) const;



	/*! Return a ConstitutionBond_sp object using the ConstitutionBond_O::atomIdMap 
	  to map Atoms to ConstitutionAtomIndex0N values */
	ConstitutionBond_sp asConstitutionBond(Atom_sp from, MapAtomsToConstitutionAtomIndex0N atomMap);

		//! For bonds that are read from XML we instruct them to join the atoms that they are supposed to join
	void		joinYourAtoms();


	string description() const;
	string describeOther(Atom_sp from) const;

	bool		invalid(Atom_sp a);
	void		failIfInvalid(Atom_sp a);

		/*! Bonds that are archived from a stream need to tell their
		 * atoms that they are bonded to each other.
		 */
	void	imposeYourself();

	Bond_O( const Bond_O& bb );


	DEFAULT_CTOR_DTOR(Bond_O);
};

extern	core::NullTerminatedEnumAssociation bondOrderKeys[];

string		bondOrderToString(BondOrder bo);
char            bondOrderToChar(BondOrder bo);
BondOrder	stringToBondOrder(string bos);

SMART(BondList);
class BondList_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,BondList_O,"BondList");
//    DECLARE_SERIALIZE();
public:
	void initialize();
public:
//	void	archiveBase(core::ArchiveP node);
private:
    gctools::Vec0<Bond_sp>	_Bonds;
//	gctools::Vec0<Bond_sp>	_Bonds;
public:

	BondList_O(const BondList_O& );


public:
//	typedef gctools::Vec0<Bond_sp>::iterator	iterator;
	typedef gctools::Vec0<Bond_sp>::iterator	iterator;
	typedef gctools::Vec0<Bond_sp>::const_iterator	const_iterator;

	iterator begin() { return this->_Bonds.begin(); };
	iterator end() { return this->_Bonds.end(); };
	Bond_sp		atIndex(int idx) { return this->_Bonds[idx];};
	void		addBond(Bond_sp b);
	void		append(Bond_sp b) { this->_Bonds.push_back(b);};
	void		removeBond(Bond_sp b);
	void		removeBondBetween(Atom_sp a,Atom_sp b);
	int		size() { return this->_Bonds.size(); };
//	BondList_sp	copy();

	string	description() const;
	string describeOthers(Atom_sp from) const;

		//! Impose all of your bonds on their atoms (serialization)
	void	imposeYourself();


	DEFAULT_CTOR_DTOR(BondList_O);
};


};

TRANSLATE(chem::Bond_O);
TRANSLATE(chem::BondList_O);


namespace kw {
    extern core::Symbol_sp _sym_singleBond;
    extern core::Symbol_sp _sym_doubleBond;
    extern core::Symbol_sp _sym_tripleBond;
    extern core::Symbol_sp _sym_aromaticBond;
    extern core::Symbol_sp _sym_hydrogenBond;
    extern core::Symbol_sp _sym_virtualBond;
    extern core::Symbol_sp _sym_dashedSingleBond;
    extern core::Symbol_sp _sym_dashedDoubleBond;
};








#endif
