#ifndef	StereoisomerAtoms_H //[
#define StereoisomerAtoms_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/adapt/symbolSet.fwd.h>
#include <cando/chem/atom.fwd.h>
#include <cando/chem/bond.fwd.h>
#include <cando/chem/residue.fwd.h>
#include <cando/chem/stereoisomerAtoms.fwd.h>
#include <cando/chem/calculatePosition.fwd.h>
#include <cando/chem/elements.fwd.h>

#include <cando/chem/chemPackage.h>


namespace chem
{





    class StereoisomerAtom_O : public core::T_O
    {
	friend class StereoisomerAtoms_O;

	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,StereoisomerAtom_O,"StereoisomerAtom");
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
	/*! Construct a StereoisomerAtom given its uniqueAtomName and its unique ConstitutionAtomIndex0N */
	static StereoisomerAtom_sp create(ConstitutionAtom_sp constitutionAtom );
    protected:
	MatterName		_AtomName;
	ConstitutionAtomIndex0N	_AtomIndex;
	double			_Charge;
        core::Symbol_sp		_AtomType;
    public:
	string __repr__() const;
	MatterName atomName() { return this->_AtomName;};
	ConstitutionAtomIndex0N atomIndex() const { return this->_AtomIndex;};
	virtual bool isVirtualAtom() { return false;};
	double getCharge() const { return this->_Charge;};
	void setCharge(double c) { this->_Charge = c;};
        core::Symbol_sp getAtomType() const { return this->_AtomType;};
	void setAtomType(core::Symbol_sp s) { this->_AtomType = s;};
	DEFAULT_CTOR_DTOR(StereoisomerAtom_O);
    };



    class StereoisomerVirtualAtom_O : public StereoisomerAtom_O
    {
	LISP_BASE1(StereoisomerAtom_O);
	LISP_CLASS(chem,ChemPkg,StereoisomerVirtualAtom_O,"StereoisomerVirtualAtom");
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
	/*! Construct a StereoisomerVirtualAtom given its uniqueAtomName and its uique ConstitutionAtomIndex0N */
	static StereoisomerVirtualAtom_sp create(MatterName uniqueAtomName, ConstitutionAtomIndex0N index,
						 CalculatePosition_sp code );
	/*! Construct a StereoisomerVirtualAtom given its uniqueAtomName and leave its ConstitutionAtomIndex0N undefined */
	static StereoisomerVirtualAtom_sp create(ConstitutionVirtualAtom_sp cva );
    public:
	virtual bool isVirtualAtom() { return true;};

	DEFAULT_CTOR_DTOR(StereoisomerVirtualAtom_O);
    };




    class StereoisomerAtoms_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,StereoisomerAtoms_O,"StereoisomerAtoms");

    public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

    private: // instance variables
	//! A list of StereoisomerAtoms
        gctools::Vec0<StereoisomerAtom_sp>	_Atoms;
    public:
	//! Create a StereoisomerAtoms object from a ConstitutionAtoms object
	static StereoisomerAtoms_sp create(ConstitutionAtoms_sp );
    public:

	/*! Return the number of atoms */
	int numberOfAtoms() const { return this->_Atoms.size();};

	/*! Return the index of the StereoisomerAtom with the given name */
	int	index(MatterName name) const;

	/*! Return the StereoisomerAtom with the give ConstitutionAtomIndex0N */
	StereoisomerAtom_sp& operator[](ConstitutionAtomIndex0N idx);

	//! Return the StereoisomerAtom that has the give name
	StereoisomerAtom_sp atomWithName(MatterName name);

	/*! Return the StereoisomerAtom at the given index */
	StereoisomerAtom_sp atomWithId(ConstitutionAtomIndex0N idx) const;


	//! Add a StereoisomerVirtualAtom to us and assign it a ConstitutionAtomIndex0N
	void addStereoisomerVirtualAtom(StereoisomerVirtualAtom_sp cva);

	//! Return a StringSet of the StereoisomerAtom names
	adapt::SymbolSet_sp atomNamesAsSymbolSet();


	/*! Create a Residue that has all the atoms/bonds and all the necessary atom/bond
	  properties set properly for this StereoisomerAtoms */
	Residue_sp makeResidue();

	StereoisomerAtoms_O( const StereoisomerAtoms_O& ss ); //!< Copy constructor
	DEFAULT_CTOR_DTOR(StereoisomerAtoms_O);
    };




};
TRANSLATE(chem::StereoisomerAtom_O);
TRANSLATE(chem::StereoisomerVirtualAtom_O);
TRANSLATE(chem::StereoisomerAtoms_O);
#endif //]
