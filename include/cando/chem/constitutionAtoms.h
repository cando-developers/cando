#ifndef	ConstitutionAtoms_H //[
#define ConstitutionAtoms_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>

#include <cando/chem/atom.fwd.h>
#include <cando/chem/bond.fwd.h>
#include <cando/chem/residue.fwd.h>
#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/chem/calculatePosition.fwd.h>
#include <cando/chem/elements.fwd.h>

#include <cando/chem/chemPackage.h>


namespace chem
{



    class ConstitutionBond_O : public core::CxxObject_O
    {
	friend class ConstitutionAtoms_O;
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,ConstitutionBond_O,"ConstitutionBond");
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
	static ConstitutionBond_sp create(BondOrder order, ConstitutionAtomIndex0N index);

    protected:
	//! The ConstitutionAtomIndex0N of the atom that is bonded to
	ConstitutionAtomIndex0N	_AtomIndex;
	//! The bond order of the bond
	BondOrder		_BondOrder;
    public:
	ConstitutionAtomIndex0N atomIndex() const { return this->_AtomIndex;};
	DEFAULT_CTOR_DTOR(ConstitutionBond_O);
    };




    class ConstitutionAtom_O : public core::CxxObject_O
    {
	friend class ConstitutionAtoms_O;
	friend class TopologyAtom_O;
	friend class TopologyVirtualAtom_O;
	friend class StereoisomerAtom_O;
	friend class StereoisomerVirtualAtom_O;

	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,ConstitutionAtom_O,"ConstitutionAtom");
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
	/*! Construct a ConstitutionAtom given its uniqueAtomName and its unique ConstitutionAtomIndex0N */
	static ConstitutionAtom_sp create(MatterName uniqueAtomName, Element element,
					  ConstitutionAtomIndex0N index );
    protected:
	MatterName		_AtomName;
	Element			_Element;
	ConstitutionAtomIndex0N	_AtomIndex;
        gctools::Vec0<ConstitutionBond_sp>	_Bonds;
    public:
	string __repr__() const;
CL_NAME("atomName");
CL_DEFMETHOD 	MatterName atomName() { return this->_AtomName;};
	ConstitutionAtomIndex0N atomIndex() const { return this->_AtomIndex;};
	virtual bool isVirtualAtom() { return false;};
	/*! Append a ConstitutionBond_sp to our list of bonds */
	void addConstitutionBond(ConstitutionBond_sp cb) {_OF(); this->_Bonds.push_back(cb);};
	DEFAULT_CTOR_DTOR(ConstitutionAtom_O);
    };



    class ConstitutionVirtualAtom_O : public ConstitutionAtom_O
    {
	LISP_BASE1(ConstitutionAtom_O);
	LISP_CLASS(chem,ChemPkg,ConstitutionVirtualAtom_O,"ConstitutionVirtualAtom");
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    private:
	CalculatePosition_sp	_CalculatePositionCode;
    public:
	/*! Construct a ConstitutionVirtualAtom given its uniqueAtomName and its uique ConstitutionAtomIndex0N */
	static ConstitutionVirtualAtom_sp create(MatterName uniqueAtomName, ConstitutionAtomIndex0N index,
						 CalculatePosition_sp code );
	/*! Construct a ConstitutionVirtualAtom given its uniqueAtomName and leave its ConstitutionAtomIndex0N undefined */
	static ConstitutionVirtualAtom_sp create(MatterName uniqueAtomName, 
						 CalculatePosition_sp code );
    public:
	virtual bool isVirtualAtom() { return true;};

	DEFAULT_CTOR_DTOR(ConstitutionVirtualAtom_O);
    };




    class ConstitutionAtoms_O : public core::CxxObject_O
    {
	friend class StereoisomerAtoms_O;
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,ConstitutionAtoms_O,"ConstitutionAtoms");

    public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

    private: // instance variables
	//! A list of ConstitutionAtoms
        gctools::Vec0<ConstitutionAtom_sp>	_Atoms;
	int				_NetCharge;
    public:
	typedef	gctools::Vec0<ConstitutionAtom_sp>::iterator	iterator;
	//! Create a ConstitutionAtoms object from a Residue
	static ConstitutionAtoms_sp create(Residue_sp residue);
    public:

	/*! Return the number of atoms */
	int numberOfAtoms() const { return this->_Atoms.size();};

	/*! Return the index of the ConstitutionAtom with the given name */
	int	index(MatterName name) const;

	/*! Return the ConstitutionAtom with the give ConstitutionAtomIndex0N */
	ConstitutionAtom_sp& operator[](ConstitutionAtomIndex0N idx);

	//! Return the ConstitutionAtom that has the give name
	ConstitutionAtom_sp atomWithName(MatterName name);

	/*! Return the ConstitutionAtom at the given index */
	ConstitutionAtom_sp atomWithId(ConstitutionAtomIndex0N idx) const;

	//! Set the net charge of this collection of atoms
CL_NAME("setNetCharge");
CL_DEFMETHOD 	void setNetCharge(int charge) { this->_NetCharge = charge;};


	//! Add a ConstitutionVirtualAtom to us and assign it a ConstitutionAtomIndex0N
	void addConstitutionVirtualAtom(ConstitutionVirtualAtom_sp cva);

	//! Return a StringSet of the ConstitutionAtom names
	adapt::SymbolSet_sp atomNamesAsSymbolSet();


	/*! Create a Residue that has all the atoms/bonds and all the necessary atom/bond
	  properties set properly for this ConstitutionAtoms */
	Residue_sp makeResidue();

	ConstitutionAtoms_O( const ConstitutionAtoms_O& ss ); //!< Copy constructor
	DEFAULT_CTOR_DTOR(ConstitutionAtoms_O);
    };




};
TRANSLATE(chem::ConstitutionBond_O);
TRANSLATE(chem::ConstitutionAtom_O);
TRANSLATE(chem::ConstitutionVirtualAtom_O);
TRANSLATE(chem::ConstitutionAtoms_O);
#endif //]
