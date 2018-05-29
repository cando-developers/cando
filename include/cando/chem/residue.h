/*
    File: residue.h
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
      
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	residue.h
 *
 *	Maintain a Residue object, a collection of bonded atoms.
 *	Residues are contained within Molecule objects and can be
 *	bonded together.
 *
 *	To facilitate rapid lookup of Atoms within the Residue,
 *	use the getHandleForAtomNamed and atomFromHandle methods.
 *
 *		getHandleForAtomNamed returns a unique AtomHandle object that
 *			can later be used to quickly lookup the Atom
 *			that cooresponds to the handle using the
 *			atomFromHandle method.
 *
 *		atomFromHandle returns a pointer to the Atom that
 *			cooresponds to the AtomHandle object that
 *			is passed to it.
 *
 */

#ifndef	RESIDUE_H
#define	RESIDUE_H
#include <iostream>
#include <string>
#include <set>
#include <clasp/core/common.h>
#include <cando/adapt/symbolSet.fwd.h>
#include <clasp/core/lispStream.fwd.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/chem/matter.h>
#include <cando/chem/residue.fwd.h>
//#include <cando/chem/candoDatabaseReference.h>
#include <cando/chem/atomIdMap.fwd.h>
#include <cando/chem/atom.h>

#include <cando/chem/chemPackage.h>


namespace chem {


  class	Loop;

  SMART(CalculatePosition);
  SMART(Restraint);
  SMART(Alias);
  SMART(Constitution);
  SMART(Residue);

  SMART(Residue );
  class Residue_O : public Matter_O
  {
    friend class Constitution_O;
    friend class Molecule_O;
    friend	class	Loop;
    LISP_CLASS(chem,ChemPkg,Residue_O,"Residue",Matter_O);
  public:
    void initialize();
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
    core::Symbol_sp    _Type;
  private:
    bool	_Selected;
    int		_NetCharge;
    int		tempInt;
    MatterName	pdbName;
    uint		_FileSequenceNumber;
//	CandoDatabaseReference_sp	_Constitution;
        /*! Stores a unique label that is set when
         * a residue is created from a monomer.
         * This is used to match residues between different
         * Molecule objects that contain the exact same
         * molecule.
         * This allows Molecules that were built from
         * Oligomers to be saved to a file and then
         * loaded back and have their atoms matched
         * up with atoms of an identical Molecule that
         * was loaded from somewhere else.
         * This also allows Molecules that are sent between
         * processes to have their coordinates put into
         * ConformationExplorers and have Scorers applied to
         * them.
         * A MoleculeMatcher takes two Molecules and matches
         * their atoms 1 to 1 and allows coordinates to
         * be copied from either molecule to the other.
         */
    gc::Nilable<MatterName>	_UniqueLabel;
    gc::Nilable<adapt::SymbolSet_sp>	_MonomerAliases;
  public:
    static Residue_sp make(MatterName name);
  private:
    void	duplicateFrom(const Residue_O* r );
  public:

    CL_DEFMETHOD core::Symbol_sp residue_type() { return this->_Type; };
    CL_DEFMETHOD core::Symbol_sp setf_residue_type(core::Symbol_sp type) { this->_Type = type; return type; };
  

    virtual	char	getClass()	{return residueId; };
    char	getMatterType() { return RESIDUE_CLASS; };

    virtual bool equal(core::T_sp obj) const;
    virtual void	transferCoordinates(Matter_sp other);


    contentIterator begin_atoms() { return this->_contents.begin(); };
    contentIterator end_atoms() {return this->_contents.end(); };

    virtual bool isResidue() { return true;};
    virtual bool applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value);


    adapt::SymbolSet_sp getMonomerAliases();
    bool recognizesMonomerAlias(core::Symbol_sp s);
    void setMonomerAliases(adapt::SymbolSet_sp s);
    void addMonomerAlias(core::Symbol_sp s);

	/*! Return true if this residue contains the atom a
	 */
    bool	containsAtom(Atom_sp a);

    void ensureAllAtomNamesAreUnique() const;
    
    CL_LISPIFY_NAME("getPdbName");
    CL_DEFMETHOD 	MatterName getPdbName() { return this->pdbName; };
    CL_LISPIFY_NAME("setPdbName");
    CL_DEFMETHOD 	void	setPdbName(MatterName p) { this->pdbName = p;};

    void	setFileSequenceNumber(uint seq) { this->_FileSequenceNumber = seq;};
    CL_LISPIFY_NAME("getFileSequenceNumber");
    CL_DEFMETHOD 	uint	getFileSequenceNumber() { return this->_FileSequenceNumber;};

    void	setUniqueLabel(MatterName str) { this->_UniqueLabel = str;};
    MatterName getUniqueLabel() { return this->_UniqueLabel; };

    CL_LISPIFY_NAME("getNetCharge");
    CL_DEFMETHOD 	int	getNetCharge() { return this->_NetCharge; };
    CL_LISPIFY_NAME("setNetCharge");
    CL_DEFMETHOD 	void	setNetCharge(int nc) { this->_NetCharge = nc; };

    int		getTempInt() {return this->tempInt;};
    void		setTempInt(int j) {this->tempInt = j;};
    adapt::SymbolSet_sp	getAtomNamesAsSymbolSet();

    void		addAtom(Atom_sp a );
    void		addAtomRetainId(Atom_sp a ) {this->addMatterRetainId(Matter_sp(a));};
    void		removeAtomsWithNames( core::List_sp atomNames );
    void		removeAtomDeleteBonds( Atom_sp a );
    void		removeAtomDontDeleteBonds( Atom_sp a );


    int	totalNetResidueCharge() { return this->_NetCharge; };

        /*! Add a virtual atom and calculate its position using
         * the provided procedure
         */
    void	addVirtualAtom(MatterName virtualAtomName, 
                               CalculatePosition_sp proc);

    VectorAtom	getAtoms();
//	gctools::Vec0<Bond_sp>	g_etUniqueIntraResidueBonds();
    BondList_sp	getOutGoingBonds();

    virtual string	subMatter() { return Atom_O::static_className(); };
    virtual string 	description() const { stringstream ss; ss << "residue("<<_rep_(this->getName())<<")@"<<std::hex<<this<<std::dec; return ss.str();};

    CL_LISPIFY_NAME("atomWithName");
    CL_DEFMETHOD 	Atom_sp		atomWithName(MatterName sName ) { return this->contentWithName(sName).as<Atom_O>(); };
    CL_LISPIFY_NAME("hasAtomWithName");
    CL_DEFMETHOD 	bool		hasAtomWithName(MatterName sName ) { return (this->hasContentWithName(sName)); };

    Atom_sp		atomWithAlias(AtomAliasName sName );
    Atom_sp		atomWithAliasOrNil(AtomAliasName sName );


    Atom_sp	aliasAtomOrNil(Alias_sp alias) {IMPLEMENT_ME();};
    Residue_sp aliasResidueOrNil(Alias_sp alias) {IMPLEMENT_ME();};

    Vector3		positionOfAtomWithName(MatterName sName );

    void setAliasesForAtoms(core::List_sp aliasAtoms, core::List_sp atomAliases);

//	bool		hasAtomWithAlias(const string& anAlias);
//	Atom_sp		atomWithAlias(const string& anAlias);

    CL_LISPIFY_NAME("atomWithId");
    CL_DEFMETHOD         Atom_sp		atomWithId( int lid ) { return this->contentWithId(lid).as<Atom_O>(); };
    CL_LISPIFY_NAME("hasAtomWithId");
    CL_DEFMETHOD 	bool		hasAtomWithId( int lid ) { return this->hasContentWithId(lid); };
    adapt::SymbolSet_sp 	getAllUniqueAtomNames();

        /*! Return true if all the atom names are unique
         * If not describe the atoms with overlapping names in (problems).
         */
    bool		testIfAllAtomNamesAreUnique(core::T_sp problemStream);

    void fillInImplicitHydrogens();

    CL_LISPIFY_NAME("firstAtom");
    CL_DEFMETHOD 	Atom_sp		firstAtom() { if ( this->contentSize() <= 0 ) { SIMPLE_ERROR(BF("This matter has no contents")); } return this->contentAt(0).as<Atom_O>(); };

    void	makeAllAtomNamesInEachResidueUnique();

    void		writeToStream(string prefix, std::ostream& out);
    CL_LISPIFY_NAME("testResidueConsistancy");
    CL_DEFMETHOD 	bool		testResidueConsistancy()	{return this->testConsistancy(Matter_sp());};

	//! For every atom, turn on the FIXED flag and copy the atom coordinates into the AnchorCoordinates
    void		useAtomCoordinatesToDefineAnchors();


    bool		isSelected() {return this->_Selected;};
    void		setSelected(bool s) {this->_Selected = s;};

    bool		invalid();

    void		failIfInvalid();

	/*! Build a map of AtomIds to Atoms */
    virtual AtomIdToAtomMap_sp buildAtomIdMap() const;

    virtual Atom_sp atomWithAtomId(AtomId_sp atomId) const;


    virtual uint	numberOfAtoms( );

  Residue_O() : Base(), _Type(_Unbound<core::Symbol_O>()) {};
    Residue_O( const Residue_O& res );


  public:
    virtual Matter_sp	copy();
  protected:
    virtual Matter_sp copyDontRedirectAtoms();
    virtual void redirectAtoms();

  };
};
#endif
