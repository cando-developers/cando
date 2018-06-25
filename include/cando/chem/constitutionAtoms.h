/*
    File: constitutionAtoms.h
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
#ifndef	ConstitutionAtoms_H //[
#define ConstitutionAtoms_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/lispList.h>
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
    CL_DOCSTRING(R"(Describes a bond from one chem:constitution-atom to another.)");
    LISP_CLASS(chem,ChemPkg,ConstitutionBond_O,"ConstitutionBond",core::CxxObject_O);
    public:
	//! The ConstitutionAtomIndex0N of the atom that is bonded to
    ConstitutionAtomIndex0N	_ToAtomIndex;
	//! The bond order of the bond
    BondOrder		_BondOrder;
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    ConstitutionAtomIndex0N toAtomIndex() const { return this->_ToAtomIndex;};
    BondOrder bondOrder() const { return this->_BondOrder; };
  public:
  ConstitutionBond_O(ConstitutionAtomIndex0N toAtomIndex, BondOrder bo) : _ToAtomIndex(toAtomIndex), _BondOrder(bo) {};
  };

  CL_DEFUN inline ConstitutionBond_sp makeConstitutionBond(chem::ConstitutionAtomIndex0N index, chem::BondOrder order) {
    // Allocate a ConstitutionBond_O
    return gctools::GC<ConstitutionBond_O>::allocate(index,order);
  }


  class ConstitutionAtom_O : public core::CxxObject_O
  {
    friend class ConstitutionAtoms_O;
    friend class TopologyAtom_O;
    friend class TopologyVirtualAtom_O;
    friend class StereoisomerAtom_O;
    friend class StereoisomerVirtualAtom_O;

    CL_DOCSTRING(R"(Maintains information to create a chem:atom within a chem:residue. 
It stores the atom name, element, properties and a vector of bonds in the form of chem:constitution-bond(s).)");
    LISP_CLASS(chem,ChemPkg,ConstitutionAtom_O,"ConstitutionAtom",core::CxxObject_O);
  public:
  public:
    MatterName                          _AtomName;
    size_t                              _Index;
    Element			        _Element;
    StereochemistryType                 _StereochemistryType;
    core::List_sp                       _Properties;
    gctools::Vec0<ConstitutionBond_sp>	_Bonds;
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    string __repr__() const;
    CL_LISPIFY_NAME("atomName");
    CL_DEFMETHOD 	MatterName atomName() { return this->_AtomName;};
    virtual bool isVirtualAtom() { return false;};
	/*! Append a ConstitutionBond_sp to our list of bonds */
    void addConstitutionBond(ConstitutionBond_sp cb) {this->_Bonds.push_back(cb);};
  ConstitutionAtom_O(MatterName atomName, Element element, size_t index, StereochemistryType stype, core::List_sp properties) : _AtomName(atomName), _Element(element), _Index(index), _StereochemistryType(stype), _Properties(core::cl__copy_list(properties)) {};

    virtual Atom_sp buildAtom() const;
  };

  // Put the namespace in front of the types - otherwise there will be problems with the wrappers
  CL_DEFUN inline ConstitutionAtom_sp makeConstitutionAtom(chem::MatterName uniqueAtomName, chem::Element element, size_t index, chem::StereochemistryType stereochemistry_type, core::List_sp properties) {
    return gctools::GC<ConstitutionAtom_O>::allocate(uniqueAtomName,element,index,stereochemistry_type,properties);
  }

  class ConstitutionVirtualAtom_O : public ConstitutionAtom_O
  {
    LISP_CLASS(chem,ChemPkg,ConstitutionVirtualAtom_O,"ConstitutionVirtualAtom",ConstitutionAtom_O);
  public:
    void fields(core::Record_sp node);
  public:
    CalculatePosition_sp	_CalculatePositionCode;
  public:
    virtual bool isVirtualAtom() { return true;};

  public:
  ConstitutionVirtualAtom_O(MatterName atomname, Element element, size_t index, StereochemistryType stype, core::List_sp properties, CalculatePosition_sp calcPos) :
    ConstitutionAtom_O(atomname,element,index,stype,properties), _CalculatePositionCode(calcPos) {};
    virtual Atom_sp buildAtom() const;
  };

  // Put the namespace in front of the types - otherwise there will be problems with the wrappers
  CL_DEFUN inline ConstitutionVirtualAtom_sp makeConstitutionVirtualAtom(core::Symbol_sp atomName, chem::Element element, size_t index, chem::StereochemistryType stereochemistry_type, core::List_sp properties, chem::CalculatePosition_sp calcPos) {
    return gctools::GC<ConstitutionVirtualAtom_O>::allocate(atomName,element,index,stereochemistry_type,properties,calcPos);
  }

  class ConstitutionAtoms_O : public core::CxxObject_O
  {
    friend class StereoisomerAtoms_O;
    CL_DOCSTRING(R"(Stores a vector of chem:constitution-atom(s) that are used to create chem:residue(s) and other things
like residues. This class can be used to create all of the atoms in the residue, with atom names and elements and connect them
up with bonds with the correct bond orders.  It doesn't connect the atoms to other residues so there will be dangling bonds
that need to be formed later to other residues. The vector of chem:constitution-atom(s) have a specific order and 
the chem:constitution-atom(s) can be accessed using integer indices.)");
    LISP_CLASS(chem,ChemPkg,ConstitutionAtoms_O,"ConstitutionAtoms",core::CxxObject_O);
  private: // instance variables
	//! A list of ConstitutionAtoms
    gctools::Vec0<ConstitutionAtom_sp>	_Atoms;
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    typedef	gctools::Vec0<ConstitutionAtom_sp>::iterator	iterator;
	//! Create a ConstitutionAtoms object from a Residue
    static ConstitutionAtoms_sp makeConstitutionAtomsFromResidue(Residue_sp residue, bool verbose=false);
  public:

	/*! Return the number of atoms */
    CL_DEFMETHOD int numberOfAtoms() const { return this->_Atoms.size();};

	/*! Return the index of the ConstitutionAtom with the given name */
    int	index(MatterName name) const;

	/*! Return the ConstitutionAtom with the give ConstitutionAtomIndex0N */
    ConstitutionAtom_sp& operator[](ConstitutionAtomIndex0N idx);

	//! Return the ConstitutionAtom that has the give name
    CL_DEFMETHOD ConstitutionAtom_sp atomWithName(MatterName name);

	/*! Return the ConstitutionAtom at the given index */
    CL_DEFMETHOD ConstitutionAtom_sp atomWithId(ConstitutionAtomIndex0N idx) const;

	//! Add a ConstitutionVirtualAtom to us and assign it a ConstitutionAtomIndex0N
    void addConstitutionVirtualAtom(ConstitutionVirtualAtom_sp cva);

    core::List_sp constitutionAtomsAsList() const;
    
	//! Return a StringSet of the ConstitutionAtom names
    adapt::SymbolSet_sp atomNamesAsSymbolSet();

    Residue_sp buildResidue();
    
    ConstitutionAtoms_O( const ConstitutionAtoms_O& ss ); //!< Copy constructor
    DEFAULT_CTOR_DTOR(ConstitutionAtoms_O);
  };

  CL_DEFUN inline ConstitutionAtoms_sp makeConstitutionAtoms() {
    return gctools::GC<ConstitutionAtoms_O>::allocate();
  }
};
#endif //]
