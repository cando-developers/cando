/*
    File: stereoisomerAtoms.h
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
  class StereoisomerAtom_O : public core::CxxObject_O
  {
    friend class StereoisomerAtoms_O;
    CL_DOCSTRING(R"(Maintain information about an atom with respect to the topology that it belongs.
A topology cares about properties that are influenced by stereochemistry and how the 
topology atoms plug into others.  So it stores things like charge and atom type.)");
    LISP_CLASS(chem,ChemPkg,StereoisomerAtom_O,"StereoisomerAtom",core::CxxObject_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
	/*! Construct a StereoisomerAtom given its uniqueAtomName and its unique ConstitutionAtomIndex0N */
    static StereoisomerAtom_sp create(ConstitutionAtom_sp constitutionAtom );
    static StereoisomerAtom_sp make(core::Symbol_sp atomName, core::Symbol_sp type, double charge, ConstitutionAtomIndex0N index);
  public:
    MatterName		        _AtomName;
    ConstitutionAtomIndex0N	_ConstitutionAtomIndex;
    double			_AtomCharge;
    core::Symbol_sp		_AtomType;
  public:
    string __repr__() const;
    CL_LISPIFY_NAME("atomName");
    CL_DEFMETHOD 	MatterName atomName() { return this->_AtomName;};
    ConstitutionAtomIndex0N constitutionAtomIndex() const { return this->_ConstitutionAtomIndex;};
    virtual bool isVirtualAtom() { return false;};
    CL_DEFMETHOD 	double getAtomCharge() const { return this->_AtomCharge;};
    CL_DEFMETHOD 	void setAtomCharge(double c) { this->_AtomCharge = c;};
    CL_DEFMETHOD         core::Symbol_sp getAtomType() const { return this->_AtomType;};
    CL_DEFMETHOD 	void setAtomType(core::Symbol_sp s) { this->_AtomType = s;};
  StereoisomerAtom_O(MatterName name, ConstitutionAtomIndex0N index, double charge, core::Symbol_sp type) : _AtomName(name), _ConstitutionAtomIndex(index), _AtomCharge(charge), _AtomType(type) {};
    StereoisomerAtom_O() {}
  };



  class StereoisomerVirtualAtom_O : public StereoisomerAtom_O
  {
    LISP_CLASS(chem,ChemPkg,StereoisomerVirtualAtom_O,"StereoisomerVirtualAtom",StereoisomerAtom_O);
  public:
  public:
	/*! Construct a StereoisomerVirtualAtom given its uniqueAtomName and its uique ConstitutionAtomIndex0N */
    static StereoisomerVirtualAtom_sp create(MatterName uniqueAtomName, ConstitutionAtomIndex0N index,
                                             CalculatePosition_sp code );
	/*! Construct a StereoisomerVirtualAtom given its uniqueAtomName and leave its ConstitutionAtomIndex0N undefined */
    static StereoisomerVirtualAtom_sp create(ConstitutionVirtualAtom_sp cva );
  public:
    virtual bool isVirtualAtom() { return true;};

  StereoisomerVirtualAtom_O(MatterName name, ConstitutionAtomIndex0N index, double charge, core::Symbol_sp type) : StereoisomerAtom_O(name,index,charge,type) {};
    StereoisomerVirtualAtom_O() {};
  };




  class StereoisomerAtoms_O : public core::CxxObject_O
  {
    CL_DOCSTRING(R"(Maintain residue information for the topology.  This
is a list of stereoisomer-atoms.)");
    LISP_CLASS(chem,ChemPkg,StereoisomerAtoms_O,"StereoisomerAtoms",core::CxxObject_O);

  public: // virtual functions inherited from Object
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
//	string	__repr__() const;

  public: // instance variables
	//! A list of StereoisomerAtoms
    core::Symbol_sp _StereoisomerName;
    gctools::Vec0<StereoisomerAtom_sp>	_Atoms;
  public:
	//! Create a StereoisomerAtoms object from a ConstitutionAtoms object
    static StereoisomerAtoms_sp create(core::Symbol_sp stereoisomerName );
	//! Create a StereoisomerAtoms object from a ConstitutionAtoms object
    static StereoisomerAtoms_sp create(core::Symbol_sp stereoisomerName, ConstitutionAtoms_sp constitutionAtoms);
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
    void addStereoisomerAtom(StereoisomerAtom_sp cva);

    CL_DEFMETHOD core::Symbol_sp getName() const;
	//! Return a StringSet of the StereoisomerAtom names
    adapt::SymbolSet_sp atomNamesAsSymbolSet();


	/*! Create a Residue that has all the atoms/bonds and all the necessary atom/bond
	  properties set properly for this StereoisomerAtoms */
    Residue_sp makeResidue();

    StereoisomerAtoms_O(core::Symbol_sp stereoisomerName);
    StereoisomerAtoms_O() {};
  };

};
#endif //]
