/*
    File: molecule.h
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





#ifndef	MOLECULE_H
#define MOLECULE_H

#include <iostream>
#include <sstream>
#include <clasp/core/common.h>
#include <cando/chem/matter.h>
#include <cando/chem/atomIdMap.fwd.h>
#include <cando/chem/residue.h>

#include <cando/chem/chemPackage.h>


namespace chem {


  class	Loop;

  SMART(Molecule);
  class Molecule_O : public Matter_O
  {
    friend class Aggregate_O;
    LISP_CLASS(chem,ChemPkg,Molecule_O,"Molecule",Matter_O);
  public:
    core::T_sp      _ForceFieldName;
    core::Symbol_sp _Type;
  public:
    void initialize();
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
    friend	class	Loop;


  public:
    static Molecule_sp make(MatterName name, core::List_sp residues = nil<core::T_O>());

  private:
    void	duplicate(const Molecule_O* m );
  public:
    CL_DEFMETHOD core::Symbol_sp molecule_type() { return this->_Type; };
    CL_DEFMETHOD core::Symbol_sp setf_molecule_type(core::Symbol_sp type) { this->_Type = type; return type; };
  
    char	getMatterType() { return MOLECULE_CLASS; };
    virtual	char	getClass() { return moleculeId;};

    Molecule_O(const Molecule_O& mol );

    virtual string __repr__() const;

    contentIterator begin_residues() { return this->_Contents.begin(); };
    contentIterator end_residues() {return this->_Contents.end(); };
    const_contentIterator begin_residues() const { return this->_Contents.begin(); };
    const_contentIterator end_residues() const {return this->_Contents.end(); };

//    virtual bool equal(core::T_sp obj) const;
    virtual void	transferCoordinates(Matter_sp other);

    // A force field name is either a symbol or a (cons symbol T/NIL)
    // The CDR of the cons is T if the given types should be used
    //                     and NIL if types should be calculated
    core::T_sp force_field_name() const;
    core::T_sp force_field_use_given_types() const;
    void setf_force_field_name(core::T_sp name);
    
    CL_LISPIFY_NAME("firstResidue");
    CL_DEFMETHOD  Residue_sp	firstResidue() {return this->contentAt(0).as<Residue_O>();};
    CL_LISPIFY_NAME("getResidue");
    CL_DEFMETHOD  Residue_sp	getResidue(int i) {return this->contentAt(i).as<Residue_O>();};
    CL_LISPIFY_NAME("residueWithId");
    CL_DEFMETHOD  Residue_sp	residueWithId( int lid ) { return this->contentWithId(lid).as<Residue_O>(); };
    CL_LISPIFY_NAME("hasResidueWithId");
    CL_DEFMETHOD 	bool		hasResidueWithId( int lid ) { return this->hasContentWithId(lid); };

    int		numberOfResiduesWithName(MatterName name);
    Residue_sp	getFirstResidueWithName(MatterName name);
    VectorResidue	getResiduesWithName(MatterName name);

    
    void addResidue( Matter_sp r );
    void addResidueRetainId( Matter_sp r );
    void removeResidue( Matter_sp a );




    virtual bool isMolecule() { return true;};

    CL_DEFMETHOD size_t residueCount() {return (this->_Contents.size());};

    void		moveAllAtomsIntoFirstResidue();

    virtual string	subMatter() { return Residue_O::static_className(); };
    virtual	string	description() const { stringstream ss; ss << "molecule("<<_rep_(this->getName())<<")@"<<std::hex<<this<<std::dec; return ss.str();};

    CL_LISPIFY_NAME("testMoleculeConsistancy");
    CL_DEFMETHOD 	bool		testMoleculeConsistancy() {IMPLEMENT_ME();};

    bool	canRender() { return true;}
#ifdef RENDER
    geom::Render_sp	rendered(core::List_sp options);
#endif

    virtual bool applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value);

    virtual uint	numberOfAtoms() const;

	/*! Build a map of AtomIds to Atoms */
    virtual AtomIdMap_sp buildAtomIdMap() const;

    virtual Atom_sp atomWithAtomId(const AtomId& atomId) const;

    Molecule_O() : Base(), _ForceFieldName(kw::_sym_default), _Type(unbound<core::Symbol_O>()) {};

  public:
    virtual Matter_sp	copy(core::T_sp new_to_old);
    
  protected:
    virtual Matter_sp copyDontRedirectAtoms(core::T_sp new_to_old);
    virtual void redirectAtoms();


  };

};
#endif
