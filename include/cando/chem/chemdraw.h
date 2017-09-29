/*
    File: chemdraw.h
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
#ifndef	ChemDraw_H //[
#define ChemDraw_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/holder.h>
#include <cando/geom/vector2.h>
#include <cando/chem/bond.fwd.h>
#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/chem/chemPackage.h>
#include <cando/adapt/quickDom.fwd.h>// chemdraw.h wants QDomNode needs quickDom.fwd.h
#include <cando/chem/atom.fwd.h>
#include <cando/chem/molecule.fwd.h>
#include <cando/chem/aggregate.fwd.h>
namespace chem
{


SMART(Bond);
SMART(Matter);
SMART(Atom);
SMART(Residue);
SMART(Binder);
SMART(QDomNode);
SMART(Lisp);


typedef	enum { singleCDBond,
		doubleCDBond,
		tripleCDBond,
		dativeCDBond,
		singleDashCDBond,
		doubleDashCDBond,
		tripleDashCDBond,
		hashCDBond,
               hollowWedgeCDBond,
               wedgeHashCDBond,
               wedgeCDBond,
                 wavyCDBond,
		unknownCDBond
    } CDBondOrder;


SMART(CDNode );
class CDNode_O : public core::CxxObject_O
{	
  friend class CDFragment_O;
  LISP_CLASS(chem,ChemPkg,CDNode_O,"CDNode",core::CxxObject_O);
 public:
  uint				_Id;
  uint                          _Color;
  geom::Vector2                 _Pos;
  std::string       		_Label;
  StereochemistryType		_StereochemistryType;
  ConfigurationEnum		_Configuration;
 private:	// generated
  gc::Nilable<Atom_sp>		_Atom;
  gc::Nilable<CDNode_sp>		_BackSpan;
  gc::Nilable<CDNode_sp>		_NextSpan;
  core::List_sp _AtomProperties;
  core::List_sp _ResidueProperties;
  core::List_sp _MoleculeProperties;
 public:
  void	initialize();
 public:
#if 0
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
#endif
 private:
  std::string	_extractLabel(adapt::QDomNode_sp node);
 public:
  uint	getId() const { return this->_Id;};
  void    setId(uint i) { this->_Id = i; };
//        core::Symbol_sp getLabel() const { return this->_Label;};
	/*! Parse the label into a name and ionization component "name/[+-] */
  void getParsedLabel(string& name, int& ionization) const;
//	void setLabel(core::Symbol_sp l) { this->_Label = l;};
  Atom_sp	getAtom() const { return this->_Atom;};
  void setAtom(Atom_sp a) { this->_Atom = a;};
  void	parseFromXml(adapt::QDomNode_sp xml, bool print);

  void	bondTo( CDNode_sp neighbor, CDBondOrder order);

  CDNode_O( const CDNode_O& ss ); //!< Copy constructor

 CDNode_O() : _Atom(_Nil<core::T_O>())
    , _BackSpan(_Nil<core::T_O>())
    , _NextSpan(_Nil<core::T_O>())
    , _AtomProperties(_Nil<core::T_O>())
    , _ResidueProperties(_Nil<core::T_O>())
    , _MoleculeProperties(_Nil<core::T_O>())
  {};
  virtual ~CDNode_O() {};
};


SMART(CDBond );
class CDBond_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,CDBond_O,"CDBond",core::CxxObject_O);
private:
	uint		_IdBegin;
	uint		_IdEnd;
	CDNode_sp	_BeginNode;
	CDNode_sp	_EndNode;
	CDBondOrder	_Order;
public:
	void	initialize();
 public:
#if 0
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
#endif
public:
	string	getOrderAsString();
	CDBondOrder	getOrder() { return this->_Order; };
	void setOrder(CDBondOrder o) { this->_Order = o; };
	BondOrder	getOrderAsBondOrder();
	uint		getIdBegin() { return this->_IdBegin;};
	void	setIdBegin(uint i) { this->_IdBegin = i;};
	uint		getIdEnd() { return this->_IdEnd;};
	void	setIdEnd(uint i) { this->_IdEnd = i;};
	void	setBeginNode(CDNode_sp n) { this->_BeginNode = n;};
	CDNode_sp getBeginNode() {_OF(); ANN(this->_BeginNode);return this->_BeginNode;};
	void	setEndNode(CDNode_sp n) { this->_EndNode = n;};
	CDNode_sp getEndNode() {_OF(); ANN(this->_EndNode);return this->_EndNode;};

	void	parseFromXml(adapt::QDomNode_sp xml, bool print);

	CDBond_O( const CDBond_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CDBond_O);
};



SMART(CDFragment );
class CDFragment_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,CDFragment_O,"CDFragment",core::CxxObject_O);
public:
    // Declare that gctools::Vec0<CDBond_sp> is managed by GC
    typedef gctools::Vec0<CDBond_sp>	CDBonds;
    typedef gctools::SmallMap<Atom_sp,CDNode_sp> AtomsToBonds; // map<Atom_sp,CDNode_sp>	AtomsToBonds;
    typedef gctools::SmallMap<int,CDNode_sp> IntMappedCDNodes; // typedef map<int,CDNode_sp>		IntMappedCDNodes;
private:
    core::Symbol_sp		_ConstitutionName;
    IntMappedCDNodes		_Nodes;
    AtomsToBonds		_AtomsToNodes;
    CDBonds			_Bonds;
    int				_LargestId;
    gc::Nilable<Molecule_sp> _Molecule;
//    core::HashTableEq_sp		_Properties;
public:
    void	initialize();
 public:
#if 0
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
#endif
private:
    bool _asKeyedObject(core::Symbol_sp label, core::Symbol_sp& keyword, core::T_sp& obj);
private:
//    Residue_sp _buildResidue(bool constitutionOnly);

public:
    void setConstitutionName(core::Symbol_sp n) { this->_ConstitutionName = n;};
CL_LISPIFY_NAME("getConstitutionName");
CL_DEFMETHOD     core::Symbol_sp getConstitutionName() { return this->_ConstitutionName;};

#if 0
    core::HashTableEq_sp getProperties() { return this->_Properties;};
    void	addProperties(core::HashTableEq_sp dict);
    //! Get the property for this ChemDraw Fragment
    core::T_sp getProperty(core::Symbol_sp key);
    core::T_sp getPropertyOrDefault(core::Symbol_sp key, core::T_sp df);
    bool hasProperty(core::Symbol_sp key);
    core::T_sp setProperty(core::Symbol_sp key, core::T_sp df);
    string describeProperties();
#endif
    Atom_sp createOneAtom(CDNode_sp n);
    void	createAtomsAndBonds();
    void	createImplicitHydrogen(CDNode_sp from, const string& name);
    void	parseFromXml(adapt::QDomNode_sp xml, bool print);
    /*! Return false if the fragment couldn't be interpreted */
    bool	interpret();

CL_LISPIFY_NAME("getMolecule");
CL_DEFMETHOD core::T_sp getMolecule() { return this->_Molecule; };


    int countNeighbors(CDNode_sp node);
    void createBonds(bool selectedAtomsOnly);
    Molecule_sp createMolecule();
    void removeAllBonds();
    void clearAtomSelected();


    Residue_sp createResidueOfSelectedAtoms();

//    Residue_sp getEntireResidue();
    ConstitutionAtoms_sp asConstitutionAtoms();

    string __repr__() const;

    CDFragment_O( const CDFragment_O& ss ); //!< Copy constructor

 CDFragment_O() : _Molecule(_Nil<core::T_O>()) {};
    virtual ~CDFragment_O() {};
};



SMART(CDText );
class CDText_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,CDText_O,"CDText",core::CxxObject_O);
 public:
    static CDText_sp make(core::HashTableEq_sp kprops);
private:
	string			_Text;
	core::HashTableEq_sp		_Properties;
public:
	void	initialize();
 public:
#if 0
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
#endif
public:

	void	parseFromXml(adapt::QDomNode_sp xml, bool print);

	bool	hasProperties();
	core::HashTableEq_sp getProperties() { return this->_Properties; };

	CDText_O( const CDText_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CDText_O);
};




SMART(ChemDraw );
class ChemDraw_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,ChemDraw_O,"ChemDraw",core::CxxObject_O);
 public:
    static ChemDraw_sp make(core::T_sp stream, bool print=false);
 public:
    static void lisp_initGlobals(core::Lisp_sp lisp);
public:
    typedef	gctools::Vec0<CDFragment_sp>	Fragments;
	typedef adapt::SymbolMap<CDFragment_O>	NamedFragments;
private:
	Fragments	_AllFragments;
	NamedFragments	_NamedFragments;
public:
//	void	archive(core::ArchiveP node);
	void	initialize();
 public:
#if 0
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
#endif
private:
	string	_getAtomName(adapt::QDomNode_sp node);
public:
        void parse(core::T_sp strm, bool print);
        void parseChild(adapt::QDomNode_sp child, bool print);

	/*! Set the properties for the named fragment.
	  @param fragmentName The name of the fragment whose properties are being set
	  @param properties is a property list (keyword symbol/object pairs) */
        void setFragmentProperties(core::Symbol_sp name
                                       , core::T_sp comment
                                       , core::T_sp chiral_centers
                                       , core::T_sp group
                                       , core::T_sp name_template
                                       , core::T_sp pdb_template
                                       , core::T_sp restraints
                                       , core::T_sp residue_charge
                                       , core::T_sp restrained_pi_bonds
                                       , core::T_sp caps);



	core::List_sp getFragments();
CL_LISPIFY_NAME("allFragmentsAsList");
CL_DEFMETHOD 	core::List_sp allFragmentsAsList() { return this->getFragments();};
	core::List_sp getSubSetOfFragments(adapt::SymbolSet_sp subsetNames );


    /*! Return the entire ChemDraw file as an Aggregate and each structure
      is a Molecule with the name given by the hashed-bond-:name attribute
    */
    Aggregate_sp asAggregate();


	ChemDraw_O( const ChemDraw_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ChemDraw_O);
};


};
#endif //]
