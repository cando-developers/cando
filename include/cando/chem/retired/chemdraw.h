#ifndef	ChemDraw_H //[
#define ChemDraw_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/holder.h>
#include <clasp/core/symbolSet.fwd.h>
#include <cando/chem/bond.fwd.h>
#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/chem/chemPackage.h>
#include "core/environment.fwd.h"// chemdraw.h wants Binder needs environment.fwd.h
#include "candoBase/quickDom.fwd.h"// chemdraw.h wants QDomNode needs quickDom.fwd.h
#include "core/lisp.fwd.h"// chemdraw.h wants Lisp needs lisp.fwd.h
#include <cando/chem/atom.fwd.h>
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
		unknownCDBond
    } CDBondOrder;


SMART(CDNode );
class CDNode_O : public core::T_O
{	
    friend class CDFragment_O;
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,CDNode_O,"CDNode");
private:
	uint					_Id;
	string					_Label;
	StereochemistryType			_StereochemistryType;
	ConfigurationEnum			_Configuration;
private:	// generated
	Atom_sp					_Atom;
		// Keep track of neighbors for spanning trees
		// 
	vector< pair<CDNode_wp,CDBondOrder> >	_Neighbors;
	CDNode_wp				_BackSpan;
	CDNode_wp				_NextSpan;
public:
	void	initialize();
private:
	string	_extractLabel(candoBase::QDomNode_sp node);
	void	unidirectionalBondTo(CDNode_sp neighbor, CDBondOrder order);
public:
	uint	getId() const { return this->_Id;};
	void    setId(uint i) { this->_Id = i; };
	string getLabel() const { return this->_Label;};
	/*! Parse the label into a name and ionization component "name/[+-] */
	void getParsedLabel(string& name, int& ionization) const;
	void setLabel(const string& l) { this->_Label = l;};
	Atom_sp	getAtom() const { return this->_Atom;};
	void setAtom(Atom_sp a) { this->_Atom = a;};
	void	parseFromXml(candoBase::QDomNode_sp xml);

	void	bondTo( CDNode_sp neighbor, CDBondOrder order);

	CDNode_O( const CDNode_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CDNode_O);
};


SMART(CDBond );
class CDBond_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,CDBond_O,"CDBond");
private:
	uint		_IdBegin;
	uint		_IdEnd;
	CDNode_wp	_BeginNode;
	CDNode_wp	_EndNode;
	CDBondOrder	_Order;
public:
	void	initialize();
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
	CDNode_sp getBeginNode() {_OF(); ANN(this->_BeginNode);return this->_BeginNode.lock();};
	void	setEndNode(CDNode_sp n) { this->_EndNode = n;};
	CDNode_sp getEndNode() {_OF(); ANN(this->_EndNode);return this->_EndNode.lock();};

	void	parseFromXml(candoBase::QDomNode_sp xml);

	CDBond_O( const CDBond_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CDBond_O);
};



SMART(CDFragment );
class CDFragment_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,CDFragment_O,"CDFragment");
public:
    typedef gctools::Vec0<CDBond_sp>	CDBonds;
    typedef gctools::SmallMap<Atom_sp,CDNode_sp> AtomsToBonds; // m a p<Atom_sp,CDNode_sp>	AtomsToBonds;
    typedef gctools::SmallMap<int,CDNode_sp> IntMappedCDNodes; // typedef m a p<int,CDNode_sp>		IntMappedCDNodes;
private:
    core::Symbol_sp		_ConstitutionName;
    IntMappedCDNodes		_Nodes;
    AtomsToBonds		_AtomsToNodes;
    CDBonds			_Bonds;
    CDNode_sp			_RootNode;
    int				_LargestId;
    core::HashTableEq_sp		_Properties;
public:
    void	initialize();
private:
    bool _asKeyedObject(const string& label, core::Symbol_sp& keyword, core::T_sp& obj);
private:
    Residue_sp _buildResidue(bool constitutionOnly);

public:
    void setConstitutionName(core::Symbol_sp n) { this->_ConstitutionName = n;};
    core::Symbol_sp getConstitutionName() { return this->_ConstitutionName;};

    core::HashTableEq_sp getProperties() { return this->_Properties;};
    void	addProperties(core::HashTableEq_sp dict);
    Atom_sp createOneAtom(CDNode_sp n);
    void	createAtoms();
    void	createImplicitHydrogen(CDNode_sp from, const string& name);
    void	parseFromXml(candoBase::QDomNode_sp xml);
    /*! Return false if the fragment couldn't be interpreted */
    bool	interpret();

    //! Get the property for this ChemDraw Fragment
    core::T_sp getProperty(core::Symbol_sp key);
    core::T_sp getPropertyOrDefault(core::Symbol_sp key, core::T_sp df);
    bool hasProperty(core::Symbol_sp key);
    core::T_sp setProperty(core::Symbol_sp key, core::T_sp df);
    string describeProperties();

    void createBonds(bool selectedAtomsOnly);
    void removeAllBonds();
    void clearAtomSelected();


    Residue_sp createResidueOfSelectedAtoms();

    Residue_sp getEntireResidue();
    ConstitutionAtoms_sp asConstitutionAtoms();

    string __repr__() const;

    CDFragment_O( const CDFragment_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(CDFragment_O);
};



SMART(CDText );
class CDText_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,CDText_O,"CDText");
#if INIT_TO_FACTORIES
 public:
    static CDText_sp make(core::HashTableEq_sp kprops);
#else
    DECLARE_INIT();
#endif
private:
	string			_Text;
	core::HashTableEq_sp		_Properties;
public:
	void	initialize();
public:

	void	parseFromXml(candoBase::QDomNode_sp xml, core::Lisp_sp);

	bool	hasProperties();
	core::HashTableEq_sp getProperties() { return this->_Properties; };

	CDText_O( const CDText_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CDText_O);
};




SMART(ChemDraw );
class ChemDraw_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ChemDraw_O,"ChemDraw");
#if INIT_TO_FACTORIES
 public:
    static ChemDraw_sp make(const string& fileName);
#else
    DECLARE_INIT();
#endif
    public:
	static void lisp_initGlobals(core::Lisp_sp lisp);
public:
    typedef	gctools::Vec0<CDFragment_sp>	Fragments;
	typedef core::SymbolMap<CDFragment_O>	NamedFragments;
private:
	Fragments	_AllFragments;
	NamedFragments	_NamedFragments;
public:
//	void	archive(core::ArchiveP node);
	void	initialize();
private:
	// instance variables
public:
    static ChemDraw_sp open(const string& fileName);
private:
	string	_getAtomName(candoBase::QDomNode_sp node);
public:
	void	parseFromFileName(const string& fileName);

	/*! Set the properties for the named fragment.
	  @param fragmentName The name of the fragment whose properties are being set
	  @param properties is a property list (keyword symbol/object pairs) */
	void setFragmentProperties( core::Symbol_sp fragment, core::Cons_sp properties);


	core::Cons_sp getFragments();
	core::Cons_sp allFragmentsAsCons() { return this->getFragments();};
	core::Cons_sp getSubSetOfFragments(core::SymbolSet_sp subsetNames );


    /*! Return the entire ChemDraw file as an Aggregate and each structure
      is a Molecule with the name given by the hashed-bond-:name attribute
    */
    Aggregate_sp asAggregate();


	ChemDraw_O( const ChemDraw_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ChemDraw_O);
};


};
TRANSLATE(chem::CDNode_O);
TRANSLATE(chem::CDBond_O);
TRANSLATE(chem::CDFragment_O);
TRANSLATE(chem::CDText_O);
TRANSLATE(chem::ChemDraw_O);
#endif //]
