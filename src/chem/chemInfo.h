       
       

//
//	chemInfo.h
//
//
#ifndef	CHEM_INFO_H
#define	CHEM_INFO_H

#include <iostream>
#include <string>

#include "core/common.h"
#include "core/stringSet.fwd.h"
#include "atom.h"
#include "core/holder.h"


#include "chemPackage.h"

namespace chem {

    SMART(BoundFrame);

    extern string	ChemInfo_XmlName;

    SMART(Root);
    SMART(AntechamberRoot);
    class	ChemInfoNode_O;
    typedef	gctools::smart_ptr<ChemInfoNode_O>	ChemInfoNode_sp;



    SMART(ChemInfoMatch);
    class ChemInfoMatch_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,ChemInfoMatch_O,"ChemInfoMatch");

    public:
	void initialize();
    public:
	void	archiveBase(core::ArchiveP node);
    private:
	bool		_Matches;
        core::HashTableEq_sp _TagLookup; // core::StringMap<Atom_O>	_TagLookup;
        core::HashTableEq_sp _ClosestMatch; // core::StringMap<Atom_O>	_ClosestMatch;
    public:
	string __repr__() const;
	bool matches();
	void setMatches(bool b) { this->_Matches = b;};
	void clearAtomTags();
	bool	recognizesAtomTag(core::Symbol_sp tag);
	void	defineAtomTag(Atom_sp a, core::Symbol_sp tag );
	bool hasAtomWithTag(core::Symbol_sp tag );
	chem::Atom_sp getAtomWithTag(core::Symbol_sp tag);
	chem::Atom_sp getAtomWithTagOrNil(core::Symbol_sp tag);
	chem::Atom_sp tag(core::Symbol_sp tag) { return this->getAtomWithTag(tag);};
	void forgetAtomTag(core::Symbol_sp tag);

	/*! Make sure all the matching atoms are unique - otherwise throw an exception */
	void throwIfInvalid();
	void describeClosestMatch();

	/*! Return a BoundFrame using atoms with tag "1","2","3" as the origin, x axis and xy plane
	 */
	BoundFrame_sp boundFrame();

	DEFAULT_CTOR_DTOR(ChemInfoMatch_O);
    };





    SMART(WildElementDict);
    class WildElementDict_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,WildElementDict_O,"WildElementDict");

    public:
	void initialize();
    public:
	void	archiveBase(core::ArchiveP node);
    private:
        core::HashTableEqual_sp _AtomWildCards; // core::StringMap<core::StringSet_O>	_AtomWildCards;
    public:
        void    addWildName(core::Symbol_sp  wildName);
        void    addWildNameMap(core::Symbol_sp wildName, core::Symbol_sp elementName );
        bool    lexWildNameMatches2Char(char wildNameChar1, char wildNameChar2);
        bool    lexWildNameMatches1Char(char wildName);
        bool    recognizesWildName(core::Symbol_sp name);
        bool    recognizesWildNameElement(core::Symbol_sp wildName, core::Symbol_sp element);
//        addon::QDomNode_sp      asXml(string name=XmlTag_WildElementDict());
//	void	parseFromXml(addon::QDomNode_sp node);

	DEFAULT_CTOR_DTOR(WildElementDict_O);
    };


    typedef	enum { noType, recursiveRoot, root, smartsRoot, antechamberRoot, chain, branch,
		       logical, tagSet, ringTest,
		       atomTest, antechamberBondTest, bondTest, afterMatchBondTest,
		       antechamberFocusAtomMatch } ChemInfoType;


    string  chemInfoTypeString(ChemInfoType type);
    ChemInfoType chemInfoTypeFromString(const string& name);


    typedef enum	{
	SABNoBond,
	SABSingleBond,
	SABSingleOrAromaticBond,
	SABDoubleOrAromaticBond,
	SABDelocalizedBond,
	SABDoubleBond,
	SABTripleBond,
	SABAromaticBond,
	SABAnyBond,
	SABDirectionalSingleUpOrUnspecified,
	SABDirectionalSingleDownOrUnspecified,
	SABDirectionalSingleUp,
	SABDirectionalSingleDown,
    } BondEnum;

    extern core::Symbol_sp _sym_STARSabBondEnumConverterSTAR;

};

ENUM_TRANSLATOR(chem::BondEnum,chem::_sym_STARSabBondEnumConverterSTAR);



namespace chem {

/*! Base class for all fragment matching
 */
    SMART(ChemInfoNode);
    class ChemInfoNode_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,ChemInfoNode_O,"ChemInfoNode");

    public:
	void	archiveBase(core::ArchiveP node);
    public:
	virtual uint depth() const;
	virtual string asSmarts() const {_OF();SUBCLASS_MUST_IMPLEMENT();};
	DEFAULT_CTOR_DTOR(ChemInfoNode_O);
    };


    SMART(ResidueList);
    class ResidueList_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,ResidueList_O,"ResidueList");

    public:
	void	archiveBase(core::ArchiveP node);
    private:
    public:
	DEFAULT_CTOR_DTOR(ResidueList_O);
    };





    SMART(RootMatchNode);
    class RootMatchNode_O : public ChemInfoNode_O
    {
	LISP_BASE1(ChemInfoNode_O);
	LISP_CLASS(chem,ChemPkg,RootMatchNode_O,"RootMatchNode");

    public:
	void	archiveBase(core::ArchiveP node);
    public:
	virtual	bool	matches(Root_sp root) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
	
	DEFAULT_CTOR_DTOR(RootMatchNode_O);
    };





    SMART(BondMatchNode);
    class BondMatchNode_O : public ChemInfoNode_O
    {
	LISP_BASE1(ChemInfoNode_O);
	LISP_CLASS(chem,ChemPkg,BondMatchNode_O,"BondMatchNode");

    public:
	void	archiveBase(core::ArchiveP node);
    public:
	virtual	bool	matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond ) {_OF();SUBCLASS_MUST_IMPLEMENT(); };

	DEFAULT_CTOR_DTOR(BondMatchNode_O);
    };






    SMART(AtomOrBondMatchNode);
    class AtomOrBondMatchNode_O : public BondMatchNode_O
    {
	LISP_BASE1(BondMatchNode_O);
	LISP_CLASS(chem,ChemPkg,AtomOrBondMatchNode_O,"AtomOrBondMatchNode");

    public:
	void	archiveBase(core::ArchiveP node);
    private:
    public:


	virtual bool    matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
	virtual	bool	matches(Root_sp root, chem::Atom_sp atom) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
	DEFAULT_CTOR_DTOR(AtomOrBondMatchNode_O);
    };



    SMART(BondListMatchNode);
    class BondListMatchNode_O : public ChemInfoNode_O
    {
	LISP_BASE1(ChemInfoNode_O);
	LISP_CLASS(chem,ChemPkg,BondListMatchNode_O,"BondListMatchNode");

    public:
	void	archiveBase(core::ArchiveP node);
    private:
    public:
	virtual	bool	matches(Root_sp root, chem::Atom_sp from, chem::BondList_sp bondList ) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
	virtual string asSmarts() const;
	DEFAULT_CTOR_DTOR(BondListMatchNode_O);
    };



    typedef	enum { logIdentity, logNot, logHighPrecedenceAnd, logOr, logLowPrecedenceAnd } LogicalOperatorType;

    extern core::Symbol_sp _sym_STARLogicalOperatorTypeConverterSTAR;
};




namespace translate {
    template <> struct to_object<chem::LogicalOperatorType>
    {								
	typedef	chem::LogicalOperatorType	GivenType;	
	static core::T_sp convert(const GivenType& val)
	{_G();
	    core::SymbolToEnumConverter_sp converter = chem::_sym_STARLogicalOperatorTypeConverterSTAR->symbolValue().as<core::SymbolToEnumConverter_O>();
	    return (converter->symbolForEnum(val));
	}
    };

    template <>
    struct from_object<chem::LogicalOperatorType>
    {								
	typedef	chem::LogicalOperatorType 	ExpectedType;
	typedef	chem::LogicalOperatorType 	DeclareType;
	DeclareType _v;
	from_object(gctools::smart_ptr<core::T_O> o) {
	    core::SymbolToEnumConverter_sp converter = chem::_sym_STARLogicalOperatorTypeConverterSTAR->symbolValue().as<core::SymbolToEnumConverter_O>();
	    _v = converter->enumForSymbol<chem::LogicalOperatorType>(o.as<core::Symbol_O>());
	}
    };


};


namespace chem {



class	AtomTest;



SMART(Logical);
class Logical_O : public AtomOrBondMatchNode_O
{
    LISP_BASE1(AtomOrBondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,Logical_O,"Logical");

public:
    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
	LogicalOperatorType	_Operator;
	AtomOrBondMatchNode_sp	_Left;
	AtomOrBondMatchNode_sp	_Right;

public:
		// Second argument can be NULL
	static Logical_sp create(LogicalOperatorType op, AtomOrBondMatchNode_sp a1, AtomOrBondMatchNode_sp a2 )
	{ _G();
	    GC_ALLOCATE(Logical_O, obj ); // RP_Create<Logical_O>(lisp);
	    obj->_Operator = op;
	    obj->_Left = a1;
	    obj->_Right = a2;
	    return obj;
	};
	static Logical_sp create( LogicalOperatorType op, AtomOrBondMatchNode_sp a1)
	{_G();
	    AtomOrBondMatchNode_sp a2;
	    a2 = _Nil<AtomOrBondMatchNode_O>();
	    return create(op,a1,a2);
	};
public:


//virtual	addon::QDomNode_sp	asXml(string name=XmlTag_Logical());
//virtual	void	parseFromXml(addon::QDomNode_sp node);
virtual	ChemInfoType	type() { return logical; };
virtual	bool		matches( Root_sp root, chem::Atom_sp atom );
virtual	bool		matches( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
virtual uint depth() const;
virtual string asSmarts() const;
DEFAULT_CTOR_DTOR(Logical_O);
};

SMART(TagSet);
    class TagSet_O : public AtomOrBondMatchNode_O
{
    LISP_BASE1(AtomOrBondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,TagSet_O,"TagSet");

public:
    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
//        BondEnum        	_Bond;
    core::Symbol_sp		_RingTag;
    AtomOrBondMatchNode_sp	_AtomTest;
public:
    static TagSet_sp create(/* BondEnum b, */ AtomOrBondMatchNode_sp at, core::Symbol_sp ri)
	{ _G();
	    LOG(BF("TagSet_sp create: ringTag = (%s)")%ri);
	    GC_ALLOCATE(TagSet_O, obj ); // RP_Create<TagSet_O>(lisp);
//	    obj->_Bond = b;
	    obj->_AtomTest = at;
	    obj->_RingTag = ri;
	    return obj;
	};
public:

virtual	ChemInfoType	type() { return tagSet;};
virtual	bool		matches(Root_sp root, chem::Atom_sp atom );
virtual string asSmarts() const;
DEFAULT_CTOR_DTOR(TagSet_O);
};



SMART(RingTest);
class RingTest_O : public AtomOrBondMatchNode_O
{
    LISP_BASE1(AtomOrBondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,RingTest_O,"RingTest");

public:
    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
        BondEnum        	_Bond;
	core::Symbol_sp 		_RingTag;
	AtomOrBondMatchNode_sp	_AtomTest;
public:
	static RingTest_sp create( BondEnum b, AtomOrBondMatchNode_sp at, core::Symbol_sp ri)
	{
	    _G();
	    LOG(BF("RingTest_sp create: ringTag = (%s)") % ri );
	    GC_ALLOCATE(RingTest_O, obj ); // RP_Create<RingTest_O>(lisp);
	    obj->_Bond = b;
	    obj->_AtomTest = at;
	    obj->_RingTag = ri;
	    return obj;
	};
public:

virtual	ChemInfoType	type() { return ringTest;};
virtual	bool		matches(Root_sp root, chem::Atom_sp atom );

DEFAULT_CTOR_DTOR(RingTest_O);
};



SMART(ResidueTest);
class ResidueTest_O : public AtomOrBondMatchNode_O
{
    LISP_BASE1(AtomOrBondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,ResidueTest_O,"ResidueTest");

public:
    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
        BondEnum        	_Bond;
	core::Symbol_sp			_RingTag;
	AtomOrBondMatchNode_sp	_AtomTest;
public:
	static ResidueTest_sp create( BondEnum b, AtomOrBondMatchNode_sp at, core::Symbol_sp ri)
	{
	    _G();
	    LOG(BF("ResidueTest_sp create: ringTag = (%s)")%ri);
	    GC_ALLOCATE(ResidueTest_O, obj ); // RP_Create<ResidueTest_O>(lisp);
	    obj->_Bond = b;
	    obj->_AtomTest = at;
	    obj->_RingTag = ri;
	    return obj;
	};
public:

virtual	ChemInfoType	type() { return ringTest;};
virtual	bool		matches(Root_sp root, chem::Atom_sp atom );
virtual string asSmarts() const;
DEFAULT_CTOR_DTOR(ResidueTest_O);
};




bool    _matchBondTypes(BondEnum be, chem::BondOrder bo);




SMART(BondTest);
class BondTest_O : public BondMatchNode_O
{
    LISP_BASE1(BondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,BondTest_O,"BondTest");

public:
	void	initialize();
	void	archiveBase(core::ArchiveP node);
private:
	BondEnum		_Bond;
	AtomOrBondMatchNode_sp	_AtomTest;

public:
	static BondTest_sp create(BondEnum b, AtomOrBondMatchNode_sp node )
	{_G();
	    GC_ALLOCATE(BondTest_O, obj ); // RP_Create<BondTest_O>(lisp);
	    obj->_Bond = b;
	    obj->_AtomTest = node;
	    return obj;
	}


public:

	BondEnum	bondType() { return this->_Bond; };

// virtual	addon::QDomNode_sp	asXml(string name=XmlTag_BondTest());
virtual	ChemInfoType	type() { return bondTest;};
virtual	bool	matches( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
virtual string asSmarts() const;

DEFAULT_CTOR_DTOR(BondTest_O);
};



SMART(AntechamberBondTest);
class AntechamberBondTest_O : public BondMatchNode_O
{
    LISP_BASE1(BondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,AntechamberBondTest_O,"AntechamberBondTest");

public:
	void initialize();
	void	archiveBase(core::ArchiveP node);
private:
    core::Symbol_sp                   _Element;
        int                     _Neighbors;
	AtomOrBondMatchNode_sp	_AtomProperties;
    core::Symbol_sp                   _Tag;
public:
	static AntechamberBondTest_sp create( core::Symbol_sp element, int neighbors,
				AtomOrBondMatchNode_sp props, core::Symbol_sp tag )
	{_G();
	    GC_ALLOCATE(AntechamberBondTest_O, obj ); // RP_Create<AntechamberBondTest_O>(lisp);
            obj->_Element = element;
            obj->_Neighbors = neighbors;
            obj->_AtomProperties = props;
            obj->_Tag = tag;
	    return obj;
	};
public:

//virtual	addon::QDomNode_sp	asXml(string name=XmlTag_AntechamberBondTest());
//virtual	void	parseFromXml(addon::QDomNode_sp node );
virtual	ChemInfoType	type() { return antechamberBondTest;};
        bool	matchBasic( AntechamberRoot_sp root, chem::Atom_sp atom );
	virtual	bool	matches( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
virtual bool    matches( Root_sp root, chem::Atom_sp atom );
virtual string asSmarts() const;


DEFAULT_CTOR_DTOR(AntechamberBondTest_O);
};





typedef	enum	{
	SAPNone		=0,
	SAPWildCard		,
	SAPDegree		,
	SAPElement		,
	SAPTotalHCount		,
	SAPImplicitHCount	,
	SAPRingTest		,
	SAPRingMembershipCount	,
	SAPRingSize		,
	SAPValence		,
	SAPConnectivity		,
	SAPNegativeCharge	,
	SAPNegativeFormalCharge	,
	SAPPositiveCharge	,
	SAPPositiveFormalCharge	,
	SAPAtomicNumber		,
	SAPChiralityAntiClockwise,
	SAPChiralityClockwise	,
	SAPAtomicMass		,
	SAPLonePair		,
	SAPTotalBondNumber	,
	SAPTotalExplicitBondNumber,
	SAPPiBondOrbital	,
	SAPAromaticPiElectron	,
	SAPHeavyAtomTotalBond	,
	SAPGroupNumber		,
	SAPElectronegativeElement,
	SAPTransitionMetal	,
        SAPBondedToPrevious	,
        SAPNotBondedToPrevious	,
        SAPInBond		,
        SAPArLevel		,
        SAPNoRing		,
	SAPResidueTest		,
	SAPAM1_BCC_x		,
	SAPAM1_BCC_y		,
	SAPLambda		,
	SAPAromaticElement	,
	SAPAliphatic		,
	SAPAromatic,
} AtomTestEnum;

extern core::Symbol_sp _sym_STARAtomTestEnumConverterSTAR;

};



namespace translate {
    template <> struct to_object<chem::AtomTestEnum>
    {								
	typedef	chem::AtomTestEnum	GivenType;	
	static core::T_sp convert(const GivenType& val)
	{_G();
	    core::SymbolToEnumConverter_sp converter = chem::_sym_STARAtomTestEnumConverterSTAR->symbolValue().as<core::SymbolToEnumConverter_O>();
	    return (converter->symbolForEnum(val));
	}
    };

    template <>
    struct from_object<chem::AtomTestEnum>
    {								
	typedef	chem::AtomTestEnum 	ExpectedType;
	typedef	chem::AtomTestEnum 	DeclareType;
	DeclareType _v;
	from_object(gctools::smart_ptr<core::T_O> o) {
	    core::SymbolToEnumConverter_sp converter = chem::_sym_STARAtomTestEnumConverterSTAR->symbolValue().as<core::SymbolToEnumConverter_O>();
	    _v = converter->enumForSymbol<chem::AtomTestEnum>(o.as<core::Symbol_O>());
	}
    };


};






namespace chem {
SMART(AtomTest);
class AtomTest_O : public AtomOrBondMatchNode_O
{
    LISP_BASE1(AtomOrBondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,AtomTest_O,"AtomTest");

public:
void initialize();
    void	archiveBase(core::ArchiveP node);
private:
	AtomTestEnum	_Test;
	int		_IntArg;
	int		_NumArg;
	string		_StringArg;
	core::Symbol_sp	_SymbolArg;
public:
	static AtomTest_sp create( AtomTestEnum t, int iArg, int num, const char* el, core::Symbol_sp sym) 
	{_G();
		LOG(BF("Creating an AtomTest"));
		GC_ALLOCATE(AtomTest_O, obj ); // RP_Create<AtomTest_O>(lisp);
		obj->_Test = t;
		obj->_IntArg = iArg;
		obj->_NumArg = num;
		if ( el != NULL ) {
		    obj->_StringArg = el;
		} else obj->_StringArg = "";
		obj->_SymbolArg = sym;
		return obj;
	};
	static AtomTest_sp create( AtomTestEnum t) 
	{_G();
	    return create( t, 0, 0, NULL, _Nil<core::Symbol_O>());
	};
	static AtomTest_sp create( AtomTestEnum t, int iArg ) 
	{_G();
	    return create(  t, iArg, 0, NULL,_Nil<core::Symbol_O>() );
	};
	static AtomTest_sp create( AtomTestEnum t, int iArg, int num ) 
	{_G();
	    return create( t, iArg, num, NULL,_Nil<core::Symbol_O>() );
	};
	static AtomTest_sp create( AtomTestEnum t, const char* el ) 
	{_G();
	    LOG(BF("Create an AtomTest with string argument: %s")%el);
	    return create( t, 0, 0, el,_Nil<core::Symbol_O>());
	};
	static AtomTest_sp create( AtomTestEnum t, core::Symbol_sp sym)
	{_G();
	    LOG(BF("Create an AtomTest with a symbol argument: %s") % sym->__repr__() );
	    return create( t, 0, 0, NULL, sym );
	}

public:
	AtomTestEnum	myType() { return this->_Test; };
	int		getIntArg() { return this->_IntArg; };
	string		testName(AtomTestEnum t) const;
	bool matchesAm1BccX(chem::Atom_sp atom) const;
	bool matchesAm1BccY(chem::Atom_sp atom) const;

virtual ChemInfoType    type() { return atomTest;};
virtual	bool	matches( Root_sp root, chem::Atom_sp atom );
virtual	bool	matches( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond);
virtual string asSmarts() const;

DEFAULT_CTOR_DTOR(AtomTest_O);
};





/*!     This is for matching the focus atom for GAFF type assignment rules
 */
SMART(AntechamberFocusAtomMatch);
class AntechamberFocusAtomMatch_O : public AtomOrBondMatchNode_O
{
    LISP_BASE1(AtomOrBondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,AntechamberFocusAtomMatch_O,"AntechamberFocusAtomMatch");

public:
	void initialize();
public:
	void	archiveBase(core::ArchiveP node);
private:
        ResidueList_sp     _ResidueNames;
        int             _AtomicNumber;
        int             _NumberOfAttachedAtoms;
        int             _NumberOfAttachedHydrogens;
        int             _NumberOfElectronWithdrawingGroups;
        AtomOrBondMatchNode_sp   _AtomicProperty;
public:
	static AntechamberFocusAtomMatch_sp create(
                        ResidueList_sp residueNames,
                        int     atomicNumber,
                        int     numberOfAttachedAtoms,
                        int     numberOfAttachedHydrogens,
                        int     numberOfElectronWithdrawingGroups,
                        AtomOrBondMatchNode_sp atomicProperty ) 
	{_G();
	    GC_ALLOCATE(AntechamberFocusAtomMatch_O, obj ); // RP_Create<AntechamberFocusAtomMatch_O>(lisp);
	    obj->_ResidueNames = residueNames;
	    obj->_AtomicNumber = atomicNumber;
	    obj->_NumberOfAttachedAtoms = numberOfAttachedAtoms;
	    obj->_NumberOfAttachedHydrogens = numberOfAttachedHydrogens;
	    obj->_NumberOfElectronWithdrawingGroups =numberOfElectronWithdrawingGroups;
	    obj->_AtomicProperty = atomicProperty;
	    LOG(BF("Created AntechamberFocusAtomMatch@%0X") % obj.get());
	    return obj;
	};
	static AntechamberFocusAtomMatch_sp create(
                        ResidueList_sp residueNames,
                        int     atomicNumber,
                        int     numberOfAttachedAtoms,
                        int     numberOfAttachedHydrogens,
                        int     numberOfElectronWithdrawingGroups )
	{_G();
	    AtomOrBondMatchNode_sp atomicProperty;
	    atomicProperty = _Nil<AtomOrBondMatchNode_O>();
	    return create( residueNames, atomicNumber, numberOfAttachedAtoms, numberOfAttachedHydrogens, numberOfElectronWithdrawingGroups, atomicProperty);
	};
public:


//virtual	addon::QDomNode_sp	asXml(string name=XmlTag_AntechamberFocusAtomMatch());
//virtual	void	parseFromXml(addon::QDomNode_sp node);
virtual	bool	matches( Root_sp root, chem::Atom_sp atom );
virtual ChemInfoType    type() { return antechamberFocusAtomMatch; };

DEFAULT_CTOR_DTOR(AntechamberFocusAtomMatch_O);
};







SMART(Chain);
class Chain_O : public BondListMatchNode_O
{
    LISP_BASE1(BondListMatchNode_O);
    LISP_CLASS(chem,ChemPkg,Chain_O,"Chain");

public:
    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
	BondMatchNode_sp		_Head;
	BondListMatchNode_sp	_Tail;
public:

	virtual uint depth() const;
virtual string asSmarts() const;
		//! Second argument can be NULL
	static Chain_sp create( BondMatchNode_sp head, BondListMatchNode_sp tail ) {
	    _G();
	    GC_ALLOCATE(Chain_O, obj ); // RP_Create<Chain_O>(lisp);
		obj->_Head = head;
		obj->_Tail = tail;
		return obj;
	}
	static Chain_sp create( BondMatchNode_sp head )
	{
	    _G();
	    BondListMatchNode_sp tail = _Nil<BondListMatchNode_O>();
	    return create(head,tail);
	}


public:


virtual	ChemInfoType	type() { return chain; };
virtual	bool		matches( Root_sp root, chem::Atom_sp from, chem::BondList_sp neighbors );

DEFAULT_CTOR_DTOR(Chain_O);
};


SMART(Branch);
class Branch_O : public BondListMatchNode_O
{
    LISP_BASE1(BondListMatchNode_O);
    LISP_CLASS(chem,ChemPkg,Branch_O,"Branch");

public:
    void	archiveBase(core::ArchiveP node);
	void	initialize();
private:
	BondListMatchNode_sp	_Left;
	BondListMatchNode_sp	_Right;
public:
	virtual uint depth() const;
virtual string asSmarts() const;
		//! Second argument can be NULL
	static Branch_sp create( BondListMatchNode_sp left, BondListMatchNode_sp right ) 
	{
	    _G();
	    GC_ALLOCATE(Branch_O, obj ); // RP_Create<Branch_O>(lisp);
		obj->_Left = left;
		obj->_Right = right;
		return obj;
	}

	static Branch_sp create( BondListMatchNode_sp left ) {
	    _G();
	    GC_ALLOCATE(Branch_O, obj ); // RP_Create<Branch_O>(lisp);
		obj->_Left = left;
		obj->_Right = _Nil<BondListMatchNode_O>();
		return obj;
	}

public:



virtual	ChemInfoType	type() { return branch; };
virtual	bool		matches( Root_sp root, chem::Atom_sp from, chem::BondList_sp neighbors );  /* CHECK CODE TO SEE IF IT HANDLES RIGHT=NULL */
//virtual	addon::QDomNode_sp	asXml(string name=XmlTag_Branch());
//virtual	void	parseFromXml(addon::QDomNode_sp node);

DEFAULT_CTOR_DTOR(Branch_O);
};

#if 0 //[
typedef	gctools::smart_ptr<O_AfterMatchTestNode>	RPAfterMatchTestNode;
class	AfterMatchTestNode : ChemInfoNode_O
{
private:
public:
	AfterMatchTestNode();
//virtual	addon::QDomNode_sp	asXml(string name) = 0;
	virtual	bool  		matches( Root_sp root ) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
	virtual	~AfterMatchTestNode();
};
#endif //]




SMART(AfterMatchBondTest);
class AfterMatchBondTest_O : public RootMatchNode_O
{
    LISP_BASE1(RootMatchNode_O);
    LISP_CLASS(chem,ChemPkg,AfterMatchBondTest_O,"AfterMatchBondTest");

public:
void initialize();
    void	archiveBase(core::ArchiveP node);
private:
    core::Symbol_sp _AtomTag1;
    core::Symbol_sp _AtomTag2;
    BondEnum	_Bond;
public:
	static AfterMatchBondTest_sp create( core::Symbol_sp tag1, core::Symbol_sp tag2, BondEnum b)
	{
	    _G();
	    GC_ALLOCATE(AfterMatchBondTest_O, obj ); // RP_Create<AfterMatchBondTest_O>(lisp);
                obj->_AtomTag1 = tag1;
                obj->_AtomTag2 = tag2;
		obj->_Bond = b;
		return obj;
	};
public:
// virtual	addon::QDomNode_sp	        asXml(string name=XmlTag_AfterMatchBondTest());
virtual	ChemInfoType            type() { return afterMatchBondTest; };
virtual	bool    		matches( Root_sp root );


DEFAULT_CTOR_DTOR(AfterMatchBondTest_O);
};







SMART(Root);
class Root_O : public AtomOrBondMatchNode_O
{
    LISP_BASE1(AtomOrBondMatchNode_O);
    LISP_CLASS(chem,ChemPkg,Root_O,"Root");

public:
	void initialize();
public:
	void	archiveBase(core::ArchiveP node);
protected:
	AtomOrBondMatchNode_sp	_FirstTest;
	BondListMatchNode_sp	_Chain;
	/*! Store tests described as lambdas that take a single atom argument and return a boolean.
	 These can be incorporated into the smarts code as <xxxx> where xxxx is the symbol
	name of the test. */
	core::HashTableEq_sp		_Tests;
protected:	// do not archive
	ChemInfoMatch_sp		_Match; //!< Match SMARTS ring ids and Antechamber Atom Tags
public:
	virtual uint depth() const;
virtual string asSmarts() const;
	static Root_sp create(AtomOrBondMatchNode_sp node, BondListMatchNode_sp chain )
	{
	    _G();
	    GC_ALLOCATE(Root_O, obj ); // RP_Create<Root_O>(lisp);
	    obj->_FirstTest = node;
	    obj->_Chain = chain;
	    ANN(obj->_Chain);
	    return obj;
	};

	void addTest(core::Symbol_sp testSymbol, core::Function_sp testCode);
	bool evaluateTest(core::Symbol_sp testSym, Atom_sp atom);

    void createNewMatch() { this->_Match = ChemInfoMatch_O::create();};
	void setMatch(ChemInfoMatch_sp match) { this->_Match = match; };
	ChemInfoMatch_sp getMatch() { return this->_Match; };

    virtual core::Symbol_sp getAssignType() { return _Nil<core::Symbol_O>();};
public:


//virtual	addon::QDomNode_sp	asXml(string name=XmlTag_Root());
//virtual	void	parseFromXml(addon::QDomNode_sp node);
virtual	ChemInfoType	type() { return root; };
virtual	bool		matches( Root_sp root, chem::Atom_sp atom );
virtual	bool		matches( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );


DEFAULT_CTOR_DTOR(Root_O);
};







SMART(SmartsRoot);
class SmartsRoot_O : public Root_O
{
    LISP_BASE1(Root_O);
    LISP_CLASS(chem,ChemPkg,SmartsRoot_O,"SmartsRoot");

public:
    void	archiveBase(core::ArchiveP node);
public:

    static SmartsRoot_sp create(AtomOrBondMatchNode_sp node,
    					BondListMatchNode_sp chain)
    { _G();
	GC_ALLOCATE(SmartsRoot_O,   obj ); // RP_Create<SmartsRoot_O>(lisp);
		obj->_FirstTest = node;
		obj->_Chain = chain;
		ANN(obj->_Chain);
		return obj;
    };
    static SmartsRoot_sp create(AtomOrBondMatchNode_sp node )
    {
	    _G();
	    BondListMatchNode_sp temp = _Nil<BondListMatchNode_O>();
		return create(node,temp);
    };

public:



//        bool	recognizesRingId(int id) { return this->_RingLookup.count(id); };
//        void	defineRingAtom( chem::Atom_sp a, int id ) { this->_RingLookup[id] = a; };
//        chem::Atom_sp	getRingAtomWithId(int id)	{ return this->_RingLookup[id]; };
//        void	forgetRingAtomWithId(int id)	{ this->_RingLookup.erase(id); };

virtual	ChemInfoType	type() { return smartsRoot; };
virtual	bool		matches( Root_sp root, chem::Atom_sp atom );
virtual	bool		matches( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );


DEFAULT_CTOR_DTOR(SmartsRoot_O);
};








SMART(AntechamberRoot);
class AntechamberRoot_O : public Root_O
{
    LISP_BASE1(Root_O);
    LISP_CLASS(chem,ChemPkg,AntechamberRoot_O,"AntechamberRoot");

public:
	void initialize();
public:
	void	archiveBase(core::ArchiveP node);
protected:
        core::Symbol_sp          _AssignType;
        RootMatchNode_sp         _AfterMatchTests;
        WildElementDict_sp       _WildElementDictionary;
public:
	static AntechamberRoot_sp create( core::Symbol_sp assignType,
                                AtomOrBondMatchNode_sp node,
                                BondListMatchNode_sp bl,
                                RootMatchNode_sp amt) 
	{_G();
	    GC_ALLOCATE(AntechamberRoot_O, obj ); // RP_Create<AntechamberRoot_O>(lisp);
	            obj->_FirstTest = node;
	            obj->_Chain = bl;
		    ANN(obj->_Chain);
		    obj->_AssignType = assignType;
		    obj->_AfterMatchTests = amt;
		    obj->_WildElementDictionary = _Nil<WildElementDict_O>();
		    return obj;
		};
	static AntechamberRoot_sp create(core::Symbol_sp assignType,
                                AtomOrBondMatchNode_sp node,
                                BondListMatchNode_sp bl )
	{_G();
	    RootMatchNode_sp amt = _Nil<RootMatchNode_O>();
		return create( assignType, node, bl, amt );
	};
	static AntechamberRoot_sp create(core::Symbol_sp assignType,
                                AtomOrBondMatchNode_sp node )
	{
	    _G();
	    BondListMatchNode_sp bl = _Nil<BondListMatchNode_O>();
	    RootMatchNode_sp amt = _Nil<RootMatchNode_O>();
		return create( assignType, node, bl, amt );
	};
	core::Symbol_sp getAssignType() { return this->_AssignType;};
public:


        void    setElementWildCardDictionary(WildElementDict_sp dict) {
                                this->_WildElementDictionary = dict; };
        WildElementDict_sp getElementWildCardDictionary() {
                                return this->_WildElementDictionary; };

//virtual	addon::QDomNode_sp	asXml(string name=XmlTag_AntechamberRoot());
//virtual	void		parseFromXml(addon::QDomNode_sp node);
virtual	ChemInfoType	type() { return antechamberRoot; };
virtual	bool		matches( Root_sp root, chem::Atom_sp atom );
virtual	bool		matches( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );


DEFAULT_CTOR_DTOR(AntechamberRoot_O);
};





SMART(ChemInfo);
class ChemInfo_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ChemInfo_O,"ChemInfo");
#if INIT_TO_FACTORIES
 public:
    static ChemInfo_sp make(core::Cons_sp tests, const string& smarts);
#else
    DECLARE_INIT();
#endif
public:
	void	initialize();
public:
	void	archiveBase(core::ArchiveP node);

private:
	string		_Code;
	Root_sp		_Root;
	string		_CompilerMessage;
public:

public: // Factory creation
//	static core::T_sp prim_ChemInfo(core::Function_sp e, core::Cons_sp args, core::Environment_sp environ, core: );

public:
	bool	compileSmarts(const string& code );
	bool	compileAntechamber(const string& code, WildElementDict_sp dict );
	string	compilerMessage() { return this->_CompilerMessage;};
	void defineTests(core::Cons_sp tests);

	bool	matches(chem::Atom_sp atom);
	ChemInfoMatch_sp getMatch();

		/*! Return the depth of the longest pattern
		 */
	uint	depth() const;


	/*! Describe local environment around an atom - for debugging pattern matching
	 Use this and an XML description of the pattern and check the matching by hand */
	string localAtomEnvironment(int maxDepth) const;
	


        bool    compileSucceeded() {_OF(); ASSERTNOTNULL(this->_Root);return this->_Root.notnilp(); };

	string asSmarts() const;

        string          getCode() { return this->_Code; };

	core::Symbol_sp getAssignType() { return this->_Root->getAssignType(); };
//	ChemInfo_O( const ChemInfo_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(ChemInfo_O);
};


#if 0
class	ChemInfoError;
extern	void	ChemInfoErrorTranslator( const ChemInfoError& e);
class	ChemInfoError: public core::CandoException_O
{
public:
    ChemInfo_sp	_code;
public:
    string message() 
    {_
        ASSERTNOTNULL(this->_code);
	return (this->_Message+this->_code->compilerMessage());
    };

    ChemInfoError(string mess, ChemInfo_sp cd)
    { this->_Message=mess;this->_code = cd; };
    inline virtual ~ChemInfoError() throw() {};
};
#endif








};
TRANSLATE(chem::ChemInfoMatch_O);
TRANSLATE(chem::WildElementDict_O);
TRANSLATE(chem::ChemInfoNode_O);
TRANSLATE(chem::ResidueList_O);
TRANSLATE(chem::RootMatchNode_O);
TRANSLATE(chem::BondMatchNode_O);
TRANSLATE(chem::AtomOrBondMatchNode_O);
TRANSLATE(chem::BondListMatchNode_O);
TRANSLATE(chem::Logical_O);
TRANSLATE(chem::TagSet_O);
TRANSLATE(chem::RingTest_O);
TRANSLATE(chem::ResidueTest_O);
TRANSLATE(chem::BondTest_O);
TRANSLATE(chem::AntechamberBondTest_O);
TRANSLATE(chem::AtomTest_O);
TRANSLATE(chem::AntechamberFocusAtomMatch_O);
TRANSLATE(chem::Chain_O);
TRANSLATE(chem::Branch_O);
TRANSLATE(chem::AfterMatchBondTest_O);
TRANSLATE(chem::Root_O);
TRANSLATE(chem::SmartsRoot_O);
TRANSLATE(chem::AntechamberRoot_O);
TRANSLATE(chem::ChemInfo_O);
#endif
