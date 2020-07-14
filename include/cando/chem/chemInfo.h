/*
    File: chemInfo.h
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
//	chemInfo.h
//
//
#ifndef	CHEM_INFO_H
#define	CHEM_INFO_H

#include <iostream>
#include <string>

#include <clasp/core/common.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/chem/atom.h>

#include <cando/chem/chemPackage.h>
#include <cando/chem/chemInfo.h>


namespace chem {

  SMART(BoundFrame);

  SMART(Root);
  SMART(AntechamberRoot);
  class	ChemInfoNode_O;
  typedef	gctools::smart_ptr<ChemInfoNode_O>	ChemInfoNode_sp;


void walk_nodes(ChemInfoNode_sp node, std::function<void(ChemInfoNode_sp)> const &fn);
size_t calculate_maxtag(ChemInfoNode_sp node);


  SMART(ChemInfoMatch);
  class ChemInfoMatch_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,ChemInfoMatch_O,"ChemInfoMatch",core::CxxObject_O);

  public:
    void initialize();
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    static ChemInfoMatch_sp make(Root_sp root,size_t maxtag, core::HashTableEql_sp ring);
  public:
    bool		_Matches;
    size_t              _MaxTagPlus1;
    Root_sp               _Root;
    core::HashTableEql_sp _RingLookup;
    core::SimpleVector_sp _TagLookup; // core::StringMap<Atom_O>	_TagLookup;
    core::List_sp         _TagHistory;
  public:
    string __repr__() const;
    bool matches();
    void setMatches(bool b) { this->_Matches = b;};
    Root_sp matchRoot() const { return this->_Root; };
    void saveTagLookup();
    void clearAtomTags();
    bool recognizesAtomTag(core::T_sp tag);
    void defineAtomTag(Atom_sp a, core::T_sp tag );
    bool hasAtomWithTag(core::T_sp tag );
    chem::Atom_sp getAtomWithTag(core::T_sp tag);
    core::T_sp getAtomWithTagOrNil(core::T_sp tag);
CL_LISPIFY_NAME("tag");
CL_DEFMETHOD     chem::Atom_sp tag(core::T_sp tag) { return this->getAtomWithTag(tag);};
    CL_DEFMETHOD core::List_sp tag_history() const { return this->_TagHistory; };
    core::HashTableEql_sp tags_as_hashtable() const;
    core::Vector_sp tags_as_vector() const;
    void forgetAtomTag(core::T_sp tag);

    void setRingTag(Atom_sp a, core::T_sp tag);
    bool matchesRingTag(Atom_sp a, core::T_sp tag);
    
	/*! Make sure all the matching atoms are unique - otherwise throw an exception */
//    void throwIfInvalid();
//    void describeClosestMatch();

	/*! Return a BoundFrame using atoms with tag "1","2","3" as the origin, x axis and xy plane
	 */
    BoundFrame_sp boundFrame();

    ChemInfoMatch_O(Root_sp root, size_t maxtagPlus1) : _Root(root), _Matches(false), _MaxTagPlus1(maxtagPlus1), _TagHistory(_Nil<core::T_O>()) {};
  };


  SMART(WildElementDict);
  class WildElementDict_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,WildElementDict_O,"WildElementDict",core::CxxObject_O);

  public:
    void initialize();
  public:
    bool fieldsp() const { return true; };
    void	fields(core::Record_sp node);
  private:
    core::HashTableEqual_sp _AtomWildCards; // core::StringMap<adapt::StringSet_O>	_AtomWildCards;
  public:
    void    addWildName(core::Symbol_sp  wildName);
    void    addWildNameMap(core::Symbol_sp wildName, core::Symbol_sp elementName );
    bool    lexWildNameMatches2Char(char wildNameChar1, char wildNameChar2);
    bool    lexWildNameMatches1Char(char wildName);
    bool    recognizesWildName(core::Symbol_sp name);
    bool    recognizesWildNameElement(core::Symbol_sp wildName, core::Symbol_sp element);
//        adapt::QDomNode_sp      asXml(string name=XmlTag_WildElementDict());
//	void	parseFromXml(adapt::QDomNode_sp node);

    DEFAULT_CTOR_DTOR(WildElementDict_O);
  };


  typedef	enum { noType, recursiveRoot, root, smartsRoot, antechamberRoot, chain, branch,
		       logical, tagSet, ringTest,
		       atomTest, antechamberBondToAtomTest, bondTest, afterMatchBondToAtomTest,
		       antechamberFocusAtomMatch } ChemInfoType;


  string  chemInfoTypeString(ChemInfoType type);
  ChemInfoType chemInfoTypeFromString(const string& name);


  typedef enum	{
      SABUseBondMatcher,
      SABNoBond,
      SABSingleBond,
      SABSingleOrAromaticBond,
      SABDoubleOrAromaticBond,
      SABTripleOrAromaticBond,
      SABDelocalizedBond,
      SABDoubleBond,
      SABTripleBond,
      SABAromaticBond,
      SABAnyBond,
      SABDirectionalSingleUpOrUnspecified,
      SABDirectionalSingleDownOrUnspecified,
      SABDirectionalSingleUp,
      SABDirectionalSingleDown,
      SABSameRingBond,
  } BondEnum;

  extern core::Symbol_sp& _sym_STARSabBondEnumConverterSTAR;

};

ENUM_TRANSLATOR(chem::BondEnum,chem::_sym_STARSabBondEnumConverterSTAR);



namespace chem {

  size_t next_ChemInfoNodeId();
  
/*! Base class for all fragment matching
 */
  SMART(ChemInfoNode);
  class ChemInfoNode_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,ChemInfoNode_O,"ChemInfoNode",core::CxxObject_O);
  public:
    size_t  _Id;
    core::T_sp _Bounds;
  public:
    bool fieldsp() const { return true; };
    void	fields(core::Record_sp node);
  public:
    CL_DEFMETHOD size_t getId() const { return this->_Id; };
    CL_DEFMETHOD core::T_sp bounds() { return this->_Bounds; };
    CL_DEFMETHOD void setf_bounds(core::T_sp bounds) { this->_Bounds = bounds; }
    string __repr__() const;
    virtual uint depth() const;
    CL_DEFMETHOD virtual string asSmarts() const {_OF();SUBCLASS_MUST_IMPLEMENT();};
    virtual bool    matches_Atom(Root_sp root, chem::Atom_sp atom) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
    virtual bool    matches_Bond(Root_sp root, chem::Atom_sp atom, chem::Bond_sp bond) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
    virtual core::T_sp children() = 0;
    
  ChemInfoNode_O() : _Id(next_ChemInfoNodeId()) {};
  };


  SMART(ResidueList);
  class ResidueList_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,ResidueList_O,"ResidueList",core::CxxObject_O);

  public:
    bool fieldsp() const { return true; };
    void	fields(core::Record_sp node);
  private:
  public:
    DEFAULT_CTOR_DTOR(ResidueList_O);
  };




  SMART(RootMatchNode);
  class RootMatchNode_O : public ChemInfoNode_O
  {
    LISP_CLASS(chem,ChemPkg,RootMatchNode_O,"RootMatchNode",ChemInfoNode_O);

  public:
    bool fieldsp() const { return true; };
    void	fields(core::Record_sp node);
  public:
    virtual	bool	matches(Root_sp root) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
	
    DEFAULT_CTOR_DTOR(RootMatchNode_O);
  };

  SMART(BondMatchNode);
  class BondMatchNode_O : public ChemInfoNode_O
  {
    LISP_CLASS(chem,ChemPkg,BondMatchNode_O,"BondMatchNode",ChemInfoNode_O);

  public:
    bool fieldsp() const { return true; };
    void	fields(core::Record_sp node);
  public:
    virtual	bool	matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond ) {_OF();SUBCLASS_MUST_IMPLEMENT(); };
    virtual core::T_sp children() = 0;
    
    DEFAULT_CTOR_DTOR(BondMatchNode_O);
  };





typedef enum {
    SARNone, SARRingSet, SARRingTest
} RingTestEnum;

};

DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::RingTestEnum,chem::_sym_STARSarRingTestEnumConverterSTAR);

namespace chem {
SMART(AtomOrBondMatchNode);
class AtomOrBondMatchNode_O : public BondMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,AtomOrBondMatchNode_O,"AtomOrBondMatchNode",BondMatchNode_O);
public:
  RingTestEnum   _RingTest;
  int            _RingId;
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
  virtual core::T_sp children() = 0;
  virtual bool    matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
  virtual	bool	matches_Atom(Root_sp root, chem::Atom_sp atom);

  void setRingTest(RingTestEnum test);
  RingTestEnum getRingTest() const;
  void setRingId(int id);
  int getRingId() const;
  string __repr__() const;


  AtomOrBondMatchNode_O() : _RingTest(SARNone), _RingId(-1) {};
};



SMART(BondListMatchNode);
class BondListMatchNode_O : public ChemInfoNode_O
{
  LISP_CLASS(chem,ChemPkg,BondListMatchNode_O,"BondListMatchNode",ChemInfoNode_O);

public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
public:
  virtual core::T_sp children() = 0;
  virtual	bool	matches_BondList(Root_sp root, chem::Atom_sp from, chem::BondList_sp bondList ) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
  virtual string asSmarts() const;
  DEFAULT_CTOR_DTOR(BondListMatchNode_O);
};



  typedef	enum { logAlwaysTrue, logIdentity, logNot, logHighPrecedenceAnd, logOr, logLowPrecedenceAnd } LogicalOperatorType;

  extern core::Symbol_sp& _sym_STARSLogicalOperatorTypeConverterSTAR;
};

DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::LogicalOperatorType,chem::_sym_STARLogicalOperatorTypeConverterSTAR);

namespace chem {

class	AtomTest;

SMART(Logical);
class Logical_O : public AtomOrBondMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,Logical_O,"Logical",AtomOrBondMatchNode_O);
public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
  void	initialize();
public:
  LogicalOperatorType	_Operator;
  gc::Nilable<ChemInfoNode_sp>	_Left;
  gc::Nilable<ChemInfoNode_sp>	_Right;

public:
		// Second argument can be NULL
  static Logical_sp create(LogicalOperatorType op,
//                             AtomOrBondMatchNode_sp a1,
//                             gc::Nilable<AtomOrBondMatchNode_sp> a2 )
                           core::T_sp a1,
                           core::T_sp a2)
  { _G();
    GC_ALLOCATE(Logical_O, obj ); // RP_Create<Logical_O>(lisp);
    obj->_Operator = op;
    if (a1.nilp()) {
      obj->_Left = _Nil<core::T_O>();
    } else {
      obj->_Left = gc::As<ChemInfoNode_sp>(a1);
    }
    if (a2.nilp()) {
      obj->_Right = _Nil<core::T_O>();
    } else {
      obj->_Right = gc::As<ChemInfoNode_sp>(a2);
    }
    return obj;
  };
  static Logical_sp create( LogicalOperatorType op, core::T_sp a1)
  {_G();
    if (a1.nilp()) {
      return create(op ,_Nil<core::T_O>(), _Nil<core::T_O>());
    } else {
      return create(op, a1, _Nil<core::T_O>());
    }

  };
public:
  virtual	ChemInfoType	type() { return logical; };
  virtual	bool		matches_Atom( Root_sp root, chem::Atom_sp atom );
  virtual	bool		matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
  CL_DEFMETHOD core::T_sp getLeft() const { return this->_Left; };
  CL_DEFMETHOD core::T_sp getRight() const { return this->_Right; };
  CL_DEFMETHOD void setLeft(ChemInfoNode_sp val) { this->_Left = val; };
  CL_DEFMETHOD void setRight(ChemInfoNode_sp val) { this->_Right = val; };
  virtual uint depth() const;
  virtual string asSmarts() const;
  virtual string __repr__() const;

  static Logical_sp create_logIdentity(core::T_sp nilOrOp);
  static Logical_sp create_logOr(core::T_sp nilOrOp1, core::T_sp nilOrOp2);
  static Logical_sp create_logNot(core::T_sp nilOrOp);
  static Logical_sp create_logLowPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2);
  static Logical_sp create_logHighPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2);
  virtual core::T_sp children();
  LogicalOperatorType logical_operator() const;
  Logical_O() : _Operator(logAlwaysTrue), _Left(_Nil<ChemInfoNode_O>()), _Right(_Nil<ChemInfoNode_O>()) {};
  virtual ~Logical_O() {};
};

SMART(ResidueTest);
class ResidueTest_O : public AtomOrBondMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,ResidueTest_O,"ResidueTest",AtomOrBondMatchNode_O);

public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
  void	initialize();
private:
  BondEnum        	_Bond;
  core::Symbol_sp			_RingTag;
  gc::Nilable<AtomOrBondMatchNode_sp>	_AtomTest;
public:
  static ResidueTest_sp create( BondEnum b, gc::Nilable<AtomOrBondMatchNode_sp> at, core::Symbol_sp ri)
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
  virtual	bool		matches_Atom(Root_sp root, chem::Atom_sp atom );
  virtual string asSmarts() const;
  virtual core::T_sp children();
  DEFAULT_CTOR_DTOR(ResidueTest_O);
};



bool    _matchInAromaticBond(Atom_sp a1, Atom_sp a2);
bool    _matchBondTypes(BondEnum be, chem::BondOrder bo, Atom_sp a1, Atom_sp a2);

FORWARD(BondMatcher);

class BondMatcher_O : public ChemInfoNode_O {
  LISP_CLASS(chem,ChemPkg,BondMatcher_O, "BondMatcher", ChemInfoNode_O );
public:
};

FORWARD(BondLogical);
class BondLogical_O : public BondMatcher_O {
  LISP_CLASS(chem,ChemPkg,BondLogical_O, "BondLogical", BondMatcher_O );
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
  static BondLogical_sp create(LogicalOperatorType op,
//                             AtomOrBondMatchNode_sp a1,
//                             gc::Nilable<AtomOrBondMatchNode_sp> a2 )
                               BondMatcher_sp a1,
                               BondMatcher_sp a2)
  { _G();
    GC_ALLOCATE(BondLogical_O, obj ); // RP_Create<Logical_O>(lisp);
    obj->_Operator = op;
    obj->_Left = a1;
    obj->_Right = a2;
    return obj;
  };
  static BondLogical_sp create(LogicalOperatorType op,
//                             AtomOrBondMatchNode_sp a1,
//                             gc::Nilable<AtomOrBondMatchNode_sp> a2 )
                               BondMatcher_sp a)
  { _G();
    GC_ALLOCATE(BondLogical_O, obj ); // RP_Create<Logical_O>(lisp);
    obj->_Operator = op;
    obj->_Left = a;
    obj->_Right = _Unbound<BondMatcher_O>();
    return obj;
  };
public:
  LogicalOperatorType _Operator;
  BondMatcher_sp  _Left;
  BondMatcher_sp  _Right;
public:
  static BondLogical_sp create_bondLogIdentity(core::T_sp nilOrOp);
  static BondLogical_sp create_bondLogOr(core::T_sp nilOrOp1, core::T_sp nilOrOp2);
  static BondLogical_sp create_bondLogNot(core::T_sp nilOrOp);
  static BondLogical_sp create_bondLogLowPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2);
  static BondLogical_sp create_bondLogHighPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2);
  virtual core::T_sp children();
public:
  LogicalOperatorType bondLogicalOperator() const;
  core::T_sp getLeft() const;
  core::T_sp getRight() const;
  void setLeft(core::T_sp val);
  void setRight(core::T_sp val);
  virtual	bool	matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond);
};

FORWARD(BondTest);
class BondTest_O : public BondMatcher_O {
  LISP_CLASS(chem,ChemPkg,BondTest_O, "BondTest", BondMatcher_O );
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
  static BondTest_sp make( BondEnum be );
public:
  BondEnum        _Bond;
public:
    virtual	bool	matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond);
  virtual core::T_sp children();
  BondEnum bondTestGetBond() const;
  BondTest_O(BondEnum be) : _Bond(be) {};
};
  

FORWARD(AtomTest);
SMART(BondToAtomTest);
class BondToAtomTest_O : public BondMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,BondToAtomTest_O,"BondToAtomTest",BondMatchNode_O);

public:
  void	initialize();
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
public:
  BondEnum		_Bond;
  BondMatcher_sp        _BondMatcher;
  gc::Nilable<AtomOrBondMatchNode_sp>	_AtomTest;

public:
#if 0
  static BondToAtomTest_sp create(BondEnum b, core::T_sp nilOrNode )
  {_G();
    GC_ALLOCATE(BondToAtomTest_O, obj ); // RP_Create<BondToAtomTest_O>(lisp);
    obj->_Bond = b;
    if (nilOrNode.nilp()) {
      obj->_AtomTest = _Nil<core::T_O>();
    } else {
      obj->_AtomTest = gc::As<AtomOrBondMatchNode_sp>(nilOrNode);
    }
    return obj;
  }
#endif
public:
  BondEnum	bondType();
  bool setfBondTypeIfOptimizable(BondEnum be);
  bool  bondMatcherBoundP() const;
  BondMatcher_sp bondMatcher() const;
  virtual	ChemInfoType	type() { return bondTest;};
  virtual	bool	matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
  virtual string asSmarts() const;

  static BondToAtomTest_sp makeBondToAtomTest(BondEnum be, core::T_sp nilOrNode);
  static BondToAtomTest_sp makeBondMatcherToAtomTest(BondMatcher_sp bm, core::T_sp nilOrNode);

  void setAtomTest(core::T_sp atomTest);
  virtual core::T_sp children();
  BondToAtomTest_O() : _Bond(SABUseBondMatcher), _BondMatcher(_Unbound<BondMatcher_O>()) {};
  BondToAtomTest_O(BondEnum be) : _Bond(be), _BondMatcher(_Unbound<BondMatcher_O>()) {};
};



SMART(AntechamberBondToAtomTest);
class AntechamberBondToAtomTest_O : public BondMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,AntechamberBondToAtomTest_O,"AntechamberBondToAtomTest",BondMatchNode_O);

public:
  void initialize();
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
private:
  core::Symbol_sp                   _Element;
  int                     _Neighbors;
  gc::Nilable<AtomOrBondMatchNode_sp>	_AtomProperties;
  core::Symbol_sp                   _Tag;
public:
  static AntechamberBondToAtomTest_sp create_args( core::Symbol_sp element, int neighbors,
                                                   AtomOrBondMatchNode_sp props, core::Symbol_sp tag );
public:

//virtual	adapt::QDomNode_sp	asXml(string name=XmlTag_AntechamberBondToAtomTest());
//virtual	void	parseFromXml(adapt::QDomNode_sp node );
  virtual	ChemInfoType	type() { return antechamberBondToAtomTest;};
  bool	matchBasic( AntechamberRoot_sp root, chem::Atom_sp atom );
  virtual	bool	matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
  virtual bool    matches_Atom( Root_sp root, chem::Atom_sp atom );
  virtual string asSmarts() const;
  virtual core::T_sp children();
  DEFAULT_CTOR_DTOR(AntechamberBondToAtomTest_O);
};





typedef	enum	{
    SAPNone		=0,
    SAPWildCard		,
    SAPDegree		,
    SAPElement		,
    SAPTotalHCount		,
    SAPImplicitHCount	,
    SAPRingConnectivity,
    SAPRingMembershipCount	,
    SAPRingTest,
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
    SAPPredicateName		,
    SAPAromaticElement	,
    SAPAliphatic		,
    SAPAromatic,
    SAPAtomMap
} AtomTestEnum;

extern core::Symbol_sp& _sym_STARAtomTestEnumConverterSTAR;
};

DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::AtomTestEnum,chem::_sym_STARAtomTestEnumConverterSTAR);



namespace chem {
SMART(AtomTest);
class AtomTest_O : public AtomOrBondMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,AtomTest_O,"AtomTest",AtomOrBondMatchNode_O);
  
public:
  void initialize();
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
private:
  AtomTestEnum	_Test;
  int		        _IntArg;
  int		        _NumArg;
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
  static AtomTest_sp create( AtomTestEnum t) { return create( t, 0, 0, NULL, _Nil<core::Symbol_O>());};
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPElectronegativeElement() { return create(SAPElectronegativeElement); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPAliphatic() { return create(SAPAliphatic); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPAromatic() { return create(SAPAromatic); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPChiralityAntiClockwise() { return create(SAPChiralityAntiClockwise); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPChiralityClockwise() { return create(SAPChiralityClockwise); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPPiBondOrbital() { return create(SAPPiBondOrbital); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPWildCard() { return create(SAPWildCard); };

  static AtomTest_sp create( AtomTestEnum t, int iArg ) { return create(  t, iArg, 0, NULL,_Nil<core::Symbol_O>() );};
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPAromaticPiElectron(int intVal) { return create(SAPAromaticPiElectron,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPAtomicMass(int intVal) { return create(SAPAtomicMass,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPAtomicNumber(int intVal) { return create(SAPAtomicNumber,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPConnectivity(int intVal) { return create(SAPConnectivity,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPDegree(int intVal) { return create(SAPDegree,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPGroupNumber(int intVal) { return create(SAPGroupNumber,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPHeavyAtomTotalBond(int intVal) { return create(SAPHeavyAtomTotalBond,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPImplicitHCount(int intVal) { return create(SAPImplicitHCount,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPLonePair(int intVal) { return create(SAPLonePair,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPNegativeCharge(int intVal) { return create(SAPNegativeCharge,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPNegativeFormalCharge(int intVal) { return create(SAPNegativeFormalCharge,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPPositiveCharge(int intVal) { return create(SAPPositiveCharge,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPPositiveFormalCharge(int intVal) { return create(SAPPositiveFormalCharge,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPResidueTest(int intVal) { return create(SAPResidueTest,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPRingMembershipCount(int intVal) { return create(SAPRingMembershipCount,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPRingSize(int intVal) { return create(SAPRingSize,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPRingTest(int intVal) { return create(SAPRingTest,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPRingConnectivity(int intVal) { return create(SAPRingConnectivity,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPTotalHCount(int intVal) { return create(SAPTotalHCount,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPTransitionMetal(int intVal) { return create(SAPTransitionMetal,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPValence(int intVal) { return create(SAPValence,intVal); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPAtomMap(int intVal) { return create(SAPAtomMap,intVal); };

  static AtomTest_sp create( AtomTestEnum t, int iArg, int num ) 
  {_G();
    return create( t, iArg, num, NULL,_Nil<core::Symbol_O>() );
  };
  static AtomTest_sp create( AtomTestEnum t, const char* el ) 
  {_G();
    LOG(BF("Create an AtomTest with string argument: %s")%el);
    return create( t, 0, 0, el,_Nil<core::Symbol_O>());
  };
  static AtomTest_sp create( AtomTestEnum t, core::Symbol_sp sym) {return create( t, 0, 0, NULL, sym );}

  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPAromaticElement(const std::string& symbolName) { return create(SAPAromaticElement,chem::chemkw_intern(symbolName)); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPElement(const std::string& symbolName) { return create(SAPElement,chem::chemkw_intern(symbolName)); };
  CL_DEF_CLASS_METHOD static AtomTest_sp create_SAPPredicateName(core::Symbol_sp name) { return create(SAPPredicateName,name); };
    
    // Special makers
  CL_LISPIFY_NAME(make-atom-test-in-bond1);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_in_bond1(BondEnum bond)
  {
    return create( SAPInBond, /*(int)*/bond, 1);
  }
  CL_LISPIFY_NAME(make-atom-test-in-bond);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_in_bond(BondEnum bond, int num)
  {
    return create( SAPInBond, /*(int)*/bond, num);
  }
  CL_LISPIFY_NAME(make-atom-test-not-bonded-to-previous);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_not_bonded_to_previous(BondEnum bond)
  {
    return create( SAPNotBondedToPrevious, bond );
  }
  CL_LISPIFY_NAME(make-atom-test-bonded-to-previous);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_bonded_to_previous(BondEnum bond)
  {
    return create( SAPBondedToPrevious, bond );
  }
  CL_LISPIFY_NAME(make-atom-test-ar-level);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_ar_level(int arLevel)
  {
    return create( SAPArLevel, arLevel );
  }
  CL_LISPIFY_NAME(make-atom-test-ring-size);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_ring_size(int size)
  {
    return create( SAPRingSize, size );
  }
  CL_LISPIFY_NAME(make-atom-test-ring-connectivity);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_ring_connectivity(int num)
  {
    return create( SAPRingConnectivity, num);
  }
  CL_LISPIFY_NAME(make-atom-test-ring-size2);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_ring_size2(int size, int num)
  {
    return create( SAPRingSize, num, size );
  }
  CL_LISPIFY_NAME(make-atom-test-no-ring-membership);
  CL_DEF_CLASS_METHOD static AtomTest_sp create_no_ring_membership()
  {
    return create( SAPNoRing);
  }
public:
  CL_DEFMETHOD AtomTestEnum	atomTestType() { return this->_Test; };
  int getIntArg();
  core::Symbol_sp getSymbolArg() const;
  
  string		testName(AtomTestEnum t) const;
  bool matchesAm1BccX(chem::Atom_sp atom) const;
  bool matchesAm1BccY(chem::Atom_sp atom) const;

  CL_DEFMETHOD void set_test(AtomTestEnum test) { this->_Test = test; };
  virtual ChemInfoType    type() { return atomTest;};
  virtual	bool	matches_Atom( Root_sp root, chem::Atom_sp atom );
  virtual	bool	matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond);
  virtual string asSmarts() const;

  virtual core::T_sp children() {
    return _Nil<core::T_O>();
  };

  AtomTest_O() : _Test(SAPNone), _IntArg(0), _NumArg(0), _StringArg(""), _SymbolArg(_Nil<core::Symbol_O>()) {};
  virtual ~AtomTest_O() {};
};





/*!     This is for matching the focus atom for GAFF type assignment rules
 */
SMART(AntechamberFocusAtomMatch);
class AntechamberFocusAtomMatch_O : public AtomOrBondMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,AntechamberFocusAtomMatch_O,"AntechamberFocusAtomMatch",AtomOrBondMatchNode_O);

public:
  void initialize();
public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
private:
  gc::Nilable<ResidueList_sp>  _ResidueNames;
  int             _AtomicNumber;
  int             _NumberOfAttachedAtoms;
  int             _NumberOfAttachedHydrogens;
  int             _NumberOfElectronWithdrawingGroups;
  gc::Nilable<AtomOrBondMatchNode_sp>   _AtomicProperty;
public:
  static AntechamberFocusAtomMatch_sp create(gc::Nilable<ResidueList_sp> residueNames,
                                             int     atomicNumber,
                                             int     numberOfAttachedAtoms,
                                             int     numberOfAttachedHydrogens,
                                             int     numberOfElectronWithdrawingGroups,
                                             gc::Nilable<AtomOrBondMatchNode_sp> atomicProperty ) 
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
                                             gc::Nilable<ResidueList_sp> residueNames,
                                             int     atomicNumber,
                                             int     numberOfAttachedAtoms,
                                             int     numberOfAttachedHydrogens,
                                             int     numberOfElectronWithdrawingGroups )
  {_G();
    gc::Nilable<AtomOrBondMatchNode_sp> atomicProperty;
    atomicProperty = _Nil<core::T_O>();
    return create( residueNames, atomicNumber, numberOfAttachedAtoms, numberOfAttachedHydrogens, numberOfElectronWithdrawingGroups, atomicProperty);
  };
public:
  string asSmarts() const;
  virtual core::T_sp children();
  virtual	bool	matches_Atom( Root_sp root, chem::Atom_sp atom );
  virtual ChemInfoType    type() { return antechamberFocusAtomMatch; };
  DEFAULT_CTOR_DTOR(AntechamberFocusAtomMatch_O);
};







SMART(Chain);
class Chain_O : public BondListMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,Chain_O,"Chain",BondListMatchNode_O);

public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
  void	initialize();
public:
  gc::Nilable<BondMatchNode_sp>	_Head;
  gc::Nilable<BondListMatchNode_sp>	_Tail;
public:

  virtual uint depth() const;
  virtual string asSmarts() const;
		//! Second argument can be NULL
  CL_LAMBDA(head tail);
  CL_LISPIFY_NAME("make-chain.head.tail");
//    CL_DEF_CLASS_METHOD static Chain_sp create_head_tail( BondMatchNode_sp head, BondListMatchNode_sp tail ) {
  CL_DEF_CLASS_METHOD static Chain_sp create_head_tail( core::T_sp head, core::T_sp tail ) {
    GC_ALLOCATE(Chain_O, obj ); // RP_Create<Chain_O>(lisp);
//      obj->_Head = head;
//      obj->_Tail = tail;
    if (head.nilp()){
      obj->_Head =  _Nil<core::T_O>();
    } else {
      obj->_Head =  gc::As<BondMatchNode_sp>(head);
    } 
    if (tail.nilp()){
      obj->_Tail =  _Nil<core::T_O>();
    } else {
      obj->_Tail =  gc::As<BondListMatchNode_sp>(tail);
    } 
    return obj;
  }
  CL_LAMBDA(head);
  CL_LISPIFY_NAME("make-chain.head");
//    CL_DEF_CLASS_METHOD static Chain_sp create_head( BondMatchNode_sp head )
  CL_DEF_CLASS_METHOD static Chain_sp create_head( core::T_sp head )
  {
    if (head.nilp()){
      return create_head_tail(_Nil<core::T_O>(),_Nil<BondListMatchNode_O>());
    } else {
      return create_head_tail(head,_Nil<BondListMatchNode_O>());
    }
  }
  CL_DEF_CLASS_METHOD static Chain_sp create_empty_chain() {
    GC_ALLOCATE(Chain_O,obj);
    return obj;
  }
public:
  virtual	ChemInfoType	type() { return chain; };
  virtual	bool		matches_Atom( Root_sp root, chem::Atom_sp from );
  virtual	bool		matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
  virtual	bool		matches_BondList( Root_sp root, chem::Atom_sp from, chem::BondList_sp neighbors );
  CL_DEFMETHOD void chain_set_head(BondMatchNode_sp head) { this->_Head = head; };
//    CL_DEFMETHOD void chain_set_tail(BondListMatchNode_sp tail) { this->_Tail = tail; };
  CL_DEFMETHOD void chain_set_tail(core::T_sp tail) {
    if (tail.nilp())
    {
      this->_Tail =   _Nil<core::T_O>();
    } else {
      this->_Tail = tail;
    }
  };

  BondMatchNode_sp chain_get_head() { return this->_Head; };
  BondListMatchNode_sp chain_get_tail();
  virtual core::T_sp children();
  Chain_O() : _Head(_Nil<BondMatchNode_O>()), _Tail(_Nil<BondListMatchNode_O>()) {};
};


SMART(Branch);
class Branch_O : public BondListMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,Branch_O,"Branch",BondListMatchNode_O);

public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
  void	initialize();
public:
  gc::Nilable<BondListMatchNode_sp>	_Left;
  gc::Nilable<BondListMatchNode_sp>	_Right;
public:
  virtual uint depth() const;
  virtual string asSmarts() const;
		//! Second argument can be NULL
  CL_LAMBDA(left right);
  CL_LISPIFY_NAME("make-branch.left.right");
//    CL_DEF_CLASS_METHOD static Branch_sp create_left_right( BondListMatchNode_sp left, BondListMatchNode_sp right )
//    {
//      _G();
//      GC_ALLOCATE(Branch_O, obj ); // RP_Create<Branch_O>(lisp);
//      obj->_Left = left;
//      obj->_Right = right;
//      return obj;
//    }
  CL_DEF_CLASS_METHOD static Branch_sp create_left_right( core::T_sp left, core::T_sp right )
  {
    _G();
    GC_ALLOCATE(Branch_O, obj ); // RP_Create<Branch_O>(lisp);
    if(left.nilp()){
      obj->_Left = _Nil<core::T_O>();
    } else {       
      obj->_Left = gc::As<BondListMatchNode_sp>(left);
    }
    if(right.nilp()){
      obj->_Right = _Nil<core::T_O>();
    } else {       
      obj->_Right = gc::As<BondListMatchNode_sp>(right);
    }
      
    return obj;
  }
  CL_LAMBDA(left);
  CL_LISPIFY_NAME("make-branch.left");
  CL_DEF_CLASS_METHOD static Branch_sp create_left( BondListMatchNode_sp left ) {
    _G();
    GC_ALLOCATE(Branch_O, obj ); // RP_Create<Branch_O>(lisp);
    obj->_Left = left;
    obj->_Right = _Nil<core::T_O>();
    return obj;
  }
public:



  virtual	ChemInfoType	type() { return branch; };
  virtual	bool		matches_BondList( Root_sp root, chem::Atom_sp from, chem::BondList_sp neighbors );  /* CHECK CODE TO SEE IF IT HANDLES RIGHT=NULL */
//virtual	adapt::QDomNode_sp	asXml(string name=XmlTag_Branch());
//virtual	void	parseFromXml(adapt::QDomNode_sp node);

  void branch_set_left(core::T_sp left);
  void branch_set_right(core::T_sp right);

  virtual core::T_sp children();
  DEFAULT_CTOR_DTOR(Branch_O);
};

SMART(AfterMatchBondToAtomTest);
class AfterMatchBondToAtomTest_O : public RootMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,AfterMatchBondToAtomTest_O,"AfterMatchBondToAtomTest",RootMatchNode_O);

public:
  void initialize();
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
private:
  core::Symbol_sp _AtomTag1;
  core::Symbol_sp _AtomTag2;
  BondEnum	_Bond;
public:
  static AfterMatchBondToAtomTest_sp create( core::Symbol_sp tag1, core::Symbol_sp tag2, BondEnum b)
  {
    _G();
    GC_ALLOCATE(AfterMatchBondToAtomTest_O, obj ); // RP_Create<AfterMatchBondToAtomTest_O>(lisp);
    obj->_AtomTag1 = tag1;
    obj->_AtomTag2 = tag2;
    obj->_Bond = b;
    return obj;
  };
public:
// virtual	adapt::QDomNode_sp	        asXml(string name=XmlTag_AfterMatchBondToAtomTest());
  virtual	ChemInfoType            type() { return afterMatchBondToAtomTest; };
  virtual	bool    		matches( Root_sp root );


  DEFAULT_CTOR_DTOR(AfterMatchBondToAtomTest_O);
};







SMART(Root);
class Root_O : public AtomOrBondMatchNode_O
{
  LISP_CLASS(chem,ChemPkg,Root_O,"Root",AtomOrBondMatchNode_O);

public:
  void initialize();
public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
public:
//  gc::Nilable<AtomOrBondMatchNode_sp>	_FirstTest;
//  gc::Nilable<BondListMatchNode_sp>	_Chain;
  gc::Nilable<ChemInfoNode_sp>          _Node;
	/*! Store tests described as lambdas that take a single atom argument and return a boolean.
	 These can be incorporated into the smarts code as <xxxx> where xxxx is the symbol
	name of the test. */
  core::T_sp		_Tests;
  size_t _MaxTag;
  std::string _Code;
//    core::HashTableEql_sp               _RingTags;
public:
  virtual uint depth() const;
  virtual string originalCode() const { return this->_Code;};
  virtual string asSmarts() const;
  static Root_sp create(const std::string& originalCode, ChemInfoNode_sp node)
  {
    size_t maxtag = calculate_maxtag(node);
    GC_ALLOCATE_VARIADIC(Root_O, obj, originalCode, node, maxtag ); // RP_Create<Root_O>(lisp);
    return obj;
  }
  core::HashTableEq_sp lazyTests();
  void setTests(core::List_sp tests);
  void addTest(core::Symbol_sp testSymbol, core::Function_sp testCode);
  bool evaluateTest(core::Symbol_sp testSym, Atom_sp atom);
  virtual core::Symbol_sp getAssignType() { return _Nil<core::Symbol_O>();};
  virtual core::T_sp children();
public:
  virtual	ChemInfoType	type() { return root; };
  virtual	bool		matches_Atom( Root_sp root, chem::Atom_sp atom );
  virtual	bool		matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );
  Root_O(const std::string& code, ChemInfoNode_sp node, size_t maxtag) : _Code(code), _Node(node), _Tests(_Nil<core::T_O>()), _MaxTag(maxtag) {}
  Root_O(const std::string& code) : _Code(code) {};
  Root_O() : _Code(""), _Node(_Nil<core::T_O>()), _Tests(_Nil<core::T_O>()), _MaxTag(0) {};
};







SMART(SmartsRoot);
class SmartsRoot_O : public Root_O
{
  LISP_CLASS(chem,ChemPkg,SmartsRoot_O,"SmartsRoot",Root_O);

public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
public:
  static SmartsRoot_sp make(const std::string& code, ChemInfoNode_sp cinode,size_t maxtag);


//        bool	recognizesRingId(int id) { return this->_RingLookup.count(id); };
//        void	defineRingAtom( chem::Atom_sp a, int id ) { this->_RingLookup[id] = a; };
//        chem::Atom_sp	getRingAtomWithId(int id)	{ return this->_RingLookup[id]; };
//        void	forgetRingAtomWithId(int id)	{ this->_RingLookup.erase(id); };

  virtual	ChemInfoType	type() { return smartsRoot; };
  virtual	bool		matches_Atom( Root_sp root, chem::Atom_sp atom );
  virtual	bool		matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );

  SmartsRoot_O(const std::string& code, ChemInfoNode_sp node, size_t maxtag) : Root_O(code,node,maxtag) {};
  SmartsRoot_O(const std::string& code) : Root_O(code) {};
  SmartsRoot_O() {};
};

FORWARD(Smirks);
class Smirks_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,Smirks_O,"Smirks",core::CxxObject_O);
public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
public:
  CL_LISPIFY_NAME("make-smirks");
  CL_DEF_CLASS_METHOD static Smirks_sp make(SmartsRoot_sp r, SmartsRoot_sp p) {
    GC_ALLOCATE_VARIADIC(Smirks_O,smirks,r,p);
    return smirks;
  }

  CL_DEFMETHOD SmartsRoot_sp reactant() { return this->_reactant; };
  CL_DEFMETHOD SmartsRoot_sp product() { return this->_product; };

  Smirks_O(SmartsRoot_sp r, SmartsRoot_sp p) : _reactant(r), _product(p) {};

  SmartsRoot_sp _reactant;
  SmartsRoot_sp _product;
  
};

};

namespace chem {
SMART(AntechamberRoot);
class AntechamberRoot_O : public Root_O
{
  LISP_CLASS(chem,ChemPkg,AntechamberRoot_O,"AntechamberRoot",Root_O);

public:
  void initialize();
public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);
public:
   gc::Nilable<RootMatchNode_sp>         _AfterMatchTests;
  gc::Nilable<WildElementDict_sp>       _WildElementDictionary;
public:
  
  void    setElementWildCardDictionary(WildElementDict_sp dict) {
    this->_WildElementDictionary = dict; };
  WildElementDict_sp getElementWildCardDictionary() {
    return this->_WildElementDictionary; };

//virtual	adapt::QDomNode_sp	asXml(string name=XmlTag_AntechamberRoot());
//virtual	void		parseFromXml(adapt::QDomNode_sp node);
  virtual	ChemInfoType	type() { return antechamberRoot; };
  virtual	bool		matches_Atom( Root_sp root, chem::Atom_sp atom );
  virtual	bool		matches_Bond( Root_sp root, chem::Atom_sp from, chem::Bond_sp bond );

  string descriptionOfContents() const;
  virtual core::T_sp children();
  AntechamberRoot_O(const std::string& code) : Root_O(code) {};
  DEFAULT_CTOR_DTOR(AntechamberRoot_O);
};
};

namespace chem {
class MoleculeGraph_O;
};

// MoleculeGraph_O objects cannot move in memory
// we pass pointers to them to boost graph C++ code.
template <>
struct gctools::GCInfo<chem::MoleculeGraph_O> {
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = collectable_immobile;
};


namespace chem {
struct MoleculeVertexData {
  MoleculeGraph_O* _moleculeGraph;
  int _AtomIndex;
  MoleculeVertexData(MoleculeGraph_O* g, int i) : _moleculeGraph(g), _AtomIndex(i) {};
  MoleculeVertexData(): _AtomIndex(-1) {};
};

#if 0
struct MoleculeEdgeData {
  BondOrder _BondOrder;
  MoleculeEdgeData(BondOrder bo) : _BondOrder(bo) {};
  bool operator<(const MoleculeEdgeData& other) const {
    return this->_BondOrder<other._BondOrder;
  }
  bool operator==(const MoleculeEdgeData& other) const {
    return this->_BondOrder==other._BondOrder;
  }
};
#endif

typedef boost::property<boost::edge_weight_t,BondOrder> BondOrderProperty;

    //create an -undirected- graph type, using vectors as the underlying containers
    //and an adjacency_list as the basic representation
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              MoleculeVertexData,
                              BondOrderProperty> MoleculeGraphType;

FORWARD(MoleculeGraph);
class MoleculeGraph_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,MoleculeGraph_O,"MoleculeGraph",core::CxxObject_O);
public:
  void initialize();
public:
  core::HashTableEq_sp     _nodes_to_index;
  core::ComplexVector_T_sp _nodes;
  MoleculeGraphType*       _moleculeGraph;
  size_t                   _num_edges;
public:

  size_t add_vertex(core::T_sp vertex);
  size_t vertex_index(core::T_sp vertex);
  core::T_sp get_vertex(size_t index);
  
  void add_edge(size_t v1, size_t v2, BondOrder bo);

  /*! Callback takes vertex */
  void walk_vertices(core::T_sp callback);
  /*! Callback takes vertex1, vertex2 and bond_order */
  void walk_edges(core::T_sp callback);
  size_t num_vertices();
  size_t num_edges();
  
  MoleculeGraph_O(Molecule_sp);
  MoleculeGraph_O();
  ~MoleculeGraph_O();
};

};

namespace chem {
class ChemInfoGraph_O;
};

template <>
struct gctools::GCInfo<chem::ChemInfoGraph_O> {
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = collectable_immobile;
};


namespace chem {



struct ChemInfoVertexData {
  ChemInfoGraph_O* _chemInfoGraph;
  int _NodeIndex;
  ChemInfoVertexData(ChemInfoGraph_O* g, int i) : _chemInfoGraph(g), _NodeIndex(i) {};
  ChemInfoVertexData() : _NodeIndex(-1) {};
};

#if 0
struct ChemInfoEdgeData {
  int _BondMatcherNodeIndex;
  BondEnum _BondEnum;
  ChemInfoEdgeData(int ni,BondEnum be) : _BondMatcherNodeIndex(ni), _BondEnum(be) {};
  bool operator<(const ChemInfoEdgeData& other) const {
    return this->_BondMatcherNodeIndex < other._BondMatcherNodeIndex;
  }
  bool operator==(const ChemInfoEdgeData& other) const {
    return this->_BondMatcherNodeIndex ==other._BondMatcherNodeIndex;
  }

};
#endif



struct RingBond {
  BondEnum _Bond;
  int _NodeIndex;
  RingBond(BondEnum b, int ni) : _Bond(b), _NodeIndex(ni) {};
};

struct RingClosers {
  bool _Active;
  int _NodeIndex;
  std::vector<RingBond> _Bonds;
  RingClosers() : _Active(false) {};
};

typedef boost::property<boost::edge_index_t,int> EdgeProperty;

    //create an -undirected- graph type, using vectors as the underlying containers
    //and an adjacency_list as the basic representation
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              ChemInfoVertexData,
                              EdgeProperty> ChemInfoGraphType;

FORWARD(ChemInfoGraph);
class ChemInfoGraph_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,ChemInfoGraph_O,"ChemInfoGraph",core::CxxObject_O);
public:
  void initialize();
public:
  Root_sp                _Root;
  core::HashTableEq_sp   _nodes_to_index;
  std::vector<size_t>               _nodeOrder;
  gctools::Vec0<ChemInfoNode_sp> _atomNodes;
  gctools::Vec0<BondToAtomTest_sp> _bondNodes;
  ChemInfoGraphType*     _chemInfoGraph;
public:

  ChemInfoGraph_O(Root_sp);
  ~ChemInfoGraph_O();
};

};









namespace chem {



core::T_mv chem__chem_info_match(Root_sp testRoot, Atom_sp atom);
SmartsRoot_sp chem__compile_smarts(const string& smarts, core::List_sp tests);
AntechamberRoot_mv chem__compile_antechamber(const string& smarts,WildElementDict_sp xpdict);
};

#endif
