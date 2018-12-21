/*
    File: chemInfo.cc
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
#define DEBUG_LEVEL_NONE
//
//
//     	chemInfo.cc
//
//
#include <clasp/core/common.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/chemInfo.h>
#include <clasp/core/hashTableEqual.h>
#include <clasp/core/hashTableEql.h>
//#include "core/archiveNode.h"
#include <clasp/core/evaluator.h>
#include <clasp/core/environment.h>
//#include "boundFrame.h"
#include <cando/adapt/symbolMap.h>
#include <cando/adapt/symbolSet.h>
#include <cando/chem/bond.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>

#if 0
#define CI_LOG(x) printf x;
#else
#define CI_LOG(x)
#endif


//extern chem::SmartsRoot_sp smarts_compile(const string &input, stringstream &errorStream);
extern chem::AntechamberRoot_sp gaff_compile(const string &input, chem::WildElementDict_sp dict, stringstream &errorStream);
extern string antechamberError();

namespace chem {

ChemInfoMatch_sp current_match() {
  ChemInfoMatch_sp match =  gc::As<ChemInfoMatch_sp>(_sym_STARcurrent_matchSTAR->symbolValue());
  return match;
}

std::atomic<size_t> global_ChemInfoNodeId;
size_t next_ChemInfoNodeId()
{
  return ++global_ChemInfoNodeId;
}

class scope_chemInfo {
public:
  scope_chemInfo(){};
  virtual ~scope_chemInfo(){};
};

#define ARGS_af_depth "(arg)"
#define DECL_af_depth ""
#define DOCS_af_depth "depth"
int af_depth(ChemInfoNode_sp n) {
  
  if (n.nilp())
    return 0;
  return n->depth();
};

string sabToString(BondEnum sabType) {
  switch (sabType) {
  case SABNoBond:
      return "";
  case SABSingleBond:
      return "-";
  case SABSingleOrAromaticBond:
      return "";
  case SABDoubleOrAromaticBond:
      return "=double-or-aromatic";
  case SABTripleOrAromaticBond:
      return "=triple-or-aromatic";
  case SABDelocalizedBond:
      return "deloc";
  case SABDoubleBond:
      return "=";
  case SABTripleBond:
      return "#";
  case SABAromaticBond:
      return ":";
  case SABAnyBond:
      return "~";
  case SABDirectionalSingleUpOrUnspecified:
      return "/?";
  case SABDirectionalSingleDownOrUnspecified:
      return "\\?";
  case SABDirectionalSingleUp:
      return "/";
  case SABDirectionalSingleDown:
      return "\\";
  }
  return "{{unknownSabType}}";
}

ChemInfoMatch_sp ChemInfoMatch_O::make(bool matches, core::HashTableEql_sp tags, core::HashTableEql_sp ringLookup) {
  GC_ALLOCATE(ChemInfoMatch_O,match);
  match->_Matches = matches;
  match->_TagLookup = tags;
  match->_RingLookup = ringLookup;
  return match;
}

void ChemInfoMatch_O::initialize() {
  this->Base::initialize();
  this->_Matches = false;
  this->_TagLookup = core::HashTableEql_O::create_default();
  this->_ClosestMatch = core::HashTableEqual_O::create_default();
}

void ChemInfoMatch_O::fields(core::Record_sp node) {
  //this->Base::fields(node); // T_O
  node->field(INTERN_(kw, matches), this->_Matches);
  node->field(INTERN_(kw, tags), this->_TagLookup);
  node->field(INTERN_(kw, closestMatch), this->_ClosestMatch);
}

CL_LISPIFY_NAME("ChemInfoMatch-matches");
CL_DEFMETHOD bool ChemInfoMatch_O::matches() {
  
  return this->_Matches;
}

string ChemInfoMatch_O::__repr__() const {
  _OF();
  stringstream ss;
  ss << "( " << this->className();
  ss << " :TagLookup '(";
  this->_TagLookup->maphash([&ss](core::T_sp key, core::T_sp val) {
      ss << (BF(":tag \"%s\" :value %s ) ") % _rep_(key) % _rep_(val) );
    });
  ss << " ))";
  return ss.str();
}

void ChemInfoMatch_O::clearAtomTags() {
  
  this->_TagLookup->clrhash();
}

bool ChemInfoMatch_O::recognizesAtomTag(core::T_sp tag) {
  
  return this->_TagLookup->contains(tag);
}

void ChemInfoMatch_O::defineAtomTag(Atom_sp a, core::T_sp tag) {
  CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
  CI_LOG(("ringtag %s atom %s \n", _rep_(tag).c_str(), a->description().c_str()));
  this->_TagLookup->setf_gethash(tag, a);
  CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
  if (this->_TagLookup->hashTableCount() > this->_ClosestMatch->hashTableCount()) {
    CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
    this->_ClosestMatch = this->_TagLookup;
  }
}

void ChemInfoMatch_O::throwIfInvalid() {
  _OF();
  gctools::SmallOrderedSet<Atom_sp> satoms;
  this->_TagLookup->mapHash([this, &satoms](core::Symbol_sp key, Atom_sp atom) {
      if ( satoms.count(atom) ) {
        SIMPLE_ERROR(BF("The ChemInfoMatch is invalid - the matching algorithm or the SMARTS pattern match gave tags with the same atoms: %s") % this->__repr__() );
      }
      satoms.insert(atom);
    });
}

bool ChemInfoMatch_O::hasAtomWithTag(core::T_sp tag) {
  
  return this->_TagLookup->gethash(tag).notnilp();
}

CL_LISPIFY_NAME("getAtomWithTagOrNil");
CL_DEFMETHOD core::T_sp ChemInfoMatch_O::getAtomWithTagOrNil(core::T_sp tag) {
  if (!this->hasAtomWithTag(tag))
    return _Nil<core::T_O>();
  return this->getAtomWithTag(tag);
}

CL_LISPIFY_NAME("getAtomWithTag");
CL_DEFMETHOD Atom_sp ChemInfoMatch_O::getAtomWithTag(core::T_sp tag) {
  
  core::T_mv tatom = this->_TagLookup->gethash(tag);
  if (tatom.nilp()) {
    SIMPLE_ERROR(BF("The ChemInfoMatch doesn't recognize tag[%s]\n"
                    "Available tags(%s)\n") %
                 _rep_(tag) % this->_TagLookup->keysAsString());
  }
  return tatom.as<Atom_O>();
}

void ChemInfoMatch_O::forgetAtomTag(core::T_sp tag) {
  this->_TagLookup->remhash(tag);
}

void ChemInfoMatch_O::setRingTag(Atom_sp atom, core::T_sp index) {
  this->_RingLookup->setf_gethash(index,atom);
}

bool ChemInfoMatch_O::matchesRingTag(Atom_sp atom, core::T_sp tag) {
  core::T_sp thing = this->_RingLookup->gethash(atom,_Nil<T_O>());
  return thing.notnilp();
}

// ------- WildElementDict_O

void WildElementDict_O::initialize() {
  this->Base::initialize();
  this->_AtomWildCards = core::HashTableEqual_O::create_default();
}

CL_LISPIFY_NAME("addWildName");
CL_DEFMETHOD void WildElementDict_O::addWildName(core::Symbol_sp name) {
  adapt::SymbolSet_sp m = adapt::SymbolSet_O::create();
  this->_AtomWildCards->setf_gethash(name, m);
}

CL_LISPIFY_NAME("addWildNameMap");
CL_DEFMETHOD void WildElementDict_O::addWildNameMap(core::Symbol_sp wildName, core::Symbol_sp elementName) {
  _OF();
  if (!this->_AtomWildCards->contains(wildName)) {
    SIMPLE_ERROR(BF("Could not find wild-card %s") % _rep_(core::Cons_O::createList(wildName)));
  }
  core::T_sp tss = this->_AtomWildCards->gethash(wildName).as<adapt::SymbolSet_O>();
  if (tss.nilp()) {
    SIMPLE_ERROR(BF("Could not find AtomWildCards for %s") % _rep_(wildName));
  }
  adapt::SymbolSet_sp ss = tss.as<adapt::SymbolSet_O>();
  ss->insert(elementName);
  this->_AtomWildCards->gethash(wildName).as<adapt::SymbolSet_O>()->insert(elementName);
}

bool WildElementDict_O::recognizesWildName(core::Symbol_sp name) {
  if (this->_AtomWildCards->gethash(name).nilp())
    return false;
  return true;
}

bool WildElementDict_O::recognizesWildNameElement(core::Symbol_sp name, core::Symbol_sp element) {
  core::T_sp symset = this->_AtomWildCards->gethash(name);
  if (symset.nilp()) {
    SIMPLE_ERROR(BF("Unrecognized wild card name(%s)") % _rep_(name));
  }
  if (symset.as<adapt::SymbolSet_O>()->contains(element) != 0)
    return true;
  return false;
}

bool WildElementDict_O::lexWildNameMatches2Char(char c1, char c2) {
  string name;
  name = c1;
  name += c2;
  core::Symbol_sp sname = chemkw_intern(name);
  return (this->_AtomWildCards->contains(sname));
}

bool WildElementDict_O::lexWildNameMatches1Char(char c1) {
  string name;
  name = c1;
  core::Symbol_sp sname = chemkw_intern(name);
  return (this->_AtomWildCards->contains(sname));
}

void WildElementDict_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, data), this->_AtomWildCards);
}

void ResidueList_O::fields(core::Record_sp node) {
  // do nothing
}

void ChemInfo_O::initialize() {
  this->Base::initialize();
  this->_Code = "";
  this->_Root = _Nil<core::T_O>();
  this->_CompilerMessage = "";
}

uint ChemInfo_O::depth() const {
  _OF();
  uint d = af_depth(this->_Root);
  LOG(BF("Returning depth=%d") % d);
  return d;
}

SYMBOL_EXPORT_SC_(ChemPkg,parse_smarts);

CL_LISPIFY_NAME("compileSmarts");
CL_DEFMETHOD bool ChemInfo_O::compileSmarts(const string &code) {
  core::SimpleBaseString_sp scode = core::SimpleBaseString_O::make(code);
  ChemInfoNode_sp node = gc::As<ChemInfoNode_sp>(core::eval::funcall(_sym_parse_smarts,scode));
  SmartsRoot_sp root = SmartsRoot_O::make(node);
  this->_Root = gc::As<Root_sp>(root);
  this->_Code = code;
#if 0
  ASSERTNOTNULL(this->_Root);
  this->_Root = _Nil<core::T_O>();
  this->_Code = code;
  stringstream errors;
  root = smarts_compile(code,errors);
  this->_Root = gc::As<Root_sp>(root);
  ANN(this->_Root);
  if (root.nilp()) {
    this->_CompilerMessage = errors.str();
    return false;
  }
  this->_CompilerMessage = "Successfully compiled SMARTS";
#endif
  return true;
}


CL_DEFMETHOD void ChemInfo_O::defineTests(core::List_sp tests) {
  ASSERTNOTNULL(this->_Root);
  ASSERTF(this->_Root.notnilp(), BF("Root must be defined to define tests"));
  for (auto cur : tests) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    core::Symbol_sp testSymbol = oCar(pair).as<core::Symbol_O>();
    core::Function_sp testCode = oCdr(pair).as<core::Function_O>();
    ASSERTF(testSymbol.notnilp(), BF("The testSymbol was nil"));
    ASSERTF(testCode.notnilp(), BF("The test code was nil"));
    this->_Root->addTest(testSymbol, testCode);
  }
}

string ChemInfo_O::descriptionOfContents() const {
  if (this->_Root.nilp()) {
    return "NIL";
  }
  stringstream ss;
  ss << this->_Root->descriptionOfContents();
  return ss.str();
}


SYMBOL_EXPORT_SC_(ChemPkg,compile_antechamber_type_rule);
CL_LISPIFY_NAME("compileAntechamber");
CL_DEFMETHOD
/*! This function invokes chem:compile-antechamber-type-rule
  chem:compile-antechamber-type-rule is a generic function with the lambda list.
  (system string &optional wild-element-dict)  */
bool ChemInfo_O::compileAntechamber(const string &code, WildElementDict_sp dict) {
  if (!chem::_sym_compile_antechamber_type_rule->fboundp()) return false;
  AntechamberRoot_sp root;
  LOG(BF("Compiling code: %s") % code.c_str());
  this->_Code = code;
  stringstream serr;
  core::SimpleBaseString_sp scode = core::SimpleBaseString_O::make(code);
  root = core::eval::funcall(chem::_sym_compile_antechamber_type_rule->symbolFunction(),_Nil<core::T_O>(),scode,dict);
  this->_Root = gc::As<Root_sp>(root);
  if (root.nilp()) {
    LOG(BF("root was null"));
    this->_CompilerMessage = serr.str();
    LOG(BF("After unsuccessful compile this->_Root = %s") % this->_Root->description().c_str());
    return false;
  }
  root->setElementWildCardDictionary(dict);
  LOG(BF("antechamber compile was successful"));
  this->_CompilerMessage = "Successfully compiled Antechamber";
  LOG(BF("After compile this->_Root = %s") % this->_Root->description().c_str());
  return true;
}

CL_LISPIFY_NAME("matches");
CL_DEFMETHOD core::T_sp ChemInfo_O::matches(chem::Atom_sp a) {
  if (this->_Root.nilp()) {
    SIMPLE_ERROR(BF("The ChemInfo root is nil!"));
  }
  LOG(BF("Starting ChemInfo::matches process with atom: %s") % a->description());
  LOG(BF("The local atom environment(depth=%d) is:\n%s") % this->depth() % a->localEnvironment(4));
  LOG(BF("The pattern to match is smarts[%s]") % this->_Code);
  core::T_mv matches_mv = chem__chem_info_match(this->_Root,a);
  if (matches_mv.notnilp()) {
    return matches_mv.second();
  }
  return _Nil<core::T_O>();
}

CL_LISPIFY_NAME("matches-atom");
CL_DEFMETHOD core::T_sp ChemInfo_O::matches_atom(chem::Atom_sp a) {
  if (this->_Root.nilp()) {
    SIMPLE_ERROR(BF("The ChemInfo root is nil!"));
  }
  LOG(BF("Starting ChemInfo::matches process with atom: %s") % a->description());
  LOG(BF("The local atom environment(depth=%d) is:\n%s") % this->depth() % a->localEnvironment(4));
  LOG(BF("The pattern to match is smarts[%s]") % this->_Code);
  core::T_mv matches_mv = chem__chem_info_match(this->_Root,a);
  return matches_mv.second();
}


void ChemInfo_O::fields(core::Record_sp node) {
  //  this->Base::fields(node);
  node->field(INTERN_(kw, root), this->_Root);
}

struct ChemInfoTypeToName {
  ChemInfoType type;
  string name;
};

ChemInfoTypeToName chemInfoTypesToName[] = {
    {root, "root"},
    {smartsRoot, "smartsRoot"},
    {antechamberRoot, "antechamberRoot"},
    {chain, "chain"},
    {branch, "branch"},
    {logical, "logical"},
    {ringTest, "ringTest"},
    {atomTest, "atomTest"},
    {antechamberBondTest, "antechamberBondTest"},
    {bondTest, "bondTest"},
    {antechamberFocusAtomMatch, "antechamberFocusAtomMatch"},
    {noType, ""}};

string chemInfoTypeString(ChemInfoType type) {
  ChemInfoTypeToName *p;
  p = chemInfoTypesToName;
  while (p->type != noType) {
    if (p->type == type)
      return p->name;
    p++;
  }
  stringstream ss;
  ss << "Unknown ChemInfoType name(" << type << ")";
  THROW_HARD_ERROR(BF("%s") % ss.str());
}

ChemInfoType chemInfoTypeFromString(const string &name) {
  ChemInfoTypeToName *p;
  p = chemInfoTypesToName;
  while (p->type != noType) {
    if (p->name == name)
      return p->type;
    p++;
  }
  stringstream ss;
  ss << "Unknown ChemInfoType name(" << name << ")";
  THROW_HARD_ERROR(BF("%s") % ss.str());
}

void ChemInfoNode_O::fields(core::Record_sp node) {
  // Nothing to do here
  node->field_if_not_default(INTERN_(kw,id),this->_Id,(size_t)0);
  node->field_if_not_nil(INTERN_(kw,bounds),this->_Bounds);
}

string ChemInfoNode_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->className() << " :id " << this->_Id << ">";
  return ss.str();
}

uint ChemInfoNode_O::depth() const {
  _OF();
  LOG(BF("Node type = %s") % this->className());
  LOG(BF("Returning 1"));
  return 1;
};

void BondListMatchNode_O::fields(core::Record_sp node) {
  this->Base::fields(node);
  // Nothing to do here
}

string BondListMatchNode_O::asSmarts() const {
  _OF();
  SUBCLASS_MUST_IMPLEMENT();
}

void AtomOrBondMatchNode_O::fields(core::Record_sp node) {
  this->Base::fields(node);
  node->field_if_not_default(INTERN_(kw,ring_test),this->_RingTest,SARNone);
  node->field_if_not_default(INTERN_(kw,ring_id),this->_RingId,-1);
}

bool AtomOrBondMatchNode_O::matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  return (this->matches(root,from) &&
          this->matches(root, bond->getOtherAtom(from)));
};

bool AtomOrBondMatchNode_O::matches(Root_sp root, chem::Atom_sp atom) {
  if (this->_RingTest==SARNone) {
    return true;
  } else if (this->_RingTest==SARRingSet) {
    current_match()->setRingTag(atom,core::make_fixnum(this->_RingId));
    return true;
  }
  // It's SARRingTest
  return current_match()->matchesRingTag(atom,core::make_fixnum(this->_RingId));
}

CL_DEFMETHOD void AtomOrBondMatchNode_O::setRingTest(RingTestEnum test) {
  this->_RingTest = test;
}
CL_DEFMETHOD RingTestEnum AtomOrBondMatchNode_O::getRingTest() const {
  return this->_RingTest;
}

CL_DEFMETHOD void AtomOrBondMatchNode_O::setRingId(int id) {
  this->_RingId = id;
}

CL_DEFMETHOD int AtomOrBondMatchNode_O::getRingId() const {
  return this->_RingId;
}


void BondMatchNode_O::fields(core::Record_sp node) {
  this->Base::fields(node);
  // Nothing to do here
}

void Logical_O::initialize() {
  this->Base::initialize();
  this->_Operator = logAlwaysTrue;
  this->_Left = _Nil<core::T_O>();
  this->_Right = _Nil<core::T_O>();
}

core::T_sp Logical_O::children() {
  ql::list result;
  if (this->_Left.notnilp()) result << this->_Left;
  if (this->_Right.notnilp()) result << this->_Right;
  return result.cons();
}

string Logical_O::asSmarts() const {
  _OF();
  stringstream ss;
  switch (this->_Operator) {
  case logAlwaysTrue:
      break;
  case logIdentity:
      ss << this->_Left->asSmarts();
      break;
  case logNot:
      ss << "!" << this->_Left->asSmarts();
      break;
  case logHighPrecedenceAnd:
      ss << this->_Left->asSmarts() << "&" << this->_Right->asSmarts();
      break;
  case logLowPrecedenceAnd:
      ss << this->_Left->asSmarts() << ";" << this->_Right->asSmarts();
      break;
  case logOr:
      ss << this->_Left->asSmarts() << "," << this->_Right->asSmarts();
      break;
  };
  return ss.str();
}

uint Logical_O::depth() const {
  _OF();
  uint d = MAX(this->_Left->depth(), this->_Right->depth());
  LOG(BF("Returning Logical depth=%d") % d);
  return d;
}

bool Logical_O::matches(Root_sp root, chem::Atom_sp atom) {
  LOG(BF("Logical pattern: %s") % this->asSmarts());
  LOG(BF("Logical match for atom: %s") % atom->description());
  switch (this->_Operator) {
  case logAlwaysTrue:
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("Always return true"));
      goto SUCCESS;
      break;
  case logIdentity:
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("Identity no-op test"));
      ASSERT(!this->_Left.nilp());
      if (this->_Left->matches(root, atom))
        goto SUCCESS;
      break;
  case logNot:
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("logNot"));
      ASSERT(!this->_Left.nilp());
      if (!(this->_Left->matches(root, atom)))
        goto SUCCESS;
      break;
  case logHighPrecedenceAnd:
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("logHighPrecedenceAnd"));
      ASSERT(!this->_Left.nilp());
      if (this->_Left->matches(root, atom) && this->_Right->matches(root, atom))
        goto SUCCESS;
      break;
  case logLowPrecedenceAnd:
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("logLowPrecedenceAnd"));
      if (this->_Left->matches(root, atom) && this->_Right->matches(root, atom))
        goto SUCCESS;
      break;
  case logOr: {
    CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
    LOG(BF("logOr"));
    ASSERT(!this->_Left.nilp());
    bool leftMatch = this->_Left->matches(root, atom);
    if (leftMatch)
      goto SUCCESS;
    LOG(BF("Left match failed"));
    ASSERT(!this->_Right.nilp());
    bool rightMatch = this->_Right->matches(root, atom);
    if (rightMatch)
      goto SUCCESS;
    LOG(BF("Right match failed"));
  } break;
  default:
      stringstream err;
      err << "Unknown logical operator(" << this->_Operator << ")";
      SIMPLE_ERROR(BF("%s") % err.str());
  }
  //FAIL:
  LOG(BF("FAIL"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!!!"));
  return this->Base::matches(root,atom);
}

bool Logical_O::matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG(BF("Logical match for bond: %s") % bond->describeOther(from));
  switch (this->_Operator) {
  case logAlwaysTrue:
      goto SUCCESS;
      break;
  case logIdentity:
      ASSERT(this->_Left.notnilp());
      if (this->_Left->matches(root, from, bond))
        goto SUCCESS;
      break;
  case logNot:
      ASSERT(this->_Left.notnilp());
      if (!(this->_Left->matches(root, from, bond)))
        goto SUCCESS;
      break;
  case logHighPrecedenceAnd:
  case logLowPrecedenceAnd:
      ASSERT(this->_Left.notnilp());
      if (this->_Left->matches(root, from, bond) && this->_Right->matches(root, from, bond))
        goto SUCCESS;
      break;
  case logOr:
      ASSERT(this->_Left.notnilp());
      if (this->_Left->matches(root, from, bond) || this->_Right->matches(root, from, bond))
        goto SUCCESS;
      break;
  default:
      stringstream err;
      err << "Unknown logical operator(" << this->_Operator << ")";
      SIMPLE_ERROR(BF("%s") % err.str());
  }
  //FAIL:
  LOG(BF("FAIL"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!!!"));
  return this->Base::matches(root,from,bond);
}

core::NullTerminatedEnumAssociation logicalEnum[] = {
    {"logAlwaysTrue",logAlwaysTrue},
    {"logIdentity", logIdentity},
    {"logNot", logNot},
    {"logHighPrecedenceAnd", logHighPrecedenceAnd},
    {"logOr", logOr},
    {"logLowPrecedenceAnd", logLowPrecedenceAnd},
    {"", -1}};

void Logical_O::fields(core::Record_sp node) {
  node->/*pod_*/field(INTERN_(kw, op), this->_Operator);
  node->field_if_not_nil(INTERN_(kw, left), this->_Left);
  node->field_if_not_nil(INTERN_(kw, right), this->_Right);
  this->Base::fields(node);
}

string Logical_O::__repr__() const {
  stringstream ss;
  ss << "#<LOGICAL ";
  ss << ":id " << this->_Id;
  ss << " :op ";
  for ( int i=0; logicalEnum[i]._Enum != -1; ++i ) {
      if (logicalEnum[i]._Enum == this->_Operator) {
          ss << logicalEnum[i]._Key;
      }
  }
      ss << ">";
            return ss.str();
}

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logIdentity(core::T_sp nilOrOp)
{
  return create(logIdentity, nilOrOp );
};

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logOr(core::T_sp nilOrOp1, core::T_sp nilOrOp2)
{
  return create(logOr, nilOrOp1, nilOrOp2 );
};

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logNot(core::T_sp nilOrOp)
{
  return create(logNot, nilOrOp);
};

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logLowPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2)
{
  return create(logLowPrecedenceAnd, nilOrOp1, nilOrOp2);
};

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logHighPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2)
{
  return create(logHighPrecedenceAnd, nilOrOp1, nilOrOp2);
};


// --- TagSet  set a tag for the atom

core::NullTerminatedEnumAssociation bondEnum[] = {
    {"SABNoBond", SABNoBond},
    {"SABSingleBond", SABSingleBond},
    {"SABSingleOrAromaticBond", SABSingleOrAromaticBond},
    {"SABDoubleOrAromaticBond", SABDoubleOrAromaticBond},
    {"SABTripleOrAromaticBond", SABTripleOrAromaticBond},
    {"SABDelocalizedBond", SABDelocalizedBond},
    {"SABDoubleBond", SABDoubleBond},
    {"SABTripleBond", SABTripleBond},
    {"SABAromaticBond", SABAromaticBond},
    {"SABAnyBond", SABAnyBond},
    {"SABDirectionalSingleUpOrUnspecified", SABDirectionalSingleUpOrUnspecified},
    {"SABDirectionalSingleDownOrUnspecified", SABDirectionalSingleDownOrUnspecified},
    {"SABDirectionalSingleUp", SABDirectionalSingleUp},
    {"SABDirectionalSingleDown", SABDirectionalSingleDown},
    {"", -1}};



// ------ ResidueTest

void ResidueTest_O::initialize() {
  this->Base::initialize();
  this->_Bond = SABAnyBond;
  this->_AtomTest = _Nil<core::T_O>();
  this->_RingTag = _Nil<core::Symbol_O>();
}

core::T_sp ResidueTest_O::children() {
  ql::list result;
  if (this->_AtomTest.notnilp()) result << this->_AtomTest;
  return result.cons();
}

string ResidueTest_O::asSmarts() const {
  _OF();
  stringstream ss;
  ss << sabToString(this->_Bond) << this->_AtomTest->asSmarts() << this->_RingTag;
  return ss.str();
}

bool ResidueTest_O::matches(Root_sp root, chem::Atom_sp atom) {
  
  Atom_sp ringAtom;
  SmartsRoot_sp smartsRoot;
  LOG(BF("ResidueTest match for atom: %s") % atom->description().c_str());
  if (root->type() != chem::smartsRoot) {
    SIMPLE_ERROR(BF("Trying to carry out a ResidueTest with a non Smarts root"));
  }
  SIMPLE_ERROR(BF("Handle atom tests in ResidueTest"));
#if 0
  // First check if the AtomTest matches
  ASSERTNOTNULL(this->_AtomTest);
  if (!this->_AtomTest->matches(root, atom)) {
    LOG(BF("The _AtomTest failed"));
    goto FAIL;
  }
  // Now handle the ring test
  smartsRoot = (root).as<SmartsRoot_O>();
  // It does recognize the ring ID so check if the current atom is bonded to it
  // and then have the root forget the atom
  if (!smartsRoot->getMatch()->recognizesAtomTag(this->_RingTag)) {
    SIMPLE_ERROR(BF("We are trying to test the atomTag (%s) but it doesn't exist!") % _rep_(this->_RingTag));
  }
  LOG(BF("      ResidueTest check close"));
  ringAtom = smartsRoot->getMatch()->getAtomWithTag(this->_RingTag);
  LOG(BF("      checking if %s matches ringAtom: %s") % atom->description().c_str() % ringAtom->description().c_str());
  if (!atom->containedByValid()) {
    SIMPLE_ERROR(BF("%s failed containedByValid()") % atom->description());
  }
  if (!ringAtom->containedByValid()) {
    SIMPLE_ERROR(BF("%s failed containedByValid()") % ringAtom->description());
  }
  LOG(BF("   atom->containedBy()->description() = %s") % atom->containedBy()->description().c_str());
  LOG(BF("   ringAtom->containedBy()->description() = %s") % ringAtom->containedBy()->description().c_str());
  if (atom->containedBy() == ringAtom->containedBy()) {
    LOG(BF("RESIDUE MATCH THEY DO!!!!"));
    goto SUCCESS;
  }
 FAIL:
  LOG(BF("FAIL"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
#endif
}

void ResidueTest_O::fields(core::Record_sp node) {
  node->/*pod_*/field( INTERN_(kw,bond), this->_Bond);
  node->field( INTERN_(kw,ringTag), this->_RingTag);
  node->field( INTERN_(kw,atomTest), this->_AtomTest);
  this->Base::fields(node);
}

// ------ BondTest

bool _matchBondTypes(BondEnum be, chem::BondOrder bo) {
  
  LOG(BF("bondOrder = %s") % bondOrderToString(bo).c_str());
  switch (be) {
  case SABSingleBond:
      LOG(BF("SMARTS BondEnum = SABSingleBond"));
      if (bo != chem::singleBond)
        goto nomatch;
      break;
  case SABSingleOrAromaticBond:
      LOG(BF("SMARTS BondEnum = SABSingleOrAromaticBond"));
      if (!(bo == chem::singleBond || bo == chem::aromaticBond))
        goto nomatch;
      break;
  case SABDoubleOrAromaticBond:
      LOG(BF("SMARTS BondEnum = SABDoubleOrAromaticBond"));
      if (!(bo == chem::doubleBond || bo == chem::aromaticBond))
        goto nomatch;
      break;
  case SABTripleOrAromaticBond:
      LOG(BF("SMARTS BondEnum = SABTripleOrAromaticBond"));
      if (!(bo == chem::tripleBond || bo == chem::aromaticBond))
        goto nomatch;
      break;
  case SABDoubleBond:
      LOG(BF("SMARTS BondEnum = SABDoubleBond"));
      if (bo != chem::doubleBond)
        goto nomatch;
      break;
  case SABTripleBond:
      LOG(BF("SMARTS BondEnum = SABTriple"));
      if (bo != chem::tripleBond)
        goto nomatch;
      break;
  case SABAromaticBond:
      LOG(BF("SMARTS BondEnum = SABAromaticBond"));
      if (bo != chem::aromaticBond)
        goto nomatch;
      break;
  case SABAnyBond:
      break;
  case SABDelocalizedBond:
      goto nomatch;
      break;
  case SABDirectionalSingleDownOrUnspecified:
      SIMPLE_ERROR(BF("Must implement SingleDirectionalDownOrUnspecified"));
  case SABDirectionalSingleUpOrUnspecified:
      SIMPLE_ERROR(BF("Must implement SingleDirectionalUpOrUnspecified"));
  case SABDirectionalSingleUp:
  case SABDirectionalSingleDown:
      SIMPLE_ERROR(BF("Must implement directional bonds"));
  default:
      goto nomatch;
  }
  LOG(BF("THEY MATCH!!"));
  return true;
 nomatch:
  LOG(BF("THEY DONT MATCH"));
  return false;
}

void BondTest_O::initialize() {
  this->Base::initialize();
  this->_Bond = SABNoBond;
  this->_AtomTest = _Nil<core::T_O>();
}

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABNoBond(core::T_sp nilOrNode) {
  return create(SABNoBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABSingleOrAromaticBond(core::T_sp nilOrNode) {
  return create(SABSingleOrAromaticBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABDoubleOrAromaticBond(core::T_sp nilOrNode) {
  return create(SABDoubleOrAromaticBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABTripleOrAromaticBond(core::T_sp nilOrNode) {
  return create(SABTripleOrAromaticBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABSingleBond(core::T_sp nilOrNode) {
  return create(SABSingleBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABDoubleBond(core::T_sp nilOrNode) {
  return create(SABDoubleBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABTripleBond(core::T_sp nilOrNode) {
  return create(SABTripleBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABAromaticBond(core::T_sp nilOrNode) {
  return create(SABAromaticBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABAnyBond(core::T_sp nilOrNode) {
  return create(SABAnyBond, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABDirectionalSingleUpOrUnspecified(core::T_sp nilOrNode) {
  return create(SABDirectionalSingleUpOrUnspecified, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABDirectionalSingleDownOrUnspecified(core::T_sp nilOrNode) {
  return create(SABDirectionalSingleDownOrUnspecified, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABDirectionalSingleUp(core::T_sp nilOrNode) {
  return create(SABDirectionalSingleUp, nilOrNode);
};

CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::create_SABDirectionalSingleDown(core::T_sp nilOrNode) {
  return create(SABDirectionalSingleDown, nilOrNode);
};

core::T_sp BondTest_O::children() {
  ql::list result;
  if (this->_AtomTest.notnilp()) result << this->_AtomTest;
  return result.cons();
}


string BondTest_O::asSmarts() const {
  _OF();
  stringstream ss;
  ss << "(" << sabToString(this->_Bond) << this->_AtomTest->asSmarts() << ")";
  return ss.str();
}

CL_DEFMETHOD void BondTest_O::setAtomTest(core::T_sp atomTest)
{
  this->_AtomTest = atomTest;
}

bool BondTest_O::matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG(BF("Checking bond: %s") % bond->describeOther(from));
  chem::BondOrder bo;
  bo = bond->getOrder();
  if (!chem::_matchBondTypes(this->_Bond, bo))
    goto FAIL;
  if (this->_AtomTest->matches(root, bond->getOtherAtom(from)))
    goto SUCCESS;
 FAIL:
  LOG(BF("FAIL!"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
}

void BondTest_O::fields(core::Record_sp node) {
  node->field( INTERN_(kw,atomTest), this->_AtomTest);
  node->/*pod_*/field( INTERN_(kw,bond), this->_Bond);
  this->Base::fields(node);
}

// ------- AtomTest

void AtomTest_O::initialize() {
  _OF();
  this->_Test = SAPNone;
  this->_IntArg = 0;
  this->_NumArg = 0;
  this->_StringArg = "";
}

CL_DEFMETHOD int AtomTest_O::getIntArg() { return this->_IntArg;};

bool AtomTest_O::matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG(BF("AtomTest_O matching pattern: %s") % this->asSmarts());
  switch (this->_Test) {
  case SAPBondedToPrevious:
      if (chem::_matchBondTypes((chem::BondEnum) this->_IntArg, bond->getOrder()))
        goto SUCCESS;
      goto FAIL;
  case SAPNotBondedToPrevious:
      if (!chem::_matchBondTypes((chem::BondEnum) this->_IntArg, bond->getOrder()))
        goto SUCCESS;
      goto FAIL;
  default:
    // Do nothing, fall through to the atom test below
      break;
  }
  // If its not a bond test then
  // test the To atom of the bond for the property
  if (this->matches(root, bond->getOtherAtom(from)))
    goto SUCCESS;
 FAIL:
  LOG(BF("FAIL"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!"));
  return this->Base::matches(root,from,bond);
}

bool AtomTest_O::matchesAm1BccX(chem::Atom_sp atom) const {
  _OF();
  Element el = atom->getElement();
  if (el == element_C) {
    // C(x3)
    return (atom->numberOfBonds() == 3 && atom->getIonization() == 0);
  }
  if (el == element_N || el == element_P) { // N(x2) || P(x2)
    if (atom->numberOfBonds() == 2 && atom->getIonization() == 0)
      return true;
    // N+(x3) || P+(x3)
    if (atom->numberOfBonds() == 3 && atom->getIonization() == 1)
      return true;
    return false;
  }
  if (el == element_O || el == element_S) {
    // O+(x2) || S+(x2)
    return (atom->numberOfBonds() == 2) && (atom->getIonization() == 1);
  }
  return false;
}

bool AtomTest_O::matchesAm1BccY(chem::Atom_sp atom) const {
  _OF();
  Element el = atom->getElement();
  if (el == element_C) {
    // C-(x2)
    return (atom->numberOfBonds() == 2 && atom->getIonization() == -1);
  }
  if (el == element_N) { // N-(x2)
    if (atom->numberOfBonds() == 2 && atom->getIonization() == -1)
      return true;
    // N(x3)
    if (atom->numberOfBonds() == 3 && atom->getIonization() == 0)
      return true;
    return false;
  }
  if (el == element_P) { // P(x3)
    if (atom->numberOfBonds() == 3 && atom->getIonization() == 0)
      return true;
    return false;
  }
  if (el == element_O || el == element_S) {
    // O(x2) || S(x2)
    return (atom->numberOfBonds() == 2) && (atom->getIonization() == 0);
  }
  return false;
}

SYMBOL_EXPORT_SC_(ChemPkg,STARcurrent_matchSTAR);

bool AtomTest_O::matches(Root_sp root, chem::Atom_sp atom) {
  _OF();
  int cnt;
  Atom_sp ringStartAtom;
  int hc = 0;
  LOG(BF("AtomTest match for atom: %s") % atom->description().c_str()); //
  switch (this->_Test) {
  case SAPWildCard:
      LOG(BF("SAPWildCard")); //
      goto SUCCESS;
  case SAPAromaticElement:
      LOG(BF("SAPAromaticElement(%s) == expecting(%s)") % _rep_(atom->getElementAsSymbol()) % _rep_(this->_SymbolArg));
      if (this->_SymbolArg == atom->getElementAsSymbol() && atom->isAromatic())
        goto SUCCESS;
      break;
  case SAPElement:
      LOG(BF("SAPElement(%s) == expecting(%s)") % _rep_(atom->getElementAsSymbol()) % _rep_(this->_SymbolArg));
      if (this->_SymbolArg == atom->getElementAsSymbol())
        goto SUCCESS;
      break;
  case SAPAliphatic:
      LOG(BF("SAPAliphatic"));
      if (!atom->isAromatic())
        goto SUCCESS;
      break;
  case SAPAromatic:
      LOG(BF("SAPAromatic"));
      if (atom->isAromatic())
        goto SUCCESS;
      break;
  case SAPPredicateName:
      LOG(BF("SAPPredicateName<%s> testing %s") % _rep_(this->_SymbolArg) % atom->description());
      if (root->evaluateTest(this->_SymbolArg, atom))
        goto SUCCESS;
      break;
  case SAPAtomicNumber:
      LOG(BF("SAPAtomicNumber(%d) == expecting(%d)") % atom->getAtomicNumber() % this->_IntArg);
      if (this->_IntArg == atom->getAtomicNumber())
        goto SUCCESS;
      break;
  case SAPTotalHCount:
      hc = atom->getBondedHydrogenCount();
      LOG(BF("SAPTotalHCount(%d) == expecting(%d)") % hc % this->_IntArg);
      LOG(BF("Checking total hydrogen count looking for(%d) got(%d)") % this->_IntArg % hc);
      if (this->_IntArg == hc) {
        LOG(BF("Match")); //
        goto SUCCESS;
      }
      LOG(BF("No match")); //
      break;
  case SAPResidueTest:
      SIMPLE_WARN(BF("What do I do with SAPResidueText"));
#if 0
      LOG(BF("SAPResidueTest looking for tag: %s") % this->_StringArg);
      if (!root->getMatch()->recognizesAtomTag(chemkw_intern(this->_StringArg))) {
        SIMPLE_ERROR(BF("We are trying to test the atomTag (" + this->_StringArg + ") but it doesn't exist!"));
      }
      LOG(BF("      SAPResidueTest check close")); //
      ringStartAtom = root->getMatch()->getAtomWithTag(chemkw_intern(this->_StringArg));
      LOG(BF("      checking if residue for %s matches residue for ringStart atom: %s") % atom->description().c_str() % ringStartAtom->description().c_str()); //
      if (atom->containedBy() == ringStartAtom->containedBy()) {
        LOG(BF("RESIDUE MATCH THEY DO!!!!"));
        goto SUCCESS;
      }
#endif
      break;
  case SAPNegativeCharge:       // next
  case SAPPositiveCharge:       // next
  case SAPNegativeFormalCharge: // next
  case SAPPositiveFormalCharge:
      LOG(BF("SAPxxxxCharge")); //
      if (this->_IntArg == atom->getIonization())
        goto SUCCESS;
      break;
  case SAPRingMembershipCount:
      LOG(BF("SAPRingMembershipCount")); //
      if (this->_IntArg == atom->getRingMembershipCount())
        goto SUCCESS;
      break;
  case SAPRingSize:
      LOG(BF("SAPRingMembershipSize")); //
      if (atom->inRingSize(this->_IntArg))
        goto SUCCESS;
      break;
  case SAPValence:
      LOG(BF("SAPRingValence"));
    // total bond valence
    // each single bond counts 1
    // each double/aromatic bond counts 2
    // each triple bond counts 3
      if (this->_IntArg == atom->getValence())
        goto SUCCESS;
      break;
  case SAPAtomMap:
    {
      // We have an atom and we have a tag in this->_IntArg
      ChemInfoMatch_sp match =  gc::As<ChemInfoMatch_sp>(_sym_STARcurrent_matchSTAR->symbolValue());
      match->defineAtomTag(atom,core::make_fixnum(this->_IntArg));
      goto SUCCESS;
    }
      break;
  case SAPConnectivity: // No implicit H's so Connectivity == Degree
  case SAPDegree:
      LOG(BF("SAPDegree testing if atom->numberOfBonds(){%d} == this->_IntArg{%d}") % atom->numberOfBonds() % this->_IntArg);
      if (this->_IntArg == atom->numberOfBonds())
        goto SUCCESS;
      break;
  case SAPInBond:
      LOG(BF("SAPInBond"));
      cnt = 0;
      for (int i = 0; i < atom->numberOfBonds(); i++) {
        if (chem::_matchBondTypes((chem::BondEnum) this->_IntArg, atom->bondedOrder(i)))
          cnt++;
      }
      if (cnt == this->_NumArg)
        goto SUCCESS;
      break;
  case SAPArLevel:
      LOG(BF("SAPArLevel"));
      switch (this->_IntArg) {
      case 1:
          if (atom->getMembershipAr1())
            goto SUCCESS;
          break;
      case 2:
          if (atom->getMembershipAr2())
            goto SUCCESS;
          break;
      case 3:
          if (atom->getMembershipAr3())
            goto SUCCESS;
          break;
      case 4:
          if (atom->getMembershipAr4())
            goto SUCCESS;
          break;
      case 5:
          if (atom->getMembershipAr5())
            goto SUCCESS;
          break;
      }
      goto FAIL;
  case SAPNoRing:
      LOG(BF("SAPNoRing"));
      if (!atom->isInRing())
        goto SUCCESS;
      break;
  case SAPElectronegativeElement: // ( O, N, F, Cl, Br )
      LOG(BF("SAPElectronegativeElement"));
      if (atom->getElement() == element_O)
        goto SUCCESS;
      if (atom->getElement() == element_N)
        goto SUCCESS;
      if (atom->getElement() == element_F)
        goto SUCCESS;
      if (atom->getElement() == element_Cl)
        goto SUCCESS;
      if (atom->getElement() == element_Br)
        goto SUCCESS;
      break;
  default:
      stringstream ss;
      ss << "You must implement a test for: " << this->testName(this->_Test);
      SIMPLE_ERROR(BF("%s") % ss.str());
  }
 FAIL:
  LOG(BF("FAIL"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!"));
  return this->Base::matches(root,atom);
}

string AtomTest_O::asSmarts() const {
  _OF();
  stringstream ss;
  switch (this->_Test) {
  case SAPWildCard:
      ss << "*";
      break;
  case SAPElement:
      ss << this->_SymbolArg->symbolName();
      break;
  case SAPAromaticElement: {
    string element = this->_SymbolArg->symbolNameAsString();
    ss << tolower(element[0]) << element.substr(1, 9999);
  } break;
  case SAPAromatic:
      ss << "a";
      break;
  case SAPAliphatic:
      ss << "A";
      break;
  case SAPPredicateName:
      ss << "<" << this->_SymbolArg->symbolNameAsString() << ">";
      break;
  case SAPAtomicNumber:
      ss << "#" << this->_IntArg;
      break;
  case SAPTotalHCount:
      ss << "H" << this->_IntArg;
      break;
  case SAPRingTest:
      ss << "?" << this->_StringArg;
      break;
  case SAPResidueTest:
      ss << "U" << this->_StringArg;
      break;
  case SAPNegativeCharge:       // next
  case SAPNegativeFormalCharge: // next
      ss << this->_IntArg;
      break;
  case SAPPositiveCharge: // next
  case SAPPositiveFormalCharge:
      ss << "+" << this->_IntArg;
      break;
  case SAPRingMembershipCount:
      ss << "R" << this->_IntArg;
      break;
  case SAPRingSize:
      ss << "r" << this->_IntArg;
      break;
  case SAPValence:
      ss << "v" << this->_IntArg;
      break;
  case SAPConnectivity: // No implicit H's so Connectivity == Degree
  case SAPDegree:
      ss << "X" << this->_IntArg;
      break;
  case SAPInBond:
      ss << "{SAPInBond}";
      break;
  case SAPArLevel:
      ss << "{SAPArLevel}";
      break;
  case SAPNoRing:
      ss << "{SAPNoRing}";
      break;
  case SAPElectronegativeElement: // ( O, N, F, Cl, Br )
      ss << "{SAPElectronegativeElement}";
      break;
  default:
      ss << "{UnknownTest-" << this->testName(this->_Test) << "}";
      break;
  }
  return ss.str();
}

core::NullTerminatedEnumAssociation testEnum[] = {
    {"SAPNone", SAPNone},
    {"SAPWildCard", SAPWildCard},
    {"SAPDegree", SAPDegree},
    {"SAPElement", SAPElement},
    {"SAPAromaticElement", SAPAromaticElement},
    {"SAPAromatic", SAPAromatic},
    {"SAPAliphatic", SAPAliphatic},
    {"SAPTotalHCount", SAPTotalHCount},
    {"SAPImplicitHCount", SAPImplicitHCount},
    {"SAPRingTest", SAPRingTest},
    {"SAPRingMembershipCount", SAPRingMembershipCount},
    {"SAPRingSize", SAPRingSize},
    {"SAPValence", SAPValence},
    {"SAPConnectivity", SAPConnectivity},
    {"SAPNegativeCharge", SAPNegativeCharge},
    {"SAPNegativeFormalCharge", SAPNegativeFormalCharge},
    {"SAPPositiveCharge", SAPPositiveCharge},
    {"SAPPositiveFormalCharge", SAPPositiveFormalCharge},
    {"SAPAtomicNumber", SAPAtomicNumber},
    {"SAPChiralityAntiClockwise", SAPChiralityAntiClockwise},
    {"SAPChiralityClockwise", SAPChiralityClockwise},
    {"SAPAtomicMass", SAPAtomicMass},
    {"SAPLonePair", SAPLonePair},
    {"SAPTotalBondNumber", SAPTotalBondNumber},
    {"SAPTotalExplicitBondNumber", SAPTotalExplicitBondNumber},
    {"SAPPiBondOrbital", SAPPiBondOrbital},
    {"SAPAromaticPiElectron", SAPAromaticPiElectron},
    {"SAPHeavyAtomTotalBond", SAPHeavyAtomTotalBond},
    {"SAPGroupNumber", SAPGroupNumber},
    {"SAPElectronegativeElement", SAPElectronegativeElement},
    {"SAPTransitionMetal", SAPTransitionMetal},
    {"SAPBondedToPrevious", SAPBondedToPrevious},
    {"SAPNotBondedToPrevious", SAPNotBondedToPrevious},
    {"SAPInBond", SAPInBond},
    {"SAPArLevel", SAPArLevel},
    {"SAPNoRing", SAPNoRing},
    {"SAPResidueTest", SAPResidueTest},
    {"SAPPredicateName", SAPPredicateName},
    {"SAPAtomMap", SAPAtomMap},
    {"", -1}};

string AtomTest_O::testName(AtomTestEnum test) const {
  
  for (uint i = 0; testEnum[i]._Enum != -1; i++) {
    if (testEnum[i]._Enum == test) {
      return testEnum[i]._Key;
    }
  }
  return "-unknownTest-";
}

void AtomTest_O::fields(core::Record_sp node) {
  node->/*pod_*/field( INTERN_(kw,test), this->_Test);
  node->/*pod_*/field_if_not_default( INTERN_(kw,int), this->_IntArg, 0);
  node->/*pod_*/field_if_not_default( INTERN_(kw,num), this->_NumArg, 0);
  node->/*pod_*/field_if_not_default( INTERN_(kw,str), this->_StringArg, std::string());
  node->field_if_not_nil( INTERN_(kw,sym), this->_SymbolArg);
  this->Base::fields(node);
}

// ------- Chain

void Chain_O::initialize() {
  this->Base::initialize();
  this->_Head = _Nil<core::T_O>();
  this->_Tail = _Nil<core::T_O>();
}

string Chain_O::asSmarts() const {
  _OF();
  stringstream ss;
  if (this->_Head.notnilp())
    ss << "(" << this->_Head->asSmarts() << ")";
  if (this->_Tail.notnilp())
    ss << "(" << this->_Tail->asSmarts() << ")";
  return ss.str();
}

core::T_sp Chain_O::children() {
  ql::list result;
  if (this->_Head.notnilp()) result << this->_Head;
  if (this->_Tail.notnilp()) result << this->_Tail;
  return result.cons();
}

uint Chain_O::depth() const {
  _OF();
  return (MAX(af_depth(this->_Head), af_depth(this->_Tail) + 1));
}

bool Chain_O::matches(Root_sp root, chem::Atom_sp from) {
  _OF();
  chem::BondList_sp bonds = from->getBondList();
  return this->matches(root,from,bonds);
}

bool Chain_O::matches(Root_sp root, chem::Atom_sp from, chem::BondList_sp neighbors) {
  _OF();
  gctools::Vec0<chem::Bond_sp>::iterator bi;
  chem::BondList_sp nextBonds, tempBondList;
  LOG(BF("Chain_O matching pattern: %s") % this->asSmarts());
  LOG(BF("There are %d neighbors bondList: %s") % neighbors->size() % neighbors->describeOthers(from));
  for (bi = neighbors->begin(); bi != neighbors->end(); bi++) {
    CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
    _BLOCK_TRACEF(BF("Checking neighbor for bond: %s") % (*bi)->describeOther(from));
    if (this->_Head->matches(root, from, *bi)) {
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("The head matches"));
      Atom_sp other = (*bi)->getOtherAtom(from);
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("The other atom: %s") % other->description());
      nextBonds = other->getBondList(); // handle new bonds
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("The others bond list before removed from: %s") % nextBonds->describeOthers(other));
      LOG(BF("Removing bond between %s and %s") % from->description() % other->description());
      nextBonds->removeBondBetween(from, other); // (*bi)->getFrom());
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      LOG(BF("The others bond list after removed from: %s") % nextBonds->describeOthers(other));
      ASSERTNOTNULL(this->_Tail);
      if (this->_Tail.notnilp()) {
        CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
        if (this->_Tail->matches(root, other, nextBonds))
          goto SUCCESS;
        LOG(BF("The tail exists but doesn't match"));
      } else {
        CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
        goto SUCCESS;
      }
    }
  }
  //FAIL:
  LOG(BF("FAIL"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
}

void Chain_O::fields(core::Record_sp node) {
  node->field_if_not_nil( INTERN_(kw,tail), this->_Tail);
  node->field_if_not_nil( INTERN_(kw,head), this->_Head);
  this->Base::fields(node);
}

CL_DEFMETHOD BondListMatchNode_sp Chain_O::chain_get_tail(){
  return this->_Tail;
}



// ------ Branch

void Branch_O::initialize() {
  this->Base::initialize();
  this->_Left = _Nil<core::T_O>();
  this->_Right = _Nil<core::T_O>();
}
string Branch_O::asSmarts() const {
  _OF();
  stringstream ss;

  if (this->_Left.notnilp())
    ss << "(" << this->_Left->asSmarts() << ")";
  if (this->_Right.notnilp())
    ss << "(" << this->_Right->asSmarts() << ")";
  return ss.str();
}

core::T_sp Branch_O::children() {
  ql::list result;
  if (this->_Left.notnilp()) result << this->_Left;
  if (this->_Right.notnilp()) result << this->_Right;
  return result.cons();
}


uint Branch_O::depth() const {
  
  return (MAX(af_depth(this->_Left), af_depth(this->_Right)));
}

bool Branch_O::matches(Root_sp root, chem::Atom_sp from, chem::BondList_sp neighbors) {
  _OF();
  gctools::Vec0<chem::Bond_sp>::iterator bi;
  chem::BondList_sp nextBonds;
  chem::BondList_sp leftBondList;
  LOG(BF("Branch_O matching pattern: %s") % this->asSmarts());
  LOG(BF("Neighbors bond list= %s") % neighbors->describeOthers(from));
  for (bi = neighbors->begin(); bi != neighbors->end(); bi++) {
    _BLOCK_TRACEF(BF("Checking neighbors focusing on %s") % (*bi)->describeOther(from));
    leftBondList = BondList_O::create();
    leftBondList->addBond(*bi);
    LOG(BF("Constructed left bond list"));
    LOG(BF("Left bond list = %s") % leftBondList->describeOthers(from));
    BondList_sp rightBondList = BondList_O::create();
    for ( auto bi = (*neighbors).begin(); bi!=(*neighbors).end(); ++bi ) {
      rightBondList->addBond(*bi);
    }
//    GC_COPY(chem::BondList_O, rightBondList, *neighbors); // = RP_Copy<BondList_O>(neighbors);
    LOG(BF("Right bond list after copy = %s") % rightBondList->describeOthers(from));
    rightBondList->removeBond(*bi);
    LOG(BF("copied neighbors into right bond list and removed left bond"));
    LOG(BF("Right bond list after remove left = %s") % rightBondList->describeOthers(from));
    LOG(BF("Checking if left matches"));
    if (this->_Left->matches(root, from, leftBondList)) {
      LOG(BF("Left matches"));
      ANN(this->_Right);
      if (this->_Right.notnilp()) {
        LOG(BF("Right branch is defined, checking if it matches"));
        if (this->_Right->matches(root, from, rightBondList)) {
          LOG(BF("Right matches does"));
          goto SUCCESS;
        } else {
          LOG(BF("Right doesn't match, keep trying"));
        }
      } else {
        LOG(BF("Left matches and there is no Right"));
        goto SUCCESS;
      }
    }
  }
  LOG(BF("No match"));
  //FAIL:
  LOG(BF("FAIL"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
}

void Branch_O::fields(core::Record_sp node) {
  node->field( INTERN_(kw,left), this->_Left);
  node->field_if_not_nil( INTERN_(kw,right), this->_Right);
  this->Base::fields(node);
}

CL_DEFMETHOD void Branch_O::branch_set_left(core::T_sp left){
  this->_Left = left;
}

CL_DEFMETHOD void Branch_O::branch_set_right(core::T_sp right){
  this->_Right = right;
}


//      AfterMatchBondNode

void RootMatchNode_O::fields(core::Record_sp node) {
  this->Base::fields(node);
// Nothing to do
}

// ------- AfterMatchBondTest

void AfterMatchBondTest_O::initialize() {
  _OF();
  this->_AtomTag1 = _Nil<core::Symbol_O>();
  this->_AtomTag2 = _Nil<core::Symbol_O>();
  this->_Bond = SABNoBond;
}

void AfterMatchBondTest_O::fields(core::Record_sp node) {
  node->field( INTERN_(kw,tag1), this->_AtomTag1);
  node->field( INTERN_(kw,tag2), this->_AtomTag2);
  node->/*pod_*/field( INTERN_(kw,bond), this->_Bond);
  this->Base::fields(node);
}

bool AfterMatchBondTest_O::matches(Root_sp root) {
  _OF();
  SIMPLE_ERROR(BF("Must implement"));
}

// -------- AntechamberFocusAtomMatch

void AntechamberFocusAtomMatch_O::initialize() {
  this->Base::initialize();
  this->_ResidueNames = _Nil<core::T_O>();
  this->_AtomicNumber = 0;
  this->_NumberOfAttachedAtoms = 0;
  this->_NumberOfAttachedHydrogens = 0;
  this->_NumberOfElectronWithdrawingGroups = 0;
  this->_AtomicProperty = _Nil<core::T_O>();
}

void AntechamberFocusAtomMatch_O::fields(core::Record_sp node) {
  node->field_if_not_nil( INTERN_(kw,residueNames), this->_ResidueNames);
  node->field_if_not_nil( INTERN_(kw,atomicProperty), this->_AtomicProperty);
  node->/*pod_*/field( INTERN_(kw,atomicNumber), this->_AtomicNumber);
  node->/*pod_*/field( INTERN_(kw,attachedAtoms), this->_NumberOfAttachedAtoms);
  node->/*pod_*/field( INTERN_(kw,attachedHs), this->_NumberOfAttachedHydrogens);
  node->/*pod_*/field( INTERN_(kw,attachedEWGs), this->_NumberOfElectronWithdrawingGroups);
  this->Base::fields(node);
}

core::T_sp AntechamberFocusAtomMatch_O::children() {
  ql::list result;
  if (this->_AtomicProperty.notnilp()) result << this->_AtomicProperty;
  return result.cons();
}

string AntechamberFocusAtomMatch_O::asSmarts() const {
  stringstream ss;
  ss << "{AntechamberFocusAtomMatch ";
  ss << _rep_(this->_ResidueNames) << "|";
  if (this->_AtomicProperty.notnilp()) {
    ss << this->_AtomicProperty->asSmarts() << "|";
  }
  ss << this->_AtomicNumber << "|";
  ss << this->_NumberOfAttachedAtoms << "|";
  ss << this->_NumberOfAttachedHydrogens << "|";
  ss << this->_NumberOfElectronWithdrawingGroups << "}";
  return ss.str();
}

bool AntechamberFocusAtomMatch_O::matches(Root_sp root, chem::Atom_sp atom) {
  
  chem::Atom_sp neighbor, nn;
  if (this->_AtomicNumber >= 0) {
    LOG(BF("Checking if atomic number(%d) == expected(%d)") % atom->getAtomicNumber() % this->_AtomicNumber);
    if (atom->getAtomicNumber() != this->_AtomicNumber)
      goto FAIL;
  }
  if (this->_NumberOfAttachedAtoms >= 0) {
    LOG(BF("Checking if number of bonds(%d) == expected(%d)") % atom->numberOfBonds() % this->_NumberOfAttachedAtoms);
    if (atom->numberOfBonds() != this->_NumberOfAttachedAtoms)
      goto FAIL;
  }
  if (this->_NumberOfAttachedHydrogens >= 0) {
    LOG(BF("Checking if number of attached hydrogens(%d) == expected(%d)") % atom->getBondedHydrogenCount() % this->_NumberOfAttachedHydrogens);
    if (atom->getBondedHydrogenCount() != this->_NumberOfAttachedHydrogens)
      goto FAIL;
  }
  if (this->_NumberOfElectronWithdrawingGroups >= 0) {
    if (atom->numberOfBonds() != 1) {
      SIMPLE_WARN(BF("Testing NumberOfElectronWithdrawingGroups must have one bond - test %s") % _rep_(this->asSmartPtr()));
      goto FAIL;
    }
    neighbor = atom->bondedNeighbor(0);
    int cnt = 0;
    for (int i = 0; i < neighbor->numberOfBonds(); i++) {
      nn = neighbor->bondedNeighbor(i);
      if (nn->getElement() == element_N)
        cnt++;
      if (nn->getElement() == element_O)
        cnt++;
      if (nn->getElement() == element_F)
        cnt++;
      if (nn->getElement() == element_Cl)
        cnt++;
      if (nn->getElement() == element_Br)
        cnt++;
    }
    LOG(BF("Checking if number of electron withdrawing groups(%d) == expected(%d)") % cnt % this->_NumberOfElectronWithdrawingGroups);
    if (cnt != this->_NumberOfElectronWithdrawingGroups)
      goto FAIL;
  }
  ANN(this->_AtomicProperty);
  if (this->_AtomicProperty.notnilp()) {
    if (!this->_AtomicProperty->matches(root, atom))
      goto FAIL;
  }
  //SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
 FAIL:
  LOG(BF("FAIL"));
  return false;
}

// ------- AntechamberBondTest

core::T_sp AntechamberBondTest_O::children() {
  ql::list result;
  if (this->_AtomProperties.notnilp()) result << this->_AtomProperties;
  return result.cons();
}

bool AntechamberBondTest_O::matchBasic(AntechamberRoot_sp root, chem::Atom_sp atom) {
  
  WildElementDict_sp dict;
  bool gotElement;
  dict = root->getElementWildCardDictionary();
  gotElement = false;
  ANN(dict);
  if (dict.notnilp()) {
    if (dict->recognizesWildName(this->_Element)) {
      LOG(BF("Matching wildcard element has(%s) == expecting (%s)") % atom->getElementAsString().c_str() % this->_Element.c_str());
      if (!dict->recognizesWildNameElement(this->_Element, atom->getElementAsSymbol()))
        goto FAIL;
      gotElement = true;
    }
  }
  if (!gotElement) {
    LOG(BF("Trying to match regular element(%s) == expected(%s)") % atom->getElementAsString().c_str() % _rep_(this->_Element));
    if (atom->getElementAsSymbol() != this->_Element)
      goto FAIL;
  }
  LOG(BF("Trying to match number of bonds(%d) == expected(%d)") % atom->numberOfBonds() % this->_Neighbors);
  if (this->_Neighbors >= 0) {
    if (atom->numberOfBonds() != this->_Neighbors)
      goto FAIL;
  }
  if (this->_Tag.notnilp()) {
    SIMPLE_WARN(BF("%s:%d What do I do with the tag %s for AntechamberBondTest %s") % __FILE__ % __LINE__ % _rep_(this->_Tag) % _rep_(this->asSmartPtr()));
//    root->getMatch()->defineAtomTag(atom, this->_Tag);
  }
  //SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
 FAIL:
  LOG(BF("FAIL"));
  return false;
}

bool AntechamberBondTest_O::matches(Root_sp root, chem::Atom_sp atom) {
  
  AntechamberRoot_sp acRoot;
  if (root->type() != antechamberRoot) {
    SIMPLE_ERROR(BF("AntechamberBondTest::matches requires an AntechamberRoot"));
  }
  acRoot = (root).as<AntechamberRoot_O>();
  if (!this->matchBasic(acRoot, atom))
    goto FAIL;
  if (this->_AtomProperties->matches(acRoot, atom))
    goto SUCCESS;
 FAIL:
  LOG(BF("FAIL"));
  return false;
 SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
}

string AntechamberBondTest_O::asSmarts() const {
  stringstream ss;
  ss << "{AntechamberBondTest";
  ss << _rep_(this->_Element) << "|";
  ss << (this->_Neighbors) <<"|";
  ss << this->_AtomProperties->asSmarts() << "}";
  return ss.str();
}

bool AntechamberBondTest_O::matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  AntechamberRoot_sp antechamberRoot;
  if (root->type() != chem::antechamberRoot) {
    stringstream ss;
    ss << "AntechamberBondTest::matches requires an AntechamberRoot "
       << "instead we have" << (unsigned long)(root.get()) << " of type: " << root->type();
    SIMPLE_ERROR(BF("%s") % ss.str());
  }
  antechamberRoot = (root).as<AntechamberRoot_O>();
  if (!this->matchBasic(antechamberRoot, bond->getOtherAtom(from)))
    goto FAIL;
  ANN(this->_AtomProperties);
  if (this->_AtomProperties.notnilp()) {
    if (!this->_AtomProperties->matches(antechamberRoot, from, bond))
      goto FAIL;
  }
  LOG(BF("Matched!!"));
  //SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
 FAIL:
  LOG(BF("FAIL"));
  return false;
}

void AntechamberBondTest_O::initialize() {
  this->Base::initialize();
  this->_Element = _Nil<core::Symbol_O>();
  this->_Neighbors = 0;
  this->_Tag = _Nil<core::Symbol_O>();
  this->_AtomProperties = _Nil<core::T_O>();
}

void AntechamberBondTest_O::fields(core::Record_sp node) {
  node->field( INTERN_(kw,element), this->_Element);
  node->/*pod_*/field( INTERN_(kw,neighbors), this->_Neighbors);
  node->field( INTERN_(kw,tag), this->_Tag);
  node->field_if_not_nil( INTERN_(kw,atomProps), this->_AtomProperties);
  this->Base::fields(node);
}

// ------- Root

uint Root_O::depth() const {
  _OF();
  uint res = MAX(af_depth(this->_FirstTest), af_depth(this->_Chain) + 1);
  LOG(BF("Returning %d") % res);
  return res;
}

core::T_sp Root_O::children() {
  ql::list result;
  if (this->_FirstTest.notnilp()) result << this->_FirstTest;
  if (this->_Chain.notnilp()) result << this->_Chain;
  return result.cons();
};

string Root_O::asSmarts() const {
  stringstream ss;
  if (this->_FirstTest.notnilp()) {
    ss << this->_FirstTest->asSmarts();
  }
  if (this->_Chain.notnilp()) {
    ss << "(" << this->_Chain->asSmarts() << ")";
  }
  return ss.str();
}

void Root_O::initialize() {
  this->Base::initialize();
  this->_FirstTest = _Nil<core::T_O>();
  this->_Chain = _Nil<core::T_O>();
  this->_Tests = core::HashTableEq_O::create_default();
}

void Root_O::addTest(core::Symbol_sp testSym, core::Function_sp testCode) {
  _OF();
  LOG(BF("Adding test<%s> with code: %s") % _rep_(testSym) % _rep_(testCode));
  this->_Tests->setf_gethash(testSym, testCode);
}

bool Root_O::evaluateTest(core::Symbol_sp testSym, Atom_sp atom) {
  ASSERTF(testSym.notnilp(), BF("The test symbol was nil! - this should never occur"));
  LOG(BF("Looking up test with symbol<%s>") % _rep_(testSym));
  core::T_mv find = this->_Tests->gethash(testSym);
  if (find.second().nilp()) {
      SIMPLE_ERROR(BF("Could not find named ChemInfo/Smarts test[%s] in Smarts object - available named tests are[%s]") % _rep_(testSym) % this->_Tests->keysAsString());
  }
  if (!gctools::IsA<core::Function_sp>(find)) {
      SIMPLE_ERROR(BF("The test ChemInfo/Smarts test[%s] must be a function - instead it is a %s") % _rep_(testSym) % _rep_(find));
  }
  core::Function_sp testCode = gctools::As_unsafe<core::Function_sp>(find);
  ASSERTF(testCode.notnilp(), BF("testCode was nil - it should never be"));
  ASSERTF(atom.notnilp(), BF("The atom arg should never be nil"));
  core::T_sp res = core::eval::funcall(testCode,atom);
  return res.isTrue();
}

void Root_O::fields(core::Record_sp node) {
  node->field_if_not_nil( INTERN_(kw,firstTest), this->_FirstTest);
  node->field_if_not_nil( INTERN_(kw,chain), this->_Chain);
  node->field_if_not_nil( INTERN_(kw,tests), this->_Tests);
#ifdef DEBUG_ON
  LOG(BF("After load chain=%s") % this->_Chain->description().c_str());
#endif
  this->Base::fields(node);
}

bool Root_O::matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  _OF();
  LOG(BF("Root_O trying to match pattern: %s") % this->asSmarts());
  return this->matches(root, bond->getOtherAtom(from));
};

bool Root_O::matches(Root_sp root, chem::Atom_sp atom) {
  chem::BondList_sp nextBonds;
  bool matches;
  matches = false;
  ANN(this->_FirstTest);
  if (this->_FirstTest.notnilp()) {
    CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
    LOG(BF("_FirstTest is notNil - testing"));
    if (!this->_FirstTest->matches(root, atom)) {
      CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
      goto FAIL;
    }
  }
  CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
  ANN(this->_Chain);
  if (this->_Chain.notnilp()) {
    LOG(BF("_Chain is not Nil - testing "));
    //	nextBonds = atom->getHeavyAtomBondList();
    nextBonds = atom->getBondList();
    if (!this->_Chain->matches(root, atom, nextBonds))
      goto FAIL;
  }
  //SUCCESS:
  LOG(BF("SUCCESS!!!"));
  return true;
 FAIL:
  LOG(BF("FAIL"));
  return false;
}

CL_LISPIFY_NAME("make-smarts-root");
CL_DEF_CLASS_METHOD SmartsRoot_sp SmartsRoot_O::make(ChemInfoNode_sp cinode)
{
  GC_ALLOCATE(SmartsRoot_O,   obj ); // RP_Create<SmartsRoot_O>(lisp);
  if (gc::IsA<Chain_sp>(cinode)) {
    Chain_sp chain = gc::As<Chain_sp>(cinode);
    obj->_FirstTest =  gc::As<AtomOrBondMatchNode_sp>(chain->chain_get_head());
    obj->_Chain = gc::As<BondListMatchNode_sp>(chain->chain_get_tail());
    return obj;
  } else if (gc::IsA<AtomOrBondMatchNode_sp>(cinode)) {
    obj->_FirstTest = gc::As<AtomOrBondMatchNode_sp>(cinode);
    obj->_Chain = _Nil<core::T_O>();
    return obj;
  }
  SIMPLE_ERROR(BF("You cannot create a smarts-root with %s") % _rep_(cinode));
};





void SmartsRoot_O::fields(core::Record_sp node) {
  this->Base::fields(node);
}

bool SmartsRoot_O::matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  return this->matches(root, bond->getOtherAtom(from));
};

bool SmartsRoot_O::matches(Root_sp root, chem::Atom_sp atom) {
  _OF();
  chem::BondList_sp nextBonds;
  bool matches;
  matches = false;
  LOG(BF("SmartsRoot matching pattern: %s") % this->asSmarts());
  LOG(BF("SmartsRoot match for atom: %s") % atom->description());
  if (!this->Root_O::matches(root, atom))
    goto FAIL;
  //SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
 FAIL:
  LOG(BF("FAIL"));
  return false;
}

void AntechamberRoot_O::initialize() {
  this->Base::initialize();
  this->_AssignType = _Nil<core::Symbol_O>();
  this->_AfterMatchTests = _Nil<core::T_O>();
  this->_WildElementDictionary = _Nil<core::T_O>();
}

void AntechamberRoot_O::fields(core::Record_sp node) {
  node->field( INTERN_(kw,assignType), this->_AssignType);
  node->field_if_not_nil( INTERN_(kw,afterMatchTests), this->_AfterMatchTests);
  node->field_if_not_nil( INTERN_(kw,wildDict), this->_WildElementDictionary);
  this->Base::fields(node);
}

core::T_sp AntechamberRoot_O::children() {
  ql::list result;
  result & this->Base::children();
  if (this->_AfterMatchTests.notnilp()) result << this->_AfterMatchTests;
  return result.cons();
};
    

string AntechamberRoot_O::descriptionOfContents() const {
  stringstream ss;
  ss << " :assign-type " << _rep_(this->_AssignType);
  return ss.str();
}

bool AntechamberRoot_O::matches(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  return this->matches(root, bond->getOtherAtom(from));
};

bool AntechamberRoot_O::matches(Root_sp root, chem::Atom_sp atom) {
  chem::BondList_sp nextBonds;
  bool matches;
  matches = false;
  if (!this->Base::matches(root, atom)) {
    goto FAIL;
  }
  ANN(this->_AfterMatchTests);
  if (this->_AfterMatchTests.notnilp()) {
    if (!this->_AfterMatchTests->matches(root)) {
      goto FAIL;
    }
  }
  // printf("%s:%d This is where I was assigning the type\n", __FILE__, __LINE__ );
//  atom->setType(this->_AssignType);
  LOG(BF("Matched and assigned type(%s)") % this->_AssignType.c_str());
  //SUCCESS:
  LOG(BF("SUCCESS!"));
  return true;
 FAIL:
  LOG(BF("FAIL"));
  return false;
}

#define ARGS_ChemInfo_O_make "(tests smarts)"
#define DECL_ChemInfo_O_make ""
#define DOCS_ChemInfo_O_make "make ChemInfo"
ChemInfo_sp ChemInfo_O::make(core::List_sp tests, const string &smarts) {
  
  GC_ALLOCATE(ChemInfo_O, me);
  me->compileSmarts(smarts);
  me->defineTests(tests);
  if (!me->compileSucceeded()) {
    SIMPLE_ERROR(BF("%s") % me->compilerMessage());
  }
  return me;
};

CL_LISPIFY_NAME("asSmarts");
CL_DEFMETHOD string ChemInfo_O::asSmarts() const {
  _OF();
  return this->_Root->asSmarts();
}



CL_DEFUN core::T_mv chem__chem_info_match(Root_sp testRoot, Atom_sp atom)
{
  core::HashTableEql_sp ringHashTable = core::HashTableEql_O::create_default();
  core::HashTableEql_sp colonOperatorHashTable = core::HashTableEql_O::create_default();
  ChemInfoMatch_sp current_match = ChemInfoMatch_O::make( false, colonOperatorHashTable, ringHashTable);
  core::DynamicScopeManager scope(_sym_STARcurrent_matchSTAR,current_match);
  bool matches = testRoot->matches(testRoot,atom);
  CI_LOG(("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ ));
  return Values(_lisp->_boolean(matches), current_match);
}

CL_DEFUN core::T_sp chem__chem_info_node_children(ChemInfoNode_sp node) {
  return node->children();
}

  

  /*! Hold nodes for the Gaff and Msmarts parsers - rewrite these in Common Lisp */
SYMBOL_EXPORT_SC_(ChemPkg, STARparserNodeHolderSTAR );

SYMBOL_EXPORT_SC_(ChemPkg,STARSarRingTestEnumConverterSTAR);
SYMBOL_EXPORT_SC_(KeywordPkg,SARNone);
SYMBOL_EXPORT_SC_(KeywordPkg,SARRingSet);
SYMBOL_EXPORT_SC_(KeywordPkg,SARRingTest);
SYMBOL_EXPORT_SC_(ChemPkg,STARSabBondEnumConverterSTAR);
CL_BEGIN_ENUM(RingTestEnum,_sym_STARSarRingTestEnumConverterSTAR,"SARRingTestEnum");
CL_VALUE_ENUM(kw::_sym_SARNone, SARNone);
CL_VALUE_ENUM(kw::_sym_SARRingSet, SARRingSet);
CL_VALUE_ENUM(kw::_sym_SARRingTest, SARRingTest);
CL_END_ENUM(_sym_STARSarRingTestEnumConverterSTAR);

SYMBOL_EXPORT_SC_(KeywordPkg,SABNoBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABSingleBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABSingleOrAromaticBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABDoubleOrAromaticBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABTripleOrAromaticBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABDelocalizedBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABDoubleBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABTripleBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABAromaticBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABAnyBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SABDirectionalSingleUpOrUnspecified);
SYMBOL_EXPORT_SC_(KeywordPkg,SABDirectionalSingleDownOrUnspecified);
SYMBOL_EXPORT_SC_(KeywordPkg,SABDirectionalSingleUp);
SYMBOL_EXPORT_SC_(KeywordPkg,SABDirectionalSingleDown);
SYMBOL_EXPORT_SC_(ChemPkg,STARSabBondEnumConverterSTAR);
CL_BEGIN_ENUM(BondEnum,_sym_STARSabBondEnumConverterSTAR,"SABBondEnum");
CL_VALUE_ENUM(kw::_sym_SABNoBond, SABNoBond );
CL_VALUE_ENUM(kw::_sym_SABSingleBond, SABSingleBond);
CL_VALUE_ENUM(kw::_sym_SABSingleOrAromaticBond, SABSingleOrAromaticBond);
CL_VALUE_ENUM(kw::_sym_SABDoubleOrAromaticBond, SABDoubleOrAromaticBond);
CL_VALUE_ENUM(kw::_sym_SABTripleOrAromaticBond, SABTripleOrAromaticBond);
CL_VALUE_ENUM(kw::_sym_SABDelocalizedBond, SABDelocalizedBond);
CL_VALUE_ENUM(kw::_sym_SABDoubleBond, SABDoubleBond);
CL_VALUE_ENUM(kw::_sym_SABTripleBond, SABTripleBond);
CL_VALUE_ENUM(kw::_sym_SABAromaticBond, SABAromaticBond);
CL_VALUE_ENUM(kw::_sym_SABAnyBond, SABAnyBond);
CL_VALUE_ENUM(kw::_sym_SABDirectionalSingleUpOrUnspecified, SABDirectionalSingleUpOrUnspecified);
CL_VALUE_ENUM(kw::_sym_SABDirectionalSingleDownOrUnspecified, SABDirectionalSingleDownOrUnspecified);
CL_VALUE_ENUM(kw::_sym_SABDirectionalSingleUp, SABDirectionalSingleUp);
CL_VALUE_ENUM(kw::_sym_SABDirectionalSingleDown, SABDirectionalSingleDown);
CL_END_ENUM(_sym_STARSabBondEnumConverterSTAR);

SYMBOL_EXPORT_SC_(KeywordPkg,logAlwaysTrue);
SYMBOL_EXPORT_SC_(KeywordPkg,logIdentity);
SYMBOL_EXPORT_SC_(KeywordPkg,logNot);
SYMBOL_EXPORT_SC_(KeywordPkg,logHighPrecedenceAnd);
SYMBOL_EXPORT_SC_(KeywordPkg,logOr);
SYMBOL_EXPORT_SC_(KeywordPkg,logLowPrecedenceAnd);
SYMBOL_EXPORT_SC_(ChemPkg, STARLogicalOperatorTypeConverterSTAR);
CL_BEGIN_ENUM(LogicalOperatorType,_sym_STARLogicalOperatorTypeConverterSTAR,"LogicalOperatorType");
CL_VALUE_ENUM(kw::_sym_logAlwaysTrue, logAlwaysTrue);
CL_VALUE_ENUM(kw::_sym_logIdentity,  logIdentity);
CL_VALUE_ENUM(kw::_sym_logNot,  logNot);
CL_VALUE_ENUM(kw::_sym_logHighPrecedenceAnd,  logHighPrecedenceAnd);
CL_VALUE_ENUM(kw::_sym_logOr,  logOr);
CL_VALUE_ENUM(kw::_sym_logLowPrecedenceAnd,  logLowPrecedenceAnd);
CL_END_ENUM(_sym_STARLogicalOperatorTypeConverterSTAR);





SYMBOL_EXPORT_SC_(KeywordPkg,SAPNone);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPWildCard);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPDegree);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPElement);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPTotalHCount);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPImplicitHCount);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPRingMembershipCount);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPRingTest);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPRingSize);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPValence);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPConnectivity);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPNegativeCharge);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPNegativeFormalCharge);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPPositiveCharge);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPPositiveFormalCharge);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPAtomicNumber);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPChiralityAntiClockwise);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPChiralityClockwise);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPAtomicMass);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPLonePair);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPTotalBondNumber);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPTotalExplicitBondNumber);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPPiBondOrbital);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPAromaticPiElectron);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPHeavyAtomTotalBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPGroupNumber);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPElectronegativeElement);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPTransitionMetal);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPBondedToPrevious);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPNotBondedToPrevious);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPInBond);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPArLevel);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPNoRing);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPResidueTest);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPPredicateName);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPAromaticElement);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPAliphatic);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPAromatic);
SYMBOL_EXPORT_SC_(KeywordPkg,SAPAtomMap);

SYMBOL_EXPORT_SC_(ChemPkg, STARAtomTestEnumConverterSTAR);
CL_BEGIN_ENUM(AtomTestEnum,_sym_STARAtomTestEnumConverterSTAR,"AtomTestEnum");
CL_VALUE_ENUM(kw::_sym_SAPNone,SAPNone);
CL_VALUE_ENUM(kw::_sym_SAPWildCard,SAPWildCard);
CL_VALUE_ENUM(kw::_sym_SAPDegree,SAPDegree);
CL_VALUE_ENUM(kw::_sym_SAPElement,SAPElement);
CL_VALUE_ENUM(kw::_sym_SAPTotalHCount,SAPTotalHCount);
CL_VALUE_ENUM(kw::_sym_SAPImplicitHCount,SAPImplicitHCount);
CL_VALUE_ENUM(kw::_sym_SAPRingMembershipCount,SAPRingMembershipCount);
CL_VALUE_ENUM(kw::_sym_SAPRingTest,SAPRingTest);
CL_VALUE_ENUM(kw::_sym_SAPRingSize,SAPRingSize);
CL_VALUE_ENUM(kw::_sym_SAPValence,SAPValence);
CL_VALUE_ENUM(kw::_sym_SAPConnectivity,SAPConnectivity);
CL_VALUE_ENUM(kw::_sym_SAPNegativeCharge,SAPNegativeCharge);
CL_VALUE_ENUM(kw::_sym_SAPNegativeFormalCharge,SAPNegativeFormalCharge);
CL_VALUE_ENUM(kw::_sym_SAPPositiveCharge,SAPPositiveCharge);
CL_VALUE_ENUM(kw::_sym_SAPPositiveFormalCharge,SAPPositiveFormalCharge);
CL_VALUE_ENUM(kw::_sym_SAPAtomicNumber,SAPAtomicNumber);
CL_VALUE_ENUM(kw::_sym_SAPChiralityAntiClockwise,SAPChiralityAntiClockwise);
CL_VALUE_ENUM(kw::_sym_SAPChiralityClockwise,SAPChiralityClockwise);
CL_VALUE_ENUM(kw::_sym_SAPAtomicMass,SAPAtomicMass);
CL_VALUE_ENUM(kw::_sym_SAPLonePair,SAPLonePair);
CL_VALUE_ENUM(kw::_sym_SAPTotalBondNumber,SAPTotalBondNumber);
CL_VALUE_ENUM(kw::_sym_SAPTotalExplicitBondNumber,SAPTotalExplicitBondNumber);
CL_VALUE_ENUM(kw::_sym_SAPPiBondOrbital,SAPPiBondOrbital);
CL_VALUE_ENUM(kw::_sym_SAPAromaticPiElectron,SAPAromaticPiElectron);
CL_VALUE_ENUM(kw::_sym_SAPHeavyAtomTotalBond,SAPHeavyAtomTotalBond);
CL_VALUE_ENUM(kw::_sym_SAPGroupNumber,SAPGroupNumber);
CL_VALUE_ENUM(kw::_sym_SAPElectronegativeElement,SAPElectronegativeElement);
CL_VALUE_ENUM(kw::_sym_SAPTransitionMetal,SAPTransitionMetal);
CL_VALUE_ENUM(kw::_sym_SAPBondedToPrevious,SAPBondedToPrevious);
CL_VALUE_ENUM(kw::_sym_SAPNotBondedToPrevious,SAPNotBondedToPrevious);
CL_VALUE_ENUM(kw::_sym_SAPInBond,SAPInBond);
CL_VALUE_ENUM(kw::_sym_SAPArLevel,SAPArLevel);
CL_VALUE_ENUM(kw::_sym_SAPNoRing,SAPNoRing);
CL_VALUE_ENUM(kw::_sym_SAPResidueTest,SAPResidueTest);
CL_VALUE_ENUM(kw::_sym_SAPPredicateName,SAPPredicateName);
CL_VALUE_ENUM(kw::_sym_SAPAromaticElement,SAPAromaticElement);
CL_VALUE_ENUM(kw::_sym_SAPAliphatic,SAPAliphatic);
CL_VALUE_ENUM(kw::_sym_SAPAromatic,SAPAromatic);
CL_VALUE_ENUM(kw::_sym_SAPAtomMap,SAPAtomMap);
CL_END_ENUM(_sym_STARAtomTestEnumConverterSTAR);





}; // namespace chem
