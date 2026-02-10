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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
//#define DEBUG_LEVEL_FULL
//
//
//     	chemInfo.cc
//
//

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/vf2_sub_graph_iso.hpp"
#include "boost/graph/mcgregor_common_subgraphs.hpp"

#include <clasp/core/common.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/designators.h>
#include <cando/chem/loop.h>
#include <cando/chem/chemInfo.h>
#include <clasp/core/hashTableEqual.h>
#include <clasp/core/hashTableEql.h>
//#include "core/archiveNode.h"
#include <clasp/core/evaluator.h>
#include <clasp/core/lispStream.h>
//#include "boundFrame.h"
#include <cando/adapt/symbolMap.h>
#include <cando/adapt/symbolSet.h>
#include <cando/chem/bond.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>

SYMBOL_EXPORT_SC_(ChemPkg, parse_smarts);
SYMBOL_EXPORT_SC_(ChemPkg, parse_smirks);
SYMBOL_EXPORT_SC_(ChemPkg, compile_antechamber_type_rule);

#if 0
#define CI_LOG(...) { __VA_ARGS__ }
#else
#define CI_LOG(...)
#endif

// extern chem::SmartsRoot_sp smarts_compile(const string &input, stringstream &errorStream);
extern chem::AntechamberRoot_sp gaff_compile(const string &input, chem::WildElementDict_sp dict, stringstream &errorStream);
extern string antechamberError();

namespace chem {

ChemInfoMatch_sp current_match() {
  ChemInfoMatch_sp match = gc::As<ChemInfoMatch_sp>(_sym_STARcurrent_matchSTAR->symbolValue());
  return match;
}

std::atomic<size_t> global_ChemInfoNodeId;
size_t next_ChemInfoNodeId() { return ++global_ChemInfoNodeId; }

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
  case SABUseBondMatcher:
    return "use-bond-matcher";
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
  case SABSameRingBond:
    return "@";
  }
  return "{{unknownSabType}}";
}

void walk_nodes(ChemInfoNode_sp node, std::function<void(ChemInfoNode_sp)> const &fn) {
  fn(node);
  core::List_sp childs = node->children();
  for (auto cur : childs) {
    core::T_sp tchild = CONS_CAR(cur);
    if (tchild.notnilp()) {
      ChemInfoNode_sp child = gc::As_unsafe<ChemInfoNode_sp>(tchild);
      walk_nodes(child, fn);
    }
  }
}

void walk_nodes_with_parent(core::T_sp parent, ChemInfoNode_sp node, std::function<bool(core::T_sp, ChemInfoNode_sp)> const &fn) {
  bool descend = fn(parent, node);
  if (descend) {
    core::List_sp childs = node->children();
    for (auto cur : childs) {
      core::T_sp tchild = CONS_CAR(cur);
      if (tchild.notnilp()) {
        ChemInfoNode_sp child = gc::As_unsafe<ChemInfoNode_sp>(tchild);
        walk_nodes_with_parent(node, gc::As<ChemInfoNode_sp>(child), fn);
      }
    }
  }
}

size_t calculate_max_tags(ChemInfoNode_sp node) {
  size_t maxtag = 0;
  walk_nodes(node, [&maxtag](ChemInfoNode_sp node) {
    if (gc::IsA<AtomTest_sp>(node)) {
      AtomTest_sp atnode = gc::As_unsafe<AtomTest_sp>(node);
      if (atnode->atomTestType() == SAPAtomMap) {
        int intarg = atnode->getIntArg();
        if (intarg > maxtag) {
          maxtag = intarg;
        }
      }
    }
  });
  return maxtag;
}

DOCGROUP(cando);
CL_DEFUN void chem__walk_chem_info_with_parent(ChemInfoNode_sp top, core::T_sp callback) {
  walk_nodes_with_parent(nil<core::T_O>(), top, [&callback](core::T_sp parentOrNil, ChemInfoNode_sp node) {
    core::eval::funcall(callback, parentOrNil, node);
    return true;
  });
};

ChemInfoMatch_sp ChemInfoMatch_O::make(Root_sp root, size_t maxTag, core::HashTable_sp ringLookup) {
  auto match = gctools::GC<ChemInfoMatch_O>::allocate(root, maxTag + 1);
  match->_TagLookup = core::SimpleVector_O::make(maxTag + 1);
  match->_RingLookup = ringLookup;
  return match;
}

void ChemInfoMatch_O::initialize() {
  this->Base::initialize();
  this->_Matches = false;
  this->_TagLookup = core::SimpleVector_O::make(this->_MaxTagPlus1);
  //  this->_ClosestMatch = core::HashTable_O::createEqual();
}

void ChemInfoMatch_O::fields(core::Record_sp node) {
  // this->Base::fields(node); // T_O
  node->field(INTERN_(kw, matches), this->_Matches);
  node->field(INTERN_(kw, tags), this->_TagLookup);
  node->field(INTERN_(kw, tag_history), this->_TagHistory);
  //  node->field(INTERN_(kw, closestMatch), this->_ClosestMatch);
}

CL_LISPIFY_NAME("ChemInfoMatch-matches");
CL_DEFMETHOD bool ChemInfoMatch_O::matches() { return this->_Matches; }

string ChemInfoMatch_O::__repr__() const {
  
  stringstream ss;
  ss << "( " << this->className();
  ss << " :TagLookup '(";
  for (size_t ii = 0; ii < this->_TagLookup->length(); ++ii) {
    core::T_sp atom = this->_TagLookup->rowMajorAref(ii);
    if (atom.notnilp()) {
      ss << fmt::format(":tag {} :value {} ) ", ii, _rep_(atom));
    };
  }
  ss << " ))";
  return ss.str();
}

void ChemInfoMatch_O::clearAtomTags() {
  for (size_t ii = 0; ii < this->_TagLookup->length(); ++ii) {
    this->_TagLookup->rowMajorAset(ii, nil<core::T_O>());
  }
}

void ChemInfoMatch_O::saveTagLookup() {
  // Only add new ones
  for (auto cur : this->_TagHistory) {
    bool same = true;
    core::SimpleVector_sp other = gc::As_unsafe<core::SimpleVector_sp>(CONS_CAR(cur));
    for (size_t ii = 0; ii < this->_TagLookup->length(); ++ii) {
      if ((*this->_TagLookup)[ii] != (*other)[ii]) {
        same = false;
        break;
      }
    }
    if (same)
      return;
  }
  core::SimpleVector_sp copy =
      core::SimpleVector_O::make(this->_TagLookup->length(), nil<core::T_O>(), false, this->_TagLookup->length(),
                                 (core::T_sp *)(this->_TagLookup->rowMajorAddressOfElement_(0)));
  this->_TagHistory = core::Cons_O::create(copy, this->_TagHistory);
}

bool ChemInfoMatch_O::recognizesAtomTag(core::T_sp tag) {
  if (tag.fixnump() && tag.unsafe_fixnum() >= 0 && tag.unsafe_fixnum() < this->_TagLookup->length()) {
    return this->_TagLookup->rowMajorAref(tag.unsafe_fixnum()).notnilp();
  }
  SIMPLE_ERROR("Illegal tag {}", _rep_(tag));
}

void ChemInfoMatch_O::defineAtomTag(Atom_sp a, core::T_sp tag) {
  if (tag.fixnump() && tag.unsafe_fixnum() >= 0 && tag.unsafe_fixnum() < this->_TagLookup->length()) {
    this->_TagLookup->rowMajorAset(tag.unsafe_fixnum(), a);
    return;
  }
  SIMPLE_ERROR("The tag {} is an illegal index into the tag vector {}", _rep_(tag), _rep_(this->_TagLookup));
}

bool ChemInfoMatch_O::hasAtomWithTag(core::T_sp tag) {
  if (tag.fixnump() && tag.unsafe_fixnum() >= 0 && tag.unsafe_fixnum() < this->_TagLookup->length()) {
    return this->_TagLookup->rowMajorAref(tag.unsafe_fixnum()).notnilp();
  }
  SIMPLE_ERROR("The tag {} is an illegal index into the tag vector {}", _rep_(tag), _rep_(this->_TagLookup));
}

CL_LISPIFY_NAME("getAtomWithTag");
CL_DEFMETHOD Atom_sp ChemInfoMatch_O::getAtomWithTag(core::T_sp tag) {
  if (oCdr(this->_TagHistory).notnilp()) {
    SIMPLE_WARN("The chem-info-match has multiple solutions - get-atom-with-tag for tag {} works when there is one - solutions {} "
                "- smarts: {}",
                _rep_(tag), _rep_(this->_TagHistory), this->_Root->originalCode());
  }
  if (tag.fixnump() && tag.unsafe_fixnum() >= 0 && tag.unsafe_fixnum() < this->_TagLookup->length()) {
    return gc::As<Atom_sp>(this->_TagLookup->rowMajorAref(tag.unsafe_fixnum()));
  }
  SIMPLE_ERROR("The tag {} is an illegal index into the tag vector {}", _rep_(tag), _rep_(this->_TagLookup));
}

CL_LISPIFY_NAME("getAtomWithTagOrNil");
CL_DEFMETHOD core::T_sp ChemInfoMatch_O::getAtomWithTagOrNil(core::T_sp tag) {
  if (!this->hasAtomWithTag(tag))
    return nil<core::T_O>();
  return this->getAtomWithTag(tag);
}

void ChemInfoMatch_O::setRingTag(Atom_sp atom, core::T_sp index) { this->_RingLookup->setf_gethash(index, atom); }

bool ChemInfoMatch_O::matchesRingTag(Atom_sp atom, core::T_sp tag) {
  Atom_sp other = gc::As<Atom_sp>(this->_RingLookup->gethash(tag, nil<T_O>()));
  return atom->isBondedTo(other);
  //  return core::cl__eq(other,atom);
}

// ------- WildElementDict_O

void WildElementDict_O::initialize() {
  this->Base::initialize();
  this->_AtomWildCards = core::HashTable_O::createEqual();
}

CL_LISPIFY_NAME("addWildName");
CL_DEFMETHOD void WildElementDict_O::addWildName(core::Symbol_sp name) {
  adapt::SymbolSet_sp m = adapt::SymbolSet_O::create();
  this->_AtomWildCards->setf_gethash(name, m);
}

CL_LISPIFY_NAME("addWildNameMap");
CL_DEFMETHOD void WildElementDict_O::addWildNameMap(core::Symbol_sp wildName, core::Symbol_sp elementName) {
  
  if (!this->_AtomWildCards->contains(wildName)) {
    SIMPLE_ERROR("Could not find wild-card {}", _rep_(core::Cons_O::createList(wildName)));
  }
  core::T_sp tss = this->_AtomWildCards->gethash(wildName).as<adapt::SymbolSet_O>();
  if (tss.nilp()) {
    SIMPLE_ERROR("Could not find AtomWildCards for {}", _rep_(wildName));
  }
  adapt::SymbolSet_sp ss = tss.as<adapt::SymbolSet_O>();
  ss->insert(elementName);
  this->_AtomWildCards->gethash(wildName).as<adapt::SymbolSet_O>()->insert(elementName);
}

bool WildElementDict_O::recognizesWildName(core::Symbol_sp name) {
  return _AtomWildCards->gethash(name).notnilp();
}

bool WildElementDict_O::recognizesWildNameElement(core::Symbol_sp name, core::Symbol_sp element) {
  core::T_sp symset = this->_AtomWildCards->gethash(name);
  if (symset.nilp()) {
    SIMPLE_ERROR("Unrecognized wild card name({})", _rep_(name));
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

void WildElementDict_O::fields(core::Record_sp node) { node->field(INTERN_(kw, data), this->_AtomWildCards); }

void ResidueList_O::fields(core::Record_sp node) {
  // do nothing
}

struct ChemInfoTypeToName {
  ChemInfoType type;
  string name;
};

ChemInfoTypeToName chemInfoTypesToName[] = {{root, "root"},
                                            {smartsRoot, "smartsRoot"},
                                            {antechamberRoot, "antechamberRoot"},
                                            {chain, "chain"},
                                            {branch, "branch"},
                                            {logical, "logical"},
                                            {ringTest, "ringTest"},
                                            {atomTest, "atomTest"},
                                            {antechamberBondToAtomTest, "antechamberBondToAtomTest"},
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
  THROW_HARD_ERROR("%s", ss.str());
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
  THROW_HARD_ERROR("%s", ss.str());
}

void ChemInfoNode_O::fields(core::Record_sp node) {
  // Nothing to do here
  node->field_if_not_default(INTERN_(kw, id), this->_Id, (size_t)0);
  node->field_if_not_nil(INTERN_(kw, bounds), this->_Bounds);
}

string ChemInfoNode_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->className();
#ifdef DEBUG_ID
  ss << " :id " << this->_Id;
#endif
  ss << ">";
  return ss.str();
}

uint ChemInfoNode_O::depth() const {
  
  LOG("Node type = {}", this->className());
  LOG("Returning 1");
  return 1;
};

void BondListMatchNode_O::fields(core::Record_sp node) {
  this->Base::fields(node);
  // Nothing to do here
}

string BondListMatchNode_O::asSmarts() const {
  
  SUBCLASS_MUST_IMPLEMENT();
}

void AtomOrBondMatchNode_O::fields(core::Record_sp node) {
  this->Base::fields(node);
  node->field_if_not_default(INTERN_(kw, ring_test), this->_RingTest, SARNone);
  node->field_if_not_default(INTERN_(kw, ring_id), this->_RingId, -1);
}

bool AtomOrBondMatchNode_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  bool matches_from = this->matches_Atom(root, from);
  bool matches_other = this->matches_Atom(root, bond->getOtherAtom(from));
  return matches_from && matches_other;
};

bool AtomOrBondMatchNode_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  if (this->_RingTest == SARNone) {
    CI_LOG(if (chem__verbose(1)) core::clasp_write_string(fmt::format("this->_RingTest == SARNone atom {}\n", _rep_(atom))););
    return true;
  } else if (this->_RingTest == SARRingSet) {
    CI_LOG(if (chem__verbose(1)) 
      core::clasp_write_string(fmt::format("this->_RingTest == SARRingSet atom = {} this->_RingId = {}\n", _rep_(atom), this->_RingId));
           );
    LOG("Setting atom {} matches _RingId {}", _rep_(atom), this->_RingId);
    current_match()->setRingTag(atom, core::make_fixnum(this->_RingId));
    return true;
  }
  // It's SARRingTest
  LOG("Checking if atom {} matches _RingId {}", _rep_(atom), this->_RingId);
  CI_LOG(if (chem__verbose(1)) {
      core::clasp_write_string(fmt::format("It's SARRingTest atom = {} this->_RingId = {}\n", _rep_(atom), this->_RingId));
    }
    );
  return current_match()->matchesRingTag(atom, core::make_fixnum(this->_RingId));
}

CL_DEFMETHOD void AtomOrBondMatchNode_O::setRingTest(RingTestEnum test) { this->_RingTest = test; }
CL_DEFMETHOD RingTestEnum AtomOrBondMatchNode_O::getRingTest() const { return this->_RingTest; }

CL_DEFMETHOD void AtomOrBondMatchNode_O::setRingId(int id) { this->_RingId = id; }

CL_DEFMETHOD int AtomOrBondMatchNode_O::getRingId() const { return this->_RingId; }

string AtomOrBondMatchNode_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->className() << " ";
#ifdef DEBUG_ID
  ss << " :id " << this->_Id;
#endif
  if (this->_RingTest != SARNone) {
    ss << " :ring-id " << this->_RingId;
  }
  ss << ">";
  return ss.str();
}

void BondMatchNode_O::fields(core::Record_sp node) {
  this->Base::fields(node);
  // Nothing to do here
}

void Logical_O::initialize() {
  this->Base::initialize();
  this->_Operator = logAlwaysTrue;
  this->_Left = nil<core::T_O>();
  this->_Right = nil<core::T_O>();
}

core::T_sp Logical_O::children() {
  ql::list result;
  if (this->_Left.notnilp())
    result << this->_Left;
  if (this->_Right.notnilp())
    result << this->_Right;
  return result.cons();
}

CL_DEFMETHOD LogicalOperatorType Logical_O::logical_operator() const { return this->_Operator; }
string Logical_O::asSmarts() const {
  stringstream ss;
  ss << "[";
  switch (this->_Operator) {
  case logAlwaysTrue:
    ss << "";
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
  ss << "]";
  return ss.str();
}

uint Logical_O::depth() const {
  
  uint d = MAX(this->_Left->depth(), this->_Right->depth());
  LOG("Returning Logical depth={}", d);
  return d;
}

bool Logical_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  switch (this->_Operator) {
  case logAlwaysTrue:
    LOG("Always return true");
    goto SUCCESS;
    break;
  case logIdentity:
    LOG("Identity no-op test");
    ASSERT(!this->_Left.nilp());
    if (this->_Left->matches_Atom(root, atom))
      goto SUCCESS;
    break;
  case logNot:
    LOG("logNot");
    ASSERT(!this->_Left.nilp());
    if (!(this->_Left->matches_Atom(root, atom)))
      goto SUCCESS;
    break;
  case logHighPrecedenceAnd:
    LOG("logHighPrecedenceAnd");
    ASSERT(!this->_Left.nilp());
    if (this->_Left->matches_Atom(root, atom) && this->_Right->matches_Atom(root, atom))
      goto SUCCESS;
    break;
  case logLowPrecedenceAnd:
    LOG("logLowPrecedenceAnd");
    if (this->_Left->matches_Atom(root, atom) && this->_Right->matches_Atom(root, atom))
      goto SUCCESS;
    break;
  case logOr: {
    LOG("logOr");
    ASSERT(!this->_Left.nilp());
    bool leftMatch = this->_Left->matches_Atom(root, atom);
    if (leftMatch)
      goto SUCCESS;
    LOG("Left match failed");
    ASSERT(!this->_Right.nilp());
    bool rightMatch = this->_Right->matches_Atom(root, atom);
    if (rightMatch)
      goto SUCCESS;
    LOG("Right match failed");
  } break;
  default:
    stringstream err;
    err << "Unknown logical operator(" << this->_Operator << ")";
    SIMPLE_ERROR("{}", err.str());
  }
  // FAIL:
  LOG("FAIL");
  return false;
SUCCESS:
  LOG("SUCCESS!!!");
  return this->Base::matches_Atom(root, atom);
}

bool Logical_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  switch (this->_Operator) {
  case logAlwaysTrue:
    goto SUCCESS;
    break;
  case logIdentity:
    ASSERT(this->_Left.notnilp());
    if (this->_Left->matches_Bond(root, from, bond))
      goto SUCCESS;
    break;
  case logNot:
    ASSERT(this->_Left.notnilp());
    if (!(this->_Left->matches_Bond(root, from, bond)))
      goto SUCCESS;
    break;
  case logHighPrecedenceAnd:
  case logLowPrecedenceAnd:
    ASSERT(this->_Left.notnilp());
    if (this->_Left->matches_Bond(root, from, bond) && this->_Right->matches_Bond(root, from, bond))
      goto SUCCESS;
    break;
  case logOr:
    ASSERT(this->_Left.notnilp());
    if (this->_Left->matches_Bond(root, from, bond) || this->_Right->matches_Bond(root, from, bond))
      goto SUCCESS;
    break;
  default:
    stringstream err;
    err << "Unknown logical operator(" << this->_Operator << ")";
    SIMPLE_ERROR("{}", err.str());
  }
  // FAIL:
  LOG("FAIL");
  return false;
SUCCESS:
  LOG("SUCCESS!!!");
  return true; // this->Base::matches_Bond(root,from,bond);
}

core::NullTerminatedEnumAssociation logicalEnum[] = {{"logAlwaysTrue", logAlwaysTrue},
                                                     {"logIdentity", logIdentity},
                                                     {"logNot", logNot},
                                                     {"logHighPrecedenceAnd", logHighPrecedenceAnd},
                                                     {"logOr", logOr},
                                                     {"logLowPrecedenceAnd", logLowPrecedenceAnd},
                                                     {"", -1}};

void Logical_O::fields(core::Record_sp node) {
  node->field_if_not_nil(INTERN_(kw, right), this->_Right);
  node->field_if_not_nil(INTERN_(kw, left), this->_Left);
  node->/*pod_*/ field(INTERN_(kw, op), this->_Operator);
  this->Base::fields(node);
}

string Logical_O::__repr__() const {
  stringstream ss;
  ss << "#<LOGICAL";
#ifdef DEBUG_ID
  ss << " :id " << this->_Id;
#endif
  ss << " :op ";
  for (int i = 0; logicalEnum[i]._Enum != -1; ++i) {
    if (logicalEnum[i]._Enum == this->_Operator) {
      ss << logicalEnum[i]._Key;
    }
  }
  if (this->_RingTest != SARNone) {
    ss << " :ring-id " << this->_RingId;
  }
  ss << "@" << (void*)this->asSmartPtr().raw_();
  ss << ">";
  return ss.str();
}

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logIdentity(core::T_sp nilOrOp) { return create(logIdentity, nilOrOp); };

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logOr(core::T_sp nilOrOp1, core::T_sp nilOrOp2) {
  return create(logOr, nilOrOp1, nilOrOp2);
};

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logNot(core::T_sp nilOrOp) { return create(logNot, nilOrOp); };

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logLowPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2) {
  return create(logLowPrecedenceAnd, nilOrOp1, nilOrOp2);
};

CL_DEF_CLASS_METHOD Logical_sp Logical_O::create_logHighPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2) {
  return create(logHighPrecedenceAnd, nilOrOp1, nilOrOp2);
};

// --- TagSet  set a tag for the atom

core::NullTerminatedEnumAssociation bondEnum[] = {{"SABNoBond", SABNoBond},
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
  this->_AtomTest = nil<core::T_O>();
  this->_RingTag = nil<core::Symbol_O>();
}

core::T_sp ResidueTest_O::children() {
  ql::list result;
  if (this->_AtomTest.notnilp())
    result << this->_AtomTest;
  return result.cons();
}

string ResidueTest_O::asSmarts() const {
  
  stringstream ss;
  ss << sabToString(this->_Bond) << this->_AtomTest->asSmarts() << this->_RingTag;
  return ss.str();
}

bool ResidueTest_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  Atom_sp ringAtom;
  SmartsRoot_sp smartsRoot;
  LOG("ResidueTest match for atom: {}", atom->description().c_str());
  if (root->type() != chem::smartsRoot) {
    SIMPLE_ERROR("Trying to carry out a ResidueTest with a non Smarts root");
  }
  SIMPLE_ERROR("Handle atom tests in ResidueTest");
#if 0
  // First check if the AtomTest matches
  ASSERTNOTNULL(this->_AtomTest);
  if (!this->_AtomTest->matches_Atom(root, atom)) {
    LOG("The _AtomTest failed");
    goto FAIL;
  }
  // Now handle the ring test
  smartsRoot = (root).as<SmartsRoot_O>();
  // It does recognize the ring ID so check if the current atom is bonded to it
  // and then have the root forget the atom
  if (!smartsRoot->getMatch()->recognizesAtomTag(this->_RingTag)) {
    SIMPLE_ERROR("We are trying to test the atomTag ({}) but it doesn't exist!" , _rep_(this->_RingTag));
  }
  LOG("      ResidueTest check close");
  ringAtom = smartsRoot->getMatch()->getAtomWithTag(this->_RingTag);
  LOG("      checking if {} matches ringAtom: {}" , atom->description().c_str() , ringAtom->description().c_str());
  if (!atom->containedByValid()) {
    SIMPLE_ERROR("{} failed containedByValid()" , atom->description());
  }
  if (!ringAtom->containedByValid()) {
    SIMPLE_ERROR("{} failed containedByValid()" , ringAtom->description());
  }
  LOG("   atom->containedBy()->description() = {}" , atom->containedBy()->description().c_str());
  LOG("   ringAtom->containedBy()->description() = {}" , ringAtom->containedBy()->description().c_str());
  if (atom->containedBy() == ringAtom->containedBy()) {
    LOG("RESIDUE MATCH THEY DO!!!!");
    goto SUCCESS;
  }
 FAIL:
  LOG("FAIL");
  return false;
 SUCCESS:
  LOG("SUCCESS!");
  return true;
#endif
}

void ResidueTest_O::fields(core::Record_sp node) {
  node->/*pod_*/ field(INTERN_(kw, bond), this->_Bond);
  node->field(INTERN_(kw, ringTag), this->_RingTag);
  node->field(INTERN_(kw, atomTest), this->_AtomTest);
  this->Base::fields(node);
}

// ------ BondToAtomTest

bool unsafe_is_aromatic(Atom_sp a1) {
  core::T_sp taromaticity_info = chem::_sym_STARcurrent_aromaticity_informationSTAR->symbolValue();
  if (taromaticity_info.unboundp()) {
    SIMPLE_ERROR("In unsafe_is_aromatic and the chem:*current-aromaticity-information* is UNBOUND and it should not be");
  }
  if (taromaticity_info.nilp()) {
    SIMPLE_ERROR("In unsafe_is_aromatic and the chem:*current-aromaticity-information* is NIL and it should not be");
  }
  core::HashTable_sp aromaticity_info = gc::As<core::HashTable_sp>(taromaticity_info);
  core::T_sp info = aromaticity_info->gethash(a1);
  CI_LOG(if (chem__verbose(1)) {
    core::clasp_write_string(fmt::format("is_aromatic {} -> {}\n", _rep_(a1), _rep_(info)));
    });
  return info.notnilp();
}

bool is_aromatic(Atom_sp a1) {
  if (chem::_sym_STARcurrent_aromaticity_informationSTAR->boundP())
    return unsafe_is_aromatic(a1);
  return false;
}

bool aromatic_information_available_p() {
  return (chem::_sym_STARcurrent_ringsSTAR->boundP() && chem::_sym_STARcurrent_aromaticity_informationSTAR->boundP());
}

bool _matchInAromaticBond(Atom_sp a1, Atom_sp a2) {
  if (a1.boundp() && a2.boundp() && chem::_sym_STARcurrent_ringsSTAR->boundP() &&
      chem::_sym_STARcurrent_aromaticity_informationSTAR->boundP()) {
    if (!(unsafe_is_aromatic(a1) && unsafe_is_aromatic(a2))) return false;
    core::List_sp rings = chem::_sym_STARcurrent_ringsSTAR->symbolValue();
    // If a1 and other are in the same ring
    CI_LOG(if (chem__verbose(1)) {
      core::clasp_write_string(fmt::format("_matchInAromaticBond a1={} a2={}\n", _rep_(a1), _rep_(a2)));
    }
      );
    bool found_a1 = false;
    bool found_a2 = false;
    for (auto cur : rings) {
      core::List_sp ring = CONS_CAR(cur);
      CI_LOG(if (chem__verbose(1)) {
        core::clasp_write_string(fmt::format("Looking at ring: {}\n", _rep_(ring)));
        });
      found_a1 = false;
      found_a2 = false;
      for (auto cur2 : ring) {
        core::T_sp one = CONS_CAR(cur2);
        if (one == a1) found_a1 = true;
        else if (one == a2) found_a2 = true;
      }
      CI_LOG(if (chem__verbose(1)) {
        core::clasp_write_string(fmt::format("found_a1 = {} found_a2 = {}\n", found_a1, found_a2));
        });
      if (found_a1 && found_a2)
        return true;
    }
    return false;
  }
  return false;
}

DOCGROUP(cando);
CL_DEFUN bool chem__in_aromatic_bond(Atom_sp a1, Atom_sp a2) {
  core::clasp_write_string(fmt::format("chem__in_aromatic_bond a1={} a2={}\n", _rep_(a1), _rep_(a2)));
  return _matchInAromaticBond(a1, a2);
}

bool _matchBondTypes(BondEnum be, chem::BondOrder bo, Atom_sp a1, Atom_sp a2) {
  LOG("bondOrder = {}", bondOrderToString(bo).c_str());
  switch (be) {
  case SABSingleBond:
    LOG("SMARTS BondEnum = SABSingleBond");
    if (!Bond_O::singleBondP(bo))
      goto nomatch;
    if (_matchInAromaticBond(a1, a2))
      goto nomatch;
    break;
  case SABSingleOrAromaticBond:
    LOG("SMARTS BondEnum = SABSingleOrAromaticBond");
    if (!(Bond_O::singleBondP(bo) || bo == chem::aromaticBond || _matchInAromaticBond(a1, a2)))
      goto nomatch;
    break;
  case SABDoubleOrAromaticBond:
    LOG("SMARTS BondEnum = SABDoubleOrAromaticBond");
    if (!(bo == chem::doubleBond || bo == chem::aromaticBond || _matchInAromaticBond(a1, a2)))
      goto nomatch;
    break;
  case SABTripleOrAromaticBond:
    LOG("SMARTS BondEnum = SABTripleOrAromaticBond");
    if (!(bo == chem::tripleBond || bo == chem::aromaticBond || _matchInAromaticBond(a1, a2)))
      goto nomatch;
    break;
  case SABDoubleBond:
    LOG("SMARTS BondEnum = SABDoubleBond");
    if (bo != chem::doubleBond)
      goto nomatch;
    if (_matchInAromaticBond(a1, a2))
      goto nomatch;
    break;
  case SABTripleBond:
    LOG("SMARTS BondEnum = SABTriple");
    if (bo != chem::tripleBond)
      goto nomatch;
    if (_matchInAromaticBond(a1, a2))
      goto nomatch;
    break;
  case SABAromaticBond:
    LOG("SMARTS BondEnum = SABAromaticBond");
    if (!(bo == chem::aromaticBond || _matchInAromaticBond(a1, a2)))
      goto nomatch;
    break;
  case SABAnyBond:
    break;
  case SABDelocalizedBond:
    goto nomatch;
    break;
  case SABDirectionalSingleDownOrUnspecified:
    SIMPLE_ERROR("Must implement SingleDirectionalDownOrUnspecified");
  case SABDirectionalSingleUpOrUnspecified:
    SIMPLE_ERROR("Must implement SingleDirectionalUpOrUnspecified");
  case SABDirectionalSingleUp:
  case SABDirectionalSingleDown:
    SIMPLE_ERROR("Must implement directional bonds");
  default:
    printf("%s:%d:%s Fell through to default with be = %s - maybe the BondMatcher should handle it\n", __FILE__, __LINE__,
           __FUNCTION__, sabToString(be).c_str());
    goto nomatch;
  }
  LOG("THEY MATCH!!");
  return true;
nomatch:
  LOG("THEY DONT MATCH");
  return false;
}

SYMBOL_EXPORT_SC_(ChemPkg, STARcurrent_ringsSTAR);
SYMBOL_EXPORT_SC_(ChemPkg, make_rings);
SYMBOL_EXPORT_SC_(ChemPkg, rings_rings);

bool _matchBondTypesWithAtoms(BondEnum be, chem::BondOrder bo, Atom_sp from, Bond_sp bond) {
  if (be == SABSameRingBond) {
    // Expect a list of lists
    Atom_sp other = bond->getOtherAtom(from);
    if (chem::_sym_STARcurrent_ringsSTAR->boundP()) {
      core::List_sp rings = chem::_sym_STARcurrent_ringsSTAR->symbolValue();
      // If from and other are in the same ring
      bool found_from = false;
      bool found_other = false;
      for (auto cur : rings) {
        core::List_sp ring = CONS_CAR(cur);
        found_from = false;
        found_other = false;
        for (auto cur2 : ring) {
          core::T_sp one = CONS_CAR(cur2);
          if (one == from)
            found_from = true;
          else if (one == other)
            found_other = true;
        }
        if (found_from && found_other)
          return true;
      }
      return false;
    }
    SIMPLE_ERROR("ring-bond detection can only be done if the chem:*current-rings* dynamic variable is bound - it is not");
  } else {
    return _matchBondTypes(be, bo, from, bond->getOtherAtom(from));
  }
}

// --------------------------------------------------
//
// BondLogical tests
//

void BondLogical_O::fields(core::Record_sp node) {
  node->field_if_not_unbound(INTERN_(kw, right), this->_Right);
  node->field_if_not_unbound(INTERN_(kw, left), this->_Left);
  node->/*pod_*/ field(INTERN_(kw, op), this->_Operator);
  this->Base::fields(node);
}

bool BondLogical_O::matches_Bond(Root_sp root, Atom_sp from, Bond_sp bond) {
  switch (this->_Operator) {
  case logAlwaysTrue:
    LOG("Always return true");
    goto SUCCESS;
    break;
  case logIdentity:
    LOG("Identity no-op test");
    ASSERT(this->_Left.boundp());
    if (this->_Left->matches_Bond(root, from, bond))
      goto SUCCESS;
    break;
  case logNot:
    LOG("logNot");
    ASSERT(this->_Left.boundp());
    if (!this->_Left->matches_Bond(root, from, bond))
      goto SUCCESS;
    break;
  case logHighPrecedenceAnd:
    LOG("logHighPrecedenceAnd");
    ASSERT(this->_Left.boundp());
    ASSERT(this->_Right.boundp());
    if (this->_Left->matches_Bond(root, from, bond) && this->_Right->matches_Bond(root, from, bond))
      goto SUCCESS;
    break;
  case logLowPrecedenceAnd:
    LOG("logLowPrecedenceAnd");
    ASSERT(this->_Left.boundp());
    ASSERT(this->_Right.boundp());
    if (this->_Left->matches_Bond(root, from, bond) && this->_Right->matches_Bond(root, from, bond))
      goto SUCCESS;
    break;
  case logOr: {
    LOG("logOr");
    ASSERT(this->_Left.boundp());
    ASSERT(this->_Right.boundp());
    if (this->_Left->matches_Bond(root, from, bond))
      goto SUCCESS;
    if (this->_Right->matches_Bond(root, from, bond))
      goto SUCCESS;
  } break;
  default:
    stringstream err;
    err << "Unknown logical operator(" << this->_Operator << ")";
    SIMPLE_ERROR("{}", err.str());
  }
  // FAIL:
  LOG("FAIL");
  return false;
SUCCESS:
  LOG("SUCCESS!!!");
  return true;
}

core::T_sp BondLogical_O::children() {
  ql::list result;
  if (this->_Left.boundp())
    result << this->_Left;
  if (this->_Right.boundp())
    result << this->_Right;
  return result.cons();
}

CL_DEFMETHOD core::T_sp BondLogical_O::getLeft() const {
  if (this->_Left.boundp()) {
    return this->_Left;
  }
  return nil<core::T_O>();
}

CL_DEFMETHOD void BondLogical_O::setLeft(core::T_sp b) {
  if (b.nilp()) {
    this->_Left = unbound<BondMatcher_O>();
  } else
    this->_Left = gc::As<BondMatcher_sp>(b);
}

CL_DEFMETHOD core::T_sp BondLogical_O::getRight() const {
  if (this->_Right.boundp()) {
    return this->_Right;
  }
  return nil<core::T_O>();
}

CL_DEFMETHOD void BondLogical_O::setRight(core::T_sp b) {
  if (b.nilp()) {
    this->_Right = unbound<BondMatcher_O>();
  } else
    this->_Right = gc::As<BondMatcher_sp>(b);
}

CL_DEFMETHOD LogicalOperatorType BondLogical_O::bondLogicalOperator() const { return this->_Operator; }

CL_DEF_CLASS_METHOD BondLogical_sp BondLogical_O::create_bondLogIdentity(core::T_sp nilOrOp1) {
  BondMatcher_sp b1(unbound<BondMatcher_O>());
  if (nilOrOp1.notnilp())
    b1 = gc::As<BondMatcher_sp>(nilOrOp1);
  return create(logIdentity, b1);
};

CL_DEF_CLASS_METHOD BondLogical_sp BondLogical_O::create_bondLogNot(core::T_sp nilOrOp1) {
  BondMatcher_sp b1(unbound<BondMatcher_O>());
  if (nilOrOp1.notnilp())
    b1 = gc::As<BondMatcher_sp>(nilOrOp1);
  return create(logNot, b1);
};

CL_DEF_CLASS_METHOD BondLogical_sp BondLogical_O::create_bondLogOr(core::T_sp nilOrOp1, core::T_sp nilOrOp2) {
  BondMatcher_sp b1(unbound<BondMatcher_O>());
  BondMatcher_sp b2(unbound<BondMatcher_O>());
  if (nilOrOp1.notnilp())
    b1 = gc::As<BondMatcher_sp>(nilOrOp1);
  if (nilOrOp2.notnilp())
    b2 = gc::As<BondMatcher_sp>(nilOrOp2);
  return create(logOr, b1, b2);
};
CL_DEF_CLASS_METHOD BondLogical_sp BondLogical_O::create_bondLogLowPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2) {
  BondMatcher_sp b1(unbound<BondMatcher_O>());
  BondMatcher_sp b2(unbound<BondMatcher_O>());
  if (nilOrOp1.notnilp())
    b1 = gc::As<BondMatcher_sp>(nilOrOp1);
  if (nilOrOp2.notnilp())
    b2 = gc::As<BondMatcher_sp>(nilOrOp2);
  return create(logLowPrecedenceAnd, b1, b2);
};

CL_DEF_CLASS_METHOD BondLogical_sp BondLogical_O::create_bondLogHighPrecedenceAnd(core::T_sp nilOrOp1, core::T_sp nilOrOp2) {
  BondMatcher_sp b1(unbound<BondMatcher_O>());
  BondMatcher_sp b2(unbound<BondMatcher_O>());
  if (nilOrOp1.notnilp())
    b1 = gc::As<BondMatcher_sp>(nilOrOp1);
  if (nilOrOp2.notnilp())
    b2 = gc::As<BondMatcher_sp>(nilOrOp2);
  return create(logHighPrecedenceAnd, b1, b2);
};

CL_LISPIFY_NAME("make-bond-test");
CL_DEF_CLASS_METHOD BondTest_sp BondTest_O::make(BondEnum be) {
  _G();
  auto obj = gctools::GC<BondTest_O>::allocate(be); // RP_Create<Logical_O>(lisp);
  return obj;
};

CL_DEFMETHOD BondEnum BondTest_O::bondTestGetBond() const { return this->_Bond; }

void BondTest_O::fields(core::Record_sp node) {
  node->/*pod_*/ field(INTERN_(kw, bond), this->_Bond);
  this->Base::fields(node);
}

bool BondTest_O::matches_Bond(Root_sp root, Atom_sp from, Bond_sp bond) {
  chem::BondOrder bo = bond->getOrderFromAtom(from);
  return _matchBondTypesWithAtoms(this->_Bond, bo, from, bond);
}

core::T_sp BondTest_O::children() { return nil<core::T_O>(); }

void BondToAtomTest_O::initialize() {
  this->Base::initialize();
  this->_AtomTest = nil<core::T_O>();
}

CL_LISPIFY_NAME(make-bond-to-atom-test);
CL_DEF_CLASS_METHOD BondToAtomTest_sp BondToAtomTest_O::makeBondToAtomTest(BondEnum be, core::T_sp nilOrNode) {
  auto bta = gctools::GC<BondToAtomTest_O>::allocate(be);
  bta->_AtomTest = nilOrNode;
  bta->_BondMatcher = unbound<BondMatcher_O>();
  return bta;
};

CL_LISPIFY_NAME(make-bond-matcher-to-atom-test);
CL_DEF_CLASS_METHOD BondToAtomTest_sp BondToAtomTest_O::makeBondMatcherToAtomTest(BondMatcher_sp bm, core::T_sp nilOrNode) {
  auto bta = gctools::GC<BondToAtomTest_O>::allocate(SABUseBondMatcher);
  bta->_BondMatcher = bm;
  bta->_AtomTest = nilOrNode;
  return bta;
};

core::T_sp BondToAtomTest_O::children() {
  ql::list result;
  if (this->_Bond == SABUseBondMatcher) {
    if (this->_BondMatcher.boundp()) {
      result << this->_BondMatcher;
    } else {
      SIMPLE_ERROR("BondToAtomTest _Bond is SABUseBondMatcher but _BondMatcher is unbound");
    }
  }
  if (this->_AtomTest.notnilp())
    result << this->_AtomTest;
  return result.cons();
}

CL_DEFMETHOD BondEnum BondToAtomTest_O::bondType() { return this->_Bond; }

CL_DEFMETHOD bool BondToAtomTest_O::bondMatcherBoundP() const { return this->_BondMatcher.boundp(); }

CL_DEFMETHOD BondMatcher_sp BondToAtomTest_O::bondMatcher() const {
  if (this->_BondMatcher.boundp()) {
    return this->_BondMatcher;
  }
  SIMPLE_ERROR("bond-matcher is unbound");
}

//! Set the BondType if it's one that is recognzed by _matchBondTypes
CL_DEFMETHOD bool BondToAtomTest_O::setfBondTypeIfOptimizable(BondEnum be) {
  switch (be) {
    // Aromatic bonds are not optimizable because there
    // is a slow test that checks if the two atoms are in the same
    // ring and if they are both aromatic.  If both of these are
    // true - then the bond between them is aromatic
  case SABSingleBond:
  case SABDoubleBond:
  case SABTripleBond:
  case SABAnyBond:
    this->_Bond = be;
    return true;
  default:
    return false;
  }
}

string BondToAtomTest_O::asSmarts() const {
  stringstream ss;
  ss << sabToString(this->_Bond) << this->_AtomTest->asSmarts();
  return ss.str();
}

CL_DEFMETHOD void BondToAtomTest_O::setAtomTest(core::T_sp atomTest) { this->_AtomTest = atomTest; }

bool BondToAtomTest_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  if (this->_Bond != SABUseBondMatcher) {
    chem::BondOrder bo;
    bo = bond->getOrderFromAtom(from);
    if (!chem::_matchBondTypesWithAtoms(this->_Bond, bo, from, bond))
      goto FAIL;
    if (this->_AtomTest->matches_Atom(root, bond->getOtherAtom(from)))
      goto SUCCESS;
  } else {
    if (this->_BondMatcher.unboundp()) {
      SIMPLE_ERROR("The BondToAtomTest BondMatcher is unbound");
    }
    if (!this->_BondMatcher->matches_Bond(root, from, bond))
      goto FAIL;
    if (this->_AtomTest->matches_Atom(root, bond->getOtherAtom(from)))
      goto SUCCESS;
  }
FAIL:
  LOG("FAIL!");
  return false;
SUCCESS:
  LOG("SUCCESS!");
  return true;
}

void BondToAtomTest_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, atomTest), this->_AtomTest);
  node->field_if_not_unbound(INTERN_(kw, bondMatcher), this->_BondMatcher);
  node->/*pod_*/ field(INTERN_(kw, bond), this->_Bond);
  this->Base::fields(node);
}

// ------- AtomTest

void AtomTest_O::initialize() {
  
  this->_Test = SAPNone;
  this->_IntArg = 0;
  this->_NumArg = 0;
  this->_StringArg = "";
}

CL_DEFMETHOD int AtomTest_O::getIntArg() { return this->_IntArg; };

CL_DEFMETHOD core::Symbol_sp AtomTest_O::getSymbolArg() const { return this->_SymbolArg; };

bool AtomTest_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  switch (this->_Test) {
  case SAPBondedToPrevious:
    if (chem::_matchBondTypesWithAtoms((chem::BondEnum)this->_IntArg, bond->getOrderFromAtom(from), from, bond))
      goto SUCCESS;
    goto FAIL;
  case SAPNotBondedToPrevious:
    if (!chem::_matchBondTypesWithAtoms((chem::BondEnum)this->_IntArg, bond->getOrderFromAtom(from), from, bond))
      goto SUCCESS;
    goto FAIL;
  default:
    // Do nothing, fall through to the atom test below
    break;
  }
  // If its not a bond test then
  // test the To atom of the bond for the property
  if (this->matches_Atom(root, bond->getOtherAtom(from)))
    goto SUCCESS;
FAIL:
  LOG("FAIL");
  return false;
SUCCESS:
  LOG("SUCCESS!");
  return true; // this->Base::matches_Bond(root,from,bond);
}

DOCGROUP(cando);
CL_DEFUN bool am1BccX(chem::Atom_sp atom) {
  
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

DOCGROUP(cando);
CL_DEFUN bool am1BccY(chem::Atom_sp atom) {
  
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

SYMBOL_EXPORT_SC_(ChemPkg, STARcurrent_matchSTAR);

bool AtomTest_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  int cnt;
  Atom_sp ringStartAtom;
  int hc = 0;
  CI_LOG(if (chem__verbose(2)) {
    core::clasp_write_string(fmt::format("AtomTest_O::matches_Atom test: {} atom: {}\n", this->testName(this->_Test), _rep_(atom)));
    });
  switch (this->_Test) {
  case SAPWildCard:
    LOG("SAPWildCard"); //
    CI_LOG(if (chem__verbose(2)) core::clasp_write_string("AtomTest_O::SAPWildCard\n"););
    goto SUCCESS;
  case SAPAromaticElement:
    LOG("SAPAromaticElement({}) == expecting({})", _rep_(atom->getElementAsSymbol()), _rep_(this->_SymbolArg));
    if (this->_SymbolArg == atom->getElementAsSymbol() && is_aromatic(atom))
      goto SUCCESS;
    break;
  case SAPElement:
    LOG("SAPElement({}) == expecting({})", _rep_(atom->getElementAsSymbol()), _rep_(this->_SymbolArg));
    if (this->_SymbolArg == atom->getElementAsSymbol())
      goto SUCCESS;
    break;
  case SAPAliphatic:
    LOG("SAPAliphatic");
    if (!is_aromatic(atom))
      goto SUCCESS;
    break;
  case SAPAromatic:
    LOG("SAPAromatic");
    if (is_aromatic(atom))
      goto SUCCESS;
    break;
  case SAPPredicateName:
    LOG("SAPPredicateName<{}> testing {}", _rep_(this->_SymbolArg), atom->description());
    if (root->evaluateTest(this->_SymbolArg, atom))
      goto SUCCESS;
    break;
  case SAPAtomicNumber:
    LOG("SAPAtomicNumber({}) == expecting({})", atom->getAtomicNumber(), this->_IntArg);
    CI_LOG(if (chem__verbose(2)) core::clasp_write_string(fmt::format("AtomTest_O::SAPAtomicNumber atom {} _IntArg {}\n", _rep_(atom), this->_IntArg)););
    if (this->_IntArg == atom->getAtomicNumber())
      goto SUCCESS;
    break;
  case SAPAtomicMass:
    LOG("SAPAtomicMass({}) == expecting({})", atom->getIntegerAtomicMass(), this->_IntArg);
    if (this->_IntArg == atom->getIntegerAtomicMass())
      goto SUCCESS;
    break;
  case SAPTotalHCount:
    hc = atom->getBondedHydrogenCount();
    LOG("SAPTotalHCount({}) == expecting({})", hc, this->_IntArg);
    LOG("Checking total hydrogen count looking for({}) got({})", this->_IntArg, hc);
    if (this->_IntArg == hc) {
      LOG("Match"); //
      goto SUCCESS;
    }
    LOG("No match"); //
    break;
  case SAPResidueTest:
    SIMPLE_WARN("What do I do with SAPResidueText");
#if 0
      LOG("SAPResidueTest looking for tag: {}" , this->_StringArg);
      if (!root->getMatch()->recognizesAtomTag(chemkw_intern(this->_StringArg))) {
        SIMPLE_ERROR(("We are trying to test the atomTag (" + this->_StringArg + ") but it doesn't exist!"));
      }
      LOG("      SAPResidueTest check close"); //
      ringStartAtom = root->getMatch()->getAtomWithTag(chemkw_intern(this->_StringArg));
      LOG("      checking if residue for {} matches residue for ringStart atom: {}" , atom->description().c_str() , ringStartAtom->description().c_str()); //
      if (atom->containedBy() == ringStartAtom->containedBy()) {
        LOG("RESIDUE MATCH THEY DO!!!!");
        goto SUCCESS;
      }
#endif
    break;
  case SAPNegativeCharge:       // next
  case SAPPositiveCharge:       // next
  case SAPNegativeFormalCharge: // next
  case SAPPositiveFormalCharge:
    LOG("SAPxxxxCharge"); //
    if (this->_IntArg == atom->getIonization())
      goto SUCCESS;
    break;
  case SAPRingConnectivity:
  case SAPRingMembershipCount: {
    int count = 0;
    if (!_sym_STARcurrent_ringsSTAR->boundP()) {
      SIMPLE_WARN("The SAPRingMembershipCount or SAPRingConnectivity test was attempted but chem:*current-rings* is not bound - so "
                  "it cannot work properly");
    } else {
      CI_LOG(if (chem__verbose(2)) core::clasp_write_string(fmt::format("SAPRingMembershipCount atom = {}\n", _rep_(atom))););
      core::List_sp rings = _sym_STARcurrent_ringsSTAR->symbolValue();
      for (auto ring_cur : rings) {
        core::List_sp ring = CONS_CAR(ring_cur);
        for (auto ring_atom_cur : ring) {
          core::T_sp ring_atom = CONS_CAR(ring_atom_cur);
          if (ring_atom == atom)
            count++;
        }
      }
    }
    if (count == this->_IntArg) {
      goto SUCCESS;
    }
  } break;
  case SAPRingSize:
    if (!_sym_STARcurrent_ringsSTAR->boundP()) {
      SIMPLE_WARN("The SAPRingSize test was attempted but chem:*current-rings* is not bound - so it cannot work properly");
    } else {
      CI_LOG(if (chem__verbose(2)) core::clasp_write_string(fmt::format("SAPRingSize _IntArg = {}\n", this->_IntArg)););
      core::List_sp rings = _sym_STARcurrent_ringsSTAR->symbolValue();
      for (auto ring_cur : rings) {
        core::T_sp ring = CONS_CAR(ring_cur);
        if (ring.consp() && (core::cl__length(ring) == this->_IntArg) &&
            gc::As_unsafe<core::Cons_sp>(ring)->memberEq(atom).notnilp()) {
          goto SUCCESS;
        }
      }
    }
    break;
  case SAPValence:
    LOG("SAPRingValence");
    // total bond valence
    // each single bond counts 1
    // each double/aromatic bond counts 2
    // each triple bond counts 3
    if (this->_IntArg == atom->getValence())
      goto SUCCESS;
    break;
  case SAPAtomMap: {
    LOG("SAPAtomMap[{}]", this->_IntArg);
    // We have an atom and we have a tag in this->_IntArg
    ChemInfoMatch_sp match = gc::As<ChemInfoMatch_sp>(_sym_STARcurrent_matchSTAR->symbolValue());
    match->defineAtomTag(atom, core::make_fixnum(this->_IntArg));
    CI_LOG(if (chem__verbose(2)) core::clasp_write_string(fmt::format("SAPAtomMap atom = {} _IntArg = {}\n", _rep_(atom), this->_IntArg)););
    goto SUCCESS;
  } break;
  case SAPConnectivity: // No implicit H's so Connectivity == Degree
  case SAPDegree:
    LOG("SAPDegree testing if atom->numberOfBonds(){{}} == this->_IntArg{{}}", atom->numberOfBonds(), this->_IntArg);
    CI_LOG(if (chem__verbose(2)) core::clasp_write_string(fmt::format("SAPConnectivity/SAPDegree _IntArg = {}\n", this->_IntArg)););
    if (this->_IntArg == atom->numberOfBonds())
      goto SUCCESS;
    break;
  case SAPInBond:
    LOG("SAPInBond");
    cnt = 0;
    for (int i = 0; i < atom->numberOfBonds(); i++) {
      if (chem::_matchBondTypes((chem::BondEnum)this->_IntArg, atom->bondedOrder(i), atom, atom->bondedNeighbor(i)))
        cnt++;
    }
    if (cnt == this->_NumArg)
      goto SUCCESS;
    break;
  case SAPArLevel:
    LOG("SAPArLevel");
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
    LOG("SAPNoRing");
    if (!atom->isInRing())
      goto SUCCESS;
    break;
  case SAPElectronegativeElement: // ( O, N, F, Cl, Br )
    LOG("SAPElectronegativeElement");
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
    SIMPLE_ERROR("{}", ss.str());
  }
FAIL:
  LOG("FAIL");
  CI_LOG(if (chem__verbose(2)) core::clasp_write_string("----> FAIL\n"););
  return false;
SUCCESS:
  LOG("SUCCESS!");
  CI_LOG(if (chem__verbose(2)) core::clasp_write_string("----> SUCCESS!!!\n"););
  return this->Base::matches_Atom(root, atom);
}

string AtomTest_O::asSmarts() const {
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
  case SAPAtomMap:
    ss << ":" << this->_IntArg;
    break;
  case SAPConnectivity: // No implicit H's so Connectivity == Degree
  case SAPDegree:
    ss << "X" << this->_IntArg;
    break;
  case SAPRingConnectivity:
    ss << "x" << this->_IntArg;
    break;
  case SAPInBond:
    ss << "<SAPInBond/" << this->_NumArg << ">";
    break;
  case SAPArLevel:
    ss << "<SAPArLevel/" << this->_IntArg << ">";
    break;
  case SAPNoRing:
    ss << "<SAPNoRing>";
    break;
  case SAPElectronegativeElement: // ( O, N, F, Cl, Br )
    ss << "<SAPElectronegativeElement>";
    break;
  case SAPBondedToPrevious:
    ss << "<SAPBondedToPrevious/" << this->_IntArg << ">";
    break;
  default:
    ss << "<UnknownTest-" << this->testName(this->_Test) << "}";
    break;
  }
  return ss.str();
}

core::NullTerminatedEnumAssociation testEnum[] = {{"SAPNone", SAPNone},
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
  node->/*pod_*/ field_if_not_default(INTERN_(kw, num), this->_NumArg, 0);
  node->/*pod_*/ field_if_not_default(INTERN_(kw, str), this->_StringArg, std::string());
  node->field_if_not_nil(INTERN_(kw, sym), this->_SymbolArg);
  node->/*pod_*/ field_if_not_default(INTERN_(kw, int), this->_IntArg, 0);
  node->/*pod_*/ field(INTERN_(kw, test), this->_Test);
  this->Base::fields(node);
}

// ------- Chain

void Chain_O::initialize() {
  this->Base::initialize();
  this->_Head = nil<core::T_O>();
  this->_Tail = nil<core::T_O>();
}

string Chain_O::asSmarts() const {
  stringstream ss;
  if (this->_Head.notnilp())
    ss << this->_Head->asSmarts();
  if (this->_Tail.notnilp())
    ss << this->_Tail->asSmarts();
  return ss.str();
}

core::T_sp Chain_O::children() {
  ql::list result;
  result << this->_Head;
  result << this->_Tail;
  return result.cons();
}

uint Chain_O::depth() const {
  
  return (MAX(af_depth(this->_Head), af_depth(this->_Tail) + 1));
}

bool Chain_O::matches_Atom(Root_sp root, chem::Atom_sp from) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(from));
#if 1
  if (this->_Head->matches_Atom(root, from)) {
    LOG("SUCCESS!\n");
    chem::BondList_sp bonds = from->getBondList();
    return this->_Tail->matches_BondList(root, from, bonds);
  }
  return false;
#else
  if (gc::IsA<AtomTest_sp>(this->_Head)) {
    AtomTest_sp atHead = gc::As_unsafe<AtomTest_sp>(this->_Head);
    if (atHead->matches_Atom(root, from)) {
      LOG("SUCCESS!\n");
      chem::BondList_sp bonds = from->getBondList();
      return this->_Tail->matches_BondList(root, from, bonds);
    }
    LOG("FAIL!\n");
    return false;
  } else if (gc::IsA<Logical_sp>(this->_Head)) {
    Logical_sp lHead = gc::As_unsafe<Logical_sp>(this->_Head);
    if (lHead->matches_Atom(root, from)) {
      LOG("SUCCESS!\n");
      chem::BondList_sp bonds = from->getBondList();
      return this->_Tail->matches_BondList(root, from, bonds);
    }
    LOG("FAIL!\n");
    return false;
  } else if (gc::IsA<AntechamberFocusAtomMatch_sp>(this->_Head)) {
    AntechamberFocusAtomMatch_sp focus = gc::As_unsafe<AntechamberFocusAtomMatch_sp>(this->_Head);
    if (focus->matches_Atom(root, from)) {
      LOG("SUCCESS!\n");
      chem::BondList_sp bonds = from->getBondList();
      return this->_Tail->matches_BondList(root, from, bonds);
    }
    LOG("FAIL!\n");
    return false;
  }
#endif
}

bool Chain_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  if (this->_Head->matches_Bond(root, from, bond)) {
    if (this->_Tail.notnilp()) {
      Atom_sp other = bond->getOtherAtom(from);
      BondList_sp nextBonds = other->getBondList();
      nextBonds->removeBondBetween(from, other);
      if (this->_Tail->matches_BondList(root, other, nextBonds)) {
        goto SUCCESS;
      }
      goto FAIL;
    }
    goto SUCCESS;
  }
FAIL:
  LOG("FAIL");
  return false;
SUCCESS:
  LOG("SUCCESS");
  return true;
}

bool Chain_O::matches_BondList(Root_sp root, chem::Atom_sp from, chem::BondList_sp neighbors) {
  
  LOG("{}\natom: {} bondList: {}", this->asSmarts(), _rep_(from), _rep_(neighbors));
  gctools::Vec0<chem::Bond_sp>::iterator bi;
  chem::BondList_sp nextBonds, tempBondList;
  LOG("Chain_O matching pattern: {}", this->asSmarts());
  LOG("There are {} neighbors bondList: {}", neighbors->size(), neighbors->describeOthers(from));
  for (bi = neighbors->begin(); bi != neighbors->end(); bi++) {
    if (this->matches_Bond(root, from, *bi)) {
      goto SUCCESS;
    }
  }
  // FAIL:
  LOG("FAIL");
  return false;
SUCCESS:
  LOG("SUCCESS!");
  return true;
}

void Chain_O::fields(core::Record_sp node) {
  node->field_if_not_nil(INTERN_(kw, tail), this->_Tail);
  node->field_if_not_nil(INTERN_(kw, head), this->_Head);
  this->Base::fields(node);
}

CL_DEFMETHOD BondListMatchNode_sp Chain_O::chain_get_tail() { return this->_Tail; }

// ------ Branch

void Branch_O::initialize() {
  this->Base::initialize();
  this->_Left = nil<core::T_O>();
  this->_Right = nil<core::T_O>();
}
string Branch_O::asSmarts() const {
  
  stringstream ss;

  if (this->_Left.notnilp())
    ss << this->_Left->asSmarts();
  if (this->_Right.notnilp())
    ss << "(" << this->_Right->asSmarts() << ")";
  return ss.str();
}

core::T_sp Branch_O::children() {
  ql::list result;
  result << this->_Left;
  result << this->_Right;
  return result.cons();
}

uint Branch_O::depth() const { return (MAX(af_depth(this->_Left), af_depth(this->_Right))); }

bool Branch_O::matches_BondList(Root_sp root, chem::Atom_sp from, chem::BondList_sp neighbors) {
  
  LOG("{}\natom: {} bondList: {}", this->asSmarts(), _rep_(from), _rep_(neighbors));
  gctools::Vec0<chem::Bond_sp>::iterator bi;
  LOG("Branch_O matching pattern: {}", this->asSmarts());
  LOG("Neighbors bond list= {}", neighbors->describeOthers(from));
  for (bi = neighbors->begin(); bi != neighbors->end(); bi++) {

    BondList_sp leftBondList = BondList_O::create();
    leftBondList->addBond(*bi);
    LOG("Constructed left bond list");
    LOG("Left bond list = {}", leftBondList->describeOthers(from));
    BondList_sp rightBondList = BondList_O::create();
    for (auto ci = (*neighbors).begin(); ci != (*neighbors).end(); ++ci) {
      if (*bi != *ci)
        rightBondList->addBond(*ci);
    }
    //    GC_COPY(chem::BondList_O, rightBondList, *neighbors); // = RP_Copy<BondList_O>(neighbors);
    LOG("Right bond list after copy without left = {}", rightBondList->describeOthers(from));
    LOG("copied neighbors into right bond list and removed left bond");
    LOG("Checking if left matches");
    if (this->_Left->matches_BondList(root, from, leftBondList)) {
      LOG("Left matches");
      if (this->_Right.notnilp()) {
        LOG("Right branch is defined, checking if it matches");
        if (this->_Right->matches_BondList(root, from, rightBondList)) {
          LOG("Right matches does");
          goto SUCCESS;
#if 0
          // A match success - save it
          ChemInfoMatch_sp match =  gc::As<ChemInfoMatch_sp>(_sym_STARcurrent_matchSTAR->symbolValue());
          match->setMatches(true);
          match->saveTagLookup();
          success = true;
          // goto SUCCESS;
#endif
        } else {
          LOG("Right doesn't match, keep trying");
        }
      } else {
        LOG("Left matches and there is no Right");
        goto SUCCESS;
#if 0
        // A match success - save it
        ChemInfoMatch_sp match =  gc::As<ChemInfoMatch_sp>(_sym_STARcurrent_matchSTAR->symbolValue());
        match->setMatches(true);
        match->saveTagLookup();
        success = true;
        // goto SUCCESS;
#endif
      }
    }
  }
  LOG("No match");
  // FAIL:
  LOG("FAIL");
  return false;
SUCCESS:
  LOG("SUCCESS!");
  return true;
}

void Branch_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, left), this->_Left);
  node->field_if_not_nil(INTERN_(kw, right), this->_Right);
  this->Base::fields(node);
}

CL_DEFMETHOD void Branch_O::branch_set_left(core::T_sp left) { this->_Left = left; }

CL_DEFMETHOD void Branch_O::branch_set_right(core::T_sp right) { this->_Right = right; }

//      AfterMatchBondNode

void RootMatchNode_O::fields(core::Record_sp node) {
  this->Base::fields(node);
  // Nothing to do
}

// ------- AfterMatchBondToAtomTest

void AfterMatchBondToAtomTest_O::initialize() {
  
  this->_AtomTag1 = nil<core::Symbol_O>();
  this->_AtomTag2 = nil<core::Symbol_O>();
  this->_Bond = SABNoBond;
}

void AfterMatchBondToAtomTest_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, tag1), this->_AtomTag1);
  node->field(INTERN_(kw, tag2), this->_AtomTag2);
  node->/*pod_*/ field(INTERN_(kw, bond), this->_Bond);
  this->Base::fields(node);
}

bool AfterMatchBondToAtomTest_O::matches(Root_sp root) {
  
  SIMPLE_ERROR("Must implement");
}

// -------- AntechamberFocusAtomMatch

void AntechamberFocusAtomMatch_O::initialize() {
  this->Base::initialize();
  this->_ResidueNames = nil<core::T_O>();
  this->_AtomicNumber = 0;
  this->_NumberOfAttachedAtoms = 0;
  this->_NumberOfAttachedHydrogens = 0;
  this->_NumberOfElectronWithdrawingGroups = 0;
  this->_AtomicProperty = nil<core::T_O>();
}

void AntechamberFocusAtomMatch_O::fields(core::Record_sp node) {
  node->field_if_not_nil(INTERN_(kw, residueNames), this->_ResidueNames);
  node->field_if_not_nil(INTERN_(kw, atomicProperty), this->_AtomicProperty);
  node->/*pod_*/ field(INTERN_(kw, atomicNumber), this->_AtomicNumber);
  node->/*pod_*/ field(INTERN_(kw, attachedAtoms), this->_NumberOfAttachedAtoms);
  node->/*pod_*/ field(INTERN_(kw, attachedHs), this->_NumberOfAttachedHydrogens);
  node->/*pod_*/ field(INTERN_(kw, attachedEWGs), this->_NumberOfElectronWithdrawingGroups);
  this->Base::fields(node);
}

core::T_sp AntechamberFocusAtomMatch_O::children() {
  ql::list result;
  if (this->_AtomicProperty.notnilp())
    result << this->_AtomicProperty;
  return result.cons();
}

string AntechamberFocusAtomMatch_O::asSmarts() const {
  stringstream ss;
  ss << "[<antechamber>";
  // ss << _rep_(this->_ResidueNames) << "|";
  ss << "#" << this->_AtomicNumber;
  if (this->_NumberOfAttachedAtoms >= 0) {
    ss << "v" << this->_NumberOfAttachedAtoms;
  }
  if (this->_NumberOfAttachedHydrogens >= 0) {
    ss << "H" << this->_NumberOfAttachedHydrogens;
  }
  if (this->_NumberOfElectronWithdrawingGroups >= 0) {
    ss << "E" << this->_NumberOfElectronWithdrawingGroups;
  }
  ss << "]";
  if (this->_AtomicProperty.notnilp()) {
    ss << this->_AtomicProperty->asSmarts();
  }

  return ss.str();
}

bool AntechamberFocusAtomMatch_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  chem::Atom_sp neighbor, nn;
  if (this->_AtomicNumber >= 0) {
    LOG("Checking if atomic number({}) == expected({})", atom->getAtomicNumber(), this->_AtomicNumber);
    if (atom->getAtomicNumber() != this->_AtomicNumber)
      goto FAIL;
  }
  if (this->_NumberOfAttachedAtoms >= 0) {
    LOG("Checking if number of bonds({}) == expected({})", atom->numberOfBonds(), this->_NumberOfAttachedAtoms);
    if (atom->numberOfBonds() != this->_NumberOfAttachedAtoms)
      goto FAIL;
  }
  if (this->_NumberOfAttachedHydrogens >= 0) {
    LOG("Checking if number of attached hydrogens({}) == expected({})", atom->getBondedHydrogenCount(),
        this->_NumberOfAttachedHydrogens);
    if (atom->getBondedHydrogenCount() != this->_NumberOfAttachedHydrogens)
      goto FAIL;
  }
  if (this->_NumberOfElectronWithdrawingGroups >= 0) {
    if (atom->numberOfBonds() != 1) {
      SIMPLE_WARN("Testing NumberOfElectronWithdrawingGroups must have one bond - test {}", _rep_(this->asSmartPtr()));
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
    LOG("Checking if number of electron withdrawing groups({}) == expected({})", cnt, this->_NumberOfElectronWithdrawingGroups);
    if (cnt != this->_NumberOfElectronWithdrawingGroups)
      goto FAIL;
  }
  ANN(this->_AtomicProperty);
  if (this->_AtomicProperty.notnilp()) {
    if (!this->_AtomicProperty->matches_Atom(root, atom))
      goto FAIL;
  }
  // SUCCESS:
  LOG("SUCCESS!");
  return true;
FAIL:
  LOG("FAIL");
  return false;
}

// ------- AntechamberBondToAtomTest
CL_LISPIFY_NAME("make-antechamber-bond-to-atom-test");
CL_LAMBDA(element neighbors props tag);
CL_DEF_CLASS_METHOD AntechamberBondToAtomTest_sp AntechamberBondToAtomTest_O::create_args(core::Symbol_sp element, int neighbors,
                                                                                          AtomOrBondMatchNode_sp props,
                                                                                          core::Symbol_sp tag) {
  _G();
  auto obj = gctools::GC<AntechamberBondToAtomTest_O>::allocate();
  obj->_Element = element;
  obj->_Neighbors = neighbors;
  obj->_AtomProperties = props;
  obj->_Tag = tag;
  return obj;
};

core::T_sp AntechamberBondToAtomTest_O::children() {
  ql::list result;
  if (this->_AtomProperties.notnilp())
    result << this->_AtomProperties;
  return result.cons();
}

bool AntechamberBondToAtomTest_O::matchBasic(AntechamberRoot_sp root, chem::Atom_sp atom) {

  WildElementDict_sp dict;
  bool gotElement;
  dict = root->getElementWildCardDictionary();
  gotElement = false;
  ANN(dict);
  if (dict.notnilp()) {
    if (dict->recognizesWildName(this->_Element)) {
      LOG("Matching wildcard element has({}) == expecting ({})", atom->getElementAsString().c_str(), _rep_(this->_Element).c_str());
      if (!dict->recognizesWildNameElement(this->_Element, atom->getElementAsSymbol()))
        goto FAIL;
      gotElement = true;
    }
  }
  if (!gotElement) {
    LOG("Trying to match regular element({}) == expected({})", atom->getElementAsString().c_str(), _rep_(this->_Element));
    if (atom->getElementAsSymbol() != this->_Element)
      goto FAIL;
  }
  LOG("Trying to match number of bonds({}) == expected({})", atom->numberOfBonds(), this->_Neighbors);
  if (this->_Neighbors >= 0) {
    if (atom->numberOfBonds() != this->_Neighbors)
      goto FAIL;
  }
  if (this->_Tag.notnilp()) {
    SIMPLE_WARN("{}:{} What do I do with the tag {} for AntechamberBondToAtomTest {}", __FILE__, __LINE__, _rep_(this->_Tag),
                _rep_(this->asSmartPtr()));
    //    root->getMatch()->defineAtomTag(atom, this->_Tag);
  }
  // SUCCESS:
  LOG("SUCCESS!");
  return true;
FAIL:
  LOG("FAIL");
  return false;
}

bool AntechamberBondToAtomTest_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  AntechamberRoot_sp acRoot;
  if (root->type() != antechamberRoot) {
    SIMPLE_ERROR("AntechamberBondToAtomTest::matches requires an AntechamberRoot");
  }
  acRoot = (root).as<AntechamberRoot_O>();
  if (!this->matchBasic(acRoot, atom))
    goto FAIL;
  if (this->_AtomProperties->matches_Atom(acRoot, atom))
    goto SUCCESS;
FAIL:
  LOG("FAIL");
  return false;
SUCCESS:
  LOG("SUCCESS!");
  return true;
}

string AntechamberBondToAtomTest_O::asSmarts() const {
  stringstream ss;
  ss << "[<AntechamberBondToAtomTest>";
  ss << "<element" << _rep_(this->_Element) << ">";
  ss << "v" << (this->_Neighbors);
  ss << this->_AtomProperties->asSmarts() << "}";
  return ss.str();
}

bool AntechamberBondToAtomTest_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  AntechamberRoot_sp antechamberRoot;
  if (root->type() != chem::antechamberRoot) {
    stringstream ss;
    ss << "AntechamberBondToAtomTest::matches requires an AntechamberRoot "
       << "instead we have" << (unsigned long)(&*root) << " of type: " << root->type();
    SIMPLE_ERROR("{}", ss.str());
  }
  antechamberRoot = (root).as<AntechamberRoot_O>();
  if (!this->matchBasic(antechamberRoot, bond->getOtherAtom(from)))
    goto FAIL;
  ANN(this->_AtomProperties);
  if (this->_AtomProperties.notnilp()) {
    if (!this->_AtomProperties->matches_Bond(antechamberRoot, from, bond))
      goto FAIL;
  }
  LOG("Matched!!");
  // SUCCESS:
  LOG("SUCCESS!");
  return true;
FAIL:
  LOG("FAIL");
  return false;
}

void AntechamberBondToAtomTest_O::initialize() {
  this->Base::initialize();
  this->_Element = nil<core::Symbol_O>();
  this->_Neighbors = 0;
  this->_Tag = nil<core::Symbol_O>();
  this->_AtomProperties = nil<core::T_O>();
}

void AntechamberBondToAtomTest_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, element), this->_Element);
  node->/*pod_*/ field(INTERN_(kw, neighbors), this->_Neighbors);
  node->field(INTERN_(kw, tag), this->_Tag);
  node->field_if_not_nil(INTERN_(kw, atomProps), this->_AtomProperties);
  this->Base::fields(node);
}

// ------- Root

core::HashTable_sp Root_O::lazyTests() {
  if (this->_Tests.nilp()) {
    this->_Tests = core::HashTable_O::createEq();
  }
  return gc::As_unsafe<core::HashTable_sp>(this->_Tests);
}

uint Root_O::depth() const {
  
  uint res = af_depth(this->_Node);
  LOG("Returning {}", res);
  return res;
}

/*! Set the tests using an alist of (name . function) */
void Root_O::setTests(core::List_sp tests) {
  for (auto cur : tests) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    core::T_sp key = oCar(pair);
    core::T_sp value = oCdr(pair);
    this->lazyTests()->setf_gethash(key, value);
  }
}

core::T_sp Root_O::children() {
  ql::list result;
  if (this->_Node.notnilp())
    result << this->_Node;
  return result.cons();
};

string Root_O::asSmarts() const {
  stringstream ss;
  if (this->_Node.notnilp()) {
    ss << this->_Node->asSmarts();
  }
  return ss.str();
}

void Root_O::initialize() { this->Base::initialize(); }

void Root_O::addTest(core::Symbol_sp testSym, core::Function_sp testCode) {
  
  LOG("Adding test<{}> with code: {}", _rep_(testSym), _rep_(testCode));
  this->lazyTests()->setf_gethash(testSym, testCode);
}

bool Root_O::evaluateTest(core::Symbol_sp testSym, Atom_sp atom) {
  ASSERTF(testSym.notnilp(), ("The test symbol was nil! - this should never occur"));
  LOG("Looking up test with symbol<{}>", _rep_(testSym));
  auto find = lazyTests()->find(testSym);
  if (!find) {
    SIMPLE_ERROR("Could not find named ChemInfo/Smarts test[{}] in Smarts object - available named tests are[{}]", _rep_(testSym),
                 this->lazyTests()->keysAsString());
  }
  core::Function_sp func = core::coerce::functionDesignator(*find);
  if (!gctools::IsA<core::Function_sp>(func)) {
    SIMPLE_ERROR("The test ChemInfo/Smarts test[{}] must be a function - instead it is a {}", _rep_(testSym), _rep_(*find));
  }
  core::Function_sp testCode = gctools::As_unsafe<core::Function_sp>(func);
  ASSERTF(testCode.notnilp(), ("testCode was nil - it should never be"));
  ASSERTF(atom.notnilp(), ("The atom arg should never be nil"));
  core::T_sp res = core::eval::funcall(testCode, atom);
  return res.isTrue();
}

void Root_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, node), this->_Node);
  node->field_if_not_nil(INTERN_(kw, tests), this->_Tests);
  this->Base::fields(node);
}

bool Root_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  return this->matches_Atom(root, bond->getOtherAtom(from));
};

bool Root_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  chem::BondList_sp nextBonds;
  bool matches;
  matches = false;
  if (this->_Node.notnilp()) {
    LOG("_Node is notNil - testing");
    if (!this->_Node->matches_Atom(root, atom)) {
      goto FAIL;
    }
  }
  // SUCCESS:
  LOG("SUCCESS!!!");
  return true;
FAIL:
  LOG("FAIL");
  return false;
}

CL_LISPIFY_NAME("make-smarts-root");
CL_DEF_CLASS_METHOD SmartsRoot_sp SmartsRoot_O::make(const std::string &code, ChemInfoNode_sp cinode, size_t maxTag) {
  auto obj = gctools::GC<SmartsRoot_O>::allocate(code, cinode, maxTag); // RP_Create<SmartsRoot_O>(lisp);
  //  printf("%s:%d:%s  cinode-> %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(cinode).c_str());
  return obj;
};

void SmartsRoot_O::fields(core::Record_sp node) { this->Base::fields(node); }

bool SmartsRoot_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  return this->matches_Atom(root, bond->getOtherAtom(from));
};

bool SmartsRoot_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  chem::BondList_sp nextBonds;
  bool matches;
  matches = false;
  if (!this->Root_O::matches_Atom(root, atom))
    goto FAIL;
  // SUCCESS:
  LOG("SUCCESS!");
  return true;
FAIL:
  LOG("FAIL");
  return false;
}

void AntechamberRoot_O::initialize() {
  this->Base::initialize();
  this->_AfterMatchTests = nil<core::T_O>();
  this->_WildElementDictionary = nil<core::T_O>();
}

void AntechamberRoot_O::fields(core::Record_sp node) {
  node->field_if_not_nil(INTERN_(kw, afterMatchTests), this->_AfterMatchTests);
  node->field_if_not_nil(INTERN_(kw, wildDict), this->_WildElementDictionary);
  this->Base::fields(node);
}

core::T_sp AntechamberRoot_O::children() {
  ql::list result;
  result &this->Base::children();
  if (this->_AfterMatchTests.notnilp())
    result << this->_AfterMatchTests;
  return result.cons();
};

string AntechamberRoot_O::descriptionOfContents() const {
  stringstream ss;
  return ss.str();
}

bool AntechamberRoot_O::matches_Bond(Root_sp root, chem::Atom_sp from, chem::Bond_sp bond) {
  
  LOG("{}\natom: {} bond: {}", this->asSmarts(), _rep_(from), _rep_(bond));
  return this->matches_Atom(root, bond->getOtherAtom(from));
};

bool AntechamberRoot_O::matches_Atom(Root_sp root, chem::Atom_sp atom) {
  
  LOG("{}\natom: {}", this->asSmarts(), _rep_(atom));
  chem::BondList_sp nextBonds;
  bool matches;
  matches = false;
  if (!this->Base::matches_Atom(root, atom)) {
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
  // SUCCESS:
  LOG("SUCCESS!");
  return true;
FAIL:
  LOG("FAIL");
  return false;
}

CL_DEFMETHOD core::Vector_sp ChemInfoMatch_O::tags_as_vector() const { return this->_TagLookup; }

CL_DEFMETHOD core::HashTable_sp ChemInfoMatch_O::tags_as_hashtable() const {
  core::HashTable_sp ht = core::HashTable_O::createEql();
  for (size_t ii = 0; ii < this->_TagLookup->length(); ++ii) {
    if (this->_TagLookup->rowMajorAref(ii).notnilp()) {
      ht->setf_gethash(core::make_fixnum(ii), this->_TagLookup->rowMajorAref(ii));
    }
  }
  return ht;
};

CL_LAMBDA(code &key tests);
DOCGROUP(cando);
CL_DEFUN SmartsRoot_sp chem__compile_smarts(const string &code, core::List_sp tests) {
  //  printf("%s:%d:%s code: %s\n", __FILE__, __LINE__, __FUNCTION__, code.c_str());
  core::SimpleBaseString_sp scode = core::SimpleBaseString_O::make(code);
  core::T_sp tnode = core::eval::funcall(_sym_parse_smarts, scode);
  if (tnode.nilp()) {
    SIMPLE_ERROR("When calling chem:parse-smarts with {} - it returned NIL", _rep_(scode));
  }
  ChemInfoNode_sp node = gc::As<ChemInfoNode_sp>(tnode);
  size_t max_tag = calculate_max_tags(node);
  SmartsRoot_sp root = SmartsRoot_O::make(code, node, max_tag);
  root->setTests(tests);
  return root;
}

CL_LAMBDA(code &key tests);
DOCGROUP(cando);
CL_DEFUN Smirks_sp chem__compile_smirks(const string &code, core::List_sp tests) {
  core::SimpleBaseString_sp scode = core::SimpleBaseString_O::make(code);
  core::T_mv nodes = core::eval::funcall(_sym_parse_smirks, scode);
  core::MultipleValues &values = core::lisp_multipleValues();
  ChemInfoNode_sp reactant = gc::As<ChemInfoNode_sp>(nodes);
  ChemInfoNode_sp product = gc::As<ChemInfoNode_sp>(values.second(nodes.number_of_values()));
  size_t max_reactant_tag = calculate_max_tags(reactant);
  size_t max_product_tag = calculate_max_tags(product);
  if (max_reactant_tag < max_product_tag) {
    SIMPLE_ERROR("The largest tag in the product must be equal to or less than the largest tag in the reactant");
  }
  SmartsRoot_sp reactant_root = SmartsRoot_O::make(code, reactant, max_reactant_tag);
  SmartsRoot_sp product_root = SmartsRoot_O::make(code, product, max_product_tag);
  reactant_root->setTests(tests);
  product_root->setTests(tests);
  return Smirks_O::make(reactant_root, product_root);
}

DOCGROUP(cando);
CL_DEFUN AntechamberRoot_mv chem__compile_antechamber(const string &code, WildElementDict_sp dict) {
  if (!chem::_sym_compile_antechamber_type_rule->fboundp()) {
    SIMPLE_ERROR("chem:compile-antechamber-type-rule is not fbound");
  }
  LOG("Compiling code: {}", code.c_str());
  core::SimpleBaseString_sp scode = core::SimpleBaseString_O::make(code);
  core::T_mv troot_mv =
      core::eval::funcall(chem::_sym_compile_antechamber_type_rule->symbolFunction(), nil<core::T_O>(), scode, dict);
  if (troot_mv.nilp()) {
    SIMPLE_ERROR("Error parsing antechamber type rule {}", code);
  }
  core::MultipleValues &values = core::lisp_multipleValues();
  AntechamberRoot_sp root = gc::As<AntechamberRoot_sp>(troot_mv);
  root->_Code = code;
  core::T_sp type = values.second(troot_mv.number_of_values());
  root->setElementWildCardDictionary(dict);
  LOG("antechamber compile was successful");
  return Values(root, type);
}

DOCGROUP(cando);
CL_DEFUN core::T_mv chem__chem_info_match(Root_sp testRoot, Atom_sp atom) {
  core::HashTable_sp ringHashTable = core::HashTable_O::createEql();
  ChemInfoMatch_sp current_match = ChemInfoMatch_O::make(testRoot, testRoot->_MaxTag, ringHashTable);
  core::DynamicScopeManager scope(_sym_STARcurrent_matchSTAR, current_match);
  bool matches = testRoot->matches_Atom(testRoot, atom);
  //  bool matches = current_match->matches();
  return Values(_lisp->_boolean(matches), current_match);
}

DOCGROUP(cando);
CL_DEFUN core::T_sp chem__matches(Root_sp testRoot, Atom_sp atom) {
  core::T_mv matches_mv = chem__chem_info_match(testRoot, atom);
  core::MultipleValues &values = core::lisp_multipleValues();
  if (matches_mv.notnilp()) {
    return values.second(matches_mv.number_of_values());
  }
  return nil<core::T_O>();
}

DOCGROUP(cando);
CL_DEFUN core::T_sp chem__chem_info_node_children(ChemInfoNode_sp node) { return node->children(); }

/*! Hold nodes for the Gaff and Msmarts parsers - rewrite these in Common Lisp */
SYMBOL_EXPORT_SC_(ChemPkg, STARparserNodeHolderSTAR);

SYMBOL_EXPORT_SC_(ChemPkg, STARSarRingTestEnumConverterSTAR);
SYMBOL_EXPORT_SC_(KeywordPkg, SARNone);
SYMBOL_EXPORT_SC_(KeywordPkg, SARRingSet);
SYMBOL_EXPORT_SC_(KeywordPkg, SARRingTest);
SYMBOL_EXPORT_SC_(ChemPkg, STARSabBondEnumConverterSTAR);
CL_BEGIN_ENUM(RingTestEnum, _sym_STARSarRingTestEnumConverterSTAR, "SARRingTestEnum");
CL_VALUE_ENUM(kw::_sym_SARNone, SARNone);
CL_VALUE_ENUM(kw::_sym_SARRingSet, SARRingSet);
CL_VALUE_ENUM(kw::_sym_SARRingTest, SARRingTest);
CL_END_ENUM(_sym_STARSarRingTestEnumConverterSTAR);

SYMBOL_EXPORT_SC_(KeywordPkg, SABUseBondMatcher);
SYMBOL_EXPORT_SC_(KeywordPkg, SABNoBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABSingleBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABSingleOrAromaticBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABDoubleOrAromaticBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABTripleOrAromaticBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABDelocalizedBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABDoubleBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABTripleBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABAromaticBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABAnyBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SABDirectionalSingleUpOrUnspecified);
SYMBOL_EXPORT_SC_(KeywordPkg, SABDirectionalSingleDownOrUnspecified);
SYMBOL_EXPORT_SC_(KeywordPkg, SABDirectionalSingleUp);
SYMBOL_EXPORT_SC_(KeywordPkg, SABDirectionalSingleDown);
SYMBOL_EXPORT_SC_(KeywordPkg, SABSameRingBond);
SYMBOL_EXPORT_SC_(ChemPkg, STARSabBondEnumConverterSTAR);
CL_BEGIN_ENUM(BondEnum, _sym_STARSabBondEnumConverterSTAR, "SABBondEnum");
CL_VALUE_ENUM(kw::_sym_SABUseBondMatcher, SABUseBondMatcher);
CL_VALUE_ENUM(kw::_sym_SABNoBond, SABNoBond);
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
CL_VALUE_ENUM(kw::_sym_SABSameRingBond, SABSameRingBond);
CL_END_ENUM(_sym_STARSabBondEnumConverterSTAR);

SYMBOL_EXPORT_SC_(KeywordPkg, logAlwaysTrue);
SYMBOL_EXPORT_SC_(KeywordPkg, logIdentity);
SYMBOL_EXPORT_SC_(KeywordPkg, logNot);
SYMBOL_EXPORT_SC_(KeywordPkg, logHighPrecedenceAnd);
SYMBOL_EXPORT_SC_(KeywordPkg, logOr);
SYMBOL_EXPORT_SC_(KeywordPkg, logLowPrecedenceAnd);
SYMBOL_EXPORT_SC_(ChemPkg, STARLogicalOperatorTypeConverterSTAR);
CL_BEGIN_ENUM(LogicalOperatorType, _sym_STARLogicalOperatorTypeConverterSTAR, "LogicalOperatorType");
CL_VALUE_ENUM(kw::_sym_logAlwaysTrue, logAlwaysTrue);
CL_VALUE_ENUM(kw::_sym_logIdentity, logIdentity);
CL_VALUE_ENUM(kw::_sym_logNot, logNot);
CL_VALUE_ENUM(kw::_sym_logHighPrecedenceAnd, logHighPrecedenceAnd);
CL_VALUE_ENUM(kw::_sym_logOr, logOr);
CL_VALUE_ENUM(kw::_sym_logLowPrecedenceAnd, logLowPrecedenceAnd);
CL_END_ENUM(_sym_STARLogicalOperatorTypeConverterSTAR);

SYMBOL_EXPORT_SC_(KeywordPkg, SAPNone);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPWildCard);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPDegree);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPElement);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPTotalHCount);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPImplicitHCount);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPRingMembershipCount);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPRingConnectivity);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPRingTest);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPRingSize);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPValence);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPConnectivity);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPNegativeCharge);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPNegativeFormalCharge);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPPositiveCharge);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPPositiveFormalCharge);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPAtomicNumber);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPChiralityAntiClockwise);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPChiralityClockwise);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPAtomicMass);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPLonePair);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPTotalBondNumber);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPTotalExplicitBondNumber);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPPiBondOrbital);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPAromaticPiElectron);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPHeavyAtomTotalBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPGroupNumber);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPElectronegativeElement);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPTransitionMetal);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPBondedToPrevious);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPNotBondedToPrevious);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPInBond);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPArLevel);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPNoRing);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPResidueTest);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPPredicateName);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPAromaticElement);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPAliphatic);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPAromatic);
SYMBOL_EXPORT_SC_(KeywordPkg, SAPAtomMap);

SYMBOL_EXPORT_SC_(ChemPkg, STARAtomTestEnumConverterSTAR);
CL_BEGIN_ENUM(AtomTestEnum, _sym_STARAtomTestEnumConverterSTAR, "AtomTestEnum");
CL_VALUE_ENUM(kw::_sym_SAPNone, SAPNone);
CL_VALUE_ENUM(kw::_sym_SAPWildCard, SAPWildCard);
CL_VALUE_ENUM(kw::_sym_SAPDegree, SAPDegree);
CL_VALUE_ENUM(kw::_sym_SAPElement, SAPElement);
CL_VALUE_ENUM(kw::_sym_SAPTotalHCount, SAPTotalHCount);
CL_VALUE_ENUM(kw::_sym_SAPImplicitHCount, SAPImplicitHCount);
CL_VALUE_ENUM(kw::_sym_SAPRingMembershipCount, SAPRingMembershipCount);
CL_VALUE_ENUM(kw::_sym_SAPRingConnectivity, SAPRingConnectivity);
CL_VALUE_ENUM(kw::_sym_SAPRingTest, SAPRingTest);
CL_VALUE_ENUM(kw::_sym_SAPRingSize, SAPRingSize);
CL_VALUE_ENUM(kw::_sym_SAPValence, SAPValence);
CL_VALUE_ENUM(kw::_sym_SAPConnectivity, SAPConnectivity);
CL_VALUE_ENUM(kw::_sym_SAPNegativeCharge, SAPNegativeCharge);
CL_VALUE_ENUM(kw::_sym_SAPNegativeFormalCharge, SAPNegativeFormalCharge);
CL_VALUE_ENUM(kw::_sym_SAPPositiveCharge, SAPPositiveCharge);
CL_VALUE_ENUM(kw::_sym_SAPPositiveFormalCharge, SAPPositiveFormalCharge);
CL_VALUE_ENUM(kw::_sym_SAPAtomicNumber, SAPAtomicNumber);
CL_VALUE_ENUM(kw::_sym_SAPChiralityAntiClockwise, SAPChiralityAntiClockwise);
CL_VALUE_ENUM(kw::_sym_SAPChiralityClockwise, SAPChiralityClockwise);
CL_VALUE_ENUM(kw::_sym_SAPAtomicMass, SAPAtomicMass);
CL_VALUE_ENUM(kw::_sym_SAPLonePair, SAPLonePair);
CL_VALUE_ENUM(kw::_sym_SAPTotalBondNumber, SAPTotalBondNumber);
CL_VALUE_ENUM(kw::_sym_SAPTotalExplicitBondNumber, SAPTotalExplicitBondNumber);
CL_VALUE_ENUM(kw::_sym_SAPPiBondOrbital, SAPPiBondOrbital);
CL_VALUE_ENUM(kw::_sym_SAPAromaticPiElectron, SAPAromaticPiElectron);
CL_VALUE_ENUM(kw::_sym_SAPHeavyAtomTotalBond, SAPHeavyAtomTotalBond);
CL_VALUE_ENUM(kw::_sym_SAPGroupNumber, SAPGroupNumber);
CL_VALUE_ENUM(kw::_sym_SAPElectronegativeElement, SAPElectronegativeElement);
CL_VALUE_ENUM(kw::_sym_SAPTransitionMetal, SAPTransitionMetal);
CL_VALUE_ENUM(kw::_sym_SAPBondedToPrevious, SAPBondedToPrevious);
CL_VALUE_ENUM(kw::_sym_SAPNotBondedToPrevious, SAPNotBondedToPrevious);
CL_VALUE_ENUM(kw::_sym_SAPInBond, SAPInBond);
CL_VALUE_ENUM(kw::_sym_SAPArLevel, SAPArLevel);
CL_VALUE_ENUM(kw::_sym_SAPNoRing, SAPNoRing);
CL_VALUE_ENUM(kw::_sym_SAPResidueTest, SAPResidueTest);
CL_VALUE_ENUM(kw::_sym_SAPPredicateName, SAPPredicateName);
CL_VALUE_ENUM(kw::_sym_SAPAromaticElement, SAPAromaticElement);
CL_VALUE_ENUM(kw::_sym_SAPAliphatic, SAPAliphatic);
CL_VALUE_ENUM(kw::_sym_SAPAromatic, SAPAromatic);
CL_VALUE_ENUM(kw::_sym_SAPAtomMap, SAPAtomMap);
CL_END_ENUM(_sym_STARAtomTestEnumConverterSTAR);

GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ChemInfoNode_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::BondToAtomTest_O>>);

MoleculeGraph_O::MoleculeGraph_O(Molecule_sp matter){};

MoleculeGraph_O::MoleculeGraph_O() : _num_edges(0){};

MoleculeGraph_O::~MoleculeGraph_O() {
  delete this->_moleculeGraph;
  this->_moleculeGraph = nullptr;
}

CL_DEFMETHOD size_t MoleculeGraph_O::num_vertices() { return this->_nodes->length(); }

CL_DEFMETHOD size_t MoleculeGraph_O::num_edges() { return this->_num_edges; }

CL_DEFMETHOD size_t MoleculeGraph_O::add_vertex(core::T_sp vertex) {
  size_t index = this->_nodes->length();
  this->_nodes_to_index->setf_gethash(vertex, core::make_fixnum(index));
  this->_nodes->vectorPushExtend(vertex);
  boost::add_vertex(MoleculeVertexData(index), *this->_moleculeGraph);
  return index;
}

CL_DEFMETHOD size_t MoleculeGraph_O::vertex_index(core::T_sp vertex) {
  core::T_sp tindex = this->_nodes_to_index->gethash(vertex);
  if (tindex.nilp()) {
    SIMPLE_ERROR("Could not find vertex {}", _rep_(vertex));
  }
  return tindex.unsafe_fixnum();
}

CL_DEFMETHOD core::T_sp MoleculeGraph_O::get_vertex(size_t index) {
  if (index >= this->_nodes->length()) {
    SIMPLE_ERROR("index {} must be less than {}", index, this->_nodes->length());
  }
  return this->_nodes->rowMajorAref(index);
}

CL_DEFMETHOD void MoleculeGraph_O::add_edge(size_t v1, size_t v2, BondOrder bo) {
  int iv1 = v1;
  int iv2 = v2;
  boost::add_edge(iv1, iv2, bo, *(this->_moleculeGraph));
  ++this->_num_edges;
}

CL_DEFMETHOD void MoleculeGraph_O::walk_vertices(core::T_sp callback) {
  for (size_t ii = 0; ii < this->_nodes->length(); ++ii) {
    core::eval::funcall(callback, this->_nodes->rowMajorAref(ii));
  }
}

CL_DEFMETHOD void MoleculeGraph_O::walk_edges(core::T_sp callback) {
  // Ok, we want to see that all our edges are now contained in the graph
  typedef boost::graph_traits<MoleculeGraphType>::edge_iterator edge_iterator;

  // Tried to make this section more clear, instead of using tie, keeping all
  // the original types so it's more clear what is going on
  boost::property_map<MoleculeGraphType, boost::edge_weight_t>::type moleculeEdge =
      boost::get(boost::edge_weight_t(), *(this->_moleculeGraph));
  std::pair<edge_iterator, edge_iterator> ei = edges(*this->_moleculeGraph);
  for (edge_iterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
    core::T_sp a1 = this->_nodes->rowMajorAref(source(*edge_iter, *this->_moleculeGraph));
    core::T_sp a2 = this->_nodes->rowMajorAref(target(*edge_iter, *this->_moleculeGraph));
    BondOrder bo = moleculeEdge[*edge_iter];
    core::eval::funcall(callback, a1, a2, translate::to_object<BondOrder>::convert(bo));
  }
}

void MoleculeGraph_O::initialize() {
  this->_nodes_to_index = core::HashTable_O::createEq();
  this->_moleculeGraph = nullptr;
  this->_nodes = core::ComplexVector_T_O::make(64, nil<core::T_O>(), core::make_fixnum(0));
}

DOCGROUP(cando);
CL_DEFUN MoleculeGraph_sp chem__make_molecule_graph() {
  auto graph = gctools::GC<MoleculeGraph_O>::allocate();
  graph->_moleculeGraph = new MoleculeGraphType();
  return graph;
}

CL_LAMBDA(matter &optional (exclude_hydrogens nil));
DOCGROUP(cando);
CL_DEFUN MoleculeGraph_sp chem__make_molecule_graph_from_molecule(Molecule_sp matter, bool exclude_hydrogens) {
  auto graph = gctools::GC<MoleculeGraph_O>::allocate(matter);
  Loop lMol;
  Loop lAtoms;
  lAtoms.loopTopGoal(matter, ATOMS);
  while (lAtoms.advanceLoopAndProcess()) {
    Atom_sp atom = lAtoms.getAtom();
    if (!exclude_hydrogens || (exclude_hydrogens && atom->getAtomicNumber() != 1)) {
      graph->_nodes_to_index->setf_gethash(atom, core::make_fixnum(graph->_nodes->length()));
      graph->_nodes->vectorPushExtend(atom);
    }
  }

  // create an -undirected- graph type, using vectors as the underlying containers
  // and an adjacency_list as the basic representation
  graph->_moleculeGraph = new MoleculeGraphType();
  for (size_t iii = 0; iii < graph->_nodes->length(); iii++) {
    add_vertex(MoleculeVertexData(iii), *graph->_moleculeGraph);
  }

  // Example uses an array, but we can easily use another container type
  // to hold our edges.
  Loop lbonds;
  lbonds.loopTopGoal(matter, BONDS);
  while (lbonds.advanceLoopAndProcess()) {
    Atom_sp a1 = lbonds.getAtom1();
    Atom_sp a2 = lbonds.getAtom2();
    if (!exclude_hydrogens || (exclude_hydrogens && (a1->getAtomicNumber() != 1) && (a2->getAtomicNumber() != 1))) {
      BondOrder bo = lbonds.getBondOrder();
      int atom1_index = graph->_nodes_to_index->gethash(a1).unsafe_fixnum();
      int atom2_index = graph->_nodes_to_index->gethash(a2).unsafe_fixnum();
      boost::add_edge(atom1_index, atom2_index, bo, *graph->_moleculeGraph);
    }
  }
  return graph;
}

DOCGROUP(cando);
CL_DEFUN void chem__molecule_graph_dump(MoleculeGraph_sp graph) {
  core::clasp_write_string(fmt::format("Number of vertices: {}\n", num_vertices(*graph->_moleculeGraph)));
  core::clasp_write_string(fmt::format("Number of edges: {}\n", num_edges(*graph->_moleculeGraph)));

  // Ok, we want to see that all our edges are now contained in the graph
  typedef boost::graph_traits<MoleculeGraphType>::edge_iterator edge_iterator;

  // Tried to make this section more clear, instead of using tie, keeping all
  // the original types so it's more clear what is going on
  std::pair<edge_iterator, edge_iterator> ei = edges(*graph->_moleculeGraph);
  for (edge_iterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
    core::T_sp a1 = graph->_nodes->rowMajorAref(source(*edge_iter, *graph->_moleculeGraph));
    core::T_sp a2 = graph->_nodes->rowMajorAref(target(*edge_iter, *graph->_moleculeGraph));
    core::clasp_write_string(fmt::format("({} - {})\n", _rep_(a1), _rep_(a2)));
  }
}

ChemInfoGraph_O::ChemInfoGraph_O(Root_sp root) : _Root(root), _chemInfoGraph(nullptr){};
ChemInfoGraph_O::~ChemInfoGraph_O() {
  if (this->_chemInfoGraph) {
    delete this->_chemInfoGraph;
    this->_chemInfoGraph = nullptr;
  }
}

void ChemInfoGraph_O::initialize() {
  this->_nodes_to_index = core::HashTable_O::createEq();
}

CL_DOCSTRING(R"dx(Make a chem:chem-info-graph from a chem:root object)dx");
DOCGROUP(cando);
CL_DEFUN ChemInfoGraph_sp chem__make_chem_info_graph(Root_sp pattern) {
  auto graph = gctools::GC<ChemInfoGraph_O>::allocate(pattern);
  graph->buildFromRoot_();
  return graph;
}

#define DEBUG_BUILD_FROM_ROOT 0
#define DEBUG_SETF_GETHASH 0

void ChemInfoGraph_O::buildFromRoot_() {
  CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
    core::clasp_write_string("Starting make-chem-info-graph\n");
    });
  this->_nodeOrder.clear();
  this->_atomNodes.clear();
  this->_bondNodes.clear();
  this->_chemInfoGraph = NULL;
  ChemInfoGraph_sp graph = this->asSmartPtr();
  Root_sp pattern = graph->_Root;
  std::vector<RingClosers> closers;
#if DEBUG_SETF_GETHASH
  std::map<void*,size_t> save_pointers;
#endif
  walk_nodes_with_parent(nil<core::T_O>(), pattern->_Node, [&graph, &closers
#if DEBUG_SETF_GETHASH
                                                            ,&save_pointers
#endif
                                                            ](core::T_sp parentOrNil, ChemInfoNode_sp node) {
    //      printf("%s:%d ChemInfoNode_sp -> %s\n", __FILE__, __LINE__, _rep_(node).c_str());
    if (gc::IsA<Chain_sp>(node)) {
      Chain_sp chain = gc::As_unsafe<Chain_sp>(node);
      ChemInfoNode_sp head = gc::As<ChemInfoNode_sp>(chain->_Head);
      CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
        core::clasp_write_string(fmt::format("Converting chain to chem-info-graph nodes head: {}\n", _rep_(head)));
        });
      if (gc::IsA<BondToAtomTest_sp>(head)) {
        head = gc::As_unsafe<BondToAtomTest_sp>(head)->_AtomTest;
      }
      size_t node_index = graph->_atomNodes.size();
      graph->_nodes_to_index->setf_gethash(head, core::make_fixnum(node_index));
#if DEBUG_SETF_GETHASH
      save_pointers[(void*)head.raw_()] = node_index;
#endif
      graph->_atomNodes.push_back(head);
      AtomOrBondMatchNode_sp ahead = gc::As<AtomOrBondMatchNode_sp>(head);
      if (ahead->_RingTest == SARRingSet) {
        CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
          core::clasp_write_string(fmt::format("Found SARRingSet _RingId = {}\n", ahead->_RingId));
          });
        if (ahead->_RingId >= closers.size()) {
          closers.resize(ahead->_RingId + 1);
        }
        closers[ahead->_RingId]._Active = true;
        closers[ahead->_RingId]._NodeIndex = node_index;
      } else if (ahead->_RingTest == SARRingTest) {
        CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
          core::clasp_write_string(fmt::format("Found SARRingTest _RingId = {}\n", ahead->_RingId));
          });
        if (ahead->_RingId >= closers.size()) {
          closers.resize(ahead->_RingId + 1);
        }
        if (!closers[ahead->_RingId]._Active) {
          SIMPLE_ERROR("Hit a SARRingTest for _RingId {} but the SARRingSet hasn't been seen", ahead->_RingId);
        }
        closers[ahead->_RingId]._Bonds.push_back(RingBond(SABSingleOrAromaticBond, node_index));
      }
      return true;
    } else if (gc::IsA<Branch_sp>(node)) {
      return true;
    } else if (parentOrNil.nilp() && gc::IsA<AtomTest_sp>(node)) {
      AtomTest_sp atomTest = gc::As_unsafe<AtomTest_sp>(node);
      CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
        core::clasp_write_string(fmt::format("Converting AtomTest to chem-info-graph nodes head: {}\n", _rep_(atomTest)));
        });
      size_t node_index = graph->_atomNodes.size();
      graph->_nodes_to_index->setf_gethash(atomTest, core::make_fixnum(node_index));
#if DEBUG_SETF_GETHASH
      save_pointers[(void*)atomTest.raw_()] = node_index;
#endif
      graph->_atomNodes.push_back(atomTest);
      return true;
    } else if (parentOrNil.nilp() && gc::IsA<Logical_sp>(node)) {
      Logical_sp logical = gc::As_unsafe<Logical_sp>(node);
      CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
        core::clasp_write_string(fmt::format("Converting logical to chem-info-graph nodes head: {}\n", _rep_(logical)));
        });
      size_t node_index = graph->_atomNodes.size();
#if DEBUG_SETF_GETHASH
      save_pointers[(void*)logical.raw_()] = node_index;
#endif
      graph->_nodes_to_index->setf_gethash(logical, core::make_fixnum(node_index));
      graph->_atomNodes.push_back(logical);
      return true;
    } else {
      return false;
    }
  });
  graph->_chemInfoGraph = new ChemInfoGraphType(); // (graph->_nodes.size());
                                                   // Our set of edges, which basically are just converted into ints (0-4)

  // Transform a ChemInfo graph into a boost graph
  // This is tricky code
  core::HashTable_sp parent_nodes = core::HashTable_O::createEq();
  walk_nodes_with_parent(
      nil<core::T_O>(), pattern->_Node, [&graph, &parent_nodes, &closers
#if DEBUG_SETF_GETHASH
                                                            ,&save_pointers
#endif
                                         ](core::T_sp parentOrNil, ChemInfoNode_sp node) {
        CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(2)) {
          core::clasp_write_string(fmt::format("Walking pattern parent: {} node: {}\n", _rep_(parentOrNil), _rep_(node)));
          });
        if (parentOrNil.nilp() && gc::IsA<AtomTest_sp>(node)) {
          AtomTest_sp atomTest = gc::As_unsafe<AtomTest_sp>(node);
          CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
            core::clasp_write_string(fmt::format("AtomTest node id({})\n", atomTest->_Id));
            });
          core::T_sp thead = atomTest;
          if (thead.nilp())
            SIMPLE_ERROR("The head of a chain is NIL");
          ChemInfoNode_sp head = gc::As<ChemInfoNode_sp>(thead);
          AtomOrBondMatchNode_sp ahead = gc::As<AtomOrBondMatchNode_sp>(head);
          core::T_sp head_index = graph->_nodes_to_index->gethash(ahead);
          if (!head_index.fixnump())
            SIMPLE_ERROR("(A) There was no index for {}", _rep_(ahead));
          if (parentOrNil.nilp()) {
            // The parent is NIL - we are at the top, create a vertex
            size_t index = head_index.unsafe_fixnum();
            CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1))
              core::clasp_write_string(fmt::format("Adding vertex: {}\n", index));
                   );
            add_vertex(ChemInfoVertexData(index), *graph->_chemInfoGraph);
            graph->_nodeOrder.push_back(index);
          } else {
            SIMPLE_ERROR("Hit an AtomTest but it wasn't the only thing in the tree - parent should be NIL but is: {}",
                         _rep_(parentOrNil));
          }
          return true;
        } else if (parentOrNil.nilp() && gc::IsA<Logical_sp>(node)) {
          Logical_sp logical = gc::As_unsafe<Logical_sp>(node);
          CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
            core::clasp_write_string(fmt::format("Logical node id({})\n", logical->_Id));
            });
          core::T_sp thead = logical;
          if (thead.nilp())
            SIMPLE_ERROR("The head of a chain is NIL");
          ChemInfoNode_sp head = gc::As<ChemInfoNode_sp>(thead);
          AtomOrBondMatchNode_sp ahead = gc::As<AtomOrBondMatchNode_sp>(head);
          core::T_sp head_index = graph->_nodes_to_index->gethash(ahead);
          if (!head_index.fixnump())
            SIMPLE_ERROR("(B) There was no index for {}", _rep_(ahead));
          if (parentOrNil.nilp()) {
            // The parent is NIL - we are at the top, create a vertex
            size_t index = head_index.unsafe_fixnum();
            CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1))
              core::clasp_write_string(fmt::format("Adding vertex: {}\n", index));
                   );
            add_vertex(ChemInfoVertexData(index), *graph->_chemInfoGraph);
            graph->_nodeOrder.push_back(index);
          } else {
            SIMPLE_ERROR("Hit an Logical but it wasn't the only thing in the tree - parent should be NIL but is: {}",
                         _rep_(parentOrNil));
          }
          return true;
        } else if (gc::IsA<Chain_sp>(node)) {
          Chain_sp chain = gc::As_unsafe<Chain_sp>(node);
          CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
            core::clasp_write_string(fmt::format("Chain node id({})\n", chain->_Id));
            core::clasp_write_string(fmt::format("Chain node     -> {}\n", _rep_(chain->_Head).c_str()));
            });
          core::T_sp thead = chain->_Head;
          if (thead.nilp())
            SIMPLE_ERROR("The head of a chain is NIL");
          ChemInfoNode_sp head = gc::As<ChemInfoNode_sp>(thead);
          BondToAtomTest_sp bond;
          if (gc::IsA<BondToAtomTest_sp>(head)) {
            bond = gc::As_unsafe<BondToAtomTest_sp>(head);
            head = gc::As_unsafe<BondToAtomTest_sp>(head)->_AtomTest;
          }
          AtomOrBondMatchNode_sp ahead = gc::As<AtomOrBondMatchNode_sp>(head);
          CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1))
                   core::clasp_write_string(fmt::format("   setting parent_node of {} to {}\n", _rep_(chain), _rep_(ahead))););
          parent_nodes->setf_gethash(chain, ahead);
          core::T_sp head_index = graph->_nodes_to_index->gethash(ahead);
          if (!head_index.fixnump()) {
            core::clasp_write_string(fmt::format("head_index must be a fixnum - it is: {}\n", _rep_(head_index) ));
            core::clasp_write_string(fmt::format("graph->_nodes_to_index->hash_table_count() = {}\n", graph->_nodes_to_index->hashTableCount()) );
#if DEBUG_SETF_GETHASH
            core::clasp_write_string(" dump of save_pointers\n");
            for ( auto i : save_pointers ) {
              core::clasp_write_string(fmt::format("   key = {}  value = {}\n", (void*)i.first, i.second ));
            }
            core::clasp_write_string(" dump of graph->_nodes_to_index\n");
            graph->_nodes_to_index->maphash([](core::T_sp key, core::T_sp value) {
              core::clasp_write_string(fmt::format("   key = {}  value = {}\n", (void*)key.raw_(), value.unsafe_fixnum() ));
            }
              );
            core::clasp_write_string(fmt::format("  dump of hash-table -> \n{} \nAborting\n", graph->_nodes_to_index->hash_table_dump() ));
            gctools::wait_for_user_signal("could not find index");
#endif
            SIMPLE_ERROR("(C) In graph @{} hash-table @{}  - there was no index for {}@{}", (void*)graph.raw_(), (void*)graph->_nodes_to_index.raw_(), _rep_(ahead), (void*)ahead.raw_());
          }
          if (parentOrNil.nilp()) {
            // The parent is NIL - we are at the top, create a vertex
            size_t index = head_index.unsafe_fixnum();
            CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1))
                     core::clasp_write_string(fmt::format("Adding vertex: {}\n", index)););
            add_vertex(ChemInfoVertexData(index), *graph->_chemInfoGraph);
            graph->_nodeOrder.push_back(index);
          } else {
            ChemInfoNode_sp parent = gc::As_unsafe<ChemInfoNode_sp>(parentOrNil);
            core::T_sp tup = parent_nodes->gethash(parent);
            if (tup.nilp()) {
              // If the parent is not NIL but has no parent node - then we are in recursive smarts
              // and we don't want to create any vertices or do anything further
              parent_nodes->setf_gethash(chain, nil<core::T_O>());
            } else {
              // If parentOrNil is not NIL and it has a parent then keep adding to the graph
              size_t index = head_index.unsafe_fixnum();
              CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1))
                       core::clasp_write_string(fmt::format("Adding vertex: {}\n", index)););
              add_vertex(ChemInfoVertexData(index), *graph->_chemInfoGraph);
              graph->_nodeOrder.push_back(index);
              ChemInfoNode_sp up = gc::As<ChemInfoNode_sp>(tup);
              core::T_sp up_index = graph->_nodes_to_index->gethash(up);
              if (!up_index.fixnump()) {
                SIMPLE_ERROR("(D) There was no index for {}", _rep_(up));
              }
              // The head must be a BondToAtomTest_sp
              BondToAtomTest_sp bondToAtomTest = gc::As<BondToAtomTest_sp>(chain->_Head);
              int edge_index = graph->_bondNodes.size();
              graph->_bondNodes.push_back(bondToAtomTest);
              CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
                core::clasp_write_string(fmt::format("Adding edge: {} {}\n", up_index.unsafe_fixnum(), head_index.unsafe_fixnum()));
                });
              add_edge(up_index.unsafe_fixnum(), head_index.unsafe_fixnum(), EdgeProperty(edge_index), *graph->_chemInfoGraph);
            }
          }
          return true;
        } else if (gc::IsA<Branch_sp>(node)) {
          if (parentOrNil.nilp()) {
            SIMPLE_ERROR("Branches should always have a parent");
          }
          Branch_sp branch = gc::As_unsafe<Branch_sp>(node);
          CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1))
                   core::clasp_write_string(fmt::format("Branch node id({})\n", branch->_Id)););
          ChemInfoNode_sp parent = gc::As_unsafe<ChemInfoNode_sp>(parentOrNil);
          core::T_sp tup = parent_nodes->gethash(parent);
          if (tup.nilp())
            SIMPLE_ERROR("The parent of {} is NIL", _rep_(parent));
          ChemInfoNode_sp up = gc::As<ChemInfoNode_sp>(tup);
          CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1))
                   core::clasp_write_string(fmt::format("   setting parent_node of {} to {}\n", _rep_(branch), _rep_(up))););
          parent_nodes->setf_gethash(branch, up);
          return true;
        } else {
          // Do nothing
          return false;
        }
      });
  // Close the rings
  for (size_t ii = 0; ii < closers.size(); ++ii) {
    if (closers[ii]._Active) {
      for (size_t jj = 0; jj < closers[ii]._Bonds.size(); ++jj) {
        auto dummyRingAtomTest = gctools::GC<BondToAtomTest_O>::allocate(closers[ii]._Bonds[jj]._Bond);
        int edge_index = graph->_bondNodes.size();
        graph->_bondNodes.push_back(dummyRingAtomTest);
        CI_LOG(if (DEBUG_BUILD_FROM_ROOT || chem__verbose(1)) {
          core::clasp_write_string(
              fmt::format("Adding ring closing edge: {} {}\n", closers[ii]._NodeIndex, closers[ii]._Bonds[jj]._NodeIndex));
          });
        add_edge(closers[ii]._NodeIndex, closers[ii]._Bonds[jj]._NodeIndex, EdgeProperty(edge_index), *graph->_chemInfoGraph);
      }
    }
  }
}

DOCGROUP(cando);
CL_DEFUN void chem__chem_info_graph_dump(ChemInfoGraph_sp graph) {
  core::clasp_write_string(fmt::format("Number of vertices: {}\n", num_vertices(*graph->_chemInfoGraph)));
  core::clasp_write_string(fmt::format("Number of edges: {}\n", num_edges(*graph->_chemInfoGraph)));
  core::clasp_write_string(fmt::format("Length of _nodeOrder -> {}\n", graph->_nodeOrder.size()));
  typedef boost::graph_traits<ChemInfoGraphType>::edge_iterator edge_iterator;
  // Tried to make this section more clear, instead of using tie, keeping all
  // the original types so it's more clear what is going on
  std::pair<edge_iterator, edge_iterator> ei = edges(*graph->_chemInfoGraph);
  for (size_t idx = 0; idx < graph->_nodeOrder.size(); ++idx) {
    ChemInfoNode_sp a1 = graph->_atomNodes[graph->_nodeOrder[idx]];
    core::clasp_write_string(fmt::format("Vertex: {}\n", _rep_(a1)));
  }
  for (edge_iterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
    ChemInfoNode_sp a1 = graph->_atomNodes[source(*edge_iter, *graph->_chemInfoGraph)];
    ChemInfoNode_sp a2 = graph->_atomNodes[target(*edge_iter, *graph->_chemInfoGraph)];
    core::clasp_write_string(fmt::format("({} - {})\n", _rep_(a1), _rep_(a2)));
  }
}

#if 0
  int vertex_index_Molecule(const MoleculeVertexData& data)
  {
    return data._AtomIndex;
  }

  int vertex_index_ChemInfo(const ChemInfoVertexData& data)
  {
    return data._NodeIndex;
  }
#endif
struct VertexComp {
  ChemInfoGraph_sp _chemInfoGraph;
  MoleculeGraph_sp _moleculeGraph;

  bool operator()(int v1, int v2) {
    ChemInfoNode_sp node = this->_chemInfoGraph->_atomNodes[this->_chemInfoGraph->_nodeOrder[v1]];
    Atom_sp atom = gc::As<Atom_sp>(this->_moleculeGraph->_nodes->rowMajorAref(v2));
    bool match = node->matches_Atom(this->_chemInfoGraph->_Root, atom);
    CI_LOG(if (chem__verbose(1)) {
      core::clasp_write_string(fmt::format("In VertexComp matching...  v1 {} to v2 {}\n", v1, v2));
      core::clasp_write_string(fmt::format(" node -> {}\n", _rep_(node)));
      core::clasp_write_string(fmt::format(" atom -> {}\n", _rep_(atom)));
      core::clasp_write_string(fmt::format(" root -> {}\n", _rep_(this->_chemInfoGraph->_Root)));
      core::clasp_write_string(fmt::format("  match -> {}\n", match));
      });
    return match;
  }
  VertexComp(ChemInfoGraph_sp cig, MoleculeGraph_sp mg) : _chemInfoGraph(cig), _moleculeGraph(mg){};
};

struct EdgeComp {
  ChemInfoGraph_sp _chemInfoGraph;
  MoleculeGraph_sp _moleculeGraph;

  bool operator()(const boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> &eci,
                  const boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> &em) {
    boost::property_map<ChemInfoGraphType, boost::edge_index_t>::type chemInfoEdge =
        boost::get(boost::edge_index_t(), *(this->_chemInfoGraph->_chemInfoGraph));
    boost::property_map<MoleculeGraphType, boost::edge_weight_t>::type moleculeEdge =
        boost::get(boost::edge_weight_t(), *(this->_moleculeGraph->_moleculeGraph));
    BondOrder bo = moleculeEdge[em];
    int ci = chemInfoEdge[eci];
    BondToAtomTest_sp bta = this->_chemInfoGraph->_bondNodes[ci];
    // If the comparison is simple then just match the bond type to the BondEnum
    if (bta->_Bond != SABUseBondMatcher) {
      if (!aromatic_information_available_p()) {
        bool match = _matchBondTypes(bta->_Bond, bo, unbound<Atom_O>(), unbound<Atom_O>());
        CI_LOG(if (chem__verbose(1)) {
          core::clasp_write_string("In EdgeComp matching with _matchBondTypes with no aromaticity info...\n");
          core::clasp_write_string(fmt::format(" bta->_Bond -> {}\n", sabToString(bta->_Bond)));
          core::clasp_write_string(fmt::format(" bo -> {}\n", (int)bo));
          core::clasp_write_string(fmt::format("  match -> {}\n", match));
          });
        return match;
      } else {
        int msource = boost::source(em, *(this->_moleculeGraph->_moleculeGraph));
        Atom_sp asource = gc::As<Atom_sp>(this->_moleculeGraph->_nodes->rowMajorAref(msource));
        int mtarget = boost::target(em, *(this->_moleculeGraph->_moleculeGraph));
        Atom_sp atarget = gc::As<Atom_sp>(this->_moleculeGraph->_nodes->rowMajorAref(mtarget));
        bool match = _matchBondTypes(bta->_Bond, bo, asource, atarget);
        CI_LOG(if (chem__verbose(1)) {
          core::clasp_write_string("In EdgeComp matching with _matchBondTypes using aromaticity info...\n");
          core::clasp_write_string(fmt::format(" bta->_Bond -> {}\n", sabToString(bta->_Bond)));
          core::clasp_write_string(fmt::format(" bo -> {}\n", (int)bo));
          core::clasp_write_string(fmt::format("  match -> {}\n", match));
          });
        return match;
      }
    }
    // The comparison is not simple - we need to use the BondMatcher
    int msource = boost::source(em, *(this->_moleculeGraph->_moleculeGraph));
    Atom_sp asource = gc::As<Atom_sp>(this->_moleculeGraph->_nodes->rowMajorAref(msource));
    int mtarget = boost::target(em, *(this->_moleculeGraph->_moleculeGraph));
    Atom_sp atarget = gc::As<Atom_sp>(this->_moleculeGraph->_nodes->rowMajorAref(mtarget));
    Bond_sp bond = asource->getBondTo(atarget);
    BondMatcher_sp bondMatcher = gc::As<BondMatcher_sp>(bta->_BondMatcher);
    bool match = bondMatcher->matches_Bond(this->_chemInfoGraph->_Root, asource, bond);
    CI_LOG(if (chem__verbose(1)) {
      core::clasp_write_string("In EdgeComp matching with matches_Bond...\n");
      core::clasp_write_string(fmt::format(" asource -> {}\n", _rep_(asource)));
      core::clasp_write_string(fmt::format(" bond -> {}\n", _rep_(bond)));
      core::clasp_write_string(fmt::format(" root -> {}\n", _rep_(this->_chemInfoGraph->_Root)));
      core::clasp_write_string(fmt::format("  match -> {}\n", match));
      });
    return match;
  }
  EdgeComp(ChemInfoGraph_sp cig, MoleculeGraph_sp mg) : _chemInfoGraph(cig), _moleculeGraph(mg){};
};

struct MatchCallback {
  ChemInfoGraph_sp _chemInfoGraph;
  MoleculeGraph_sp _moleculeGraph;
  ChemInfoMatch_sp _currentMatch;
  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1 g) const {
    CI_LOG(if (chem__verbose(2)) {
      core::clasp_write_string("vf2 found a match\n");
      });
    core::SimpleVector_sp copy = core::SimpleVector_O::make(
        this->_currentMatch->_TagLookup->length(), nil<core::T_O>(), false, this->_currentMatch->_TagLookup->length(),
        (core::T_sp *)(this->_currentMatch->_TagLookup->rowMajorAddressOfElement_(0)));
    this->_currentMatch->_TagHistory = core::Cons_O::create(copy, this->_currentMatch->_TagHistory);
    // printf("%s:%d In match\n", __FILE__, __LINE__);
    //  printf("    -> %s\n", _rep_(this->_currentMatch->_TagLookup).c_str());
    return true; // chem::_matchBondTypes(chemInfoEdge._BondEnum,moleculeEdge._BondOrder);
  }
  MatchCallback(ChemInfoGraph_sp cig, MoleculeGraph_sp mg, ChemInfoMatch_sp cm)
      : _chemInfoGraph(cig), _moleculeGraph(mg), _currentMatch(cm){};
};

#if 0
// Leading example? https://git.skewed.de/count0/graph-tool/blob/master/src/graph/topology/graph_subgraph_isomorphism.cc
// http://sophie.aero.jussieu.fr/distrib/Fedora/20/i386/media/updates/by-pkgid/875ec90536d12b055fdbc9a697ef0042/files/942



  template <typename T, typename U>
    bool user_callback(T t, U u) {
    printf("%s:%d in user_callback\n");
  }

#endif

DOCGROUP(cando);
CL_DEFUN core::List_sp chem__boost_graph_vf2(ChemInfoGraph_sp chemInfoGraph, MoleculeGraph_sp moleculeGraph) {
  if (boost::num_vertices(*chemInfoGraph->_chemInfoGraph) > boost::num_vertices(*moleculeGraph->_moleculeGraph)) {
    return nil<core::T_O>();
  }
  core::HashTable_sp ringHashTable = core::HashTable_O::createEql();
  ChemInfoMatch_sp current_match = ChemInfoMatch_O::make(chemInfoGraph->_Root, chemInfoGraph->_Root->_MaxTag, ringHashTable);
  core::DynamicScopeManager scope(_sym_STARcurrent_matchSTAR, current_match);
  EdgeComp edge_comp(chemInfoGraph, moleculeGraph);
  VertexComp vertex_comp(chemInfoGraph, moleculeGraph);
  boost::vf2_print_callback<ChemInfoGraphType, MoleculeGraphType> callback(*chemInfoGraph->_chemInfoGraph,
                                                                           *moleculeGraph->_moleculeGraph);
  MatchCallback matchCallback(chemInfoGraph, moleculeGraph, current_match);
  CI_LOG(if (chem__verbose(1)) {
    core::clasp_write_string("About to run boost::vf2_subgraph_iso\n");
    });
  boost::vf2_subgraph_iso(*chemInfoGraph->_chemInfoGraph, *moleculeGraph->_moleculeGraph, matchCallback, chemInfoGraph->_nodeOrder,
                          boost::edges_equivalent(edge_comp).vertices_equivalent(vertex_comp));
  return current_match->_TagHistory;
}
//#endif

struct CommonVertexComp {
  MoleculeGraph_sp _moleculeGraph1;
  MoleculeGraph_sp _moleculeGraph2;
  core::T_sp _matcher;
  bool operator()(int v1, int v2) {
    core::T_sp node1 = this->_moleculeGraph1->_nodes->rowMajorAref(v1);
    core::T_sp node2 = this->_moleculeGraph2->_nodes->rowMajorAref(v2);
    core::T_sp matched = core::eval::funcall(this->_matcher, node1, node2);
    bool match = matched.notnilp();
    return match;
  }
  CommonVertexComp(MoleculeGraph_sp mg1, MoleculeGraph_sp mg2, core::T_sp matcher)
      : _moleculeGraph1(mg1), _moleculeGraph2(mg2), _matcher(matcher){};
};

struct CommonEdgeComp {
  MoleculeGraph_sp _moleculeGraph1;
  MoleculeGraph_sp _moleculeGraph2;

  bool operator()(const boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> &em1,
                  const boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> &em2) {
    boost::property_map<MoleculeGraphType, boost::edge_weight_t>::type moleculeEdge1 =
        boost::get(boost::edge_weight_t(), *(this->_moleculeGraph1->_moleculeGraph));
    boost::property_map<MoleculeGraphType, boost::edge_weight_t>::type moleculeEdge2 =
        boost::get(boost::edge_weight_t(), *(this->_moleculeGraph2->_moleculeGraph));
    BondOrder bo1 = moleculeEdge1[em1];
    BondOrder bo2 = moleculeEdge2[em2];
    bool match = (bo1 == bo2);
    CI_LOG(if (chem__verbose(1)) {
      core::clasp_write_string("In EdgeComp matching with _matchBondTypes with no aromaticity info...\n");
      core::clasp_write_string(fmt::format(" bo1 -> {}\n", (int)bo1));
      core::clasp_write_string(fmt::format(" bo2 -> {}\n", (int)bo2));
      core::clasp_write_string(fmt::format("  match -> {}\n", match));
      });
    return match;
  }
  CommonEdgeComp(MoleculeGraph_sp mg1, MoleculeGraph_sp mg2) : _moleculeGraph1(mg1), _moleculeGraph2(mg2){};
};

struct print_callback {
public:
  const MoleculeGraph_sp m_graph1;
  const MoleculeGraph_sp m_graph2;
  core::T_sp m_callback;
  core::ComplexVector_T_sp _results;
  size_t _max_callbacks;
  size_t *_count_callbacks;

public:
  print_callback(MoleculeGraph_sp graph1, MoleculeGraph_sp graph2, core::T_sp callback, size_t max_callbacks)
      : m_graph1(graph1), m_graph2(graph2), m_callback(callback), _max_callbacks(max_callbacks), _count_callbacks(NULL) {
    this->_results = core::ComplexVector_T_O::make(32, nil<core::T_O>(), core::make_fixnum(0));
  }

  template <typename CorrespondenceMapFirstToSecond, typename CorrespondenceMapSecondToFirst>
  bool operator()(CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
                  CorrespondenceMapSecondToFirst correspondence_map_2_to_1,
                  typename boost::graph_traits<MoleculeGraphType>::vertices_size_type subgraph_size) {

    if (subgraph_size > this->_results->fillPointer()) {
      // Only ever pass the longest solution found so far.
      this->_results->fillPointerSet(0);
      // Print out correspondences between vertices
      BGL_FORALL_VERTICES_T(vertex1, *m_graph1->_moleculeGraph, MoleculeGraphType) {
        // Skip unmapped vertices
        if (get(correspondence_map_1_to_2, vertex1) != boost::graph_traits<MoleculeGraphType>::null_vertex()) {
          this->_results->vectorPushExtend(core::make_fixnum(vertex1));
          this->_results->vectorPushExtend(core::make_fixnum(get(correspondence_map_1_to_2, vertex1)));
        }
      }
      [[maybe_unused]]core::T_sp bres = core::eval::funcall(this->m_callback, this->_results, core::make_fixnum((*this->_count_callbacks)),
                                            this->m_graph1, this->m_graph2);
    }
    if (this->_count_callbacks != NULL) {
      (*this->_count_callbacks)++;
      return ((*this->_count_callbacks) < this->_max_callbacks);
    }
    return true;
  }
};

CL_DOCSTRING(R"dx(A wrapper for boost::graph::mcgregor_common_subgraphs)dx");
CL_DOCSTRING_LONG(R"dx(The **maximum-callbacks** argument sets the maximum number
of times the mcgregor algorithm can call its internal match callback function (which calls **match-callback**) before this function returns.
The **match-callback** will only be invoked if the match is longer than any previous result.
The **vertex-matcher** argument is a lambda that takes two atoms and returns T if they are considered equivalent.
The **match-callback** argument is a lambda that takes three arguments (RESULTS GRAPH1 GRAPH2) and its return value is ignored.
Edges are matched using bond orders.)dx")
CL_LAMBDA(molecule-graph1 molecule-graph2 only-connected-subgraphs vertex-matcher match-callback &optional (maximum-callbacks 1000));
DOCGROUP(cando);
CL_DEFUN core::List_sp chem__boost_graph_mcgregor_common_subgraphs(MoleculeGraph_sp moleculeGraph1, MoleculeGraph_sp moleculeGraph2,
                                                                   bool only_connected_subgraphs, core::T_sp vertex_matcher,
                                                                   core::T_sp match_callback, size_t maximum_callbacks) {
  CommonEdgeComp edge_comp(moleculeGraph1, moleculeGraph2);
  CommonVertexComp vertex_comp(moleculeGraph1, moleculeGraph2, vertex_matcher);
  print_callback my_callback(moleculeGraph1, moleculeGraph2, match_callback, maximum_callbacks);
  my_callback._count_callbacks = new size_t(0);
  boost::mcgregor_common_subgraphs(*moleculeGraph1->_moleculeGraph, *moleculeGraph2->_moleculeGraph, only_connected_subgraphs,
                                   my_callback, boost::edges_equivalent(edge_comp).vertices_equivalent(vertex_comp));
  delete my_callback._count_callbacks;
  return nil<core::T_O>();
}

void Smirks_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, reactant), this->_reactant);
  node->field(INTERN_(kw, product), this->_product);
}

#if 0
DOCGROUP(cando);
CL_DEFUN core::List_sp chem__boost_graph_mcgregor_common_subgraphs_maximum(MoleculeGraph_sp moleculeGraph1, MoleculeGraph_sp moleculeGraph2,
                                                                   bool only_connected_subgraphs) {
  CommonEdgeComp edge_comp(moleculeGraph1,moleculeGraph2);
  CommonVertexComp vertex_comp(moleculeGraph1,moleculeGraph2);
  print_callback<MoleculeGraphType,MoleculeGraphType> my_callback(*moleculeGraph1->_moleculeGraph,*moleculeGraph2->_moleculeGraph);
  boost::mcgregor_common_subgraphs_maximum(
                          *moleculeGraph1->_moleculeGraph,
                          *moleculeGraph2->_moleculeGraph,
                          only_connected_subgraphs,
                          my_callback,
                          boost::edges_equivalent(edge_comp).vertices_equivalent(vertex_comp));
  return _Nil<core::T_O>();
}

DOCGROUP(cando);
CL_DEFUN core::List_sp chem__boost_graph_mcgregor_common_subgraphs_maximum_unique(MoleculeGraph_sp moleculeGraph1, MoleculeGraph_sp moleculeGraph2,
                                                                   bool only_connected_subgraphs) {
  CommonEdgeComp edge_comp(moleculeGraph1,moleculeGraph2);
  CommonVertexComp vertex_comp(moleculeGraph1,moleculeGraph2);
  print_callback<MoleculeGraphType,MoleculeGraphType> my_callback(*moleculeGraph1->_moleculeGraph,*moleculeGraph2->_moleculeGraph);
  boost::mcgregor_common_subgraphs_maximum_unique(
                          *moleculeGraph1->_moleculeGraph,
                          *moleculeGraph2->_moleculeGraph,
                          only_connected_subgraphs,
                          my_callback,
                          boost::edges_equivalent(edge_comp).vertices_equivalent(vertex_comp));
  return _Nil<core::T_O>();
}
#endif

//#endif

}; // namespace chem
