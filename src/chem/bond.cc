/*
    File: bond.cc
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
#define	DEBUG_LEVEL_NONE
//
// (C) 2004 Christian E. Schafmeister
//



// Comment
#include <clasp/core/common.h>
#include <clasp/core/exceptions.h>
#include <clasp/core/lispList.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/constitutionAtoms.h>
#include <clasp/core/wrappers.h>



namespace chem 
{
SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_bondOrderToSymbolConverter_PLUS_);
SYMBOL_EXPORT_SC_(ChemKwPkg,_singleBond);
SYMBOL_EXPORT_SC_(ChemKwPkg,_doubleBond);
SYMBOL_EXPORT_SC_(ChemKwPkg,_tripleBond);
SYMBOL_EXPORT_SC_(ChemKwPkg,_aromaticBond);
SYMBOL_EXPORT_SC_(ChemKwPkg,_hydrogenBond);
SYMBOL_EXPORT_SC_(ChemKwPkg,_virtualBond);
SYMBOL_EXPORT_SC_(ChemKwPkg,_dashedSingleBond);
SYMBOL_EXPORT_SC_(ChemKwPkg,_dashedDoubleBond);
SYMBOL_EXPORT_SC_(ChemKwPkg,singleWedgeBegin);
SYMBOL_EXPORT_SC_(ChemKwPkg,singleHashBegin);
SYMBOL_EXPORT_SC_(ChemKwPkg,singleWedgeEnd);
SYMBOL_EXPORT_SC_(ChemKwPkg,singleHashEnd);


string	XmlName_Bond = "bond";

DOCGROUP(cando);
CL_DEFUN bool chem__single_bond_p(chem::BondOrder o)
{
  return Bond_O::singleBondP(o);
}


Bond_sp	Bond_O::create(Atom_sp from, Atom_sp to, BondOrder o)
{
  auto  bond  = gctools::GC<Bond_O>::allocate_with_default_constructor();
  bond->_Atom1 = from;
  bond->_Atom2 = to;
  bond->setOrder(o);
#if 0
  if (from->getName()->symbolNameAsString() == "C8" &&
      to->getName()->symbolNameAsString() == "C7") {
    printf("%s:%d Created bond %s -> %s %s\n", __FILE__, __LINE__, _rep_(from).c_str(), _rep_(to).c_str(), bondOrderToString(o).c_str());
  }
#endif
  LOG("created bond from={} to={}"
      , bond->_Atom1->description()
      , bond->_Atom2->description() );
  return bond;
}


void Bond_O::initialize()
{
  this->Base::initialize();
  this->_DirectionalOrder = singleBond;
  LOG("Setting fromAtom and toAtom to nil" ); //
  this->_Atom1 = nil<Atom_O>();
  this->_Atom2 = nil<Atom_O>();
}

Bond_O::Bond_O(const Bond_O& bb)  : core::CxxObject_O(bb)
{
  LOG("Copying Bond_O" ); //
  this->_DirectionalOrder = bb._DirectionalOrder;
  this->_Atom1 = bb._Atom1;
  this->_Atom2 = bb._Atom2;
//  printf("%s:%d  Bond_O copy ctor %s -> %s %s\n", __FILE__, __LINE__, _rep_(bb._Atom1).c_str(), _rep_(bb._Atom2).c_str(), bondOrderToString(bb.order).c_str());
  this->_Properties = bb._Properties; // You can't call allocators from ctors core::cl__copy_list(bb._Properties);
  LOG("copy _Atom1={} _Atom2={}"
      , this->_Atom1->description()
      , this->_Atom2->description() );
}


bool Bond_O::isAtom1(Atom_sp a) const
{
  ASSERTNOTNULL(this->_Atom1);
  if (this->_Atom1==a) return true;
  ASSERTNOTNULL(this->_Atom2);
  if ( this->_Atom2==a) return false;
  SIMPLE_ERROR("Atom[{}] is neither atom1[{}] or atom2[{}] of bond" , _rep_(a) , _rep_(this->_Atom1) , _rep_(this->_Atom2) );
}




Bond_sp Bond_O::copyDontRedirectAtoms()
{
  auto  bondNew  = gctools::GC<Bond_O>::copy( *this); // = RP_Copy<Bond_O>(this);
  return bondNew;
}

void Bond_O::addYourselfToCopiedAtoms()
{
//  printf("%s:%d Redirecting bond from atoms %s@%p - %s@%p\n", __FILE__, __LINE__, _rep_(this->_Atom1).c_str(), this->_Atom1.raw_(), _rep_(this->_Atom2).c_str(), this->_Atom2.raw_() );
//  printf("%s:%d                    to atoms %s@%p - %s@%p\n", __FILE__, __LINE__, _rep_(this->_Atom1->getCopyAtom()).c_str(), this->_Atom1->getCopyAtom().raw_(), _rep_(this->_Atom2->getCopyAtom()).c_str(), this->_Atom2->getCopyAtom().raw_());
  ASSERTNOTNULL(this->_Atom1);
  this->_Atom1 = this->_Atom1->getCopyAtom();
  this->_Atom1->addBond(this->asSmartPtr());
  ASSERTNOTNULL(this->_Atom2);
  this->_Atom2 = this->_Atom2->getCopyAtom();
  ASSERTNOTNULL(this->_Atom2);
  this->_Atom2->addBond(this->asSmartPtr());
}




//
//	isInterResidueBond
//
//	Return true if this is an inter-residue bond
//
bool	Bond_O::isInterResidueBond(core::HashTable_sp atomToResidue)
{
  core::T_sp res1 = atomToResidue->gethash(this->_Atom1);
  core::T_sp res2 = atomToResidue->gethash(this->_Atom2);
  return (res1!=res2);
}

core::NullTerminatedEnumAssociation bondOrderKeys[] = {
    { "noBond", noBond },
    { "single", singleBond },
    { "double", doubleBond },
    { "triple", tripleBond },
    { "aromatic", aromaticBond },
    { "hbond", hydrogenBond },
    { "virtual", virtualBond },
    { "dashedSingle", dashedSingleBond },
    { "dashedDouble", dashedDoubleBond },
    { "singleWedgeBegin", singleWedgeBegin },
    { "singleHashBegin", singleHashBegin },
    { "singleWedgeEnd", singleWedgeEnd },
    { "singleHashEnd", singleHashEnd },
    { "", -1 }
};


/*!
 * Saving: The Archive getCurrentNode() will be a Node that
 * 	we will fill the node with attributes and children.
 * 	If it is then we just create a node that links to the existing one.
 * Loading: The Archive _CurrentNode will be pointing at the deep node for this object.
 * 	This function will extract the data.
 *
 */
void Bond_O::fields(core::Record_sp node)
{
  LOG("archive direction = {}" , (node->loading()?"loading":"saving") ); //
  node->field(INTERN_(kw,order), this->_DirectionalOrder ); // attributeSymbolEnumHiddenConverter("order",this->order,_sym__PLUS_bondOrderToSymbolConverter_PLUS_);
  node->field_if_not_nil(INTERN_(kw,properties),this->_Properties);
  node->field( INTERN_(kw,a1), this->_Atom1 );
  node->field( INTERN_(kw,a2), this->_Atom2 );
};

CL_LISPIFY_NAME(Bond/properties);
CL_DEFMETHOD
core::List_sp Bond_O::properties() const {
  return this->_Properties;
}

CL_LISPIFY_NAME(Bond/setProperties);
CL_DEFMETHOD
void Bond_O::setProperties(core::List_sp p) {
  this->_Properties = p;
}

CL_LISPIFY_NAME(Bond/clearProperty);
CL_DEFMETHOD
void Bond_O::clearProperty(core::Symbol_sp prop)
{
  this->_Properties = core::core__rem_f(this->_Properties,prop);
}

CL_DEFMETHOD BondOrder Bond_O::getOrderFromAtom(Atom_sp firstAtom) {
  return this->maybeFlipOrder(firstAtom);
}


CL_LISPIFY_NAME(Bond/setProperty);
CL_DEFMETHOD
void Bond_O::setProperty(core::Symbol_sp prop, core::T_sp val)
    {
      this->_Properties = core::core__put_f(this->_Properties,val,prop);
    }


CL_LISPIFY_NAME(Bond/getProperty);
CL_DEFMETHOD core::T_sp Bond_O::getProperty(core::Symbol_sp prop, core::T_sp defval)
{
  core::T_sp res = core::cl__getf(this->_Properties,prop,unbound<core::T_O>());
  if (res.unboundp()) {
    return defval;
  }
  return res;
}

CL_LISPIFY_NAME("Bond/hasProperty");
CL_DEFMETHOD bool Bond_O::hasProperty(core::Symbol_sp prop)
{
  return !core::cl__getf(this->_Properties,prop,unbound<core::T_O>()).unboundp();
}


CL_DEFUN float chem__bond_length(Bond_sp bond) {
  Vector3 p1 = bond->_Atom1->getPosition();
  Vector3 p2 = bond->_Atom2->getPosition();
  Vector3 diff = p1-p2;
  return diff.length();
}


CL_LISPIFY_NAME("Bond/getOtherAtom");
CL_DEFMETHOD Atom_sp Bond_O::getOtherAtom(Atom_sp atom) const
{
  if ( atom == this->_Atom1 ) return this->_Atom2;
  return this->_Atom1;
}


bool Bond_O::isBondBetween(Atom_sp a, Atom_sp b) const
{
  if ( this->_Atom1 == a )
  {
    if ( this->_Atom2 == b ) return true;
  }
  if ( this->_Atom2 == a )
  {
    if ( this->_Atom1 == b ) return true;
  }
  return false;
}


string  Bond_O::description() const
{
  stringstream    ss;
  ss << "(Bond ";
  Atom_sp wa1 = this->_Atom1;
  ss << _rep_(wa1);
  ss << bondOrderToChar(this->_DirectionalOrder);
  Atom_sp wa2 = this->_Atom2;
  ss << _rep_(wa2);
  ss << " " << bondOrderToString(this->_DirectionalOrder);
  ss << " @"<<std::hex<<this<<std::dec<<")";
  return ss.str();
}

string  Bond_O::describeOther(Atom_sp from) const
{
  ASSERTF(from==this->_Atom1||from==this->_Atom2,"describeFrom failed because from atom[{}] is not part of bond: {}" , from->description() , this->description() );
  stringstream    ss;
  ss << "Bond(";
  ss << from->description();
  ss << bondOrderToChar(this->_DirectionalOrder);
  Atom_sp a2 = this->getOtherAtom(from);
  ss << a2->description();
  ss << ")";
  return ss.str();
}



bool	Bond_O::invalid(Atom_sp a)
{
  ASSERTNOTNULL(this->_Atom1);
  ASSERTNOTNULL(this->_Atom2);
  if ( this->_Atom1.nilp() ) return true;
  if ( this->_Atom2.nilp() ) return true;
  return false;
}

void	Bond_O::failIfInvalid(Atom_sp a)
{
  if ( this->invalid(a) ) {
    SIMPLE_ERROR("INVALID {}" , this->description());
  }
}



void Bond_O::redirectToAtomCopies()
{
  LOG("Redirecting bond@{}" , this ); //
  Atom_sp fa =  this->_Atom1;
  Atom_sp ta = this->_Atom2;
#ifdef	DEBUG_ON
  if ( !fa.objectp() ) {
    SIMPLE_ERROR("redirectToAtomCopies _Atom1 is NULL");
  }
  LOG("  original from atom@{}" , fa ); //
  if ( !ta.objectp() ) {
    SIMPLE_ERROR("redirectToAtomCopies _Atom2 is NULL");
  }
  LOG("  original   to atom@{}" , ta ); //
#endif
  Atom_sp fc = fa->getCopyAtom();
  ASSERTNOTNULL(fc);
  LOG("    new from {}" , fc->description() );
  Atom_sp tc = ta->getCopyAtom();
  ASSERTNOTNULL(tc);
  LOG("    new   to {}" , tc->description() );
  this->_Atom1 = fc;
  this->_Atom2 = tc;
}

void	Bond_O::imposeYourself()
{
  Atom_sp a1 = this->_Atom1;
  Atom_sp a2 = this->_Atom2;
  a1->_addExistingBond(this->sharedThis<Bond_O>());
  a2->_addExistingBond(this->sharedThis<Bond_O>());
}

/*! If the from atom matches this->_Atom1 then return the order,
otherwise if it's a directional stereochemical bond then flip
the direction */
BondOrder Bond_O::maybeFlipOrder(Atom_sp from)
{
  if (from == this->_Atom1) return this->_DirectionalOrder;
  switch (this->_DirectionalOrder) {
  case singleWedgeBegin: return singleWedgeEnd;
  case singleWedgeEnd: return singleWedgeBegin;
  case singleHashBegin: return singleHashEnd;
  case singleHashEnd: return singleHashBegin;
  default:
      return this->_DirectionalOrder;
  }
}

void Bond_O::canonicalizeBondOrder(Atom_sp& a1, Atom_sp& a2, BondOrder& order) {
  switch (order) {
  case singleWedgeEnd: {
    order = singleWedgeBegin;
    Atom_sp temp = a2;
    a2 = a1;
    a1 = temp;
    return;
  }
      break;
  case singleHashEnd: {
    order = singleHashBegin;
    Atom_sp temp = a2;
    a2 = a1;
    a1 = temp;
    return;
  }
      break;
  default: return;
  }
}

#ifdef USE_TOPOLOGY
ConstitutionBond_sp Bond_O::asConstitutionBond(Atom_sp from, const MapAtomsToConstitutionAtomIndex0N& atomMap)
{
  Atom_sp to = this->getOtherAtom(from);
  MapAtomsToConstitutionAtomIndex0N::const_iterator it = atomMap.find(to);
  if ( it == atomMap.end() )
  {
    SIMPLE_ERROR("Could not find atom[{}] in atomMap" , _rep_(to));
  }
  ConstitutionAtomIndex0N index = it->second;
  BondOrder order = this->maybeFlipOrder(from);
  ConstitutionBond_sp cb = makeConstitutionBond(index, order);
  return cb;
}
#endif




CL_LISPIFY_NAME("getOrderFromAtomAsString");
CL_DEFMETHOD string	Bond_O::getOrderFromAtomAsString(Atom_sp fromAtom)
{
  return bondOrderToString(this->getOrderFromAtom(fromAtom));
}

void	Bond_O::joinYourAtoms()
{
  IMPLEMENT_ME(); // handle new way of storing bonds
#if 0
  ASSERTNOTNULLP(this->fromAtom,"Bond_O::joinYourAtoms from atom is undefined!");
  ASSERTNOTNULLP(this->toAtom,"Bond_O::joinYourAtoms to atom is undefined!");
  this->fromAtom->bondTo(this->toAtom,this->order);
#endif
}




//typedef	enum 	BondOrderEnum {
//		noBond=NoBond,
//		singleBond=SingleBond,
//    		doubleBond = DoubleBond,
//		tripleBond = TripleBond,
//		aromaticBond = AromaticBond } BondOrder;

string	bondOrderToString(BondOrder bo) {

  switch (bo) {
  case noBond:
      return "none";
      break;
  case singleBond:
      return "single";
      break;
  case singleWedgeBegin:
      return "singleWedgeBegin";
      break;
  case singleHashBegin:
      return "singleHashBegin";
      break;
  case singleWedgeEnd:
      return "singleWedgeEnd";
      break;
  case singleHashEnd:
      return "singleHashEnd";
      break;
  case doubleBond:
      return "double";
      break;
  case tripleBond:
      return "triple";
      break;
  case aromaticBond:
      return "aromatic";
      break;
  case hydrogenBond:
      return "HBond";
      break;
  default:
      return "-unknown-";
  }
  printf( "Illegal bond order: %d\n", (int)(bo));
  exit(1);
}



BondOrder	stringToBondOrder(string bos) {
  if ( bos=="none" ) return noBond;
  if ( bos=="single" ) return singleBond;
  if ( bos=="singleWedgeBegin" ) return singleWedgeBegin;
  if ( bos=="singleHashBegin" ) return singleHashBegin;
  if ( bos=="singleWedgeEnd" ) return singleWedgeEnd;
  if ( bos=="singleHashEnd" ) return singleHashEnd;
  if ( bos=="double" ) return doubleBond;
  if ( bos=="triple" ) return tripleBond;
  if ( bos=="aromatic" ) return aromaticBond;
  if ( bos=="HBond" ) return hydrogenBond;
  printf( "Illegal bond order: %s\n", bos.c_str() );
  exit(1);
}


char bondOrderToChar(BondOrder bo) 
{
  switch (bo) {
  case noBond:
      return '|';
      break;
  case singleBond:
      return '-';
      break;
  case singleWedgeBegin:
      return '/';
      break;
  case singleHashBegin:
      return '\\';
      break;
  case singleWedgeEnd:
      return '\'';
      break;
  case singleHashEnd:
      return '"';
      break;
  case doubleBond:
      return '=';
      break;
  case tripleBond:
      return '#';
      break;
  case aromaticBond:
      return '~';
      break;
  case hydrogenBond:
      return '&';
      break;
  case virtualBond:
      return '^';
      break;
  case dashedSingleBond:
      return '!';
      break;
  case dashedDoubleBond:
      return '%';
      break;
  case unknownOrderBond:
      return '?';
      break;
  default:
      return '/';
      break;
  }
  printf( "Illegal bond order: %d\n", (int)(bo));
  exit(1);
}


BondList_O::BondList_O(const BondList_O& orig) : Base(orig), _Bonds(orig._Bonds)
{
}



void	BondList_O::initialize()
{
  this->Base::initialize();
  this->_Bonds.clear();
};




void	BondList_O::imposeYourself()
{ 
  gctools::Vec0<Bond_sp>::iterator	bi;
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
    (*bi)->imposeYourself();
  }
}



void	BondList_O::removeBond(Bond_sp b)
{ 
  BondList_O::iterator	bi;
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
    if ( (*bi) == b ) {
      this->_Bonds.erase(bi);
      return;
    }
  }
	// If the bond wasn't found then just return
}

void	BondList_O::removeBondBetween(Atom_sp a,Atom_sp b)
{
  BondList_O::iterator	bi;
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    if ( (*bi)->isBondBetween(a,b) )
    {
      this->_Bonds.erase(bi);
      return;
    }
  }
  SIMPLE_ERROR("Could not find bond between {} and {} in bondList: {}" , _rep_(a) , _rep_(b) , _rep_(this->asSmartPtr()));
	// If the bond wasn't found then just return
}


void	BondList_O::addBond(Bond_sp b)
{
  this->_Bonds.push_back(b);
};

bool Bond_O::equalp(core::T_sp b) const {
  if ( Bond_sp bo = b.asOrNull<Bond_O>() ) {
    if ( (this->_Atom1 == bo->_Atom1) && (this->_Atom2 == bo->_Atom2) ) return true;
    if ( (this->_Atom1 == bo->_Atom2) && (this->_Atom2 == bo->_Atom1) ) return true;
  }
  return false;
}

/*! deep(ish) copy of a bond list.
Create a new BondList_O and create a new copy of the Vec0 in the BondList_O but
use the same Bond_O's */
BondList_sp	BondList_O::deepishCopy() const
{
  BondList_O::const_iterator bi;
  auto  bl  = gctools::GC<BondList_O>::allocate_with_default_constructor();
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) bl->addBond(*bi);
  return bl;
}

string	BondList_O::description() const
{
  BondList_O::const_iterator bi;
  stringstream ss;
  ss << "(BondList numBonds=" << this->_Bonds.size() << std::endl;
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
    ss << (*bi)->description() << " " << std::endl;
  }
  ss << " )" << std::endl;
  return ss.str();
}

string	BondList_O::describeOthers(Atom_sp from) const
{
  BondList_O::const_iterator bi;
  stringstream ss;
  ss << "BondList[[[" << std::endl;
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
    ss << (*bi)->describeOther(from) << "," << std::endl;
  }
  ss << "]]]" << std::endl;
  return ss.str();
}

void BondList_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,bonds),this->_Bonds);
}

#ifdef XML_ARCHIVE
void	BondList_O::archiveBase(core::ArchiveP node)
{
  BondList_O::iterator	ib;
  if ( node->saving() ) {
    LOG("About to save BondList with {} members" , this->size()  ); //
    for ( ib=this->_Bonds.begin(); ib!=this->_Bonds.end(); ib++ ) {
      
      node->archiveSaveObjectAsChildAssignAutoUniqueId<Bond_O>(*ib);
      LOG("After archiveObject just print" ); //
    }
  } else {
    LOG("About to load BondList" ); //
    this->_Bonds.clear();
    core::VectorNodes::iterator	ic;
    Bond_sp			b,bEmpty;
    for ( ic=node->begin_Children(); ic!=node->end_Children(); ic++ ) {
      (*ic)->archiveLoadObjectDirectly<Bond_O>(b);
      ASSERTNOTNULL(b);
      this->_Bonds.push_back(b);
    }
  }
}
#endif


  SYMBOL_EXPORT_SC_(ChemKwPkg,noBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,singleBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,singleWedgeBegin);
  SYMBOL_EXPORT_SC_(ChemKwPkg,singleHashBegin);
  SYMBOL_EXPORT_SC_(ChemKwPkg,singleWedgeEnd);
  SYMBOL_EXPORT_SC_(ChemKwPkg,singleHashEnd);
  SYMBOL_EXPORT_SC_(ChemKwPkg,doubleBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,tripleBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,aromaticBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,hydrogenBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,virtualBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,dashedSingleBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,dashedDoubleBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,unknownOrderBond);
  CL_BEGIN_ENUM(BondOrder,_sym__PLUS_bondOrderToSymbolConverter_PLUS_,"BondOrder");
  CL_VALUE_ENUM(chemkw::_sym_noBond, noBond );
  CL_VALUE_ENUM(chemkw::_sym_singleBond, singleBond );
  CL_VALUE_ENUM(chemkw::_sym_singleWedgeBegin, singleWedgeBegin );
  CL_VALUE_ENUM(chemkw::_sym_singleHashBegin, singleHashBegin );
  CL_VALUE_ENUM(chemkw::_sym_singleWedgeEnd, singleWedgeEnd );
  CL_VALUE_ENUM(chemkw::_sym_singleHashEnd, singleHashEnd );
  CL_VALUE_ENUM(chemkw::_sym_doubleBond, doubleBond );
  CL_VALUE_ENUM(chemkw::_sym_tripleBond, tripleBond );
  CL_VALUE_ENUM(chemkw::_sym_aromaticBond, aromaticBond );
  CL_VALUE_ENUM(chemkw::_sym_hydrogenBond, hydrogenBond );
  CL_VALUE_ENUM(chemkw::_sym_virtualBond, virtualBond);
  CL_VALUE_ENUM(chemkw::_sym_dashedSingleBond,dashedSingleBond);
  CL_VALUE_ENUM(chemkw::_sym_dashedDoubleBond,dashedDoubleBond);
  CL_VALUE_ENUM(chemkw::_sym_unknownOrderBond,unknownOrderBond);
  CL_END_ENUM(_sym__PLUS_bondOrderToSymbolConverter_PLUS_);













}; // namespace chem
