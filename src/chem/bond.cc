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


string	XmlName_Bond = "bond";


Bond_sp	Bond_O::create(Atom_sp from, Atom_sp to, BondOrder o)
{
  GC_ALLOCATE(Bond_O, bond );
  bond->_Atom1 = from;
  bond->_Atom2 = to;
  bond->setOrder(o);
  LOG(BF("created bond from=%s to=%s")
      % bond->_Atom1->description()
      % bond->_Atom2->description() );
  return bond;
}


void Bond_O::initialize()
{_OF();
  this->Base::initialize();
  this->order = singleBond;
  LOG(BF("Setting fromAtom and toAtom to nil") ); //
  this->_Atom1 = _Nil<Atom_O>();
  this->_Atom2 = _Nil<Atom_O>();
}

Bond_O::Bond_O(const Bond_O& bb)  : core::CxxObject_O(bb)
{
  LOG(BF("Copying Bond_O") ); //
  this->order = bb.order;
  this->_Atom1 = bb._Atom1;
  this->_Atom2 = bb._Atom2;
  this->_Properties = core::cl__copy_list(bb._Properties);
  LOG(BF("copy _Atom1=%s _Atom2=%s")
      % this->_Atom1->description()
      % this->_Atom2->description() );
}


bool Bond_O::isAtom1(Atom_sp a) const
{_OF();
  ASSERTNOTNULL(this->_Atom1);
  if (this->_Atom1==a) return true;
  ASSERTNOTNULL(this->_Atom2);
  if ( this->_Atom2==a) return false;
  SIMPLE_ERROR(BF("Atom[%s] is neither atom1[%s] or atom2[%s] of bond")
               % _rep_(a) % _rep_(this->_Atom1) % _rep_(this->_Atom2) );
}




Bond_sp Bond_O::copyDontRedirectAtoms()
{_OF();
  GC_COPY(Bond_O, bondNew , *this); // = RP_Copy<Bond_O>(this);
  return bondNew;
}

void Bond_O::addYourselfToCopiedAtoms()
{_OF();
//  printf("%s:%d Redirecting bond from atoms %p - %p\n", __FILE__, __LINE__, this->_Atom1.raw_(), this->_Atom2.raw_());
//  printf("%s:%d                    to atoms %p - %p\n", __FILE__, __LINE__, this->_Atom1->getCopyAtom().raw_(), this->_Atom2->getCopyAtom().raw_());
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
bool	Bond_O::isInterResidueBond()
{
  Matter_wp	wfromCont, wtoCont;
  Matter_sp	fromCont, toCont;
  Atom_sp a1 = this->_Atom1;
  ASSERTP(a1.objectp(), "Bond_O::isInterResidueBond from atom is undefined!");
  LOG(BF("from atom = %s") % a1->description() );
  wfromCont = a1->getResidueContainedBy();
  ASSERTNOTNULL(wfromCont);
  if ( wfromCont.nilp() ) SIMPLE_ERROR(BF("From atom isnt in a residue"));
  Atom_sp a2 = this->_Atom2;
  ASSERTNOTNULLP(a2, "Bond_O::isInterResidueBond to atom is undefined!");
  wtoCont = a2->getResidueContainedBy();
  ASSERTNOTNULL(wtoCont);
  if ( wtoCont.nilp() ) SIMPLE_ERROR(BF("To atom isnt in a residue"));
  fromCont = wfromCont;
  toCont = wtoCont;
  return ( fromCont != toCont );
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
  LOG(BF("archive direction = %s") % (node->loading()?"loading":"saving") ); //
  node->field(INTERN_(kw,order), this->order ); // attributeSymbolEnumHiddenConverter("order",this->order,_sym__PLUS_bondOrderToSymbolConverter_PLUS_);
  node->field_if_not_nil(INTERN_(kw,properties),this->_Properties);
  node->field( INTERN_(kw,a1), this->_Atom1 );
  node->field( INTERN_(kw,a2), this->_Atom2 );
};

string Bond_O::propertiesAsString() const
{_OF();
  IMPLEMENT_ME();
}

void Bond_O::clearProperty(core::Symbol_sp prop)
{
  this->_Properties = core::core__rem_f(this->_Properties,prop);
}

CL_LISPIFY_NAME("bond-setProperty");
CL_DEFMETHOD void Bond_O::setProperty(core::Symbol_sp prop, core::T_sp val)
    {
      this->_Properties = core::core__put_f(this->_Properties,val,prop);
    }


CL_LISPIFY_NAME("bond-getProperty");
CL_DEFMETHOD core::T_sp Bond_O::getProperty(core::Symbol_sp prop, core::T_sp defval)
{
  core::T_sp res = core::cl__getf(this->_Properties,prop,_Unbound<core::T_O>());
  if (res.unboundp()) {
    return defval;
  }
  return res;
}

CL_LISPIFY_NAME("bond-hasProperty");
CL_DEFMETHOD bool Bond_O::hasProperty(core::Symbol_sp prop)
{
  return !core::cl__getf(this->_Properties,prop,_Unbound<core::T_O>()).unboundp();
}



CL_LISPIFY_NAME("getOtherAtom");
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
{_OF();
  stringstream    ss;
  ss << "(Bond ";
  Atom_sp wa1 = this->_Atom1;
  ss << _rep_(wa1);
  ss << bondOrderToChar(this->order);
  Atom_sp wa2 = this->_Atom2;
  ss << _rep_(wa2);
  ss << " @"<<std::hex<<this<<std::dec<<")";
  return ss.str();
}

string  Bond_O::describeOther(Atom_sp from) const
{_OF();
  ASSERTF(from==this->_Atom1||from==this->_Atom2,BF("describeFrom failed because from atom[%s] is not part of bond: %s") % from->description() % this->description() );
  stringstream    ss;
  ss << "Bond(";
  ss << from->description();
  ss << bondOrderToChar(this->order);
  Atom_sp a2 = this->getOtherAtom(from);
  ss << a2->description();
  ss << ")";
  return ss.str();
}



bool	Bond_O::invalid(Atom_sp a)
{_OF();
  ASSERTNOTNULL(this->_Atom1);
  ASSERTNOTNULL(this->_Atom2);
  if ( this->_Atom1.nilp() ) return true;
  if ( this->_Atom2.nilp() ) return true;
  return false;
}

void	Bond_O::failIfInvalid(Atom_sp a)
{
  if ( this->invalid(a) ) {
    SIMPLE_ERROR(BF("INVALID %s")%this->description());
  }
}



void Bond_O::redirectToAtomCopies()
{
  LOG(BF("Redirecting bond@%p") % this ); //
  Atom_sp fa =  this->_Atom1;
  Atom_sp ta = this->_Atom2;
#ifdef	DEBUG_ON
  if ( !fa.objectp() ) {
    SIMPLE_ERROR(BF("redirectToAtomCopies _Atom1 is NULL"));
  }
  LOG(BF("  original from atom@%p") % fa ); //
  if ( !ta.objectp() ) {
    SIMPLE_ERROR(BF("redirectToAtomCopies _Atom2 is NULL"));
  }
  LOG(BF("  original   to atom@%p") % ta ); //
#endif
  Atom_wp fc = fa->getCopyAtom();
  ASSERTNOTNULL(fc);
  LOG(BF("    new from %s") % fc->description() );
  Atom_wp tc = ta->getCopyAtom();
  ASSERTNOTNULL(tc);
  LOG(BF("    new   to %s") % tc->description() );
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



ConstitutionBond_sp Bond_O::asConstitutionBond(Atom_sp from, MapAtomsToConstitutionAtomIndex0N atomMap)
{_OF();
  Atom_sp to = this->getOtherAtom(from);
  MapAtomsToConstitutionAtomIndex0N::iterator it = atomMap.find(to);
  if ( it == atomMap.end() )
  {
    SIMPLE_ERROR(BF("Could not find atom[%s] in atomMap") % _rep_(to));
  }
  ConstitutionAtomIndex0N index = it->second;
  ConstitutionBond_sp cb = makeConstitutionBond(index, this->getOrder());
  return cb;
}




CL_LISPIFY_NAME("getOrderAsString");
CL_DEFMETHOD string	Bond_O::getOrderAsString()
{
  return bondOrderToString(this->order);
}

void	Bond_O::joinYourAtoms()
{_OF();
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
      return '?';
      break;
  case singleBond:
      return '-';
      break;
  case doubleBond:
      return '=';
      break;
  case tripleBond:
      return '#';
      break;
  case aromaticBond:
      return '@';
      break;
  case dashedSingleBond:
      return '!';
  case dashedDoubleBond:
      return '%';
      break;
  default:
      THROW_HARD_ERROR(BF("No char available for bond"));
  }
  printf( "Illegal bond order: %d\n", (int)(bo));
  exit(1);
}


BondList_O::BondList_O(const BondList_O& orig) : Base(orig), _Bonds(orig._Bonds)
{_OF();
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
{_OF();
  BondList_O::iterator	bi;
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    if ( (*bi)->isBondBetween(a,b) )
    {
      this->_Bonds.erase(bi);
      return;
    }
  }
  SIMPLE_ERROR(BF("Could not find bond between %s and %s in bondList: %s") % a->description() % b->description() % this->description() );
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

#if 0
/*! Shallow copy of a bond list */
BondList_sp	BondList_O::copy()
{
  BondList_O::iterator	bi;
  GC_ALLOCATE(BondList_O, bl );
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    bl->addBond(*bi);
  }
  return bl;
}
#endif

string	BondList_O::description() const
{
  BondList_O::const_iterator	bi;
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
  BondList_O::const_iterator	bi;
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
    LOG(BF("About to save BondList with %d members") % this->size()  ); //
    for ( ib=this->_Bonds.begin(); ib!=this->_Bonds.end(); ib++ ) {
      _BLOCK_TRACE("Archiving one bond of bondList");
      node->archiveSaveObjectAsChildAssignAutoUniqueId<Bond_O>(*ib);
      LOG(BF("After archiveObject just print") ); //
    }
  } else {
    LOG(BF("About to load BondList") ); //
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
  SYMBOL_EXPORT_SC_(ChemKwPkg,doubleBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,tripleBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,aromaticBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,hydrogenBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,virtualBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,dashedSingleBond);
  SYMBOL_EXPORT_SC_(ChemKwPkg,dashedDoubleBond);
  CL_BEGIN_ENUM(BondOrder,_sym__PLUS_bondOrderToSymbolConverter_PLUS_,"BondOrder");
  CL_VALUE_ENUM(chemkw::_sym_noBond, noBond );
  CL_VALUE_ENUM(chemkw::_sym_singleBond, singleBond );
  CL_VALUE_ENUM(chemkw::_sym_doubleBond, doubleBond );
  CL_VALUE_ENUM(chemkw::_sym_tripleBond, tripleBond );
  CL_VALUE_ENUM(chemkw::_sym_aromaticBond, aromaticBond );
  CL_VALUE_ENUM(chemkw::_sym_hydrogenBond, hydrogenBond );
  CL_VALUE_ENUM(chemkw::_sym_virtualBond, virtualBond);
  CL_VALUE_ENUM(chemkw::_sym_dashedSingleBond,dashedSingleBond);
  CL_VALUE_ENUM(chemkw::_sym_dashedDoubleBond,dashedDoubleBond);
  CL_END_ENUM(_sym__PLUS_bondOrderToSymbolConverter_PLUS_);













}; // namespace chem
