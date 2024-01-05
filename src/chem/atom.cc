/*
    File: atom.cc
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
#include <string.h>
#include <iomanip>
#include <clasp/core/foundation.h>
#include <clasp/core/common.h>
#include <clasp/core/record.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/array.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/matter.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <cando/chem/molecule.h>
#include <cando/chem/restraint.h>
#include <clasp/core/numerics.h>
#include <clasp/core/hashTable.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/sort.h>
#include <cando/geom/omatrix.h>
#include <clasp/core/numbers.h>
#include <cando/chem/elements.h>
#include <clasp/core/translators.h>
#include <cando/chem/constitutionAtoms.h>

#include <clasp/core/wrappers.h>


namespace chem
{

SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_configurationEnumConverter_PLUS_);
SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_stereochemistryTypeConverter_PLUS_);
SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_atomFlagSymbolConverter_PLUS_);

/* Important properties
 * used by RenderMatter to render atoms
 *
 * :colorByElement true/false - if true overrides color; 
 * 				if not set and color is set then use color
 * 				if not set and color not set then colorByElement
 *
 * :color type(GrColor)		- defines color  (default byelement)
 * :renderStyle			- one of :none, :line, :ballAndStick or :cpk (default :none if it has bonds and :ball if it doesn't)
 * :showLabel			- true or false (default false)
 * :label			- Label to show (default "")
 */


std::atomic<size_t> globalUniqueAtomOrder;


core::NullTerminatedEnumAssociation configurationEnum[] = {
    { "UNDEFINED", undefinedConfiguration },
    { "S", S_Configuration },
    { "R", R_Configuration },
    { "RightHanded", RightHanded_Configuration },
    { "LeftHanded", LeftHanded_Configuration },
    { "", -1 }
};

core::NullTerminatedEnumAssociation stereochemistryTypes[] = {
    { "UNDEFINED", undefinedCenter },
    { "CHIRAL", chiralCenter },
    { "PROCHIRAL", prochiralCenter },
    { "", -1 }
};


struct	VdwRadiiPairs
{
  string	_Element;
  double	_VdwRadius;
};



CL_LISPIFY_NAME(make-atom);
DOCGROUP(cando);
CL_DEFUN Atom_sp Atom_O::make(core::Symbol_sp name, chem::Element element)
{
  auto atom = gctools::GC<Atom_O>::allocate_with_default_constructor();
  atom->setName(name);
  atom->setElement(element);
  return atom;
};

CL_LISPIFY_NAME("randomizeAtomPosition");
CL_LAMBDA((core::a chem:atom) &optional (width 10.0));
CL_DEFMETHOD void Atom_O::randomizeAtomPosition(double width)
{
  double				x,y,z;
  Vector3				v;
  double half_width = width/2.0;
  x = core::randomNumber01()*width-half_width;
  y = core::randomNumber01()*width-half_width;
  z = core::randomNumber01()*width-half_width;
  v.set(x,y,z);
  this->setPosition(v);
}



//
//	perturbAtomPosition
//
//	Set the coordinates of all atoms to random values
//		between 0 and 100 angstroms
//
//
CL_LISPIFY_NAME("perturbAtomPosition");
CL_DEFMETHOD void	Atom_O::perturbAtomPosition(double dist)
{
  double				xd,yd,zd;
  Vector3				v,vd;
  xd = (core::randomNumber01()*2.0-1.0)*dist;
  yd = (core::randomNumber01()*2.0-1.0)*dist;
  zd = (core::randomNumber01()*2.0-1.0)*dist;
  v = this->getPosition();
  vd.set(xd,yd,zd);
  v = v+vd;
  this->setPosition(v);
}


CL_NAME("GET-POSITION!");
CL_DEFMETHOD void Atom_O::getPosition_BANG_(Vector3& pos)
{
  pos = this->_Position;
}

/*! Return +1 if priority(a)>priority(b)
 *  return -1 if priority(a)<priority(b)
 *  return 0 if priority(a) == priority(b)
 *  Priority is decided by comparing RelativePriority and if they are the same
 *  by comparing names.
 */
int Atom_O::priorityOrder(Atom_sp a, Atom_sp b, core::HashTable_sp cip_priority)
{
  if ( a.nilp() && b.nilp() ) return 0;
  if ( b.nilp() ) return 1;
  if ( a.nilp() ) return -1;
  Fixnum apriority = core::clasp_to_fixnum(cip_priority->gethash(a));
  Fixnum bpriority = core::clasp_to_fixnum(cip_priority->gethash(b));
  if ( apriority>bpriority ) return 1;
  if ( apriority<bpriority ) return -1;
  if ( a->getName()->symbolNameAsString() > b->getName()->symbolNameAsString() ) return 1;
  if ( a->getName()->symbolNameAsString() < b->getName()->symbolNameAsString() ) return -1;
  return 0;
}

CL_DOCSTRING(R"dx(Calculate the stereochemical configuration of the atom.
This requires that relative cip priorities are defined using CipPrioritizer_O::assignCahnIngoldPrelogPriorityToAtomsRelativePriority.)dx")
CL_LISPIFY_NAME("calculateStereochemicalConfiguration");
CL_DEFMETHOD     ConfigurationEnum Atom_O::calculateStereochemicalConfiguration(core::HashTable_sp cip_priority)
{
  if ( this->numberOfBonds() != 4 ) return undefinedConfiguration;
  core::List_sp neighborsByPriority = this->getNeighborsByRelativePriority(cip_priority);
  Atom_sp a1 = core::oCar(neighborsByPriority).as<Atom_O>();
  Atom_sp a2 = core::oCadr(neighborsByPriority).as<Atom_O>();
  Atom_sp a3 = core::oCaddr(neighborsByPriority).as<Atom_O>();
  Atom_sp a4 = core::oCadddr(neighborsByPriority).as<Atom_O>();
  Vector3 vme = this->getPosition();
  Vector3 v1 = a1->getPosition().sub(vme);
  Vector3 v2 = a2->getPosition().sub(vme);
  Vector3 v3 = a3->getPosition().sub(vme);
  Vector3 v4 = a4->getPosition().sub(vme);
  Vector3 v43cross = v4.crossProduct(v3);
  double dir1 = v1.dotProduct(v43cross); // For R this should be >0
  double dir2 = v2.dotProduct(v43cross); // For R this should be <0
  if ( dir1 > 0.0 ) {
    if ( dir2 < 0.0 ) {
      return R_Configuration;
    }
  } else if ( dir1 < 0.0 ) {
    if ( dir2 > 0.0 ) {
      return S_Configuration;
    }
  }
  return undefinedConfiguration;
}


string Atom_O::calculateStereochemicalConfigurationAsString(core::HashTable_sp cip)
{
  string s;
  ConfigurationEnum config = this->calculateStereochemicalConfiguration(cip);
  switch ( config )
  {
  case S_Configuration:
      s = "S";
      break;
  case R_Configuration:
      s = "R";
      break;
  case RightHanded_Configuration:
      s = "RightHanded";
      break;
  case LeftHanded_Configuration:
      s = "LeftHanded";
      break;
  default:
      s = "undefinedConfiguration";
      break;
  }
  return s;
}

CL_DEFMETHOD bool Atom_O::Atom_equal(core::T_sp obj) const
{
  if ( this->eq(obj) ) goto T;
  if ( obj.isA<Atom_O>() )
  {
    Atom_sp other = obj.as<Atom_O>();
    if ( other->getName() != this->getName() )
    {
      LOG("Atom names[ this->getName()={}  other->getName()={} ]  don't match!" , this->getName() , other->getName() );
      goto F;
    }
    if ( other->numberOfBonds() != this->numberOfBonds() )
    {
      LOG("Numbers of bonds don't match this->numberOfBonds()={}  other->numberOfBonds()={}!" , this->numberOfBonds() , other->numberOfBonds() );
      goto F;
    }
    goto T;
  }
  LOG("The other object is not an ATOM!!!!");
 F:
  LOG("atoms this({}) =equal= other({})  equal==false" , this->description() , obj->description() );
  return false;
 T:
  LOG("atoms this({}) =equal= other({})  equal==true" , this->description() , obj->description() );
  return true;
}

void Atom_O::transferCoordinates(Matter_sp obj)
{
  if ( !this->Atom_equal(obj) )
  {
    SIMPLE_ERROR("This atom({}) is not equal to {} so you cannot transfer coordinates" , this->description() , obj->description() );
  }
  Atom_sp other = obj.as<Atom_O>();
  this->setPosition(other->getPosition());
}

bool Atom_O::atomWithinAngstroms(Atom_sp other, float angstroms) const
{
  Vector3 delta = this->_Position-other->_Position;
  if (fabs(delta.getX())<angstroms &&
      fabs(delta.getY())<angstroms &&
      fabs(delta.getZ())<angstroms) {
    if (delta.length()<angstroms) return true;
  }
  return false;
}

//
// Copy constructor
//
Atom_O::Atom_O(const Atom_O& ss) :Matter_O(ss)
{
  this->_UniqueAtomOrder = nextUniqueAtomOrder();
  this->_Element = ss._Element;
  this->_Alias = ss._Alias;
  this->_Hybridization = ss._Hybridization;
  this->_Flags = ss._Flags;
  this->_Mask = ss._Mask;
  this->_Position = ss._Position;
  this->_Charge = ss._Charge;
  this->_RingMembershipCount = ss._RingMembershipCount;
  this->_VdwRadius = ss._VdwRadius;
  this->_CovalentRadius = ss._CovalentRadius;
  this->_CopyAtom = Atom_sp();
  this->_TempInt = ss._TempInt;
//  this->moeIndex = ss.moeIndex;
//  this->moeType = ss.moeType;
  this->_Ionization = ss._Ionization;
  this->_Configuration = ss._Configuration;
  this->_StereochemistryType = ss._StereochemistryType;
}

//
// Destructor
//
//	Remove all bonds to this atom
//


CL_LISPIFY_NAME("setConfiguration");
CL_DEFMETHOD     void Atom_O::setConfiguration(ConfigurationEnum conf)
{
  this->_Configuration = conf;
  if (chem__verbose(1)) {
    core::clasp_write_string(fmt::format("Changed configuration of atom {} to {}\n" , _rep_(this->asSmartPtr()) , this->getConfigurationAsString()));
  }
  LOG("Changed configuration of atom[{}] to [{}]"
      , this->__repr__() , this->getConfigurationAsString() );
}

void Atom_O::setElementFromAtomName()
{
  string nameStr = this->getName()->symbolName()->get_std_string();
  Element element = elementFromAtomNameString(nameStr);
  LOG(" Resulting element= |{}|" %element );
  this->setElement(element);
}

void Atom_O::setElementFromString(const string& str)
{
  Element element = elementFromAtomNameStringCaseInsensitive(str);
  LOG(" Resulting element= |{}|" %element );
  this->setElement(element);
}


void Atom_O::setHybridizationFromString(const string& h)
{
  core::SimpleBaseString_sp sbs = core::SimpleBaseString_O::make(h);
  core::SimpleString_sp usbs = core::cl__string_upcase(sbs);
  core::Symbol_sp sym = _lisp->internKeyword(usbs->get_std_string());
  this->_Hybridization = hybridizationForSymbol(sym);
}

void Atom_O::_addHydrogenWithName(Residue_sp residueContainedBy, MatterName name)
{
  Atom_sp h = Atom_O::create();
  h->setf_needs_build(true);
  h->setName(name);
  h->setElement(element_H);
  residueContainedBy->addAtom(h);
  this->bondTo(h,singleBond);
}

Bond_sp Atom_O::bondAtIndex(int i)
{
  return this->_Bonds[i];
}

CL_LISPIFY_NAME("bondedNeighbor");
CL_DEFMETHOD     Atom_sp Atom_O::bondedNeighbor(int i)
{
  return this->_Bonds[i]->getOtherAtom(this->sharedThis<Atom_O>());
}

CL_LISPIFY_NAME("bondedOrder");
CL_DEFMETHOD     BondOrder Atom_O::bondedOrder(int i) 
{
  return this->_Bonds[i]->getOrderFromAtom(this->asSmartPtr());
};

CL_LISPIFY_NAME("numberOfBonds");
CL_DEFMETHOD     int	Atom_O::numberOfBonds() const
{
  return this->_Bonds.size();
};

CL_LISPIFY_NAME("coordination");
CL_DEFMETHOD     int Atom_O::coordination() 
{
  return this->_Bonds.size();
};

/*! If this atom is a carbon then count the number of bonds and add
 * additional hydrogens to get the number of bonds to 4.
 * Only do this on neutral carbons.
 * Give the hydrogens names based on this atoms name.
 * eg:
 * This atom name :  Number of Hydrogens to add :  Hydrogen names
 *       CA                     1                       HA
 *       CB                     2                       HB1, HB2
 *       CC1                    3                       HC11, HC12, HC13 
 * Also handle N and O
 */ 
CL_LISPIFY_NAME("numberOfOpenValence");
CL_DEFMETHOD uint Atom_O::numberOfOpenValence()
{
  int maxHydrogens = 0;
  switch (this->getElement()) {
  case element_C:
      if ( this->getIonization() != 0 )
      {
        core::Warn(core::Str_O::create("Add support for ~d ionization of carbon"),
                   core::Cons_O::create(core::clasp_make_fixnum(this->getIonization()),nil<core::T_O>()));
      }
      maxHydrogens = 4;
      break;
  case element_N:
      maxHydrogens = 3+this->getIonization();
      break;
  case element_O:
      maxHydrogens = 2+this->getIonization();
      break;
  default:
      //printf("%s:%d createImplicitHydrogenNames skipping atom: %s element: %s\n", __FILE__, __LINE__, _rep_(this->getName()).c_str(), this->getElementAsString().c_str());
      return 0;
  }
  uint totalBondOrder = this->totalBondOrder();
  uint addHydrogens = 0;
  if ( maxHydrogens >= totalBondOrder ) {
    addHydrogens = maxHydrogens - totalBondOrder;
  }
  if ( addHydrogens > 10 ) {
    SIMPLE_ERROR("There are too many hydrogens to be added for element {} totalBondOrder: {} addHydrogens: {}" , _rep_(symbolFromElement(this->getElement())) , totalBondOrder , addHydrogens );
  }
  return addHydrogens;
}

/*! If this atom is a carbon then count the number of bonds and add
 * additional hydrogens to get the number of bonds to 4.
 * Only do this on neutral carbons.
 * Give the hydrogens names based on this atoms name.
 * eg:
 * This atom name :  Number of Hydrogens to add :  Hydrogen names
 *       CA                     1                       HA
 *       CB                     2                       HB1, HB2
 *       CC1                    3                       HC11, HC12, HC13 
 * Also handle N and O
 */ 
CL_LISPIFY_NAME("createImplicitHydrogenNames");
CL_DEFMETHOD core::List_sp Atom_O::createImplicitHydrogenNames()
{
  uint addHydrogens = this->numberOfOpenValence();
  string nameSuffix = this->getName()->symbolName()->get_std_string().substr(1,9999);
  core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
  core::Cons_sp cons = first;
  if ( addHydrogens == 1 )
  {
    MatterName hname = chemkw_intern("H"+nameSuffix);
    core::Cons_sp one = core::Cons_O::create(hname,nil<core::T_O>());
    cons->setCdr(one);
    cons = one;
    return oCdr(first);
  }
  for ( uint i=0; i<addHydrogens; i++ )
  {
    stringstream ss;
    ss << "H" << nameSuffix << (i+1);
    MatterName hname = chemkw_intern(ss.str());
    core::Cons_sp one = core::Cons_O::create(hname,nil<core::T_O>());
    cons->setCdr(one);
    cons = one;
  }
  return oCdr(first);
}



CL_DEFMETHOD
size_t Atom_O::fillInImplicitHydrogensWithResidue(Residue_sp residue)
{
  core::List_sp names = this->createImplicitHydrogenNames();
  if ( names.nilp() ) return 0;
  size_t num = 0;
  for ( auto cur : names ) {
    this->_addHydrogenWithName(residue,oCar(cur).as<MatterName::Type>());
    ++num;
  }
  return num;  
}

size_t Atom_O::fillInImplicitHydrogens()
{
  SIMPLE_ERROR("You cannot invoke chem:fill-in-implicit with an atom - any new hydrogens must be added to the residue that contains this atom");
}

//
//	randomizeAtomPosition
//
//	Set the coordinates of all atoms to random values
//		between 0 and 100 angstroms
//
//



gc::Nilable<Atom_sp> Atom_O::highestPriorityNeighborThatIsnt(gc::Nilable<Atom_sp> avoid, core::HashTable_sp cip)
{
  Atom_sp bestAtom = nil<Atom_O>();
  Atom_sp atom;
  VectorBond::iterator	b;
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
    atom = (*b)->getOtherAtom(this->sharedThis<Atom_O>());
    if ( avoid.notnilp() && atom == avoid ) continue;
    if ( priorityOrder(atom,bestAtom,cip)>0 )
    {
      bestAtom = atom;
    }
  }
  return bestAtom;
}
    

gc::Nilable<Atom_sp> Atom_O::lowestPriorityNeighborThatIsnt(gc::Nilable<Atom_sp> avoid, core::HashTable_sp cip)
{
  gc::Nilable<Atom_sp> bestAtom = nil<core::T_O>();
  Atom_sp atom;
  VectorBond::iterator	b;
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
    atom = (*b)->getOtherAtom(this->sharedThis<Atom_O>());
    if ( avoid.notnilp() && atom == avoid ) continue;
    if ( bestAtom.nilp() ) {
      bestAtom = atom;
      continue;
    }
    if ( priorityOrder(atom,bestAtom,cip)<0 ) bestAtom = atom;
  }
  return bestAtom;
}



CL_LISPIFY_NAME("getBondTo");
CL_DEFMETHOD Bond_sp Atom_O::getBondTo(Atom_sp a)
{
  VectorBond::iterator	b;
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ )
  {
    if ( a == (*b)->getOtherAtom(this->sharedThis<Atom_O>()) ) return *b;
  }
  SIMPLE_ERROR("Could not find bond");
}

void Atom_O::_addExistingBond(Bond_sp const& bond)
{
  this->_Bonds.push_back(bond);
}

void Atom_O::addBond(Bond_sp bond)
{
  this->_Bonds.push_back(bond);
}

CL_DOCSTRING(R"doc(Create a bond between two atoms with the desired bond-order.
If error-if-exists is T and the bond already exists then signal an error.
If error-if-exists is NIL and the bond already exists then remove the old bond and create the new one.)doc");
CL_LISPIFY_NAME("bondTo");
CL_LAMBDA((from chem:atom) to bond_order &optional (error_if_exists t) (error_if_exceed_valence t));
CL_DEFMETHOD     Bond_sp Atom_O::bondTo( Atom_sp to, BondOrder o, bool error_if_exists, bool error_if_exceed_valence )
{
	// Check if there is already a bond to this atom and
	// throw an exception if there is
  Atom_sp from = this->sharedThis<Atom_O>();
  VectorBond::iterator	b;
  VectorBond::iterator  begin(this->_Bonds.begin());;
  VectorBond::iterator  end(this->_Bonds.end());
  for ( b=begin;b!=end ; b++ ) {
    if ( (*b)->getOtherAtom(from) == to ) {
      if (error_if_exists) {
        SIMPLE_ERROR("You tried to form a bond from[{}]-to[{}] but there is already one there!!" , this->__repr__() , _rep_(to) );
      } else {
        this->_Bonds.erase(b);
        break;
      }
    }
  }
  Bond_sp bn = Bond_O::create(from,to,o);
  this->_Bonds.push_back(bn);
  to->_Bonds.push_back(bn);
  if (error_if_exceed_valence) {
    if (this->_Element == element_C && this->_Bonds.size()>4) {
      SIMPLE_ERROR("More than four bonds to carbon will be made to {}" , _rep_(this->asSmartPtr()));
    }
    if (to->_Element == element_C && to->_Bonds.size()>4 ) {
      SIMPLE_ERROR("More than four bonds to carbon will be made to {}" , _rep_(to));
    }
  }
  return bn;
}

CL_LISPIFY_NAME("bondToOrderInt");
CL_DEFMETHOD     Bond_sp Atom_O::bondToOrderInt( Atom_sp to, int o )
{
  return this->bondTo(to,(BondOrder)(o));
};





CL_LISPIFY_NAME("bondToSingle");
CL_DEFMETHOD     Bond_sp Atom_O::bondToSingle(  Atom_sp a )
{
  return this->bondTo(a,singleBond);
}



CL_LISPIFY_NAME("getNameIndex");
CL_DEFMETHOD     string	Atom_O::getNameIndex()
{
  FIX_ME();
#if 0
  string	s;
  string	sIndex;
  char	caIndex[255];

  s = this->getName()+":";
  sprintf( caIndex, "%d", this->moeIndex );
  s = s + caIndex;
  return(s);
#endif
}


class	OrderByPriorityAndName
{
  core::HashTable_sp _priority;
public:
  bool operator() ( Atom_sp x, Atom_sp y )
  {
    Fixnum xpriority = core::clasp_to_fixnum(this->_priority->gethash(x));
    Fixnum ypriority = core::clasp_to_fixnum(this->_priority->gethash(y));
    if ( xpriority<ypriority ) return true;
    if ( xpriority>ypriority ) return false;
    if ( x->getName()->symbolNameAsString()<y->getName()->symbolNameAsString() ) return true;
    return false;
  }
  OrderByPriorityAndName(core::HashTable_sp cip_priority) : _priority(cip_priority) {};
};

CL_DEFMETHOD int Atom_O::getRelativePriority(core::HashTable_sp cip) const {
    return core::clasp_to_fixnum(cip->gethash(this->asSmartPtr()));
}

CL_LISPIFY_NAME(Atom/getNeighborsForAbsoluteConfiguration);
CL_DEFMETHOD
core::List_sp	Atom_O::getNeighborsForAbsoluteConfiguration()
{
  ql::list l;
  for ( auto b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
    l << (*b)->getOtherAtom(this->sharedThis<Atom_O>());
  }
  return l.cons();
}

CL_DOCSTRING(R"dx(Return the neighbors sorted from highest CIP priority to lowest)dx");
CL_DEFMETHOD
core::List_sp	Atom_O::getNeighborsByRelativePriority(core::HashTable_sp cip_priority)
{
  LOG("Ordering neighbors around({}) by priority and name" , this->getName()  );
  VectorAtom	reversedNeighbors;
  VectorBond::iterator	b;
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
    reversedNeighbors.push_back((*b)->getOtherAtom(this->sharedThis<Atom_O>()));
  }
  OrderByPriorityAndName orderer(cip_priority);
  core::List_sp ncons = nil<core::T_O>();
  sort::quickSortVec0(reversedNeighbors, 0, reversedNeighbors.size(), orderer);
	// At this point the sorted in reverse order
	// 
	// Now copy them into a Cons list backwards to get the
	// proper order
  LOG("The sorted neighbors in REVERSE order is:" );
  for ( VectorAtom::iterator ni=reversedNeighbors.begin();
        ni!=reversedNeighbors.end(); ni++ )
  {
    LOG("    neighbor priority({}) name({})" , (*ni)->getRelativePriority() , (*ni)->getName()  );
    ncons = core::Cons_O::create(*ni,ncons);
  }
  return ncons;
}





void	Atom_O::modifyFlags(int op, ATOM_FLAGS fl ) {
  if ( op == ATOMFLAG_OR ) {
    this->_Flags |= fl;
  } else if ( op==ATOMFLAG_AND ) {
    this->_Flags &= fl;
  } else if ( op==ATOMFLAG_ON ) {
    this->_Flags |= fl;
  } else if ( op==ATOMFLAG_OFF ) {
    this->_Flags &= (~fl);
  }
}




CL_LISPIFY_NAME("isBondedToAtomNamed");
CL_DEFMETHOD     bool	Atom_O::isBondedToAtomNamed(core::Symbol_sp name)
{
  VectorBond::iterator	b;

  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
    if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>())->getName() == name ) {
      return true;
    }
  }
  return false;
}

CL_LISPIFY_NAME("bondedNeighborWithName");
CL_DEFMETHOD     Atom_sp	Atom_O::bondedNeighborWithName(MatterName name)
{
  VectorBond::iterator	b;
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
    if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>())->getName() == name ) {
      return (*b)->getOtherAtom(this->sharedThis<Atom_O>());
    }
  }
  Atom_sp a = nil<Atom_O>();
  return a;
}


CL_LISPIFY_NAME("isBondedToElementOrder");
CL_DEFMETHOD     bool	Atom_O::isBondedToElementOrder(Element el, BondOrder o)
{
  VectorBond::iterator	b;
  if (Bond_O::singleBondP(o)) {
    for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
      if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>())->getElement() == el ) {
        if ( Bond_O::singleBondP((*b)->getRawOrder())) return true;
      }
    }
  } else {
    for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
      if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>())->getElement() == el ) {
        if ( (*b)->getRawOrder() == o ) return true;
      }
    }
  }
  return false;
}


CL_LISPIFY_NAME("hasBondWithOrder");
CL_DEFMETHOD     bool Atom_O::hasBondWithOrder(BondOrder o) const
{
  VectorBond::const_iterator	b;
  if (Bond_O::singleBondP(o)) {
    for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
      if ( Bond_O::singleBondP((*b)->getRawOrder())) return true;
    }
  } else {
    for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
      if ( (*b)->getRawOrder() == o) return true;
    }
  }
  return false;
}


CL_LISPIFY_NAME("isBondedToElementHybridization");
CL_DEFMETHOD     bool	Atom_O::isBondedToElementHybridization(Element el, Hybridization hy)
{
  VectorBond::const_iterator	b;

  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) 
  {
    if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>())->getElement() == el ) 
    {
      if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>())->getHybridization() == hy ) return true;
    }
  }
  return false;
}




CL_LISPIFY_NAME("isBondedToElementHybridizationElementHybridization");
CL_DEFMETHOD     bool	Atom_O::isBondedToElementHybridizationElementHybridization(Element el1, Hybridization hy1, Element el2, Hybridization hy2)
{
  IMPLEMENT_ME(); // Handle new way of handling bonds
#if 0
  VectorBond::const_iterator	b;
  VectorBond::const_iterator	b2;
  Atom_sp	a1;
  Atom_sp	a2;
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) 
  {
    a1 = Atom_sp((*b)->getTo());
    if ( a1->getElement() == el1 ) 
    {
      if ( a1->getHybridization() == hy1 ) 
      {
        for ( b2=a1->_Bonds.begin();b2!=a1->_Bonds.end();b2++) 
        {
          a2 = Atom_sp((*b2)->getTo());
          if ( a2 == this->sharedThis<Atom_O>() ) continue;
          if ( a2->getElement() == el2 ) 
          {
            if ( (*b2)->getTo()->getHybridization() == hy2 ) return true;
          }
        }
      }
    }
  }
  return false;
#endif
}

float Atom_O::distanceSquaredToAtom(Atom_sp other)
{
  Vector3 delta = this->_Position-other->_Position;
  return delta.getX()*delta.getX()+delta.getY()*delta.getY()+delta.getZ()*delta.getZ();
}

CL_DOCSTRING(R"dx(Return the distance squared between two atoms)dx");
CL_DOCSTRING_LONG(R"dx(This provides a way to compare atom positions without allocating memory.)dx");
DOCGROUP(cando);
CL_DEFUN
float chem__distance_squared_between_two_atoms(Atom_sp atom1, Atom_sp atom2)
{
  return atom1->distanceSquaredToAtom(atom2);
}


void Atom_O::makeAllAtomNamesInEachResidueUnique()
{
  SIMPLE_ERROR("This should never be called for an Atom_sp");
}


void Atom_O::basicRemoveBondTo(Atom_sp a)
{
  VectorBond::iterator	b;
  Bond_sp				bb;
  Atom_sp				aa;
  Atom_sp				wpa;
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ )
  {
    if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>()) == a )
    {
      this->_Bonds.erase(b);
      return;
    }
  }
  stringstream ss;
  ss << "Trying to remove bond from atom(" << this->description()
     << ") but I can't find the to atom(" << a->description() << ")";
  SIMPLE_ERROR("{}" , ss.str());
}

CL_LISPIFY_NAME("removeBondTo");
CL_DEFMETHOD     void	Atom_O::removeBondTo(Atom_sp a)
{
  Atom_sp	atemp;
  LOG("Bond_O::removeBondTo" );
  this->basicRemoveBondTo(a);
  atemp = this->sharedThis<Atom_O>();
  a->basicRemoveBondTo(atemp);
}

CL_LISPIFY_NAME("removeAllBonds");
CL_DEFMETHOD     void	Atom_O::removeAllBonds()
{
  VectorBond::const_iterator	b;
  Atom_sp				atemp;
	// Remove back bond
  LOG("Bond_O::removeAllBonds for {}" , this->description() );
  atemp = this->sharedThis<Atom_O>();
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ )
  {
    Atom_sp ato = (*b)->getOtherAtom(atemp);
    ato->basicRemoveBondTo(atemp);
  }
  this->_Bonds.clear();
}



//
//	applyTransformToAtoms
//
CL_LISPIFY_NAME("applyTransformToAtoms");
CL_DEFMETHOD     void	Atom_O::applyTransformToAtoms( const Matrix& m )
{
  this->_Position = (m)*this->_Position;
  if ( this->_RestraintList.notnilp() ) {
    this->Base::applyTransformToRestraints(m);
  }
}





string Atom_O::getHybridizationAsString()
{
  return hybridizationSymbolFromHybridization(this->getHybridization())->symbolName()->get_std_string();
}

CL_LISPIFY_NAME("getElementAsString");
CL_DEFMETHOD     string Atom_O::getElementAsString() const
{
  return atomicSymbolFromElement(this->getElement())->symbolName()->get_std_string();
}

CL_DEFMETHOD core::Symbol_sp Atom_O::getElementAsSymbol() const
{
  return symbolFromElement(this->getElement());
}

SYMBOL_EXPORT_SC_(KeywordPkg,given_atom_type);

CL_LISPIFY_NAME("atomType");
CL_DEFMETHOD AtomType Atom_O::atomType() const {
  return this->getProperty(kw::_sym_given_atom_type);
};
CL_LISPIFY_NAME("setAtomType");
CL_DEFMETHOD 	void	Atom_O::setAtomType(AtomType o) {
  this->setProperty(kw::_sym_given_atom_type,o);
}

void	Atom_O::fields(core::Record_sp node)
{
  node->field_if_not_nil( INTERN_(kw,alias), this->_Alias );
  node->/*pod_*/field_if_not_default<size_t>( INTERN_(kw,flags), this->_Flags, 0 );
  node->/*pod_*/field( INTERN_(kw,element), this->_Element);
  node->/*pod_*/field_if_not_default( INTERN_(kw,hybridization), this->_Hybridization,hybridization_sp3 );
  node->/*pod_*/field_if_not_default( INTERN_(kw,chg), this->_Charge, 0.0 );
  node->/*pod_*/field_if_not_default<short>( INTERN_(kw,ion), this->_Ionization, 0 );
  node->/*pod_*/field_if_not_default<ushort>( INTERN_(kw,rings), this->_RingMembershipCount, 0 );
  node->/*pod_*/field_if_not_default( INTERN_(kw,mask), this->_Mask, (uint)(0) );
  node->/*pod_*/field_if_not_default( INTERN_(kw,configuration), this->_Configuration, undefinedConfiguration  );
  node->/*pod_*/field_if_not_default( INTERN_(kw,stereochemistryType), this->_StereochemistryType, undefinedCenter );
  node->/*pod_*/field_if_not_default( INTERN_(kw,pos), this->_Position, Vector3());
//  node->field_if_not_empty(INTERN_(kw,bonds),this->_Bonds);
  this->Base::fields(node);
}

bool Atom_O::applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value ) {
  if ( this->Base::applyPropertyToSlot(prop,value) ) return true;
  return false;
}
CL_LISPIFY_NAME("setPositionInNanometers");
CL_DEFMETHOD     void Atom_O::setPositionInNanometers(Vector3 o)
{
  Vector3 angpos = o.multiplyByScalar(10.0);
  this->_Position = angpos;
}


CL_LISPIFY_NAME("getConfigurationAsString");
CL_DEFMETHOD     string	Atom_O::getConfigurationAsString() const
{
  string config, stereo;
  if ( this->_Configuration == undefinedConfiguration ) {
    config = "UndefinedConfiguration";
  } else if ( this->_Configuration == S_Configuration ) {
    config = "(S)";
  } else if ( this->_Configuration == R_Configuration ) {
    config = "(R)";
  } else if ( this->_Configuration == RightHanded_Configuration ) {
    config = "(rh)";
  } else if ( this->_Configuration == LeftHanded_Configuration ) {
    config = "(lh)";
  } else
  {
    config = "???";
  }
  if ( this->_StereochemistryType == undefinedCenter )
  {
    stereo = "UndefinedCenter";
  } else if ( this->_StereochemistryType == chiralCenter )
  {
    stereo = "ChiralCenter";
  } else if ( this->_StereochemistryType == prochiralCenter )
  {
    stereo = "ProchiralCenter";
  } else
  {
    stereo = "???";
  }
  return config+"/"+stereo;
}

string	Atom_O::__repr__() const
{
  stringstream ss;
  ss << "#<" << this->className() << " " << this->_Name << "/" << _rep_(symbolFromElement(this->_Element)) << " :id " << this->_Id << " 0x" << std::setbase(16) << gctools::lisp_general_badge(this->asSmartPtr());
  if (this->_StereochemistryType!=undefinedCenter) {
    ss << this->getConfigurationAsString();
  }
  ss << ">";
  return ss.str();
}

string	Atom_O::description() const
{
  stringstream ss;
  ss << this->className() << "("<<this->getName();
  ss << " :element " << _rep_(this->getElementAsSymbol());
  ss << " bonds[" << _rep_(this->bondsAsList()) << "]";
  ss <<")";
  return ss.str();
}


/*!
 * Add all unique intra residue bonds to the BondList_sp
 */
void	Atom_O::addUniqueIntraResidueBondCopiesToBondList(core::HashTable_sp atomToResidue,BondList_sp list)
{
  VectorBond			bonds;
  VectorBond::iterator	b;
  ASSERTNOTNULL(list);
  bonds = this->getBonds();
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=bonds.begin();b!=bonds.end() ; b++ ) {
    if ( me->atLowerUniqueAtomOrderThan( (*b)->getOtherAtom(me))) {
      if ( !(*b)->isInterResidueBond(atomToResidue) ) {
        LOG("Original bond: {}" , (*b)->description() );
        auto bondCopy = gctools::GC<Bond_O>::copy(**b);
        LOG("Copy bond: {}" , bondCopy->description() );
        list->addBond(bondCopy);
      }
    }
  }
}

/*!
 * Add all unique inter residue bonds to the BondList_sp
 */
void Atom_O::addUniqueInterResidueBondCopiesToBondList(core::HashTable_sp atomToResidue, BondList_sp list)
{
  VectorBond			bonds;
  VectorBond::iterator	b;
  ASSERTNOTNULL(list);
  bonds = this->getBonds();
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=bonds.begin();b!=bonds.end() ; b++ ) {
    if ( me->atLowerUniqueAtomOrderThan( (*b)->getOtherAtom(me))) {
      if ( (*b)->isInterResidueBond(atomToResidue) ) {
        auto  bondCopy = gctools::GC<Bond_O>::copy( **b ); // = RP_Copy<Bond_O>(*b);
        list->addBond(bondCopy);
      }
    }
  }
}

void Atom_O::addInterResidueBondsToBondList(core::HashTable_sp atomToResidue, BondList_sp bondlist)
{
  VectorBond::iterator	b;
  this->_Bonds = this->getBonds();
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
    if ( (*b)->isInterResidueBond(atomToResidue) ) {
      bondlist->append(*b);
    }
  }
}


//
// isBondedTo
//
// Return true if we are bonded to the atom described by the required
// bond order
CL_LISPIFY_NAME("isBondedToWithBondOrder");
CL_DEFMETHOD     bool Atom_O::isBondedToWithBondOrder( Atom_sp aTarget, BondOrder o )
{
  VectorBond::iterator	b;
  Atom_sp				a2;
  Atom_sp me = this->sharedThis<Atom_O>();
  if (Bond_O::singleBondP(o)) {
    for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
      a2 = (*b)->getOtherAtom(me);
      if ( a2 == aTarget ) {
        if ( Bond_O::singleBondP((*b)->getRawOrder())) return true;
      }
    }
  } else {
    for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
      a2 = (*b)->getOtherAtom(me);
      if ( a2 == aTarget ) {
        if ( o == (*b)->getRawOrder() ) return true;
      }
    }
  }
  return false;
}

//
// isBondedTo
//
// Return true if we are bonded to the atom described by the required
// bond order
CL_LISPIFY_NAME("isBondedTo");
CL_DEFMETHOD     bool	Atom_O::isBondedTo( Atom_sp aTarget)
{
  VectorBond::iterator	b;
  Atom_sp				a2;
  LOG("Looking at atom({}) for bonded atom({})"
      , this->getName()
      , aTarget->getName() );
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ )
  {
    a2 = (*b)->getOtherAtom(me);
    LOG("  looking at atom({})" , a2->getName() );
    if ( a2 == aTarget )
    {
      LOG("Found match" );
      return true;
    }
  }
  LOG("No match" );
  return false;
}
//
// bondOrderTo
//
// Return the bond order by which we are bonded to the atom.
// Return noBond if there is no bond.
CL_DOCSTRING(R"dx(Return the bond-order to the OTHER atom.)dx");
CL_DEFMETHOD
BondOrder	Atom_O::bondOrderTo( Atom_sp other )
{
  VectorBond::iterator	b;
  Atom_sp				a2;
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ ) {
    a2 = (*b)->getOtherAtom(me);
    if ( a2 == other ) {
      return (*b)->getOrderFromAtom(me);
    }
  }
  return noBond;
}

CL_DOCSTRING(R"dx(Reorder the bonds.)dx");
CL_LISPIFY_NAME(Atom/reorderBonds);
CL_DEFMETHOD
void Atom_O::reorderBonds(core::List_sp atoms) {
  VectorBond newBonds;
  size_t numAtoms = core::cl__length(atoms);
  Atom_sp me = this->sharedThis<Atom_O>();
  if ( numAtoms != this->_Bonds.size()) {
    SIMPLE_ERROR("The number of atoms {} must match the number of bonds {}", numAtoms, this->_Bonds.size() );
  }
  for ( auto cur : atoms ) {
    Atom_sp atom = gc::As<Atom_sp>(CONS_CAR(cur));
    VectorBond::iterator b;
    for ( b = this->_Bonds.begin(); b != this->_Bonds.end(); b++ ) {
      if ( (*b)->getOtherAtom(me) == atom ) {
        newBonds.push_back(*b);
        break;
      }
    }
    if (b == this->_Bonds.end()) {
      SIMPLE_ERROR("For atom {} not find bond with other atom {} bonded-atoms: {}", _rep_(me), _rep_(atom), _rep_(this->bondedAtomsAsList()));
    }
  }
  this->_Bonds = newBonds;
}

//
// testConsistancy
//
//	Make sure that all of the bonds are consistant
CL_LISPIFY_NAME("testConsistancy");
CL_DEFMETHOD     bool Atom_O::testConsistancy(Matter_sp parentShouldBe)
{
  VectorBond::iterator	b;
  Atom_sp				a2,atemp;
  LOG("bonds array start = 0x%08x" , ((void*)&(this->_Bonds[0]))  );
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=this->_Bonds.begin();b!=this->_Bonds.end() ; b++ )
  {
    a2 = (*b)->getOtherAtom(me);
    if ( !(a2->isBondedToWithBondOrder(me,(*b)->getRawOrder())) )
    {
      core::clasp_write_string("Atom_O::testConsistancy failed\n");
      return false;
    }
  }
  return true;
}



Matter_sp Atom_O::copy(core::T_sp new_to_old)
{
  auto atm = gc::GC<Atom_O>::copy(*this);
  atm->_Bonds.clear();
  return atm;
}




//
//	copy
//
//	Allocate a copy of this atom and return it.
//	Keep track of the new atom so that we can redirect pointers to the copy.
//
Matter_sp Atom_O::copyDontRedirectAtoms(core::T_sp new_to_old)
{
  LOG("Copying atom @{}" , this );
  auto  aNew  = gctools::GC<Atom_O>::copy( *this); // = RP_Copy<Atom_O>(this);
  aNew->_Bonds.clear();
  this->_CopyAtom = aNew;
  if (this->_Properties.notnilp()) {
    aNew->_Properties = core::cl__copy_seq(this->_Properties);
  }
  aNew->copyRestraintsDontRedirectAtoms(this->asSmartPtr());
  LOG("    copy atom== {}" , aNew->description());
  return(aNew);
}



//
//	redirectAtom
//
//	Redirect the bonds of this atom to their copies.
void Atom_O::redirectAtoms()
{
  this->redirectRestraintAtoms();
}




VectorAtom Atom_O::getBondedAtoms()
{
  VectorBond::iterator	b;
  VectorAtom			vAtoms;
  Atom_sp me = this->sharedThis<Atom_O>();
  for (b=this->_Bonds.begin();b!=this->_Bonds.end(); b++ )
  {
    vAtoms.push_back((*b)->getOtherAtom(me));
  }
  return vAtoms;
}

CL_LISPIFY_NAME("Atom/bondedAtomsAsList");
CL_DEFMETHOD     core::List_sp Atom_O::bondedAtomsAsList() const
{
  ql::list ll;
  VectorBond::const_iterator	b;
  Atom_sp me(const_cast<Atom_O*>(this));
  for (b=this->_Bonds.begin();b!=this->_Bonds.end(); b++ ) {
    ll << (*b)->getOtherAtom(me);
  }
  return ll.result();
}


CL_LISPIFY_NAME("getBondList");
CL_DEFMETHOD     BondList_sp	Atom_O::getBondList()
{
  VectorBond::const_iterator	b;
  auto  bl  = gctools::GC<BondList_O>::allocate_with_default_constructor();
  for (b=this->_Bonds.begin();b!=this->_Bonds.end(); b++ )
  {
    bl->addBond(*b);
  }
  return bl;
}

uint Atom_O::maxTotalBondOrder()
{
  return chem::maxTotalBondOrderForElement(this->getElement());
}


CL_LISPIFY_NAME("totalBondOrder");
CL_DEFMETHOD     uint Atom_O::totalBondOrder()
{
  // If the hybridization is defined - then it was calculated or
  // set and the bond orders are not reliable - in this case return
  // what the bond order should be based on hybridization and element.
  if (this->getElement() == element_H) return 1;
  if (this->getHybridization() == hybridization_sp3) {
    if (this->_Element==element_C) return 4;
    if (this->_Element==element_N) return 3;
    if (this->_Element==element_O) return 2;
    if (this->_Element==element_F) return 1;
    if (this->_Element==element_Cl) return 1;
    if (this->_Element==element_Br) return 1;
    if (this->_Element==element_I) return 1;
    return 4;
  }
  if (this->getHybridization() == hybridization_sp2) {
    return 3;
  }
  if (this->getHybridization() == hybridization_sp) {
    return 2;
  }
      // Uncertain hybridization - count bond order
  VectorBond::iterator	b;
  core::List_sp		list;
  uint twice = 0;
  for (b=this->_Bonds.begin();b!=this->_Bonds.end(); b++ ) 
  {
    BondOrder bo = (*b)->getRawOrder();
    if ( Bond_O::singleBondP(bo)) twice += 2;
    else if (bo == doubleBond ) twice += 4;
    else if (bo == aromaticBond ) twice += 3;
    else if (bo == tripleBond ) twice += 6;
    else if (bo == dashedSingleBond ) twice += 2;
    else if (bo == dashedDoubleBond ) twice += 4;
    else twice += 2;
  }
  if ( (twice & 1) != 0 )
  {
    SIMPLE_ERROR("The total bond order for {} will not be a whole number", this->description());
  }
  int bondOrder = twice/2;
  if ( bondOrder < 0 || bondOrder > 10 ) {
    SIMPLE_ERROR("The total bond order for {} is an unreasonable value {}" , this->description() , bondOrder );
  }
  return bondOrder;
}


CL_LISPIFY_NAME("bondsAsList");
CL_DEFMETHOD     core::List_sp Atom_O::bondsAsList() const
{
  ql::list ll;
  VectorBond::const_iterator	b;
  for (b=this->_Bonds.begin();b!=this->_Bonds.end(); b++ ) {
    ll << *b;
  }
  return ll.result();
}


bool	Atom_O::isHeavyAtom()
{
  return ((this->_Element!=element_H) && (this->_Element!=element_D) && (this->_Element!=element_T));
}


CL_LISPIFY_NAME("getHeavyAtomBondList");
CL_DEFMETHOD     BondList_sp	Atom_O::getHeavyAtomBondList()
{
  VectorBond::const_iterator	b;
  auto  bl  = gctools::GC<BondList_O>::allocate_with_default_constructor();
  Atom_sp me = this->sharedThis<Atom_O>();
  for (b=this->_Bonds.begin();b!=this->_Bonds.end(); b++ )
  {
    if ( (*b)->getOtherAtom(me)->getElement() != element_H )
    {
      bl->addBond(*b);
    }
  }
  return bl;
}




CL_LISPIFY_NAME("getAtomicNumber");
CL_DEFMETHOD     int     Atom_O::getAtomicNumber()
{
  return atomicNumberForElement(this->getElement());

}

CL_LISPIFY_NAME("getIntegerAtomicMass");
CL_DEFMETHOD     int     Atom_O::getIntegerAtomicMass()
{
  return integerAtomicMassForElement(this->getElement());

}

double	Atom_O::getAtomicWeight()
{
  return atomicWeightForElement(this->getElement());
}


CL_LISPIFY_NAME("getBondedHydrogenCount");
CL_DEFMETHOD     int     Atom_O::getBondedHydrogenCount()
{
  VectorBond::iterator	b;
  int                             hs;
  Atom_sp				weakAtomTo;
  hs=0;
  Atom_sp me = this->sharedThis<Atom_O>();
  for (b=this->_Bonds.begin();b!=this->_Bonds.end(); b++ )
  {
    Atom_sp other = (*b)->getOtherAtom(me);
    if ( other->getElement() == element_H ) hs++;
  }
  return hs;
}


CL_LISPIFY_NAME("inRingSize");
CL_DEFMETHOD     bool    Atom_O::inRingSize(int r) const
{
  bool result = false;
  switch (r)
  {
  case 3:
      result = this->testAnyFlags(in3MemberRing);
      break;
  case 4:
      result = this->testAnyFlags(in4MemberRing);
      break;
  case 5:
      result = this->testAnyFlags(in5MemberRing);
      break;
  case 6:
      result = this->testAnyFlags(in6MemberRing);
      break;
  case 7:
      result = this->testAnyFlags(in7MemberRing);
      break;
  case 8:
      result = this->testAnyFlags(in8MemberRing);
      break;
  default:
      result = false;
      break;
  }
  LOG("Tested ring membership[{}]={} of {}" , r , result , this->description() );
  return result;
}


void Atom_O::clearAllRingMembershipFlags()
{
  this->turnOffFlags(InRing);
}


void Atom_O::setInRingOfSize(int r)
{
  LOG("Setting {} in ring of size[{}]" , this->description() , r );
  switch (r) {
  case 3:
      this->turnOnFlags(in3MemberRing);
      break;
  case 4:
      this->turnOnFlags(in4MemberRing);
      break;
  case 5:
      this->turnOnFlags(in5MemberRing);
      break;
  case 6:
      this->turnOnFlags(in6MemberRing);
      break;
  case 7:
      this->turnOnFlags(in7MemberRing);
      break;
  case 8:
      this->turnOnFlags(in8MemberRing);
      break;
  default:
	    // We ignore larger rings
      break;
  }
}



CL_LISPIFY_NAME("getValence");
CL_DEFMETHOD     int     Atom_O::getValence()
{
  VectorBond::iterator	b;
  BondOrder       bo;
  int valence2 = 0;
  for (b=this->_Bonds.begin();b!=this->_Bonds.end(); b++ ) {
    bo = (*b)->getRawOrder();
    if (Bond_O::singleBondP(bo)) {
      valence2 += 2;
    } else {
      switch ( bo ) {
      case doubleBond:
        valence2 += 4;
        break;
      case aromaticBond:
        valence2 += 3;
        break;
      case tripleBond:
        valence2 += 6;
          break;
      default:
        valence2 += 0;
      }
    }
  }
  return valence2/2;
}




core::List_sp Atom_O::_expandLocalSpanningTree(Atom_sp avoidAtom, Bond_sp incomingBond, uint depth)
{
  core::Cons_sp localTree = core::Cons_O::create(incomingBond,nil<core::T_O>());
  if ( depth <= 0 ) return localTree;
  core::List_sp tail = localTree;
  for ( VectorBond::iterator bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    Atom_sp neighborAtom = (*bi)->getOtherAtom(this->sharedThis<Atom_O>());
    if ( neighborAtom == avoidAtom ) continue;
    core::List_sp one = neighborAtom->_expandLocalSpanningTree(this->sharedThis<Atom_O>(),*bi,depth-1);
    tail.asCons()->setCdr(one);
    tail = one;
  }
  return localTree;
}

CL_LISPIFY_NAME("localSpanningTree");
CL_DEFMETHOD     core::List_sp Atom_O::localSpanningTree(uint depth)
{
  core::Cons_sp localTree = core::Cons_O::create(this->sharedThis<Atom_O>(),nil<core::T_O>());
  if ( depth <= 0 ) return localTree;
  core::List_sp tail = localTree;
  for ( VectorBond::iterator bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    Atom_sp neighborAtom = (*bi)->getOtherAtom(this->sharedThis<Atom_O>());
    core::List_sp one = neighborAtom->_expandLocalSpanningTree(this->sharedThis<Atom_O>(),*bi,depth-1);
    tail.asCons()->setCdr(one);
    tail = one;
  }
  return localTree;
}


CL_LISPIFY_NAME("isConfigurable");
CL_DEFMETHOD     bool Atom_O::isConfigurable()
{
  return this->_StereochemistryType != undefinedCenter;
}

CL_DEFMETHOD void Atom_O::setAbsoluteConfiguration(ConfigurationEnum config, Atom_sp n1, Atom_sp n2, Atom_sp n3)
{
  Bond_sp b1(unbound<Bond_O>());
  Bond_sp b2(unbound<Bond_O>());
  Bond_sp b3(unbound<Bond_O>());
  Bond_sp b4(unbound<Bond_O>());
  if (!(this->_Bonds.size()==4||this->_Bonds.size()==3)) {
    SIMPLE_ERROR("You cannot define the absolute configuration of an atom that has {} bonds - it must have 3 or 4" , this->_Bonds.size());
  }
  for ( int bi=0; bi<this->_Bonds.size(); bi++ ) {
    Bond_sp bc = this->_Bonds[bi];
    Atom_sp other = bc->getOtherAtom(this->asSmartPtr());
    if (other==n1) b1 = bc;
    else if (other==n2) b2 = bc;
    else if (other==n3) b3 = bc;
    else b4 = bc;
  }
  this->_Bonds[0] = b1;
  this->_Bonds[1] = b2;
  this->_Bonds[2] = b3;
  if (!b4.unboundp()) {
    this->_Bonds[3] = b4;
  }
  this->_StereochemistryType = chiralCenter;
  this->setConfiguration(config);
}



void Atom_O::invertStructureAndRestraints()
{
	/*! Invert the Z coordinate */
  Vector3& pos=this->_Position;
  this->_Position.set(pos.getX(),pos.getY(),pos.getZ()*-1.0);
  this->invertStereochemistryOfRestraints();
  ConfigurationEnum config = this->_Configuration;
  ConfigurationEnum newconfig = config;
  if ( config == S_Configuration ) {
    newconfig = R_Configuration;
  } else if ( config == R_Configuration ) {
    newconfig = S_Configuration;
  } else if ( config == RightHanded_Configuration ) {
    newconfig = LeftHanded_Configuration;
  } else if ( config == LeftHanded_Configuration ) {
    newconfig = RightHanded_Configuration;
  }
  this->setConfiguration(newconfig);
}




#if 0
ConfigurationEnum Atom_O::getConfiguration()
{
  return this->_Configuration;
}
void	Atom_O::setConfigurationPriorityHighToLow(string highest, string high, string low, string lowest )
{
  ASSERTP(this->isConfigurable(), "You cannot assign priorities to an unconfigurable atom");
  this->_ConfigurationInfo._HighestPriority = highest;
  this->_ConfigurationInfo._HighPriority = high;
  this->_ConfigurationInfo._LowPriority = low;
  this->_ConfigurationInfo._LowestPriority = lowest;
}

string	Atom_O::getConfigurationPriorityHighest()
{
  ASSERTP(this->isConfigurable(), "You cannot assign priorities to an unconfigurable atom");
  return this->_ConfigurationInfo._HighestPriority;
}

string	Atom_O::getConfigurationPriorityHigh()
{
  ASSERTP(this->isConfigurable(), "You cannot assign priorities to an unconfigurable atom");
  return this->_ConfigurationInfo._HighPriority;
}

string	Atom_O::getConfigurationPriorityLow()
{
  ASSERTP(this->isConfigurable(), "You cannot assign priorities to an unconfigurable atom");
  return this->_ConfigurationInfo._LowPriority;
}

string	Atom_O::getConfigurationPriorityLowest()
{
  ASSERTP(this->isConfigurable(), "You cannot assign priorities to an unconfigurable atom");
  return this->_ConfigurationInfo._LowestPriority;
}
#endif


void Atom_O::_describeAtomRecursively(string prefix, Atom_sp parent, BondOrder order, int maxDepth, stringstream& ss) const
{
  if ( maxDepth == 0 ) return;
  ss << prefix;
  if ( order != noBond )
  {
    ss << bondOrderToString(order);
  } else
  {
    ss << " ";
  }
  ss << this->description() << std::endl;
  for ( VectorBond::const_iterator bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    Atom_sp other = (*bi)->getOtherAtom(this->const_sharedThis<Atom_O>());
    if ( other!=parent)
    {
      other->_describeAtomRecursively(prefix+"  ",this->const_sharedThis<Atom_O>(),(*bi)->getRawOrder(),maxDepth-1,ss);
    }
  }
}


string Atom_O::localEnvironment(int maxDepth) const
{
  stringstream ss;
  this->_describeAtomRecursively("",nil<Atom_O>(),noBond,maxDepth,ss);
  return ss.str();
}









bool	Atom_O::invalid()
{
  IMPLEMENT_ME();
}

CL_LISPIFY_NAME("getType");
CL_DEFMETHOD core::T_sp Atom_O::getType(core::HashTable_sp atomTypes) const {
  return atomTypes->gethash(this->asSmartPtr());
}



AtomIdMap_sp Atom_O::buildAtomIdMap() const
{
  SIMPLE_ERROR("Atom should never buildAtomIdMap");
};


Atom_sp Atom_O::atomWithAtomId(const AtomId& atomId) const
{
  SIMPLE_ERROR("Atom should never return atomWithAtomId");
};




void	Atom_O::bumpPosition(double dist)
{
  double		x, y, z;
  x = core::randomNumber01()*dist*2.0-dist+this->_Position.getX();
  y = core::randomNumber01()*dist*2.0-dist+this->_Position.getY();
  z = core::randomNumber01()*dist*2.0-dist+this->_Position.getZ();
  this->_Position.set(x,y,z);
}


void	Atom_O::failIfInvalid()
{
  if ( this->invalid() ) {
    SIMPLE_ERROR("INVALID {}" , this->description());
  }
}



SYMBOL_EXPORT_SC_(ChemPkg,atomNeedsMinimizer);
SYMBOL_EXPORT_SC_(ChemPkg,selected);
SYMBOL_EXPORT_SC_(ChemPkg,onGrid);
SYMBOL_EXPORT_SC_(ChemPkg,atomFixed);
SYMBOL_EXPORT_SC_(ChemPkg,in3MemberRing);
SYMBOL_EXPORT_SC_(ChemPkg,in4MemberRing);
SYMBOL_EXPORT_SC_(ChemPkg,in5MemberRing);
SYMBOL_EXPORT_SC_(ChemPkg,in6MemberRing);
SYMBOL_EXPORT_SC_(ChemPkg,in7MemberRing);
SYMBOL_EXPORT_SC_(ChemPkg,in8MemberRing);
SYMBOL_EXPORT_SC_(ChemPkg,inRing);
SYMBOL_EXPORT_SC_(ChemPkg,needs_build);
CL_BEGIN_ENUM(AtomFlagEnum,_sym__PLUS_atomFlagSymbolConverter_PLUS_,"AtomFlags");
CL_VALUE_ENUM(_sym_atomNeedsMinimizer, atomNeedsMinimizer);
CL_VALUE_ENUM(_sym_selected, selected);
CL_VALUE_ENUM(_sym_onGrid, onGrid);
CL_VALUE_ENUM(_sym_atomFixed, atomFixed);
CL_VALUE_ENUM(_sym_in3MemberRing, in3MemberRing);
CL_VALUE_ENUM(_sym_in4MemberRing, in4MemberRing);
CL_VALUE_ENUM(_sym_in5MemberRing, in5MemberRing);
CL_VALUE_ENUM(_sym_in6MemberRing, in6MemberRing);
CL_VALUE_ENUM(_sym_in7MemberRing, in7MemberRing);
CL_VALUE_ENUM(_sym_in8MemberRing, in8MemberRing);
CL_VALUE_ENUM(_sym_inRing, inRing);
CL_VALUE_ENUM(_sym_needs_build, needsBuild);
CL_END_ENUM(_sym__PLUS_atomFlagSymbolConverter_PLUS_);


#ifdef USE_TOPOLOGY
ConstitutionAtom_sp Atom_O::asConstitutionAtom(ConstitutionAtomIndex0N index)
{
  ConstitutionAtom_sp ca = makeConstitutionAtom(this->getName(),this->_Element,/*this->_Type,*/ index,this->_StereochemistryType, this->_Properties);
  return ca;
}

void	Atom_O::defineConstitutionAtomBonding(ConstitutionAtom_sp consAtom, MapAtomsToConstitutionAtomIndex0N atomMap)
{
  Atom_O::bondIterator bi;
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    Atom_sp otherAtom = (*bi)->getOtherAtom(me);
    if (atomMap.find(otherAtom) != atomMap.end()) {
    // Only use atoms in the same residue
      ConstitutionBond_sp cb = (*bi)->asConstitutionBond(me,atomMap);
      consAtom->addConstitutionBond(cb);
    }
  }
}
#endif


CL_DEFUN double chem__calculate_dihedral_for_atoms(Atom_sp a, Atom_sp b, Atom_sp c, Atom_sp d) {
  Vector3 va = a->getPosition();
  Vector3 vb = b->getPosition();
  Vector3 vc = c->getPosition();
  Vector3 vd = d->getPosition();
  return geom::calculateDihedral(va,vb,vc,vd);
}

CL_DOCSTRING(R"doc(Fill the x,y,z position of the atom into the nvector starting at index3.
Return the index3 for the next atom)doc");
CL_DEFUN size_t chem__write_position(Atom_sp atom, NVector_sp pos, size_t index3)
{
  if (index3+2<pos->length()) {
    const Vector3& vec = atom->getPosition();
    (*pos)[index3] = vec.getX();
    (*pos)[index3+1] = vec.getY();
    (*pos)[index3+2] = vec.getZ();
    return index3+3;
  }
  SIMPLE_ERROR("The index {} will write the Vector3 out of bounds in nvector of length {}", index3, pos->length());
}


}; // namespace chem
