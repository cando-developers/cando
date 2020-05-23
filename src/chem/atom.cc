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
#include <clasp/core/foundation.h>
#include <clasp/core/common.h>
#include <clasp/core/record.h>
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


size_t globalUniqueAtomOrder = 1;

long	__AtomBuildCounter = 1;







core::NullTerminatedEnumAssociation configurationEnum[] = {
    { "UNDEFINED", undefinedConfiguration },
    { "S", S_Configuration },
    { "R", R_Configuration },
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



#if ATOMIC_ANCHOR
void AnchorRestraint::archive( core::ArchiveP node )
{
  node->archivePlainObject<Vector3>("anchorPos","AnchorPos",this->_AnchorPos);
}
#endif



CL_LISPIFY_NAME(make-atom);
CL_DEFUN Atom_sp Atom_O::make(core::Symbol_sp name, chem::Element element)
{
  GC_ALLOCATE(Atom_O,atom);
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

/*! Return +1 if priority(a)>priority(b)
 *  return -1 if priority(a)<priority(b)
 *  return 0 if priority(a) == priority(b)
 *  Priority is decided by comparing RelativePriority and if they are the same
 *  by comparing names.
 */
int Atom_O::priorityOrder(Atom_sp a, Atom_sp b)
{
  if ( a.nilp() && b.nilp() ) return 0;
  if ( b.nilp() ) return 1;
  if ( a.nilp() ) return -1;
  if ( a->getRelativePriority()>b->getRelativePriority() ) return 1;
  if ( a->getRelativePriority()<b->getRelativePriority() ) return -1;
  if ( a->getName()->symbolNameAsString() > b->getName()->symbolNameAsString() ) return 1;
  if ( a->getName()->symbolNameAsString() < b->getName()->symbolNameAsString() ) return -1;
  return 0;
}

CL_DOCSTRING(R"doc(Calculate the stereochemical configuration of the atom.
This requires that relative cip priorities are defined using CipPrioritizer_O::assignCahnIngoldPrelogPriorityToAtomsRelativePriority.)doc");
CL_LISPIFY_NAME("calculateStereochemicalConfiguration");
CL_DEFMETHOD     ConfigurationEnum Atom_O::calculateStereochemicalConfiguration()
{
  if ( this->numberOfBonds() != 4 ) return undefinedConfiguration;
  core::List_sp neighborsByPriority = this->getNeighborsByRelativePriority();
  Atom_sp a4 = core::oCar(neighborsByPriority).as<Atom_O>();
  Atom_sp a3 = core::oCadr(neighborsByPriority).as<Atom_O>();
  Atom_sp a2 = core::oCaddr(neighborsByPriority).as<Atom_O>();
  Atom_sp a1 = core::oCadddr(neighborsByPriority).as<Atom_O>();
  Vector3 vme = this->getPosition();
  Vector3 v4 = a4->getPosition().sub(vme);
  Vector3 v3 = a3->getPosition().sub(vme);
  Vector3 v2 = a2->getPosition().sub(vme);
  Vector3 v1 = a1->getPosition().sub(vme);
  double dir2 = v2.dotProduct(v4.crossProduct(v3)); // For R this should be >0
  double dir1 = v1.dotProduct(v4.crossProduct(v3)); // For R this should be <0
  if ( dir1 > 0.0 ) {
    if ( dir2 < 0.0 ) {
      return S_Configuration;
    }
  } else if ( dir1 < 0.0 ) {
    if ( dir2 > 0.0 ) {
      return R_Configuration;
    }
  }
  return undefinedConfiguration;
}


string Atom_O::calculateStereochemicalConfigurationAsString()
{
  string s;
  ConfigurationEnum config = this->calculateStereochemicalConfiguration();
  switch ( config )
  {
  case S_Configuration:
      s = "S";
      break;
  case R_Configuration:
      s = "R";
      break;
  default:
      s = "undefinedConfiguration";
      break;
  }
  return s;
}



//
// Constructor
//

#if 0
void	Atom_O::initialize()
{
  this->Base::initialize();
  this->copyAtom = _Nil<Atom_O>();
  this->flags = ATOM_NEEDS_MINIMIZER;
  this->_Mask = 0;
  this->tempInt = 0;
  this->_Element = element_Undefined;
  this->_Hybridization = hybridization_undefined;
  this->_Alias = _Nil<core::Symbol_O>();
  this->_Ionization = 0; // neutral
#if ATOMIC_ANCHOR
  this->_AnchorRestraint.setIsDefined(false);
  this->_AnchorRestraint._AnchorPos.set(0.0,0.0,0.0);
#endif
  this->charge= 0;
  this->_RingMembershipCount = 0;
  this->type= _Nil<core::Symbol_O>();
  this->invalidateBackSpan();
  this->invalidateNextSpan();
  this->_RelativePriority = 0;
  this->_Configuration = undefinedConfiguration;
  this->_StereochemistryType = undefinedCenter;
}
#endif

CL_DEFMETHOD bool Atom_O::Atom_equal(core::T_sp obj) const
{
  if ( this->eq(obj) ) goto T;
  if ( obj.isA<Atom_O>() )
  {
    Atom_sp other = obj.as<Atom_O>();
    if ( other->getName() != this->getName() )
    {
      LOG(BF("Atom names[ this->getName()=%s  other->getName()=%s ]  don't match!") % this->getName() % other->getName() );
      goto F;
    }
    if ( other->numberOfBonds() != this->numberOfBonds() )
    {
      LOG(BF("Numbers of bonds don't match this->numberOfBonds()=%d  other->numberOfBonds()=%d!") % this->numberOfBonds() % other->numberOfBonds() );
      goto F;
    }
    goto T;
  }
  LOG(BF("The other object is not an ATOM!!!!"));
 F:
  LOG(BF("atoms this(%s) =equal= other(%s)  equal==false") % this->description() % obj->description() );
  return false;
 T:
  LOG(BF("atoms this(%s) =equal= other(%s)  equal==true") % this->description() % obj->description() );
  return true;
}

void Atom_O::transferCoordinates(Matter_sp obj)
{_OF();
  if ( !this->Atom_equal(obj) )
  {
    SIMPLE_ERROR(BF("This atom(%s) is not equal to %s so you cannot transfer coordinates")%this->description() % obj->description() );
  }
  Atom_sp other = obj.as<Atom_O>();
  this->setPosition(other->getPosition());
}

bool Atom_O::atomWithinAngstroms(Atom_sp other, float angstroms) const
{
  Vector3 delta = this->position-other->position;
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
  this->type = ss.type;
  this->_Hybridization = ss._Hybridization;
  this->flags = ss.flags;
  this->_Mask = ss._Mask;
  this->position = ss.position;
#if ATOMIC_ANCHOR
  this->_AnchorRestraint = ss._AnchorRestraint;
#endif
  this->charge = ss.charge;
  this->_RingMembershipCount = ss._RingMembershipCount;
  this->occupancy = ss.occupancy;
  this->tempFactor = ss.tempFactor;
  this->vdwRadius = ss.vdwRadius;
  this->covalentRadius = ss.covalentRadius;
  this->copyAtom = Atom_sp();
  this->tempInt = ss.tempInt;
//  this->moeIndex = ss.moeIndex;
//  this->moeType = ss.moeType;
  this->_Ionization = ss._Ionization;
  this->_RelativePriority = ss._RelativePriority;
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
{_OF();
  this->_Configuration = conf;
  LOG(BF("Changed configuration of atom[%s] to [%s]")
      % this->__repr__() % this->getConfigurationAsString() );
}

void Atom_O::setElementFromAtomName()
{
  string nameStr = this->getName()->symbolName()->get_std_string();
  Element element = elementFromAtomNameString(nameStr);
  LOG(BF(" Resulting element= |%d|") %element );
  this->setElement(element);
}

void Atom_O::setElementFromString(const string& str)
{
  Element element = elementFromAtomNameStringCaseInsensitive(str);
  LOG(BF(" Resulting element= |%d|") %element );
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
  return this->bonds[i];
}

CL_LISPIFY_NAME("bondedNeighbor");
CL_DEFMETHOD     Atom_sp Atom_O::bondedNeighbor(int i)
{
  return this->bonds[i]->getOtherAtom(this->sharedThis<Atom_O>());
}

CL_LISPIFY_NAME("bondedOrder");
CL_DEFMETHOD     BondOrder Atom_O::bondedOrder(int i) 
{
  return this->bonds[i]->getOrder();
};

CL_LISPIFY_NAME("numberOfBonds");
CL_DEFMETHOD     int	Atom_O::numberOfBonds() const
{
  return this->bonds.size();
};

CL_LISPIFY_NAME("coordination");
CL_DEFMETHOD     int Atom_O::coordination() 
{
  return this->bonds.size();
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
                   core::Cons_O::create(core::clasp_make_fixnum(this->getIonization()),_Nil<core::T_O>()));
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
    SIMPLE_ERROR(BF("There are too many hydrogens to be added for element %s totalBondOrder: %u addHydrogens: %u") % _rep_(symbolFromElement(this->getElement())) % totalBondOrder % addHydrogens );
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
  core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::T_O>());
  core::Cons_sp cons = first;
  if ( addHydrogens == 1 )
  {
    MatterName hname = chemkw_intern("H"+nameSuffix);
    core::Cons_sp one = core::Cons_O::create(hname,_Nil<core::T_O>());
    cons->setCdr(one);
    cons = one;
    return oCdr(first);
  }
  for ( uint i=0; i<addHydrogens; i++ )
  {
    stringstream ss;
    ss << "H" << nameSuffix << (i+1);
    MatterName hname = chemkw_intern(ss.str());
    core::Cons_sp one = core::Cons_O::create(hname,_Nil<core::T_O>());
    cons->setCdr(one);
    cons = one;
  }
  return oCdr(first);
}




void Atom_O::fillInImplicitHydrogensWithResidue(Residue_sp residue)
{
  core::List_sp names = this->createImplicitHydrogenNames();
  if ( names.nilp() ) return;
  for ( auto cur : names ) {
    this->_addHydrogenWithName(residue,oCar(cur).as<MatterName::Type>());
  }
}

void Atom_O::fillInImplicitHydrogens()
{
  SIMPLE_ERROR(BF("You cannot invoke chem:fill-in-implicit with an atom - any new hydrogens must be added to the residue that contains this atom"));
}

//
//	randomizeAtomPosition
//
//	Set the coordinates of all atoms to random values
//		between 0 and 100 angstroms
//
//



gc::Nilable<Atom_sp> Atom_O::highestPriorityNeighborThatIsnt(gc::Nilable<Atom_sp> avoid)
{
  Atom_sp bestAtom = _Nil<Atom_O>();
  Atom_sp atom;
  VectorBond::iterator	b;
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) {
    atom = (*b)->getOtherAtom(this->sharedThis<Atom_O>());
    if ( avoid.notnilp() && atom == avoid ) continue;
    if ( priorityOrder(atom,bestAtom)>0 )
    {
      bestAtom = atom;
    }
  }
  return bestAtom;
}
    

gc::Nilable<Atom_sp> Atom_O::lowestPriorityNeighborThatIsnt(gc::Nilable<Atom_sp> avoid)
{
  gc::Nilable<Atom_sp> bestAtom = _Nil<core::T_O>();
  Atom_sp atom;
  VectorBond::iterator	b;
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) {
    atom = (*b)->getOtherAtom(this->sharedThis<Atom_O>());
    if ( avoid.notnilp() && atom == avoid ) continue;
    if ( bestAtom.nilp() ) {
      bestAtom = atom;
      continue;
    }
    if ( priorityOrder(atom,bestAtom)<0 ) bestAtom = atom;
  }
  return bestAtom;
}



CL_LISPIFY_NAME("getBondTo")
CL_DEFMETHOD Bond_sp Atom_O::getBondTo(Atom_sp a)
{_OF();
  VectorBond::iterator	b;
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ )
  {
    if ( a == (*b)->getOtherAtom(this->sharedThis<Atom_O>()) ) return *b;
  }
  SIMPLE_ERROR(BF("Could not find bond"));
}

void Atom_O::_addExistingBond(Bond_sp const& bond)
{_OF();
  this->bonds.push_back(bond);
}

void Atom_O::addBond(Bond_sp bond)
{_OF();
  this->bonds.push_back(bond);
}

CL_LISPIFY_NAME("bondTo");
CL_DEFMETHOD     Bond_sp Atom_O::bondTo( Atom_sp to, BondOrder o )
{_OF();
	// Check if there is already a bond to this atom and
	// throw an exception if there is
  Atom_sp from = this->sharedThis<Atom_O>();
  VectorBond::iterator	b;
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ )
  {
    if ( (*b)->getOtherAtom(from) == to )
    {
      SIMPLE_ERROR(BF("You tried to form a bond from[%s]-to[%s] but there is already one there!!") 
                   % this->__repr__() % _rep_(to) );
    }
  }
  Bond_sp bn = Bond_O::create(from,to,o);
  this->bonds.push_back(bn);
  to->bonds.push_back(bn);
  if (this->_Element == element_C && this->bonds.size()>4) {
    SIMPLE_WARN(BF("More than four bonds to carbon will be made to %s") % _rep_(this->asSmartPtr()));
  }
  if (to->_Element == element_C && to->bonds.size()>4 ) {
    SIMPLE_WARN(BF("More than four bonds to carbon will be made to %s") % _rep_(to));
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
{_OF();
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
public:
  bool operator() ( Atom_sp x, Atom_sp y )
  {
    if ( x->getRelativePriority()<y->getRelativePriority() ) return true;
    if ( x->getRelativePriority()>y->getRelativePriority() ) return false;
    if ( x->getName()->symbolNameAsString()<y->getName()->symbolNameAsString() ) return true;
    return false;
  }
};

core::List_sp	Atom_O::getNeighborsByRelativePriority()
{
  LOG(BF("Ordering neighbors around(%s) by priority and name") % this->getName()  );
  VectorAtom	reversedNeighbors;
  VectorBond::iterator	b;
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) {
    reversedNeighbors.push_back((*b)->getOtherAtom(this->sharedThis<Atom_O>()));
  }
  OrderByPriorityAndName orderer;
  core::List_sp ncons = _Nil<core::T_O>();
  sort::quickSort(reversedNeighbors.begin(),reversedNeighbors.end(),orderer);
	// At this point the sorted in reverse order
	// 
	// Now copy them into a Cons list backwards to get the
	// proper order
  LOG(BF("The sorted neighbors in REVERSE order is:") );
  for ( VectorAtom::iterator ni=reversedNeighbors.begin();
        ni!=reversedNeighbors.end(); ni++ )
  {
    LOG(BF("    neighbor priority(%u) name(%s)") % (*ni)->getRelativePriority() % (*ni)->getName()  );
    ncons = core::Cons_O::create(*ni,ncons);
  }
  return ncons;
}





void	Atom_O::modifyFlags(int op, ATOM_FLAGS fl ) {
  if ( op == ATOMFLAG_OR ) {
    this->flags |= fl;
  } else if ( op==ATOMFLAG_AND ) {
    this->flags &= fl;
  } else if ( op==ATOMFLAG_ON ) {
    this->flags |= fl;
  } else if ( op==ATOMFLAG_OFF ) {
    this->flags &= (~fl);
  }
}




CL_LISPIFY_NAME("isBondedToAtomNamed");
CL_DEFMETHOD     bool	Atom_O::isBondedToAtomNamed(core::Symbol_sp name)
{
  VectorBond::iterator	b;

  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) {
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
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) {
    if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>())->getName() == name ) {
      return (*b)->getOtherAtom(this->sharedThis<Atom_O>());
    }
  }
  Atom_sp a = _Nil<Atom_O>();
  return a;
}


CL_LISPIFY_NAME("isBondedToElementOrder");
CL_DEFMETHOD     bool	Atom_O::isBondedToElementOrder(Element el, BondOrder o)
{
  VectorBond::iterator	b;

  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) 
  {
    if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>())->getElement() == el ) 
    {
      if ( (*b)->getOrder() == o ) return true;
    }
  }
  return false;
}


CL_LISPIFY_NAME("hasBondWithOrder");
CL_DEFMETHOD     bool Atom_O::hasBondWithOrder(BondOrder o) const
{
  VectorBond::const_iterator	b;
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) 
  {
    if ( (*b)->getOrder() == o) return true;
  }
  return false;
}



CL_LISPIFY_NAME("isBondedToElementHybridization");
CL_DEFMETHOD     bool	Atom_O::isBondedToElementHybridization(Element el, Hybridization hy)
{
  VectorBond::iterator	b;

  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) 
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
{_OF();
  IMPLEMENT_ME(); // Handle new way of handling bonds
#if 0
  VectorBond::iterator	b;
  VectorBond::iterator	b2;
  Atom_sp	a1;
  Atom_sp	a2;
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ ) 
  {
    a1 = Atom_sp((*b)->getTo());
    if ( a1->getElement() == el1 ) 
    {
      if ( a1->getHybridization() == hy1 ) 
      {
        for ( b2=a1->bonds.begin();b2!=a1->bonds.end();b2++) 
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
  Vector3 delta = this->position-other->position;
  return delta.getX()*delta.getX()+delta.getY()*delta.getY()+delta.getZ()*delta.getZ();
}

CL_DOCSTRING(R"doc(Return the distance squared between twoatoms. This provides a way to compare atom positions without allocating memory.)doc");
CL_DEFUN
float chem__distance_squared_between_two_atoms(Atom_sp atom1, Atom_sp atom2)
{
  return atom1->distanceSquaredToAtom(atom2);
}


void Atom_O::makeAllAtomNamesInEachResidueUnique()
{
  SIMPLE_ERROR(BF("This should never be called for an Atom_sp"));
}


void Atom_O::basicRemoveBondTo(Atom_sp a)
{_OF();
  VectorBond::iterator	b;
  Bond_sp				bb;
  Atom_sp				aa;
  Atom_sp				wpa;
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ )
  {
    if ( (*b)->getOtherAtom(this->sharedThis<Atom_O>()) == a )
    {
      this->bonds.erase(b);
      return;
    }
  }
  stringstream ss;
  ss << "Trying to remove bond from atom(" << this->description()
     << ") but I can't find the to atom(" << a->description() << ")";
  SIMPLE_ERROR(BF("%s")%ss.str());
}

CL_LISPIFY_NAME("removeBondTo");
CL_DEFMETHOD     void	Atom_O::removeBondTo(Atom_sp a)
{_OF();
  Atom_sp	atemp;
  LOG(BF("Bond_O::removeBondTo") );
  this->basicRemoveBondTo(a);
  atemp = this->sharedThis<Atom_O>();
  a->basicRemoveBondTo(atemp);
}

CL_LISPIFY_NAME("removeAllBonds");
CL_DEFMETHOD     void	Atom_O::removeAllBonds()
{_OF();
  VectorBond::iterator	b;
  Atom_sp				atemp;
	// Remove back bond
  LOG(BF("Bond_O::removeAllBonds for %s") % this->description() );
  atemp = this->sharedThis<Atom_O>();
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ )
  {
    Atom_sp ato = (*b)->getOtherAtom(atemp);
    ato->basicRemoveBondTo(atemp);
  }
  this->bonds.clear();
}



//
//	applyTransformToAtoms
//
CL_LISPIFY_NAME("applyTransformToAtoms");
CL_DEFMETHOD     void	Atom_O::applyTransformToAtoms( const Matrix& m )
{
  this->position = (m)*this->position;
  if ( this->_RestraintList.notnilp() ) {
    this->Base::applyTransformToRestraints(m);
  }
}





string Atom_O::getHybridizationAsString()
{
  return hybridizationSymbolFromHybridization(this->getHybridization())->symbolName()->get_std_string();
}

CL_LISPIFY_NAME("getElementAsString");
CL_DEFMETHOD     string Atom_O::getElementAsString()
{
  return atomicSymbolFromElement(this->getElement())->symbolName()->get_std_string();
}

CL_DEFMETHOD core::Symbol_sp Atom_O::getElementAsSymbol() const
{
  return symbolFromElement(this->getElement());
}


#ifdef OLD_SERIALIZE
void	Atom_O::serialize(serialize::SNode snode)
{_OF();
  this->Base::serialize(snode);
  snode->attributeIfNotDefault<ATOM_FLAGS>( "flags", this->flags, (unsigned long)0 );
  snode->attributeSymbolEnumHiddenConverter( "element", this->_Element, _sym__PLUS_elementToSymbolConverter_PLUS_);
  snode->attributeSymbolEnumHiddenConverterIfNotDefault( "hybrid", this->_Hybridization, _sym__PLUS_hybridizationToSymbolConverter_PLUS_,hybridization_sp3 );
  snode->attributeIfNotNil( "alias", this->_Alias );
  snode->attributeIfNotDefault<uint>( "priority", this->_RelativePriority, 0 );
  snode->attributeIfNotDefault<double>( "chg", this->charge, 0.0 );
  snode->attributeSymbolEnumHiddenConverterIfNotDefault( "configuration", this->_Configuration, _sym__PLUS_configurationEnumConverter_PLUS_, undefinedConfiguration  );
  snode->attributeSymbolEnumHiddenConverterIfNotDefault( "stereochemistryType", this->_StereochemistryType, _sym__PLUS_stereochemistryTypeConverter_PLUS_, undefinedCenter );
  snode->attributeIfNotDefault<int>( "ion", this->_Ionization, 0 );
  snode->attributeIfNotDefault<int>( "rings", this->_RingMembershipCount, 0 );
  snode->attributeIfNotDefault<int>( "tempInt", this->tempInt, 0 );
  snode->attributeIfNotDefault<string>( "type", this->typeString, "" );
  snode->attributeIfNotDefault<uint>( "mask", this->_Mask, (unsigned int)(0) );
  snode->archivePlainObjectIfDefined<Vector3>( "pos","Vector3",
                                               this->position.isDefined(), this->position );
  LOG(BF("After pos archived Atom position@%p = %s") % &(this->position) % this->position.asString()  );
  LOG(BF("Atom position = %s") % this->position.asString() );
}
#endif





void	Atom_O::fields(core::Record_sp node)
{
  node->field_if_not_nil( INTERN_(kw,alias), this->_Alias );
  node->/*pod_*/field_if_not_default<size_t>( INTERN_(kw,flags), this->flags, 0 );
  node->/*pod_*/field( INTERN_(kw,element), this->_Element);
  node->/*pod_*/field_if_not_default( INTERN_(kw,hybridization), this->_Hybridization,hybridization_sp3 );
  node->/*pod_*/field_if_not_default( INTERN_(kw,priority), this->_RelativePriority, 0 );
  node->/*pod_*/field_if_not_default( INTERN_(kw,chg), this->charge, 0.0 );
  node->/*pod_*/field_if_not_default<short>( INTERN_(kw,ion), this->_Ionization, 0 );
  node->field_if_not_nil( INTERN_(kw,type), this->type);
  node->/*pod_*/field_if_not_default<ushort>( INTERN_(kw,rings), this->_RingMembershipCount, 0 );
  node->/*pod_*/field_if_not_default( INTERN_(kw,mask), this->_Mask, (uint)(0) );
  node->/*pod_*/field_if_not_default( INTERN_(kw,configuration), this->_Configuration, undefinedConfiguration  );
  node->/*pod_*/field_if_not_default( INTERN_(kw,stereochemistryType), this->_StereochemistryType, undefinedCenter );
  node->/*pod_*/field_if_not_default( INTERN_(kw,pos), this->position, Vector3());
  node->field_if_not_empty(INTERN_(kw,bonds),this->bonds);
  this->Base::fields(node);
}

bool Atom_O::applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value ) {
  if ( this->Base::applyPropertyToSlot(prop,value) ) return true;
  return false;
}
CL_LISPIFY_NAME("setPositionInNanometers");
CL_DEFMETHOD     void Atom_O::setPositionInNanometers(Vector3 o)
{_OF();
  Vector3 angpos = o.multiplyByScalar(10.0);
  this->position = angpos;
}


CL_LISPIFY_NAME("getConfigurationAsString");
CL_DEFMETHOD     string	Atom_O::getConfigurationAsString()
{
  string config, stereo;
  if ( this->_Configuration == undefinedConfiguration )
  {
    config = "UndefinedConfiguration";
  } else if ( this->_Configuration == S_Configuration )
  {
    config = "(S)";
  } else if ( this->_Configuration == R_Configuration )
  {
    config = "(R)";
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






//
// writeUniqueBondsToStream
//
// Stream the atoms bonds description to the output
//
#if 0
void Atom_O::writeUniqueBondsToStream(string prefix, std::ostream& out )
{
  VectorBond	bonds;
  VectorBond::iterator	b;
  bonds = this->getBonds();
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=bonds.begin();b!=bonds.end() ; b++ ) {
    if ( me->atLowerUniqueAtomOrderThan((*b)->getOtherAtom(me)))
    {
      out << prefix << "bond ";
      switch ( (*b)->getOrder() ) {
      case singleBond:
          out << "single";
          break;
      case doubleBond:
          out << "double";
          break;
      case tripleBond:
          out << "triple";
          break;
      case aromaticBond:
          out << "aromatic";
          break;
      default:
          out << "--unknown--";
      }
      out << " " << me->getId() << " ";
      out << (*b)->getOtherAtom(me)->getId() << " ;" << std::endl;
    }
  }
}
#endif



string	Atom_O::__repr__() const
{
  stringstream ss;
  ss << "#<" << this->className() << " " << this->name << "/" << _rep_(symbolFromElement(this->_Element)) << " " << this->getId() << ">";
  return ss.str();
}

string	Atom_O::description() const
{
  stringstream ss;
  ss << this->className() << "("<<this->getName();
  if ( this->containerContainedBy.notnilp() ) {
    ss <<"["<<this->getResidueContainedBy_const()->description()<<"]";
  } else {
    ss << "[residue-NULL-]";
  }
  ss << "bonds[" << this->numberOfBonds() << "]";
  ss <<")";
  return ss.str();
}


/*!
 * Add all unique intra residue bonds to the BondList_sp
 */
void	Atom_O::addUniqueIntraResidueBondCopiesToBondList(BondList_sp list)
{_OF();
  VectorBond			bonds;
  VectorBond::iterator	b;
  ASSERTNOTNULL(list);
  bonds = this->getBonds();
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=bonds.begin();b!=bonds.end() ; b++ ) {
    if ( me->atLowerUniqueAtomOrderThan( (*b)->getOtherAtom(me))) {
      if ( !(*b)->isInterResidueBond() ) {
        LOG(BF("Original bond: %s") % (*b)->description() );
        GC_COPY(Bond_O,bondCopy,*(b->get()));
        LOG(BF("Copy bond: %s") % bondCopy->description() );
        list->addBond(bondCopy);
      }
    }
  }
}

/*!
 * Add all unique inter residue bonds to the BondList_sp
 */
void Atom_O::addUniqueInterResidueBondCopiesToBondList(BondList_sp list)
{_OF();
  VectorBond			bonds;
  VectorBond::iterator	b;
  ASSERTNOTNULL(list);
  bonds = this->getBonds();
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=bonds.begin();b!=bonds.end() ; b++ ) {
    if ( me->atLowerUniqueAtomOrderThan( (*b)->getOtherAtom(me))) {
      if ( (*b)->isInterResidueBond() ) {
        GC_COPY(Bond_O, bondCopy, *(b->get()) ); // = RP_Copy<Bond_O>(*b);
        list->addBond(bondCopy);
      }
    }
  }
}

void Atom_O::addInterResidueBondsToBondList(BondList_sp bondlist)
{
  VectorBond::iterator	b;
  bonds = this->getBonds();
  for ( b=bonds.begin();b!=bonds.end() ; b++ ) {
    if ( (*b)->isInterResidueBond() ) {
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
{_OF();
  VectorBond::iterator	b;
  Atom_sp				a2;
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ )
  {
    a2 = (*b)->getOtherAtom(me);
    if ( a2 == aTarget ) {
      if ( o == (*b)->getOrder() ) {
        return true;
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
  LOG(BF("Looking at atom(%s) for bonded atom(%s)")
      % this->getName()
      % aTarget->getName() );
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ )
  {
    a2 = (*b)->getOtherAtom(me);
    LOG(BF("  looking at atom(%s)") % a2->getName() );
    if ( a2 == aTarget )
    {
      LOG(BF("Found match") );
      return true;
    }
  }
  LOG(BF("No match") );
  return false;
}
//
// bondOrderTo
//
// Return the bond order by which we are bonded to the atom.
// Return noBond if there is no bond.
CL_DOCSTRING("Return the bond-order to the OTHER atom.");
CL_DEFMETHOD
BondOrder	Atom_O::bondOrderTo( Atom_sp other )
{
  VectorBond::iterator	b;
  Atom_sp				a2;
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ )
  {
    a2 = (*b)->getOtherAtom(me);
    if ( a2 == other )
    {
      return (*b)->getOrder();
    }
  }
  return noBond;
}


//
// testConsistancy
//
//	Make sure that all of the bonds are consistant
CL_LISPIFY_NAME("testConsistancy");
CL_DEFMETHOD     bool Atom_O::testConsistancy(Matter_sp parentShouldBe)
{_OF();
  VectorBond::iterator	b;
  Atom_sp				a2,atemp;
  LOG(BF("bonds array start = 0x%08x") % ((void*)&(this->bonds[0]))  );
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( b=this->bonds.begin();b!=this->bonds.end() ; b++ )
  {
    a2 = (*b)->getOtherAtom(me);
    if ( !(a2->isBondedToWithBondOrder(me,(*b)->getOrder())) )
    {
      _lisp->print(BF( "Atom_O::testConsistancy failed" ));
      return false;
    }
  }
  return true;
}



Matter_sp Atom_O::copy(core::T_sp new_to_old)
{_OF();
  SIMPLE_ERROR(BF("Don't copy single atoms - bonds will be messed up"));
}




//
//	copy
//
//	Allocate a copy of this atom and return it.
//	Keep track of the new atom so that we can redirect pointers to the copy.
//
Matter_sp Atom_O::copyDontRedirectAtoms(core::T_sp new_to_old)
{_OF();
  LOG(BF("Copying atom @%p") % this );
  Atom_sp myself = this->sharedThis<Atom_O>();
  GC_COPY(Atom_O, aNew , *this); // = RP_Copy<Atom_O>(this);
#if 1
  aNew->bonds.clear();
#else
  printf("%s:%d copyDontRedirectAtom original: %p  copy %p\n", __FILE__, __LINE__, myself.raw_(), aNew.raw_());
  for ( VectorBond::const_iterator b=bonds.begin();b!=bonds.end() ; b++ ) {
	    /* Only copy Bond object if we are the Atom1
	     This is to ensure that bonds aren't duplicated*/
    if ( (*b)->isAtom1(myself) ) {
      Bond_sp bNew = (*b)->copyDontRedirectAtoms();
      printf("%s:%d copyDontRedirectAtom newBond between atoms: %p - %p\n", __FILE__, __LINE__, bNew->_Atom1.raw_(), bNew->_Atom2.raw_());
      aNew->bonds.push_back(bNew);
    }
  }
#endif
  this->copyAtom = aNew;
  aNew->copyRestraintsDontRedirectAtoms(this->asSmartPtr());
  LOG(BF("    copy atom== %s") % aNew->description());
  return(aNew);
}



//
//	redirectAtom
//
//	Redirect the bonds of this atom to their copies.
void Atom_O::redirectAtoms()
{
#if 0
  VectorBond::iterator	b;
  LOG(BF("Redirecting bonds for %s") % this->description());
  printf("%s:%d redirecting atoms for %zu bonds\n", __FILE__, __LINE__, bonds.size());
  Atom_sp myself = this->sharedThis<Atom_O>();
  for ( b=bonds.begin();b!=bonds.end() ; b++ ) {
    if ( (*b)->isAtom1(myself->getCopyAtom()) ) {
      printf("%s:%d Redirecting atoms for bond: %zu\n", __FILE__, __LINE__, (b-bonds.begin()));
      (*b)->redirectAtoms();
    }
  }
#endif
  this->redirectRestraintAtoms();
}




VectorAtom	Atom_O::getBondedAtoms()
{
  VectorBond::iterator	b;
  VectorAtom			vAtoms;
  Atom_sp me = this->sharedThis<Atom_O>();
  for (b=this->bonds.begin();b!=this->bonds.end(); b++ )
  {
    vAtoms.push_back((*b)->getOtherAtom(me));
  }
  return vAtoms;
}

CL_LISPIFY_NAME("bondedAtomsAsList");
CL_DEFMETHOD     core::List_sp Atom_O::bondedAtomsAsList()
{
  core::List_sp list = _Nil<core::T_O>();
  VectorBond::iterator	b;
  Atom_sp me = this->sharedThis<Atom_O>();
  for (b=this->bonds.begin();b!=this->bonds.end(); b++ )
  {
    list = core::Cons_O::create((*b)->getOtherAtom(me),list);
  }
  return list;
}

void	advanceAtomBuildCounter()
{
  __AtomBuildCounter++;
  if (__AtomBuildCounter==0) __AtomBuildCounter++;
};

long	getAtomBuildCounter()
{
  return __AtomBuildCounter;
};


CL_LISPIFY_NAME("getBondList");
CL_DEFMETHOD     BondList_sp	Atom_O::getBondList()
{
  VectorBond::iterator	b;
  GC_ALLOCATE(BondList_O, bl );
  for (b=this->bonds.begin();b!=this->bonds.end(); b++ )
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
  for (b=this->bonds.begin();b!=this->bonds.end(); b++ ) 
  {
    if ( (*b)->getOrder() == singleBond ) twice += 2;
    else if ((*b)->getOrder() == doubleBond ) twice += 4;
    else if ((*b)->getOrder() == aromaticBond ) twice += 3;
    else if ((*b)->getOrder() == tripleBond ) twice += 6;
    else if ((*b)->getOrder() == dashedSingleBond ) twice += 2;
    else if ((*b)->getOrder() == dashedDoubleBond ) twice += 4;
    else twice += 2;
  }
  if ( (twice & 1) != 0 )
  {
    SIMPLE_ERROR(BF("The total bond order for "+this->description()+" will not be a whole number"));
  }
  int bondOrder = twice/2;
  if ( bondOrder < 0 || bondOrder > 10 ) {
    SIMPLE_ERROR(BF("The total bond order for %s is an unreasonable value %d") % this->description() % bondOrder );
  }
  return bondOrder;
}


CL_LISPIFY_NAME("bondsAsList");
CL_DEFMETHOD     core::List_sp Atom_O::bondsAsList()
{
  VectorBond::iterator	b;
  core::List_sp		list;
  list = _Nil<core::T_O>();
  for (b=this->bonds.begin();b!=this->bonds.end(); b++ ) {
    list = core::Cons_O::create(*b,list);
  }
  return list;
}

bool	Atom_O::isHeavyAtom()
{
  return ((this->_Element!=element_H) && (this->_Element!=element_D) && (this->_Element!=element_T));
}


CL_LISPIFY_NAME("getHeavyAtomBondList");
CL_DEFMETHOD     BondList_sp	Atom_O::getHeavyAtomBondList()
{
  VectorBond::iterator	b;
  GC_ALLOCATE(BondList_O, bl );
  Atom_sp me = this->sharedThis<Atom_O>();
  for (b=this->bonds.begin();b!=this->bonds.end(); b++ )
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
{_OF();
  VectorBond::iterator	b;
  int                             hs;
  Atom_sp				weakAtomTo;
  hs=0;
  Atom_sp me = this->sharedThis<Atom_O>();
  for (b=this->bonds.begin();b!=this->bonds.end(); b++ )
  {
    Atom_sp other = (*b)->getOtherAtom(me);
    if ( other->getElement() == element_H ) hs++;
  }
  return hs;
}


CL_LISPIFY_NAME("inRingSize");
CL_DEFMETHOD     bool    Atom_O::inRingSize(int r) const
{_OF();
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
  LOG(BF("Tested ring membership[%s]=%d of %s") % r % result % this->description() );
  return result;
}


void Atom_O::clearAllRingMembershipFlags()
{
  this->turnOffFlags(InRing);
}


void Atom_O::setInRingOfSize(int r)
{_OF();
  LOG(BF("Setting %s in ring of size[%d]") % this->description() % r );
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
  int             valence;
  BondOrder       bo;
  valence = 0;
  for (b=this->bonds.begin();b!=this->bonds.end(); b++ )
  {
    bo = (*b)->getOrder();
    switch ( bo ) {
    case singleBond:
        valence += 1;
        break;
    case doubleBond:
    case aromaticBond:
        valence += 2;
        break;
    case tripleBond:
        valence += 3;
        break;
    default:
        valence += 0;
    }
  }
  return valence;
}




core::List_sp Atom_O::_expandLocalSpanningTree(Atom_sp avoidAtom, Bond_sp incomingBond, uint depth)
{
  core::Cons_sp localTree = core::Cons_O::create(incomingBond,_Nil<core::T_O>());
  if ( depth <= 0 ) return localTree;
  core::List_sp tail = localTree;
  for ( VectorBond::iterator bi=this->bonds.begin(); bi!=this->bonds.end(); bi++ )
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
  core::Cons_sp localTree = core::Cons_O::create(this->sharedThis<Atom_O>(),_Nil<core::T_O>());
  if ( depth <= 0 ) return localTree;
  core::List_sp tail = localTree;
  for ( VectorBond::iterator bi=this->bonds.begin(); bi!=this->bonds.end(); bi++ )
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


void Atom_O::invertStructureAndRestraints()
{
	/*! Invert the Z coordinate */
  Vector3& pos=this->position;
  this->position.set(pos.getX(),pos.getY(),pos.getZ()*-1.0);
#if ATOMIC_ANCHOR
  Vector3& anchor = this->_AnchorRestraint._AnchorPos;
  anchor.set(anchor.getX(),anchor.getY(),anchor.getZ()*-1.0);
#else
  this->invertStereochemistryOfRestraints();
#endif
	
  ConfigurationEnum config = this->_Configuration;
  if ( config == S_Configuration )
  {
    this->_Configuration = R_Configuration;
  } else if ( config == R_Configuration )
  {
    this->_Configuration = S_Configuration;
  }
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

Matter_sp	Atom_O::getMatterContainedBy()
{_OF();
  Matter_sp	containedBy;
  Residue_sp	wpresidue;
  Residue_sp	residue;
  ASSERTNOTNULL(this->containedBy());
  containedBy = this->containedBy();
  if ( containedBy.nilp() )
  {
    SIMPLE_ERROR(BF("This atom isnt contained by anything"));
  }
  return containedBy;

//    wpresidue = downcast <Residue_O> (containedBy);
//    residue = wpresidue;
//    return residue;
}


CL_LISPIFY_NAME("getResidueContainedBy");
CL_DEFMETHOD     Residue_sp	Atom_O::getResidueContainedBy() 
{
  return this->containedBy().as<Residue_O>();
}

Residue_sp	Atom_O::getResidueContainedBy_const() const
{
  return this->containedBy().as<Residue_O>();
}







void Atom_O::_describeAtomRecursively(string prefix, Atom_sp parent, BondOrder order, int maxDepth, stringstream& ss) const
{_OF();
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
  for ( VectorBond::const_iterator bi=this->bonds.begin(); bi!=this->bonds.end(); bi++ )
  {
    Atom_sp other = (*bi)->getOtherAtom(this->const_sharedThis<Atom_O>());
    if ( other!=parent)
    {
      other->_describeAtomRecursively(prefix+"  ",this->const_sharedThis<Atom_O>(),(*bi)->getOrder(),maxDepth-1,ss);
    }
  }
}


string Atom_O::localEnvironment(int maxDepth) const
{_OF();
  stringstream ss;
  this->_describeAtomRecursively("",_Nil<Atom_O>(),noBond,maxDepth,ss);
  return ss.str();
}









bool	Atom_O::invalid()
{
  IMPLEMENT_ME();
}



AtomIdToAtomMap_sp Atom_O::buildAtomIdMap() const
{_OF();
  SIMPLE_ERROR(BF("Atom should never buildAtomIdMap"));
};


Atom_sp Atom_O::atomWithAtomId(const AtomId& atomId) const
{_OF();
  SIMPLE_ERROR(BF("Atom should never return atomWithAtomId"));
};




void	Atom_O::bumpPosition(double dist)
{_OF();
  double		x, y, z;
  x = core::randomNumber01()*dist*2.0-dist+this->position.getX();
  y = core::randomNumber01()*dist*2.0-dist+this->position.getY();
  z = core::randomNumber01()*dist*2.0-dist+this->position.getZ();
  this->position.set(x,y,z);
}


void	Atom_O::failIfInvalid()
{_OF();
  if ( this->invalid() ) {
    SIMPLE_ERROR(BF("INVALID %s")% this->description());
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


ConstitutionAtom_sp Atom_O::asConstitutionAtom(ConstitutionAtomIndex0N index)
{_OF();
  ConstitutionAtom_sp ca = makeConstitutionAtom(this->getName(),this->_Element,this->type, index,this->_StereochemistryType, this->_Properties);
  return ca;
}

void	Atom_O::defineConstitutionAtomBonding(ConstitutionAtom_sp consAtom, MapAtomsToConstitutionAtomIndex0N atomMap)
{_OF();
  Atom_O::bondIterator bi;
  Atom_sp me = this->sharedThis<Atom_O>();
  for ( bi=this->bonds.begin(); bi!=this->bonds.end(); bi++ )
  {
    Atom_sp otherAtom = (*bi)->getOtherAtom(me);
    if (atomMap.find(otherAtom) != atomMap.end()) {
    // Only use atoms in the same residue
      ConstitutionBond_sp cb = (*bi)->asConstitutionBond(me,atomMap);
      consAtom->addConstitutionBond(cb);
    }
  }
}


}; // namespace chem
