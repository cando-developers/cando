/*
    File: atom.h
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
// (C) 2004 Christian E. Schafmeister
//

/*
 *	atom.h
 *
 *	Maintain an Atom
 *
 */

#ifndef	ATOM_H
#define	ATOM_H
#include <clasp/core/common.h>

#include <string>
#include <vector>
#include <cando/chem/matter.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atom.fwd.h>
#include <cando/chem/elements.h>
#include <cando/geom/matrix.h>

#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/chem/bond.fwd.h>
#include <cando/chem/atom.fwd.h>
#include <cando/chem/atomIdMap.fwd.h>
#include <clasp/core/serialize.h>
#include <cando/chem/atomId.fwd.h>
#include <cando/chem/chemPackage.h>

#if !BUILT_IN_RESTRAINT_LIST
#include <cando/chem/restraint.fwd.h>
#endif

namespace chem
{
  SMART(Conformation);
  SMART(Alias);

#define	ATOM_FLAGS		size_t
#define	ATOMFLAG_OR		1
#define	ATOMFLAG_AND		2
#define	ATOMFLAG_ON		3
#define	ATOMFLAG_OFF		4

  extern size_t globalUniqueAtomOrder;

  inline size_t nextUniqueAtomOrder() {
    ++globalUniqueAtomOrder;
    return globalUniqueAtomOrder;
  }

#include <cando/chem/elements.h>


#if ATOMIC_ANCHOR
  struct	AnchorRestraint {
    bool	_Defined;
    Vector3	_AnchorPos;

	//! If the restraint is defined return true
    bool isDefined() { return this->_Defined; };
	//! Set if the restraint is defined
    void setIsDefined(bool d) { this->_Defined = d; };
	//! Serialize the object
  };
#endif
  inline	string	XmlTag_Atom() { return "Atom";};

    // Atom flags
#define	ALL_OFF			0x00000000
#define	ATOM_NEEDS_MINIMIZER	0x00000001
#define	SELECTED		0x00000002
#define	ON_GRID			0x00000004
#define	ATOM_FIXED		0x00000008
#define In3MemberRing           0x00000010
#define In4MemberRing           0x00000020
#define In5MemberRing           0x00000040
#define In6MemberRing           0x00000080
#define In7MemberRing           0x00000100
#define In8MemberRing           0x00000200
#define	Aromatic		0x00000400
#define InRing                  In3MemberRing|In4MemberRing|In5MemberRing|In6MemberRing|In7MemberRing|In8MemberRing
#define	TEMP_IN_RING		0x00000400
#define NEEDS_BUILD             0x00000800
#define MEMBERSHIP_AR1          0x00001000
#define MEMBERSHIP_AR2          0x00002000
#define MEMBERSHIP_AR3          0x00004000
#define MEMBERSHIP_AR4          0x00008000
#define MEMBERSHIP_AR5          0x00010000
#define HINT_LP                 0x00020000
#define LAST_FLAG               0x80000000

  typedef enum {  flagOr = ATOMFLAG_OR,
                  flagAnd = ATOMFLAG_AND,
                  flagOn = ATOMFLAG_ON,
                  flagOff = ATOMFLAG_OFF
  } AtomFlagOperation;

  typedef enum {  allOff = ALL_OFF,
                  atomNeedsMinimizer = ATOM_NEEDS_MINIMIZER,
                  selected = SELECTED,
                  onGrid = ON_GRID,
                  atomFixed = ATOM_FIXED,
                  in3MemberRing = In3MemberRing,
                  in4MemberRing = In4MemberRing,
                  in5MemberRing = In5MemberRing,
                  in6MemberRing = In6MemberRing,
                  in7MemberRing = In7MemberRing,
                  in8MemberRing = In8MemberRing,
                  inRing = InRing,
                  tempInRing = TEMP_IN_RING,
                  needsBuild = NEEDS_BUILD,
                  MembershipAr1 = MEMBERSHIP_AR1,
                  MembershipAr2 = MEMBERSHIP_AR2,
                  MembershipAr3 = MEMBERSHIP_AR3,
                  MembershipAr4 = MEMBERSHIP_AR4,
                  MembershipAr5 = MEMBERSHIP_AR5,
                  HintLp        = HINT_LP,
                  LastFlag      = LAST_FLAG
  } AtomFlagEnum;


  struct PertFields {
    Element _Element;
    core::Symbol_sp _Name;
    AtomType _Type;
    double _Charge;
    double _Polar;
  };

  extern	long	__AtomBuildCounter;

  extern	void	advanceAtomBuildCounter();
  extern	long	getAtomBuildCounter();


  SMART(Bond);
  SMART(BondList);

  typedef	gctools::Vec0_uncopyable<Bond_sp>	VectorBond;

};

template <>
struct gctools::GCInfo<chem::Atom_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {
  // Every atom gets a unique atom 64-bit index
  // 64 bits means we should never run out of unique
  SMART(Atom);
  class Atom_O : public Matter_O {
    LISP_CLASS(chem,ChemPkg,Atom_O,"Atom",Matter_O);	
    friend	class	Loop;
    friend	class	Bond_O;
    friend	class	Residue_O;
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    typedef VectorBond::iterator bondIterator;
  public:
    size_t              _UniqueAtomOrder;
    core::Symbol_sp	_Alias;	//!< alias name
    ATOM_FLAGS	        flags;
    Vector3		position;
//    Vector3		_Force;
#if ATOMIC_ANCHOR
    AnchorRestraint	_AnchorRestraint;
#endif
    AtomType	        type;
    VectorBond	        bonds;
    double		charge;
	// copy atoms
    Atom_sp		copyAtom;
	// Specific to MOE
//    MoeType 	moeType;	// type string as read from MOE
//    int		moeIndex;
    int		backCount;
    int		tempInt;
    double      dAtomTemp; // a temporary double value - terrible idea
    int		seenId;
    float		occupancy;
    float		tempFactor;
    float		vdwRadius;
    float		covalentRadius;
    StereochemistryType _StereochemistryType;
    ConfigurationEnum   _Configuration;
    int		        _RelativePriority;
	// Selection mask
    unsigned int	_Mask;
    Element		_Element;
    Hybridization 	_Hybridization;
	// Ring membership
    ushort       _RingMembershipCount;
    short        _Ionization;
  public:
  public:
    static Atom_sp make(MatterName name, Element element);
  public:
    static int priorityOrder(Atom_sp a, Atom_sp b);
//	static double vdwRadiusForElement(Element element);
  public:
    uint	getMask() { return this->_Mask; };
    void	setMask(uint m) { this->_Mask = m; };

//	void	duplicateFrom(Atom_O* a);
    Atom_sp	getCopyAtom()			{return this->copyAtom;};


    void	_addHydrogenWithName(Residue_sp residueContainedBy, core::Symbol_sp name);
    void _addExistingBond(Bond_sp const& bond);
  public:

    inline size_t atLowerUniqueAtomOrderThan(Atom_sp a) const {
      return this->_UniqueAtomOrder < a->_UniqueAtomOrder;
    };
    CL_LISPIFY_NAME("isAromatic");
    CL_DEFMETHOD 	bool isAromatic() const { return this->testAnyFlags(Aromatic);};
    CL_LISPIFY_NAME("setIsAromatic");
    CL_DEFMETHOD 	void setIsAromatic(bool b) { if (b) turnOnFlags(Aromatic); else turnOffFlags(Aromatic);};

    virtual	char	getClass()	{return atomId;};

    CL_LISPIFY_NAME("atomName");
    CL_DEFMETHOD         core::Symbol_sp atomName() const { return this->getName();};
    bool	Atom_equal(core::T_sp obj) const;
    virtual void	transferCoordinates(Matter_sp other);

    int	totalNetResidueCharge() { THROW_HARD_ERROR(BF("Atoms can't calculate totalNetResidueCharge"));};
    VectorBond::iterator bonds_begin() { return this->bonds.begin(); };
    VectorBond::iterator bonds_end() {return this->bonds.end(); };

    char	getMatterType()	{ return ATOM_CLASS; };

//    CL_LISPIFY_NAME("getAtomId");
//    CL_DEFMETHOD 	int	getAtomId() { return this->getId(); };
    int	getTempInt() { return this->tempInt; };
    void	setTempInt(int o) { this->tempInt = o;};

    CL_LISPIFY_NAME("getAlias");
    CL_DEFMETHOD 	core::Symbol_sp	getAlias() { return this->_Alias; };
    CL_LISPIFY_NAME("setAlias");
    CL_DEFMETHOD 	void setAlias(core::Symbol_sp alias) { this->_Alias = alias;};


    CL_LISPIFY_NAME("getFlags");
    CL_DEFMETHOD 	ATOM_FLAGS	getFlags() { return this->flags; };
    CL_LISPIFY_NAME("resetFlags");
    CL_DEFMETHOD 	void	resetFlags() { this->flags = 0; };
    CL_LISPIFY_NAME("turnOnFlags");
    CL_DEFMETHOD 	void	turnOnFlags(ATOM_FLAGS o) { this->flags |= o; };
    CL_LISPIFY_NAME("turnOffFlags");
    CL_DEFMETHOD 	void	turnOffFlags(ATOM_FLAGS o) { this->flags &= ~o; };
    void	modifyFlags(int op, ATOM_FLAGS o);
    CL_LISPIFY_NAME("testAllFlags");
    CL_DEFMETHOD 	bool    testAllFlags( ATOM_FLAGS o ) const { return (this->flags&o)==o;};
    CL_LISPIFY_NAME("testAnyFlags");
    CL_DEFMETHOD 	bool    testAnyFlags( ATOM_FLAGS o ) const { return (this->flags&o)!=0;};
    void	clearAllRingMembershipFlags();
    CL_LISPIFY_NAME("isInRing");
    CL_DEFMETHOD 	bool    isInRing() const { return (this->flags&inRing)!= 0; };
    CL_DEFMETHOD bool needs_build() const { return (this->flags&needsBuild)!=0;};
    CL_DEFMETHOD void setf_needs_build(bool val) {
      if (val) this->turnOnFlags(needsBuild);
      else this->turnOffFlags(needsBuild);
    }
    bool    inRingSize(int s) const;
    void	setInRingOfSize(int s);
    void	incrementRingMembershipCount() { this->_RingMembershipCount++;};
    CL_LISPIFY_NAME("getRingMembershipCount");
    CL_DEFMETHOD 	int     getRingMembershipCount() { return this->_RingMembershipCount;};
    CL_LISPIFY_NAME("setRingMembershipCount");
    CL_DEFMETHOD 	void    setRingMembershipCount(int r) {this->_RingMembershipCount = r;};
    string	getNameIndex();
    void	setElementFromAtomName();
    void	setElementFromString(const string& elementSymbol);


    virtual bool isAtom() { return true;};
    gc::Nilable<Atom_sp> highestPriorityNeighborThatIsnt(gc::Nilable<Atom_sp> a);
    gc::Nilable<Atom_sp> lowestPriorityNeighborThatIsnt(gc::Nilable<Atom_sp> a);

    bool	hasAtomNameAttribute(char c);

    CL_LISPIFY_NAME("getType");
    CL_DEFMETHOD         AtomType getType() { return this->type; };
    CL_LISPIFY_NAME("setType");
    CL_DEFMETHOD 	void	setType(AtomType o) { this->type = o; };
    CL_LISPIFY_NAME("getHybridization");
    CL_DEFMETHOD         Hybridization getHybridization() { return this->_Hybridization; };
    string getHybridizationAsString();
    CL_LISPIFY_NAME("setHybridization");
    CL_DEFMETHOD 	void	setHybridization(Hybridization o) { this->_Hybridization = (o); };
    void	setHybridizationFromString(const string& h);
    bool	getHintLP() { return this->testAnyFlags(HintLp); };
    void	setHintLP(bool o) { if (o) turnOnFlags(HintLp); else turnOffFlags(HintLp); };
    CL_LISPIFY_NAME("GET-ELEMENT");
    CL_DEFMETHOD 	Element getElement() const { return this->_Element; };
    string	getElementAsString() const;
    core::Symbol_sp	getElementAsSymbol() const;
    CL_LISPIFY_NAME("SET-ELEMENT");
    CL_DEFMETHOD 	void	setElement(Element o) { this->_Element= o; };
    int     getAtomicNumber();
    int     getIntegerAtomicMass();
    double	getAtomicWeight();
    int     getValence();
    virtual bool isVirtual() const { return false;};
#if 0
    CL_LISPIFY_NAME("getSelected");
    CL_DEFMETHOD 	bool	getSelected() { return this->_selected; };
    CL_LISPIFY_NAME("setSelected");
    CL_DEFMETHOD 	void	setSelected(bool o) { this->_selected = o; };
#endif
    CL_LISPIFY_NAME("getPosition");
    CL_DEFMETHOD 	Vector3 getPosition() { return this->position; };
    CL_DEFMETHOD bool atomWithinAngstroms(Atom_sp other, float angstroms) const;
    CL_LISPIFY_NAME("getPositionInNanometers");
    CL_DEFMETHOD 	Vector3 getPositionInNanometers() { return this->position.inNanometers(); };
    Vector3& getPositionRef() { return this->position; };
    CL_LISPIFY_NAME("setPosition");
    CL_DEFMETHOD 	void	setPosition(const Vector3& o) { this->position= o; };
    CL_DEFMETHOD 	void	setPositionXYZ(double x, double y, double z) {
      this->position.getX() = x;
      this->position.getY() = y;
      this->position.getZ() = z;
    }
    void	setPositionInNanometers(Vector3 o);
#if 0
    Vector3 getForce() { return this->_Force; };
    void	setForce(Vector3 o) { this->_Force= o; };
#endif
#if ATOMIC_ANCHOR
    void	setAnchorRestraintOn() { this->_AnchorRestraint.setIsDefined(true); };
    void	setAnchorRestraintOff() { this->_AnchorRestraint.setIsDefined(false); };
    bool	isAnchorRestraintOn() { return this->_AnchorRestraint.isDefined(); };
    Vector3 getAnchorPos() { return this->_AnchorRestraint._AnchorPos; };
    void	setAnchorPos(Vector3 o) { this->_AnchorRestraint._AnchorPos = o; };
#endif
    float distanceSquaredToAtom(Atom_sp other);
	//! Atoms should throw an exception
    void	makeAllAtomNamesInEachResidueUnique();
#if 0
    CL_LISPIFY_NAME("getTouched");
    CL_DEFMETHOD 	bool	getTouched() { return this->touched; };
    CL_LISPIFY_NAME("setTouched");
    CL_DEFMETHOD 	void	setTouched(bool o) { this->touched = o; };
#endif
    CL_DEFMETHOD 	void setMembershipAr1(bool b) { if (b) turnOnFlags(MembershipAr1); else turnOffFlags(MembershipAr1);};
    CL_DEFMETHOD 	void setMembershipAr2(bool b) { if (b) turnOnFlags(MembershipAr2); else turnOffFlags(MembershipAr2);};
    CL_DEFMETHOD 	void setMembershipAr3(bool b) { if (b) turnOnFlags(MembershipAr3); else turnOffFlags(MembershipAr3);};
    CL_DEFMETHOD 	void setMembershipAr4(bool b) { if (b) turnOnFlags(MembershipAr4); else turnOffFlags(MembershipAr4);};
    CL_DEFMETHOD 	void setMembershipAr5(bool b) { if (b) turnOnFlags(MembershipAr5); else turnOffFlags(MembershipAr5);};
    CL_DEFMETHOD 	bool     getMembershipAr1()	{ return testAllFlags(MembershipAr1); };
    CL_DEFMETHOD 	bool     getMembershipAr2()	{ return testAllFlags(MembershipAr2); };
    CL_DEFMETHOD 	bool     getMembershipAr3()	{ return testAllFlags(MembershipAr3); };
    CL_DEFMETHOD 	bool     getMembershipAr4()	{ return testAllFlags(MembershipAr4); };
    CL_DEFMETHOD 	bool     getMembershipAr5()	{ return testAllFlags(MembershipAr5); };
    CL_LISPIFY_NAME("getIonization");
    CL_DEFMETHOD 	int     getIonization()	{ return this->_Ionization; };
    CL_LISPIFY_NAME("setIonization");
    CL_DEFMETHOD 	void	setIonization(int c)	{ this->_Ionization = c; };
    CL_LISPIFY_NAME("getCharge");
    CL_DEFMETHOD 	double	getCharge()	{ return this->charge; };
    CL_LISPIFY_NAME("setCharge");
    CL_DEFMETHOD 	void	setCharge(double c)	{ this->charge = c; };
    CL_LISPIFY_NAME("getOccupancy");
    CL_DEFMETHOD 	double	getOccupancy()	{ return this->occupancy; };
    CL_LISPIFY_NAME("setOccupancy");
    CL_DEFMETHOD 	void	setOccupancy(double c)	{ this->occupancy = c; };
    CL_LISPIFY_NAME("getTempFactor");
    CL_DEFMETHOD 	double	getTempFactor()	{ return this->tempFactor; };
    CL_LISPIFY_NAME("setTempFactor");
    CL_DEFMETHOD 	void	setTempFactor(double c)	{ this->tempFactor = c; };
    CL_LISPIFY_NAME("getVdwRadius");
    CL_DEFMETHOD 	double	getVdwRadius()	{ return this->vdwRadius; };
    CL_LISPIFY_NAME("setVdwRadius");
    CL_DEFMETHOD 	void	setVdwRadius(double c)	{ this->vdwRadius = c; };
    CL_LISPIFY_NAME("getCovalentRadius");
    CL_DEFMETHOD 	double	getCovalentRadius()	{ return this->covalentRadius; };
    CL_LISPIFY_NAME("setCovalentRadius");
    CL_DEFMETHOD 	void	setCovalentRadius(double c)	{ this->covalentRadius = c; };
    Matter_sp	getMatterContainedBy();
    Residue_sp	getResidueContainedBy();
    Residue_sp	getResidueContainedBy_const() const;

#if 0
    CL_LISPIFY_NAME("getMoeIndex");
    CL_DEFMETHOD 	int	getMoeIndex() { return this->moeIndex; };
    CL_LISPIFY_NAME("setMoeIndex");
    CL_DEFMETHOD 	void	setMoeIndex(int o) { this->moeIndex = o; };
    CL_LISPIFY_NAME("setMoeType");
    CL_DEFMETHOD 	void	setMoeType(MoeType type) {this->moeType=type;};
    CL_LISPIFY_NAME("getMoeType");
    CL_DEFMETHOD 	MoeType	getMoeType() {return this->moeType;};
#endif
	/*! Return a ConstitutionAtom for this atom and give it the ConstitutionAtomIndex0N (index)
	  @param index The ConstitutionAtomIndex0N that will be assigned to the new ConstitutionAtom */
    ConstitutionAtom_sp asConstitutionAtom(ConstitutionAtomIndex0N index);
	/*! For every bond that this atom makes, define a ConstitutionBond and add it to the (consAtom).
	  Use the (atomMap) to map bonded Atoms to ConstiutionAtomIndex0N indexes.
	  @param consAtom  The ConstitutionAtom to add the ConstitutionBonds to
	  @param atomMap   A map of Atoms to ConstitutionAtomIndex0N values */
    void defineConstitutionAtomBonding(ConstitutionAtom_sp consAtom, MapAtomsToConstitutionAtomIndex0N atomMap);

  private:
    core::List_sp 	_expandLocalSpanningTree(Atom_sp avoid, Bond_sp bond, uint depth);
  public:
    core::List_sp	localSpanningTree(uint depth);


    Atom_sp	aliasAtomOrNil(Alias_sp alias) {IMPLEMENT_ME();};
    Residue_sp aliasResidueOrNil(Alias_sp alias) {IMPLEMENT_ME();};

    CL_LISPIFY_NAME("flagsSet");
    CL_DEFMETHOD 	bool	flagsSet(int f)	{return (bool)((this->flags&f)!=0);};
    void	applyTransformToAtoms(const Matrix& m);
    VectorBond& getBonds()		{return this->bonds; };
    VectorAtom getBondedAtoms();
    core::List_sp	bondedAtomsAsList();

    Bond_sp	basicBondTo(  Atom_sp a, BondOrder o );
    Bond_sp	bondTo(  Atom_sp a, BondOrder o );
    Bond_sp	bondToSingle(  Atom_sp a );
    Bond_sp	bondToOrderInt(  Atom_sp a, int o );
    Bond_sp	getBondTo(Atom_sp a);
    bool 	hasBondWithOrder(BondOrder o) const;
    bool	isBondedToAtomNamed( MatterName name );
    bool	isBondedToElementOrder( Element element, BondOrder o );
    bool	isBondedToElementHybridization( Element element, Hybridization hybrid );
    bool	isBondedToElementHybridizationElementHybridization(Element element1 , Hybridization hybrid1, Element element2 , Hybridization hybrid2 );
    bool	isBondedToWithBondOrder( Atom_sp a1, BondOrder o );
    bool	isBondedTo( Atom_sp a1 );
    void	basicRemoveBondTo( Atom_sp a );
    void	removeBondTo(Atom_sp a);
    void	removeAllBonds();
    int	coordination();
    Bond_sp bondAtIndex(int i);
    Atom_sp	bondedNeighbor(int i);
    Atom_sp	bondedNeighborWithName(MatterName name);
    BondOrder bondOrderTo( Atom_sp aTarget );
    int	numberOfBonds() const;
    BondOrder	bondedOrder(int i);
	/*! Total all of the bond orders.  */
    uint totalBondOrder();
    uint maxTotalBondOrder();
    uint numberOfOpenValence();
    core::List_sp createImplicitHydrogenNames();
    void fillInImplicitHydrogensWithResidue(Residue_sp residue);
    
    void fillInImplicitHydrogens(); // signals error - use fillInImplicitHydrogensWithResidue

    void randomizeAtomPosition(double width);
    void perturbAtomPosition(double dist);

	/*! Calculate the stereochemical configuration "R" or "S"
	 * based on the positions of this atom and the four attached atoms
	 */
    ConfigurationEnum calculateStereochemicalConfiguration();
    string calculateStereochemicalConfigurationAsString();

    core::List_sp	bondsAsList();
    BondList_sp	getBondList();
    bool		isHeavyAtom();
    BondList_sp	getHeavyAtomBondList();
    int             getBondedHydrogenCount();

    virtual	bool testConsistancy(Matter_sp c);
    void invertStructureAndRestraints();
	//
	// Configuration
	//
	//! If the atom can have a stereochemical configuration return true
    core::List_sp	getNeighborsByRelativePriority();
    core::List_sp	getNeighborsForAbsoluteConfiguration();
    CL_LISPIFY_NAME("getRelativePriority");
    CL_DEFMETHOD 	int	getRelativePriority() { return this->_RelativePriority;};
    void	setRelativePriority(int p) { this->_RelativePriority = p;};
    bool	isConfigurable();
//	void	setConfigurationInfo(ConfigurationInfo& conf);
    void	setConfiguration(ConfigurationEnum conf);
    CL_LISPIFY_NAME("setStereochemistryType");
    CL_DEFMETHOD 	void	setStereochemistryType(StereochemistryType conf) { this->_StereochemistryType = conf;};

    virtual bool applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value);

    CL_LISPIFY_NAME("getStereochemistryType");
    CL_DEFMETHOD 	StereochemistryType getStereochemistryType() { return this->_StereochemistryType; };
    CL_LISPIFY_NAME("getConfiguration");
    CL_DEFMETHOD 	ConfigurationEnum getConfiguration() { return this->_Configuration;};
    string	getConfigurationAsString();

    void setAbsoluteConfiguration(ConfigurationEnum config, Atom_sp n1, Atom_sp n2, Atom_sp n3);

    void	addUniqueIntraResidueBondCopiesToBondList(BondList_sp list);
    void	addUniqueInterResidueBondCopiesToBondList(BondList_sp list);

	//! Downgrade this
//	void	addUniqueIntraResidueBondsToQDomNodeAsChildren(adapt::QDomNode_sp node);
	//! Downgrade this
    void	addUniqueIntraResidueBondsToVectorBonds(VectorBond& bonds);
	//! Downgrade this
    void	addInterResidueBondsToBondList(BondList_sp bonds);

    void	failIfInvalid();
    bool	invalid();

	//! Add a random perturbation to the atom position
    void	bumpPosition(double dist);


	//
	// storage id handling
	//
//	void	setAtomStorageId(int sid) { this->setStorageId(sid);};
//	int	getAtomStorageId() { return this->getStorageId();};

    virtual string	subMatter() { return "Undefined"; };
    virtual	string	description() const;
    string __repr__() const;

    void setConformation(Conformation_sp b);
    Conformation_sp getWeakConformation();

    void setAtomHolderIndex(uint i);
    uint getAtomHolderIndex();

  private:
    void _describeAtomRecursively(string prefix, Atom_sp parent, BondOrder order, int maxDepth, stringstream& ss) const;
  public:
    string localEnvironment(int maxDepth) const;
//	virtual Atom_sp	copyDropExternalResidueBonds();
	//! Copy constructor
    Atom_O( const Atom_O& ss );

	/*! Build a map of AtomIds to Atoms */
    virtual AtomIdToAtomMap_sp buildAtomIdMap() const;

    virtual Atom_sp atomWithAtomId(const AtomId& atomId) const;



  public:
    virtual Matter_sp	copy(core::T_sp new_to_old);
  protected:
    virtual Matter_sp copyDontRedirectAtoms(core::T_sp new_to_old);
    virtual void redirectAtoms();

	/*! Used by Bond to connect a copied atom to a bond */
  public:
    void addBond(Bond_sp bond);

/*! There is so much crap in this class that should be removed.
    Use alists for some of these things like aromatic ring membership
*/
  Atom_O() :
    _UniqueAtomOrder(nextUniqueAtomOrder()),
    _Alias(_Nil<core::Symbol_O>()),
      _Element(element_Undefined),
      _Hybridization(hybridization_undefined),
      type(_Nil<core::Symbol_O>()),
      flags(0),
      _StereochemistryType(undefinedCenter),
      _Configuration(undefinedConfiguration),
      _RelativePriority(0),
      _Ionization(0),
      charge(0.0),
      occupancy(0.0),
      tempFactor(1.0),
      vdwRadius(0.0),
      covalentRadius(0.0),
      _Mask(0),
      copyAtom(_Unbound<chem::Atom_O>()),
      _RingMembershipCount(0),
      tempInt(0)
//      moeIndex(0),
//      moeType(_Nil<core::Symbol_O>())
      {};
    virtual ~Atom_O() {};
  };

  extern core::Symbol_sp& _sym__PLUS_atomFlagSymbolConverter_PLUS_;

};

DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::AtomFlagEnum,chem::_sym__PLUS_atomFlagSymbolConverter_PLUS_);

#endif
