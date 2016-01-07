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

#define	ATOM_FLAGS		uint
#define	ATOMFLAG_OR		1
#define	ATOMFLAG_AND		2
#define	ATOMFLAG_ON		3
#define	ATOMFLAG_OFF		4

    extern uint nextUniqueAtomOrder();

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
	void archive(core::ArchiveP node );

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
		    tempInRing = TEMP_IN_RING
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

    typedef	gctools::Vec0<Bond_sp>	VectorBond;

    SMART(Atom);
    class Atom_O : public Matter_O
    {
	LISP_CLASS(chem,ChemPkg,Atom_O,"Atom",Matter_O);	
    public:
	static void lisp_initGlobals(core::Lisp_sp lisp);
	void initialize();
    public:
//	void	serialize(serialize::SNode node);
        bool fieldsp() const { return true; };
        void	fields(core::Record_sp node);

	friend	class	Loop;
	friend	class	Bond_O;
	friend	class	Residue_O;
    public:
	typedef VectorBond::iterator bondIterator;
    protected:
	core::Symbol_sp		_Alias;	//!< alias name
        Element		_Element;
        Hybridization 	_Hybridization;
        AtomType	type;
	bool		_HintLP;
	bool		_selected;
	ATOM_FLAGS	flags;
	StereochemistryType _StereochemistryType;
	ConfigurationEnum _Configuration;
	uint		_RelativePriority;
	Vector3		position;
	Vector3		_Force;
#if ATOMIC_ANCHOR
	AnchorRestraint	_AnchorRestraint;
#endif
	bool		touched;
	VectorBond	bonds;
	int             _Ionization;
	double		charge;
	double		occupancy;
	double		tempFactor;
	double		vdwRadius;
	double		covalentRadius;
        PertFields      _Pert;
	// Selection mask
	unsigned int	_Mask;
	// copy atoms
	Atom_sp		copyAtom;
	// Ring membership
	int             _RingMembershipCount;
	// Gaff aromaticity flags
	int             _MembershipAr1;
	int             _MembershipAr2;
	int             _MembershipAr3;
	int             _MembershipAr4;
	int             _MembershipAr5;
	// Spanning tree stuff
	int		seenId;
	Atom_sp		_BackSpan;
	Atom_sp		_NextSpan;
	int		backCount;
	int		tempInt;

	// Specific to MOE
	int		moeIndex;
	MoeType 	moeType;	// type string as read from MOE
    public:
#ifdef XML_ARCHIVE
	void archiveBase(core::SNode_sp node);
#endif
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


	void	_addHydrogenWithName(core::Symbol_sp name);
	void _addExistingBond(Bond_sp const& bond);
    public:

CL_NAME("isAromatic");
CL_DEFMETHOD 	bool isAromatic() const { return this->testAnyFlags(Aromatic);};
CL_NAME("setIsAromatic");
CL_DEFMETHOD 	void setIsAromatic(bool b) { if (b) turnOnFlags(Aromatic); else turnOffFlags(Aromatic);};

	virtual	char	getClass()	{return atomId;};

CL_NAME("atomName");
CL_DEFMETHOD         core::Symbol_sp atomName() const { return this->getName();};
	virtual bool	equal(core::T_sp obj) const;
	virtual void	transferCoordinates(Matter_sp other);

	int	totalNetResidueCharge() { THROW_HARD_ERROR(BF("Atoms can't calculate totalNetResidueCharge"));};
	VectorBond::iterator bonds_begin() { return this->bonds.begin(); };
	VectorBond::iterator bonds_end() {return this->bonds.end(); };

	char	getMatterType()	{ return ATOM_CLASS; };

	bool	atLowerAddressThan(Atom_sp a2);
//	uint	getOrder() { return (uint)(this) };
CL_NAME("getAtomId");
CL_DEFMETHOD 	int	getAtomId() { return this->getId(); };
	int	getTempInt();
	void	setTempInt(int o);

CL_NAME("getAlias");
CL_DEFMETHOD 	core::Symbol_sp	getAlias() { return this->_Alias; };
CL_NAME("setAlias");
CL_DEFMETHOD 	void setAlias(core::Symbol_sp alias) { this->_Alias = alias;};


CL_NAME("getFlags");
CL_DEFMETHOD 	ATOM_FLAGS	getFlags() { return this->flags; };
CL_NAME("resetFlags");
CL_DEFMETHOD 	void	resetFlags() { this->flags = 0; };
CL_NAME("turnOnFlags");
CL_DEFMETHOD 	void	turnOnFlags(ATOM_FLAGS o) { this->flags |= o; };
CL_NAME("turnOffFlags");
CL_DEFMETHOD 	void	turnOffFlags(ATOM_FLAGS o) { this->flags &= ~o; };
	void	modifyFlags(int op, ATOM_FLAGS o);
CL_NAME("testAllFlags");
CL_DEFMETHOD 	bool    testAllFlags( ATOM_FLAGS o ) const { return (this->flags&o)==o;};
CL_NAME("testAnyFlags");
CL_DEFMETHOD 	bool    testAnyFlags( ATOM_FLAGS o ) const { return (this->flags&o)!=0;};
	void	clearAllRingMembershipFlags();
CL_NAME("isInRing");
CL_DEFMETHOD 	bool    isInRing() const { return (this->flags&inRing)!= 0; };
	bool    inRingSize(int s) const;
	void	setInRingOfSize(int s);
	void	incrementRingMembershipCount() { this->_RingMembershipCount++;};
CL_NAME("getRingMembershipCount");
CL_DEFMETHOD 	int     getRingMembershipCount() { return this->_RingMembershipCount;};
CL_NAME("setRingMembershipCount");
CL_DEFMETHOD 	void    setRingMembershipCount(int r) {this->_RingMembershipCount = r;};
	string	getNameIndex();
	void	setElementFromAtomName();
	void	setElementFromString(const string& elementSymbol);


	virtual bool isAtom() { return true;};
        gc::Nilable<Atom_sp> highestPriorityNeighborThatIsnt(gc::Nilable<Atom_sp> a);
        gc::Nilable<Atom_sp> lowestPriorityNeighborThatIsnt(gc::Nilable<Atom_sp> a);

	bool	hasAtomNameAttribute(char c);

CL_NAME("getType");
CL_DEFMETHOD         AtomType getType() { return this->type; };
CL_NAME("setType");
CL_DEFMETHOD 	void	setType(AtomType o) { this->type = o; };
CL_NAME("getHybridization");
CL_DEFMETHOD         Hybridization getHybridization() { return this->_Hybridization; };
        string getHybridizationAsString();
CL_NAME("setHybridization");
CL_DEFMETHOD 	void	setHybridization(Hybridization o) { this->_Hybridization = (o); };
	void	setHybridizationFromString(const string& h);
	bool	getHintLP() { return this->_HintLP; };
	void	setHintLP(bool o) { this->_HintLP = o; };
CL_NAME("getElement");
CL_DEFMETHOD 	Element getElement() const { return this->_Element; };
	string	getElementAsString();
        core::Symbol_sp	getElementAsSymbol() const;
CL_NAME("setElement");
CL_DEFMETHOD 	void	setElement(Element o) { this->_Element= o; };
	int     getAtomicNumber();
	double	getAtomicWeight();
	int     getValence();
	virtual bool isVirtual() const { return false;};
CL_NAME("getSelected");
CL_DEFMETHOD 	bool	getSelected() { return this->_selected; };
CL_NAME("setSelected");
CL_DEFMETHOD 	void	setSelected(bool o) { this->_selected = o; };
CL_NAME("getPosition");
CL_DEFMETHOD 	Vector3 getPosition() { return this->position; };
CL_NAME("getPositionInNanometers");
CL_DEFMETHOD 	Vector3 getPositionInNanometers() { return this->position.inNanometers(); };
	Vector3& getPositionRef() { return this->position; };
CL_NAME("setPosition");
CL_DEFMETHOD 	void	setPosition(const Vector3& o) { this->position= o; };
	void	setPositionInNanometers(Vector3 o);
	Vector3 getForce() { return this->_Force; };
	void	setForce(Vector3 o) { this->_Force= o; };
#if ATOMIC_ANCHOR
	void	setAnchorRestraintOn() { this->_AnchorRestraint.setIsDefined(true); };
	void	setAnchorRestraintOff() { this->_AnchorRestraint.setIsDefined(false); };
	bool	isAnchorRestraintOn() { return this->_AnchorRestraint.isDefined(); };
	Vector3 getAnchorPos() { return this->_AnchorRestraint._AnchorPos; };
	void	setAnchorPos(Vector3 o) { this->_AnchorRestraint._AnchorPos = o; };
#endif
	//! Atoms should throw an exception
	void	makeAllAtomNamesInEachResidueUnique();
CL_NAME("getTouched");
CL_DEFMETHOD 	bool	getTouched() { return this->touched; };
CL_NAME("setTouched");
CL_DEFMETHOD 	void	setTouched(bool o) { this->touched = o; };
CL_NAME("getMembershipAr1");
CL_DEFMETHOD 	int     getMembershipAr1()	{ return this->_MembershipAr1; };
CL_NAME("setMembershipAr1");
CL_DEFMETHOD 	void	setMembershipAr1(int c) { this->_MembershipAr1 = c; };
CL_NAME("getMembershipAr2");
CL_DEFMETHOD 	int     getMembershipAr2()	{ return this->_MembershipAr2; };
CL_NAME("setMembershipAr2");
CL_DEFMETHOD 	void	setMembershipAr2(int c) { this->_MembershipAr2 = c; };
CL_NAME("getMembershipAr3");
CL_DEFMETHOD 	int     getMembershipAr3()	{ return this->_MembershipAr3; };
CL_NAME("setMembershipAr3");
CL_DEFMETHOD 	void	setMembershipAr3(int c) { this->_MembershipAr3 = c; };
CL_NAME("getMembershipAr4");
CL_DEFMETHOD 	int     getMembershipAr4()	{ return this->_MembershipAr4; };
CL_NAME("setMembershipAr4");
CL_DEFMETHOD 	void	setMembershipAr4(int c) { this->_MembershipAr4 = c; };
CL_NAME("getMembershipAr5");
CL_DEFMETHOD 	int     getMembershipAr5()	{ return this->_MembershipAr5; };
CL_NAME("setMembershipAr5");
CL_DEFMETHOD 	void	setMembershipAr5(int c) { this->_MembershipAr5 = c; };
CL_NAME("getIonization");
CL_DEFMETHOD 	int     getIonization()	{ return this->_Ionization; };
CL_NAME("setIonization");
CL_DEFMETHOD 	void	setIonization(int c)	{ this->_Ionization = c; };
CL_NAME("getCharge");
CL_DEFMETHOD 	double	getCharge()	{ return this->charge; };
CL_NAME("setCharge");
CL_DEFMETHOD 	void	setCharge(double c)	{ this->charge = c; };
CL_NAME("getOccupancy");
CL_DEFMETHOD 	double	getOccupancy()	{ return this->occupancy; };
CL_NAME("setOccupancy");
CL_DEFMETHOD 	void	setOccupancy(double c)	{ this->occupancy = c; };
CL_NAME("getTempFactor");
CL_DEFMETHOD 	double	getTempFactor()	{ return this->tempFactor; };
CL_NAME("setTempFactor");
CL_DEFMETHOD 	void	setTempFactor(double c)	{ this->tempFactor = c; };
CL_NAME("getVdwRadius");
CL_DEFMETHOD 	double	getVdwRadius()	{ return this->vdwRadius; };
CL_NAME("setVdwRadius");
CL_DEFMETHOD 	void	setVdwRadius(double c)	{ this->vdwRadius = c; };
CL_NAME("getCovalentRadius");
CL_DEFMETHOD 	double	getCovalentRadius()	{ return this->covalentRadius; };
CL_NAME("setCovalentRadius");
CL_DEFMETHOD 	void	setCovalentRadius(double c)	{ this->covalentRadius = c; };
	Matter_wp	getMatterContainedBy();
	Residue_sp	getResidueContainedBy();
	Residue_sp	getResidueContainedBy_const() const;

CL_NAME("getMoeIndex");
CL_DEFMETHOD 	int	getMoeIndex() { return this->moeIndex; };
CL_NAME("setMoeIndex");
CL_DEFMETHOD 	void	setMoeIndex(int o) { this->moeIndex = o; };
CL_NAME("setMoeType");
CL_DEFMETHOD 	void	setMoeType(MoeType t) {this->moeType=t;};
CL_NAME("getMoeType");
CL_DEFMETHOD 	MoeType	getMoeType() {return this->moeType;};

	/*! Return a ConstitutionAtom for this atom and give it the ConstitutionAtomIndex0N (index)
	  @param index The ConstitutionAtomIndex0N that will be assigned to the new ConstitutionAtom */
	ConstitutionAtom_sp asConstitutionAtom(ConstitutionAtomIndex0N index);
	/*! For every bond that this atom makes, define a ConstitutionBond and add it to the (consAtom).
	  Use the (atomMap) to map bonded Atoms to ConstiutionAtomIndex0N indexes.
	  @param consAtom  The ConstitutionAtom to add the ConstitutionBonds to
	  @param atomMap   A map of Atoms to ConstitutionAtomIndex0N values */
	void defineConstitutionAtomBonding(ConstitutionAtom_sp consAtom, MapAtomsToConstitutionAtomIndex0N atomMap);


	// Spanning tree stuff

	/*! When a spanning tree has crawled over this atom assign the atom a unique value (seenId) to indicate this */
CL_NAME("setSeenId");
CL_DEFMETHOD 	void	setSeenId(int i) {this->seenId=i;};

CL_NAME("getSeenId");
CL_DEFMETHOD 	int	getSeenId() {return this->seenId;};
CL_NAME("setBackCount");
CL_DEFMETHOD 	void	setBackCount(int i) {this->backCount=i;};
CL_NAME("getBackCount");
CL_DEFMETHOD 	int	getBackCount() {return this->backCount;};
	void	invalidateBackSpan() {this->_BackSpan = _Nil<Atom_O>(); };// = false;
	bool	isBackSpanValid() {_OF();ANN(this->_BackSpan);return this->_BackSpan.notnilp(); };
CL_NAME("setBackSpan");
CL_DEFMETHOD 	void	setBackSpan(Atom_sp a) {this->_BackSpan=a;}//this->validBackSpan=true;
	Atom_sp	getBackSpan();
	void	invalidateNextSpan() {this->_NextSpan = _Nil<Atom_O>(); };
	bool	isNextSpanValid() { return this->_NextSpan.notnilp(); };
CL_NAME("setNextSpan");
CL_DEFMETHOD 	void	setNextSpan(Atom_sp a) {this->_NextSpan=a;};//;this->validNextSpan=true;
	Atom_sp	getNextSpan();

	/*! Return a local spanning tree as a list of lists of atoms
	 * It will have the form ( rootAtom ( neighbor1 (neighbor1-1) ... ) (neighbor2 (neighbor2-1)...) )
	 * The depth describes how deep to go 
	 *   - 0 means return yourself
	 *   - 1 means return yourself and immediatly attached neighbors
	 */
    private:
	core::List_sp 	_expandLocalSpanningTree(Atom_sp avoid, Bond_sp bond, uint depth);
    public:
	core::List_sp	localSpanningTree(uint depth);


	Atom_sp	aliasAtomOrNil(Alias_sp alias) {IMPLEMENT_ME();};
	Residue_sp aliasResidueOrNil(Alias_sp alias) {IMPLEMENT_ME();};

CL_NAME("flagsSet");
CL_DEFMETHOD 	bool	flagsSet(int f)	{return (bool)((this->flags&f)!=0);};
	void	applyTransformToAtoms(const Matrix& m);
	VectorBond& getBonds()		{return this->bonds; };
        VectorAtom getBondedAtoms();
	core::List_sp	bondedAtomsAsList();

	Bond_sp	basicBondTo(  Atom_sp a, BondOrder o );
	Bond_sp	bondTo(  Atom_sp a, BondOrder o );
	Bond_sp	bondToSingle(  Atom_sp a );
	Bond_sp	bondToBondOrderInt(  Atom_sp a, int o )
	{
	    return this->bondTo(a,(BondOrder)(o));
	};
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
	void fillInImplicitHydrogens();

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
CL_NAME("getRelativePriority");
CL_DEFMETHOD 	uint	getRelativePriority() { return this->_RelativePriority;};
	void	setRelativePriority(uint p) { this->_RelativePriority = p;};
	bool	isConfigurable();
//	void	setConfigurationInfo(ConfigurationInfo& conf);
	void	setConfiguration(ConfigurationEnum conf);
CL_NAME("setStereochemistryType");
CL_DEFMETHOD 	void	setStereochemistryType(StereochemistryType conf) { this->_StereochemistryType = conf;};

        virtual bool applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value);

CL_NAME("getStereochemistryType");
CL_DEFMETHOD 	StereochemistryType getStereochemistryType() { return this->_StereochemistryType; };
CL_NAME("getConfiguration");
CL_DEFMETHOD 	ConfigurationEnum getConfiguration() { return this->_Configuration;};
	string	getConfigurationAsString();
#if 0
	void	setConfigurationPriorityHighToLow(string highest, string high, string low, string lowest);
	string	getConfigurationPriorityHighest();
	string	getConfigurationPriorityHigh();
	string	getConfigurationPriorityLow();
	string	getConfigurationPriorityLowest();
#endif
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
	Conformation_wp getWeakConformation();

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

	virtual Atom_sp atomWithAtomId(AtomId_sp atomId) const;



    public:
	virtual Matter_sp	copy();
    protected:
	virtual Matter_sp copyDontRedirectAtoms();
	virtual void redirectAtoms();

	/*! Used by Bond to connect a copied atom to a bond */
    public:
	void addBond(Bond_sp bond);

        Atom_O() {};
        virtual ~Atom_O() {};
    };

#if 0
    inline	Atom_sp	Atom()	{
	Atom_O*	da;
	Atom_sp	ra;
	ra = Atom_sp(da);
	return ra;
    }
#endif

/*! new_Atom_sp is the constructor for atoms.  It allocates the atom and
 * saves a weak_ptr to that atom so that we don't run into problems with shared_from_this
 * with python
 */
//extern	Atom_sp new_Atom_sp();
//extern	Atom_sp _copy_Atom_sp(Atom_O*);	// Used when you want to create a copy




    extern core::Symbol_sp& _sym__PLUS_atomFlagSymbolConverter_PLUS_;

};

TRANSLATE(chem::Atom_O);


DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::AtomFlagEnum,chem::_sym__PLUS_atomFlagSymbolConverter_PLUS_);

#endif
