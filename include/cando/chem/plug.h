       
       
#ifndef	Plug_H //[
#define Plug_H



#include <clasp/core/common.h>


#include <cando/chem/entityNameSet.h>
#include <clasp/core/holder.h>




#include <cando/chem/chemPackage.h>


namespace chem
{


    typedef	uint	PlugId;

class	CandoDatabase_O;



SMART(Mate);
class Mate_O : public EntityNameSetWithCap_O
{
    LISP_BASE1(EntityNameSetWithCap_O);
    LISP_CLASS(chem,ChemPkg,Mate_O,"Mate");
#if INIT_TO_FACTORIES
 public:
    static Mate_sp make(core::Symbol_sp cap);
#else
    DECLARE_INIT();
#endif
public:
    void	initialize();
//    void	archiveBase(core::ArchiveP node);

private:
    /*! The cap is the name of a monomer that will cap
     * the current plug.  It should be a member of
     * one of the Mates but if OverrideCapMateRequirement
     * then it can be anything.  Make sure its chemically
     * reasonable
     */
    core::Symbol_sp	_Cap;

public:
//	static  Mate_sp create(core::Lisp_sp e,CandoDatabase_sp db);
public:

    core::Symbol_sp getName() const;
    void	setCap(core::Symbol_sp cc ) {this->_Cap = cc; };
    core::Symbol_sp	getCap() { return this->_Cap; };

    string descriptionOfContents() const;

    string __repr__() const;

    Mate_O( const Mate_O& ss ); //!< Copy constructor


    DEFAULT_CTOR_DTOR(Mate_O);
};





SMART(RingClosingMate);
class RingClosingMate_O : public EntityNameSetWithCap_O
{
    LISP_BASE1(EntityNameSetWithCap_O);
    LISP_CLASS(chem,ChemPkg,RingClosingMate_O,"RingClosingMate");
public:
//    void	archiveBase(core::ArchiveP node);

public:
//	static  RingClosingMate_sp create(core::Lisp_sp e,CandoDatabase_sp db);
public:

	core::Symbol_sp getName() const;

	string descriptionOfContents() const;

	RingClosingMate_O( const RingClosingMate_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(RingClosingMate_O);
};






SMART(Plug);
class Plug_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,Plug_O,"Plug");
#if INIT_TO_FACTORIES
 public:
    static Plug_sp make(core::Symbol_sp name);
#else
    DECLARE_INIT();
#endif
public:
//    void	archiveBase(core::ArchiveP node);
    void	initialize();
public:
    typedef	gctools::Vec0<Mate_sp>	Mates;
private:
    Constitution_wp		_WeakConstitution;
    /*! Name of plug, prefix of '-' means its an incoming plug
     * prefix of '+' means its an outgoing plug
     */
    core::Symbol_sp			_Name;
public:
    Constitution_sp	getConstitution();
    void	setConstitution(Constitution_sp con);

    string descriptionOfContents() const;

    void	setName(core::Symbol_sp s) { this->_Name = s; };
    core::Symbol_sp getName() const { return this->_Name; };

    virtual bool getIsIn() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
    virtual bool getIsRingClosing() { return false; };

    Plug_sp getPlug() { return this->sharedThis<Plug_O>();};

    virtual core::Symbol_sp otherSidePlugName();
    virtual bool hasMates() { return false;};

    virtual void	addMate(Mate_sp z ) {_OF();SUBCLASS_MUST_IMPLEMENT();};

    virtual Mates::iterator begin_Mates();
    virtual Mates::iterator end_Mates();

    virtual int		numberOfMates() { return 0; };
//    virtual bool	recognizesMateNameOrPdb(const string& name) { return false;};

    virtual core::List_sp	matesAsCons() { return _Nil<core::T_O>(); };

    Plug_O( const Plug_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(Plug_O);
};





SMART(PlugWithMates);
class PlugWithMates_O : public Plug_O
{
    LISP_BASE1(Plug_O);
    LISP_CLASS(chem,ChemPkg,PlugWithMates_O,"PlugWithMates");
#if INIT_TO_FACTORIES
 public:
    static PlugWithMates_sp make(core::Symbol_sp bond0, core::Symbol_sp bond1, core::List_sp mates);
#else
    DECLARE_INIT();
#endif
public:
//    void	archiveBase(core::ArchiveP node);
    void	initialize();

protected:
    //! Name of atom for first bond
    core::Symbol_sp			_B0;
    //! Name of atom for second bond
    core::Symbol_sp			_B1;

    //! RepresentedEntityNameSets that this plug can plug into
    Mates		_Mates;
public:
    string descriptionOfContents() const;

    bool getIsIn() {_OF(); SUBCLASS_MUST_IMPLEMENT();};

    PlugWithMates_sp getPlugWithMates() { return this->sharedThis<PlugWithMates_O>();};

    void	setB0(core::Symbol_sp s) {this->_B0 = s;};
    core::Symbol_sp	getB0() const	{ return this->_B0; };

    void	setB1(core::Symbol_sp s) {this->_B1 = s;};
    core::Symbol_sp getB1() const	{ return this->_B1; };
    bool	hasB1() const { return this->_B1.notnilp();};

    virtual bool hasMates() { return true;};

    string __repr__() const;


//	string otherSidePlugWithMatesName();

    void	addMate(Mate_sp z ) {this->_Mates.push_back(z);};

    Mates::iterator begin_Mates() { return this->_Mates.begin(); };
    Mates::iterator end_Mates() { return this->_Mates.end(); };

    int		numberOfMates() { return this->_Mates.size(); };
    bool	recognizesMateNameOrPdb(core::Symbol_sp name);

    core::List_sp	matesAsCons() { return core::Cons_O::createFromVec0(this->_Mates); };

    PlugWithMates_O( const PlugWithMates_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(PlugWithMates_O);
};



SMART(OutPlug);
class OutPlug_O : public PlugWithMates_O
{
    LISP_BASE1(PlugWithMates_O);
    LISP_CLASS(chem,ChemPkg,OutPlug_O,"OutPlug");
#if INIT_TO_FACTORIES
 public:
    static OutPlug_sp make(MatterName stubPivotAtom);
#else
    DECLARE_INIT();
#endif
public:
//    void	archiveBase(core::ArchiveP node);
	void	initialize();

private:
	/*! This contains the name of the atom that we will force to be the third atom
	  that defines the Stub of the Bond0 atom */
    MatterName	_StubPivotAtom;
public:

    bool hasStubPivotAtom() { return this->_StubPivotAtom.notnilp(); };
	MatterName getStubPivotAtom() { return this->_StubPivotAtom;};

	string descriptionOfContents() const;

	bool getIsIn() { return false;};

	OutPlug_O( const OutPlug_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(OutPlug_O);
};







SMART(InPlug);
class InPlug_O : public PlugWithMates_O
{
    LISP_BASE1(PlugWithMates_O);
    LISP_CLASS(chem,ChemPkg,InPlug_O,"InPlug");
public:
	void initialize();
private:
public:
	/*! Regular InPlugs use the first bonded atom as the root */
    virtual MatterName rootAtomName() const {return this->_B0;};


	string descriptionOfContents() const;

	bool getIsIn() { return true;};

	InPlug_O( const InPlug_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(InPlug_O);
};





SMART(JumpPlug);

/*! @class Define a Plug that is connected to the AtomTree through a Jump */

class JumpPlug_O : public Plug_O
{
    LISP_BASE1(Plug_O);
    LISP_CLASS(chem,ChemPkg,JumpPlug_O,"JumpPlug");
#if INIT_TO_FACTORIES
 public:
    static JumpPlug_sp make(MatterName jumpAtomName);
#else
    DECLARE_INIT();
#endif
public:
//    void	archiveBase(core::ArchiveP node);
    void	initialize();

private:
    MatterName		_JumpAtomName;
public:

    /*! Return the name of the root atom */
    virtual MatterName rootAtomName() const { return this->_JumpAtomName;};

	string descriptionOfContents() const;

		//! JumpPlugs are a kind of InPlugs
	bool getIsIn() { return true;};

	JumpPlug_O( const JumpPlug_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(JumpPlug_O);
};




SMART(RingClosingPlug);
class RingClosingPlug_O : public OutPlug_O
{
    LISP_BASE1(OutPlug_O);
    LISP_CLASS(chem,ChemPkg,RingClosingPlug_O,"RingClosingPlug");
#if INIT_TO_FACTORIES
 public:
    static RingClosingPlug_sp make(core::List_sp ringClosingMates);
#else
    DECLARE_INIT();
#endif
public:
//    void	archiveBase(core::ArchiveP node);
	void	initialize();

private:
    gctools::Vec0<RingClosingMate_sp>	_RingClosingMates;
public:
	string descriptionOfContents() const;
	bool getIsIn() { return false;};
	virtual bool getIsRingClosing() { return true; };
	bool recognizesRingClosingMate(core::Symbol_sp mateName);
	core::List_sp ringClosingMatesAsCons();
	RingClosingPlug_O( const RingClosingPlug_O& ss ); //!< Copy constructor
	DEFAULT_CTOR_DTOR(RingClosingPlug_O);
};








};

TRANSLATE(chem::Mate_O);
TRANSLATE(chem::RingClosingMate_O);
TRANSLATE(chem::Plug_O);
TRANSLATE(chem::PlugWithMates_O);
TRANSLATE(chem::OutPlug_O);
TRANSLATE(chem::InPlug_O);
TRANSLATE(chem::JumpPlug_O);
TRANSLATE(chem::RingClosingPlug_O);
#endif //]
