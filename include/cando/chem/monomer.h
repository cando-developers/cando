       
       
//
// (C) 2004 Christian E. Schafmeister
//




/*
 *	monomer.h
 *
 *	A Constitution maintains a Residue and the names of atoms and links
 *	that allow Monomers to be linked together into an Oligomer.
 *
 *	A Constitution can't be actually connected together, rather it can be
 *	used to create a Monomer that can be connected to other Monomers.
 *
 *	Member function "monomerWithName" returns a new monomer
 *		stamped from this constitution
 *		all atom pointers redirected into the new monomer.
 *		When connecting monomers together use only "copy"d monomers and
 *		leave the original alone.
 *
 */

#ifndef Monomer_H  //[
#define Monomer_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/holder.h>
#include <cando/geom/vector2.h>
#include <cando/geom/vector3.h>
#include <cando/chem/entity.h>

#include <cando/chem/chemPackage.h>
#include <cando/adapt/symbolList.fwd.h>
#include <cando/adapt/symbolSet.fwd.h>
#include <cando/adapt/objectSet.fwd.h>

#include <cando/chem/residue.fwd.h>


namespace chem {

    typedef core::Symbol_sp MonomerName;

    SMART(RepresentedEntityNameSet);
    SMART(MonomerContext);
    SMART(CandoDatabase);
    SMART(AtomIndexer);
    SMART(ObjectSet);
    SMART(Alias);
    SMART(Plug);
    SMART(Coupling);
    SMART(DirectionalCoupling);
    SMART(RingCoupling);
    SMART(Oligomer);
    SMART(Topology);
    SMART(RingClosingPlug);
    SMART(StatusTracker);
    SMART(Constitution);
    SMART(SymbolSet);

/*!
  \class Monomer_O
  \brief Monomer_O stores everything that a monomer needs to store
  Published signals:
  changed2D - sent if Position2D or Selected changes
  changed - sent if contents change
*/



    SMART(Monomer);
    class Monomer_O : public Entity_O
    {
        LISP_BASE1(Entity_O);
        LISP_CLASS(chem,ChemPkg,Monomer_O,"Monomer");

	friend class Topology_O;
	friend class Constitution_O;
	friend class MonomerContext_O;
	friend class AlchemistState_O;
	friend class AlchemistFingerprint_O;
    public:
        typedef gctools::SmallMultimap<core::Symbol_sp,Coupling_sp,core::SymbolComparer> Couplings;
    public:
//	void	archiveBase(core::ArchiveP node);
    private:
	Oligomer_sp			_Oligomer;
    protected: // archive
        /*! When a monomer is defined using "monomer" and "link..monomer" commands the Id is
         * used to uniquely refer to the monomer
         */
	core::Symbol_sp			_Id;
	string				_Comment; //!< this will help identify monomers that have errors
	uint				_SequenceNumber;
    public: // archive
	Couplings			_Couplings;
	Residue_wp			_TempResidue;
	int				_TemporaryInt;
    protected: // archive
        /*! Aliases identify monomers for searches.
         * Multiple monomers can share aliases as long as the atoms that are to be
         * accessed through those aliases are have unique names across those monomers.
         */
	adapt::SymbolSet_sp			_Aliases; //!< this will identify monomers for searches
	geom::Vector2				_Position2D;
	bool				_Selected;
    public: // Do not archive/in error information
        /*! If there is a problem with anythin register errors here */
	bool				_Verbose;
	StatusTracker_sp			_Status;
    protected:
	friend class Oligomer_O;
	friend class Coupling_O;
	virtual Residue_sp	createResidue() {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
	void setTempResidue(Residue_sp res)
	{
	    this->_TempResidue = res;
	};
	adapt::SymbolSet_sp getMonomerAliases()	{_OF(); ANN(this->_Aliases);return this->_Aliases; };
	virtual uint getMonomerIndex() {return 0;};
	virtual void setMonomerIndex(uint i) {_OF(); SUBCLASS_MUST_IMPLEMENT();};
        //! Return true if inc worked and false if we carry
	virtual	bool incrementMonomerIndex() {_OF(); SUBCLASS_MUST_IMPLEMENT();};

    public:

	/*! Return true if we have a RingClosingOutPlug */
	bool hasRingClosingOutPlug() const;

	void	setId(core::Symbol_sp id) { this->_Id = id; };
	core::Symbol_sp getId() { return this->_Id; };
	void		setVerbose(bool b) {this->_Verbose = b;};
	bool		isVerbose() { return this->_Verbose;};
	void		setSelected(bool b) { this->_Selected = b; };
	bool		isSelected() { return this->_Selected; };
	Topology_sp	getTopology();
	RingClosingPlug_sp getMissingRingClosingPlug(Monomer_sp mate);
	Residue_sp	getTemporaryResidue()
	{ _G();
	    ASSERTNOTNULL(this->_TempResidue);
	    return this->_TempResidue;
	}
	uint	getTemporaryInt() { return this->_TemporaryInt; }
	void	setTemporaryInt(uint i) { this->_TemporaryInt = i; }
	bool hasTemporaryResidue();
        virtual adapt::SymbolList_sp allAtomAliases() {_OF();SUBCLASS_MUST_IMPLEMENT();};

        //! Return true if we recognize the alias with the form "![monomerAlias]@[atomAlias]"
	virtual bool recognizesAlias(Alias_sp alias) { return false; };
    public:
	void	initialize();

    public:
	core::List_sp plugNamesAndCouplingsAsCons();

	adapt::SymbolSet_sp	plugNames() const;

	void	setComment(const string& s) { this->_Comment = s;};
	string	getComment() { return this->_Comment; };

	void	setAliasesFromSymbolList(adapt::SymbolList_sp aliases);
	void	setAliasesFromCons(core::List_sp aliases);
//	void	setAliasesFromString(const string& s);
	string	getAliasesAsString();

    public:

	Oligomer_sp	getOligomer() {_OF();
	    return this->_Oligomer;
	};

	void	setPosition2D_xy(double x, double y) { this->_Position2D.set(x,y);}; // this->notify(Monomer_changed2D);;
	void	setPosition2D(const geom::Vector2& pos ) { this->_Position2D = pos;}; //  this->notify(Monomer_changed2D);;
	geom::Vector2 getPosition2D() { return this->_Position2D; };
	void	setSequenceNumber(uint idx) { this->_SequenceNumber = idx;};
	uint	getSequenceNumber() { return this->_SequenceNumber; };

//	bool	sameAs(Monomer_sp mon) { return this == mon.get(); };
	string	getAddress();

	virtual void updateGroupName() { return; };


        uint	numberOfCouplings() const;

	void	addCoupling(core::Symbol_sp plugName, Coupling_sp out);
        /*! Generate the in-coupling name and add this
         * coupling by that name.
         * Make sure that no other in couplings are defined.
         */
	void	resetInCoupling();
	void	setInCoupling(Coupling_sp out);
        /*! Generate the out-coupling name and add this
         * coupling by that name.
         * Make sure that no other out couplings with this
         * name are defined.
         */
	void	addOutCoupling(Coupling_sp out);
	void	removeCouplingToMonomer(Monomer_sp mon);

        /*! Fix the plug name key for the coupling.
         * This has to be done when a coupling name is changed.
         */
	void	fixPlugNameForCoupling(Coupling_sp coup);

	void	removeCoupling(Coupling_sp coup);
	bool	hasCouplingWithPlugName(core::Symbol_sp s);
	Coupling_sp getCouplingWithPlugName(core::Symbol_sp s);
        /*! Return true if this monomer has a plug name that points to the coupling
         */
	bool	hasMatchingPlugNameAndCoupling(core::Symbol_sp plugName, Coupling_sp coup);
	DirectionalCoupling_sp getInCoupling();
	core::Symbol_sp getInCouplingName();
	bool	hasInCoupling();
	bool	hasOutCouplings();
//	bool	hasOutCouplingWithName(const string& name);

        gctools::Vec0<Coupling_sp> getOutCouplings();


        /*! Return true if the current monomer context is a valid one
         */
	bool	isMonomerContextValid();

        /*! Return the current monomer context, throw error if its invalid
         * Use the specific names of the monomers
         */
	MonomerContext_sp getSpecificMonomerContext();

        /*! Return the current monomer context, throw error if its invalid
         * Use the general/group names of the monomers if defined
         */
	MonomerContext_sp getGeneralMonomerContext();

	void	throwIfBadConnections();
	bool	checkForBadConnections();

	virtual void checkForErrorsAndUnknownContexts(CandoDatabase_sp bdb);
	void	clearError();
	bool	getHasError();
	string	getStatusMessage();
	void	addErrorMessage(const string& s);
	void	addStatusMessage(const string& s);


	StatusTracker_sp getStatusTracker();

        virtual uint	numberOfPossibleMonomers() { return 1; };
        virtual void randomizeMonomer() { /*Do nothing, subclasses with multiple monomers will do something */ };

	void	setOligomer(Oligomer_sp o );
	virtual AtomIndexer_sp getInterestingAtomIndexer() {_OF(); SUBCLASS_MUST_IMPLEMENT();};

        virtual adapt::ObjectSet_sp getAllAliases();


    public:
        virtual core::Symbol_sp getName() const {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
        virtual core::Symbol_sp getGroupName() const {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
        virtual void setGroupName(core::Symbol_sp name) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
        Constitution_sp getConstitution();
        Plug_sp getPlugNamed(core::Symbol_sp s);
    public:


        DEFAULT_CTOR_DTOR(Monomer_O);
    };


    SMART(OneMonomer);
    class OneMonomer_O : public core::CxxObject_O
    {
        LISP_BASE1(core::CxxObject_O);
        LISP_CLASS(chem,ChemPkg,OneMonomer_O,"OneMonomer");
        friend class AlchemistFingerprint_O;
    public:
//    void	archiveBase(core::ArchiveP node);
    private: // archive
        core::Symbol_sp	_Name;
    public:
    public:
        core::Symbol_sp	getName() const { return this->_Name;};
        void	setName(core::Symbol_sp s) { this->_Name = s; };

	DEFAULT_CTOR_DTOR(OneMonomer_O);
    };




    FORWARD(MultiMonomer);
    class MultiMonomer_O : public Monomer_O
    {
        LISP_BASE1(Monomer_O);
        LISP_CLASS(chem,ChemPkg,MultiMonomer_O,"MultiMonomer");
#if INIT_TO_FACTORIES
    public:
        static MultiMonomer_sp make(core::Symbol_sp name, const string& comment);
#else
        DECLARE_INIT();
#endif
    public:
	void initialize();
    public:
//	void	archiveBase(core::ArchiveP n);
	friend class	AlchemistState_O;
    public:
        typedef gctools::Vec0<OneMonomer_sp>	Monomers;
    private:
	mutable uint			_CurrentMonomerIndex;
	core::Symbol_sp			_GroupName;
	Monomers			_Monomers;
        //! When a residue gets created, temporarily point to it

    protected:	// accessible to AlchemistState
	virtual uint getMonomerIndex() { return this->_CurrentMonomerIndex; };
	virtual void setMonomerIndex(uint i) { this->_CurrentMonomerIndex = i; };
	void	_expandGroupName();


    public:
//	void setCandoDatabase(CandoDatabase_sp hold);

        //!< Add another monomer name
	void addMonomerName(core::Symbol_sp name);
	Residue_sp createResidue();

//	void	updateGroupName() { this->_expandGroupName(); };


    public:
        virtual void checkForErrorsAndUnknownContexts(CandoDatabase_sp bdb);

        virtual core::Symbol_sp getName() const;
        core::Symbol_sp getGroupName() const { return this->_GroupName; };
        virtual AtomIndexer_sp getInterestingAtomIndexer();
        void setGroupName(core::Symbol_sp groupName);
        OneMonomer_sp	getOneMonomer() const;
        OneMonomer_sp	getOneMonomerWithName(core::Symbol_sp monomerName) const;
        string description() const;
        uint numberOfPossibleMonomers() { return this->_Monomers.size(); };
        void randomizeMonomer();
        bool incrementMonomerIndex();


        virtual adapt::SymbolList_sp allAtomAliases();

        //! Return true if we recognize the alias with the form ";[monomerAlias]@[atomAlias]"
        virtual bool recognizesAlias(Alias_sp alias);


    public:


    public:

        DEFAULT_CTOR_DTOR(MultiMonomer_O);
    };





#if 0
    extern MultiMonomer_sp create_MultiMonomer(CandoDatabase_sp bdb);

    inline MultiMonomer_sp create_MultiMonomer(CandoDatabase_sp bdb, const string& name, const string& comment )
    {
        MultiMonomer_sp mon;
        mon = create_MultiMonomer(bdb);
        mon->setGroupName(name);
        mon->setComment(comment);
        return mon;
    }
    // Thin wrappers
    inline MultiMonomer_sp create_MultiMonomer1(CandoDatabase_sp bdb)
    {
        return create_MultiMonomer(bdb);
    }
    // Thin wrappers
    inline MultiMonomer_sp create_MultiMonomer2(CandoDatabase_sp bdb, const string& name, const string& comment)
    {
        return create_MultiMonomer(bdb,name,comment);
    }
#endif




};



TRANSLATE(chem::Monomer_O);
TRANSLATE(chem::OneMonomer_O);
TRANSLATE(chem::MultiMonomer_O);
#endif //]
