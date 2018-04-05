/*
    File: monomer.h
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

  FORWARD(RepresentedEntityNameSet);
  FORWARD(MonomerContext);
  FORWARD(CandoDatabase);
  FORWARD(AtomIndexer);
  FORWARD(ObjectSet);
  FORWARD(Alias);
  FORWARD(Plug);
  FORWARD(Coupling);
  FORWARD(DirectionalCoupling);
  FORWARD(RingCoupling);
  FORWARD(Oligomer);
  FORWARD(Topology);
  FORWARD(RingClosingPlug);
  FORWARD(StatusTracker);
  FORWARD(Constitution);
  FORWARD(SymbolSet);
  FORWARD(Monomer);
  FORWARD(MultiMonomer);
  FORWARD(MonoMonomer);
/*!
  \class Monomer_O
  \brief Monomer_O stores everything that a monomer needs to store
  Published signals:
  changed2D - sent if Position2D or Selected changes
  changed - sent if contents change
*/



  class Monomer_O : public Entity_O
  {
    LISP_CLASS(chem,ChemPkg,Monomer_O,"Monomer",Entity_O);

    friend class Topology_O;
    friend class Constitution_O;
    friend class MonomerContext_O;
    friend class AlchemistState_O;
    friend class AlchemistFingerprint_O;
  public:
    typedef gctools::SmallMultimap<core::Symbol_sp,Coupling_sp,core::SymbolComparer> Couplings;
  public:
//	void	archiveBase(core::ArchiveP node);
  protected: // archive
        /*! When a monomer is defined using "monomer" and "link..monomer" commands the Id is
         * used to uniquely refer to the monomer
         */
    core::Symbol_sp			_Id;
    uint				_SequenceNumber;
  public: // archive
    Couplings			_Couplings;
    Residue_sp			_TempResidue;
    int				_TemporaryInt;
  protected: // archive
        /*! Aliases identify monomers for searches.
         * Multiple monomers can share aliases as long as the atoms that are to be
         * accessed through those aliases are have unique names across those monomers.
         */
    adapt::SymbolSet_sp			_Aliases; //!< this will identify monomers for searches
  protected:
    friend class Oligomer_O;
    friend class Coupling_O;
  public:
    CL_DEF_CLASS_METHOD static Monomer_sp makeMonomer(core::Symbol_sp id) {
      GC_ALLOCATE(Monomer_O,mon);
      mon->setId(id);
      return mon;
    }
  public:
//    virtual Residue_sp	createResidue();
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

    CL_LISPIFY_NAME("setId");
    CL_DEFMETHOD 	void	setId(core::Symbol_sp id) { this->_Id = id; };
    CL_LISPIFY_NAME("getId");
    CL_DEFMETHOD 	core::Symbol_sp getId() { return this->_Id; };
//    Topology_sp	getTopology();
    RingClosingPlug_sp getMissingRingClosingPlug(Monomer_sp mate);
    CL_LISPIFY_NAME("getTemporaryResidue");
    CL_DEFMETHOD 	Residue_sp	getTemporaryResidue()
    { _G();
      ASSERTNOTNULL(this->_TempResidue);
      return this->_TempResidue;
    }
    uint	getTemporaryInt() { return this->_TemporaryInt; }
    void	setTemporaryInt(uint i) { this->_TemporaryInt = i; }
    bool hasTemporaryResidue();
    virtual core::List_sp allAtomAliases() {_OF();SUBCLASS_MUST_IMPLEMENT();};

        //! Return true if we recognize the alias with the form "![monomerAlias]@[atomAlias]"
    CL_LISPIFY_NAME("recognizesAlias");
    CL_DEFMETHOD 	virtual bool recognizesAlias(Alias_sp alias) { return false; };
  public:
    void	initialize();

  public:
    core::List_sp plugNamesAndCouplingsAsList();

    adapt::SymbolSet_sp	plugNames() const;

    void	setAliasesFromSymbolList(adapt::SymbolList_sp aliases);
    void	setAliasesFromCons(core::List_sp aliases);
//	void	setAliasesFromString(const string& s);
    string	getAliasesAsString();

  public:

    void	setSequenceNumber(uint idx) { this->_SequenceNumber = idx;};
    CL_LISPIFY_NAME("getSequenceNumber");
    CL_DEFMETHOD 	uint	getSequenceNumber() { return this->_SequenceNumber; };

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

    CL_DEFMETHOD void	removeCoupling(Coupling_sp coup);
    CL_DEFMETHOD bool	hasCouplingWithPlugName(core::Symbol_sp s);
    CL_DEFMETHOD Coupling_sp getCouplingWithPlugName(core::Symbol_sp s);
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

    virtual uint	numberOfPossibleMonomers() { return 1; };
    virtual void randomizeMonomer() { /*Do nothing, subclasses with multiple monomers will do something */ };

    virtual AtomIndexer_sp getInterestingAtomIndexer() {_OF(); SUBCLASS_MUST_IMPLEMENT();};

    virtual adapt::ObjectSet_sp getAllAliases();


  public:
    CL_LISPIFY_NAME("getName");
    CL_DEFMETHOD         virtual core::Symbol_sp getName() const {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
    virtual core::Symbol_sp getGroupName() const {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
    virtual void setGroupName(core::Symbol_sp name) {_OF(); SUBCLASS_MUST_IMPLEMENT(); };
    Constitution_sp getConstitution();
    Plug_sp getPlugNamed(core::Symbol_sp s);
  public:
    Monomer_O() {};
  };

  FORWARD(MonoMonomer);
  class MonoMonomer_O : public Monomer_O
  {
    LISP_CLASS(chem,ChemPkg,MonoMonomer_O,"MonoMonomer",Monomer_O);
  public:
    static MonoMonomer_sp make(core::Symbol_sp name);
  public:
//	void	archiveBase(core::ArchiveP n);
    friend class	AlchemistState_O;
  private:
    core::Symbol_sp                 _Name;
    core::Symbol_sp			_GroupName;
        //! When a residue gets created, temporarily point to it
  protected:	// accessible to AlchemistState
    virtual uint getMonomerIndex() { return 0; };
    void	_expandGroupName();
  public:
    virtual void checkForErrorsAndUnknownContexts(CandoDatabase_sp bdb);

    CL_DEFMETHOD void setName(core::Symbol_sp name) { this->_Name = name; };
    virtual core::Symbol_sp getName() const { return this->_Name;};
    core::Symbol_sp getGroupName() const { return this->_GroupName; };
    virtual AtomIndexer_sp getInterestingAtomIndexer() {IMPLEMENT_ME();};
    void setGroupName(core::Symbol_sp groupName) { this->_GroupName = groupName; };
    core::Symbol_sp	getOneMonomer() const { return this->_Name;};
    uint numberOfPossibleMonomers() { return 1; };
    virtual core::List_sp allAtomAliases() { IMPLEMENT_ME(); };

        //! Return true if we recognize the alias with the form ";[monomerAlias]@[atomAlias]"
    virtual bool recognizesAlias(Alias_sp alias) {IMPLEMENT_ME(); };
  public:
    DEFAULT_CTOR_DTOR(MonoMonomer_O);
  };

};


template <>
struct gctools::GCInfo<chem::MultiMonomer_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {

    FORWARD(MultiMonomer);
    class MultiMonomer_O : public Monomer_O
    {
        LISP_CLASS(chem,ChemPkg,MultiMonomer_O,"MultiMonomer",Monomer_O);
#if INIT_TO_FACTORIES
    public:
        static MultiMonomer_sp make(core::Symbol_sp name);
#else
        DECLARE_INIT();
#endif
    public:
	void initialize();
    public:
//	void	archiveBase(core::ArchiveP n);
	friend class	AlchemistState_O;
    public:
        typedef gctools::Vec0<core::Symbol_sp>	Monomers;
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
        //!< Add another monomer name
	void addMonomerName(core::Symbol_sp name);
    public:
        virtual void checkForErrorsAndUnknownContexts(CandoDatabase_sp bdb);

        virtual core::Symbol_sp getName() const;
        core::Symbol_sp getGroupName() const { return this->_GroupName; };
        virtual AtomIndexer_sp getInterestingAtomIndexer();
        void setGroupName(core::Symbol_sp groupName);
        core::Symbol_sp	getOneMonomer() const;
        string description() const;
        uint numberOfPossibleMonomers() { return this->_Monomers.size(); };
        void randomizeMonomer();
        bool incrementMonomerIndex();

        virtual core::List_sp allAtomAliases();

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

#endif //]
