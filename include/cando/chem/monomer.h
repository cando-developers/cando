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
/*!
  \class Monomer_O
  \brief Monomer_O stores everything that a monomer needs to store
  Published signals:
  changed2D - sent if Position2D or Selected changes
  changed - sent if contents change
*/
  FORWARD(Isoname);
  class Isoname_O : public core::CxxObject_O {
    LISP_CLASS(chem,ChemPkg,Isoname_O,"Isoname",core::CxxObject_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    core::Symbol_sp  _Name;
    Fixnum           _Isomer;
    Topology_sp      _Topology;
  public:
    CL_LAMBDA(name topology &optional (isomer 0));
    CL_DEF_CLASS_METHOD static Isoname_sp make_isoname(core::Symbol_sp name, Topology_sp topology, Fixnum isomer) {
      GC_ALLOCATE_VARIADIC(Isoname_O,isoname,name,topology,isomer);
      return isoname;
    }
    CL_DEFMETHOD core::Symbol_sp name() const { return this->_Name; };
    CL_DEFMETHOD Fixnum isomer() const { return this->_Isomer; };
    CL_DEFMETHOD Topology_sp topology() const { return this->_Topology; };
  Isoname_O(core::Symbol_sp name, Topology_sp topology, Fixnum isomer ) : _Name(name), _Topology(topology), _Isomer(isomer) {};
  };
};



#if 0
namespace chem {
  class MonomerBase_O : public Entity_O
  {
    LISP_CLASS(chem,ChemPkg,MonomerBase_O,"MonomerBase",Entity_O);

    friend class Topology_O;
    friend class Constitution_O;
    friend class MonomerContext_O;
    friend class AlchemistState_O;
    friend class AlchemistFingerprint_O;
  public:
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
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
    CL_DEF_CLASS_METHOD static MonomerBase_sp makeMonomer() {
      GC_ALLOCATE(MonomerBase_O,mon);
      mon->setId(id);
      return mon;
    }
  public:
  public:
    MonomerBase_O() {};
  };
};
#endif


template <>
struct gctools::GCInfo<chem::Monomer_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {

  FORWARD(Monomer);
  class Monomer_O : public Entity_O
  {
    LISP_CLASS(chem,ChemPkg,Monomer_O,"Monomer",Entity_O);
  public:
    CL_DEF_CLASS_METHOD static Monomer_sp makeMonomer();
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    void initialize();
  public:
//	void	archiveBase(core::ArchiveP n);
    friend class	AlchemistState_O;
  public:
    typedef gctools::SmallMultimap<core::Symbol_sp,Coupling_sp,core::SymbolComparer> Couplings;
    typedef gctools::Vec0<Isoname_sp>	Monomers;
  public:
    core::Symbol_sp		_Id;
    uint			_SequenceNumber;
    Couplings			_Couplings;
    mutable uint		_CurrentMonomerIndex;
//    core::Symbol_sp		_GroupName;
    Monomers			_Monomers;
    Residue_sp			_TempResidue;
    int				_TemporaryInt;
//    adapt::SymbolSet_sp			_Aliases; //!< this will identify monomers for searches
  public:
//    virtual Residue_sp	createResidue();
    void setTempResidue(Residue_sp res)
    {
      this->_TempResidue = res;
    };
//    adapt::SymbolSet_sp getMonomerAliases()	{_OF(); ANN(this->_Aliases);return this->_Aliases; };
        //! Return true if inc worked and false if we carry
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
  public:
    core::List_sp plugNamesAndCouplingsAsList();

    adapt::SymbolSet_sp	plugNames() const;

#if 0
    void	setAliasesFromSymbolList(adapt::SymbolList_sp aliases);
    void	setAliasesFromCons(core::List_sp aliases);
//	void	setAliasesFromString(const string& s);
//    string	getAliasesAsString();
#endif
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
    Coupling_sp getInCoupling();
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
//    virtual adapt::ObjectSet_sp getAllAliases();
  public:
    Constitution_sp getConstitution();
    Plug_sp getPlugNamed(core::Symbol_sp s);

//! When a residue gets created, temporarily point to it
  public:	// accessible to AlchemistState
    virtual uint getMonomerIndex() { return this->_CurrentMonomerIndex; };
    virtual void setMonomerIndex(uint i) { this->_CurrentMonomerIndex = i; };
//    void	_expandGroupName();
  public:
        //!< Add another monomer name
    void addIsoname(Isoname_sp name);
  public:
    string __repr__() const;
    CL_DEFMETHOD virtual core::Symbol_sp monomerName() const;
    CL_DEFMETHOD virtual Fixnum getIsomer() const;
    //core::Symbol_sp getGroupName() const { return this->_GroupName; };
//    virtual AtomIndexer_sp getInterestingAtomIndexer();
//    void setGroupName(core::Symbol_sp groupName);
    CL_DEFMETHOD Isoname_sp	getOneMonomer() const;
    string description() const;
    uint numberOfPossibleMonomers() { return this->_Monomers.size(); };
    void randomizeMonomer();
    bool incrementMonomerIndex();
//    virtual core::List_sp allAtomAliases();
        //! Return true if we recognize the alias with the form ";[monomerAlias]@[atomAlias]"
//    CL_DEFMETHOD virtual bool recognizesAlias(Alias_sp alias);
  public:
  Monomer_O() : _CurrentMonomerIndex(0) {};
  };
};

#endif //]



