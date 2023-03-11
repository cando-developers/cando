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
#include <cando/chem/plug.fwd.h>
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
      auto isoname = gctools::GC<Isoname_O>::allocate(name,topology,isomer);
      return isoname;
    }
    string __repr__() const;
    CL_DEFMETHOD core::Symbol_sp name() const { return this->_Name; };
    CL_DEFMETHOD Fixnum isomer() const { return this->_Isomer; };
    CL_DEFMETHOD Topology_sp topology() const { return this->_Topology; };
    Isoname_O(core::Symbol_sp name, Topology_sp topology, Fixnum isomer ) : _Name(name), _Isomer(isomer), _Topology(topology) {};
  };
};



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
  static Monomer_sp makeMonomer(core::List_sp topology_list);
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
  void initialize();
public:
  friend class	AlchemistState_O;
public:
  typedef gctools::SmallMultimap_uncopyable<core::Symbol_sp,Coupling_sp,core::SymbolComparer> Couplings;
public:
  core::Symbol_sp		_Id;
  size_t			_SequenceNumber;
  Couplings			_Couplings;
  mutable size_t		_CurrentMonomerIndex;
  mutable size_t                _CurrentStereoisomerOffset;
  gctools::Vec0_uncopyable<core::Symbol_sp>    _Monomers;
//    adapt::SymbolSet_sp			_Aliases; //!< this will identify monomers for searches
public:
//    adapt::SymbolSet_sp getMonomerAliases()	{_OF(); ANN(this->_Aliases);return this->_Aliases; };
        //! Return true if inc worked and false if we carry
public:
	/*! Return true if we have a RingClosingOutPlug */
  bool hasRingClosingOutPlug() const;

  CL_LISPIFY_NAME("setId");
  CL_DEFMETHOD 	void	setId(core::Symbol_sp id) { this->_Id = id; };
  CL_LISPIFY_NAME("getId");
  CL_DEFMETHOD 	core::Symbol_sp getId() { return this->_Id; };
  Topology_sp	currentTopology();
  core::T_sp getMissingRingClosingPlug(Monomer_sp mate);
public:
  core::List_sp plugNamesAndCouplingsAsList();

  CL_DOCSTRING("Return a list of all of the plug names that this monomer has.");
  CL_DEFMETHOD core::List_sp plugNamesAsList() const;

public:

  void	setSequenceNumber(uint idx) { this->_SequenceNumber = idx;};
  CL_LISPIFY_NAME("getSequenceNumber");
  CL_DEFMETHOD 	uint	getSequenceNumber() { return this->_SequenceNumber; };

//	bool	sameAs(Monomer_sp mon) { return this == mon.get(); };
  string	getAddress();

  virtual void updateGroupName() { return; };


  uint	numberOfCouplings() const;

  Monomer_sp copyAndRemoveCouplings() const;
  
  void eraseCouplings();
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
  size_t numberOfStereoisomers() const;
//  void randomizeMonomer();
  bool incrementMonomerIndex();
  virtual size_t getMonomerIndex() const;
  virtual void setMonomerIndex(size_t i);
//    void	_expandGroupName();
public:
        //!< Add another monomer name
  void addTopologyName(core::Symbol_sp name);
public:
  string __repr__() const;
//  CL_DEFMETHOD virtual core::Symbol_sp monomerName() const;
  core::Symbol_sp currentStereoisomerName() const;
//  core::Symbol_sp stereoisomerName() const;
  CL_DEFMETHOD virtual Fixnum getStereoisomerOffset() const { return this->_CurrentStereoisomerOffset;}; // getIsomer
    //core::Symbol_sp getGroupName() const { return this->_GroupName; };
//    virtual AtomIndexer_sp getInterestingAtomIndexer();
//    void setGroupName(core::Symbol_sp groupName);
//  CL_DEFMETHOD Isoname_sp	getIsoname() const;
  string description() const;
//    virtual core::List_sp allAtomAliases();
        //! Return true if we recognize the alias with the form ";[monomerAlias]@[atomAlias]"
//    CL_DEFMETHOD virtual bool recognizesAlias(Alias_sp alias);
public:
  Monomer_O(const Monomer_O& other);
  Monomer_O() : _Id(nil<core::T_O>()), _CurrentMonomerIndex(0), _CurrentStereoisomerOffset(0) {};
};
};

#endif //



