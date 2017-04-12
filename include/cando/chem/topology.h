/*
    File: topology.h
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





#ifndef Topology_H
#define Topology_H
#include <clasp/core/common.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/adapt/symbolSet.fwd.h>
#include <cando/chem/bond.h>
#include <cando/chem/candoDatabase.h>
#include <cando/geom/vector2.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/stringSet.fwd.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/stereoisomerAtoms.fwd.h>
//#include "kinematics/atomTemplate.fwd.h"
//#include "kinematics/chiDihedrals.fwd.h"
#include <cando/chem/topology.fwd.h>
#include <cando/chem/constitutionAtoms.fwd.h>


#include <cando/chem/chemPackage.h>


namespace chem
{

  SMART(MonomerContext);
  SMART(Constitution);
  SMART(ExtractFragment);
  SMART(Plug);
  SMART(CandoDatabase);
  SMART(ExtractScaffold);
  SMART(Atom);
  SMART(Residue);
  SMART(Monomer);
  SMART(RingClosingPlug);
  SMART(ExtractCoreFragment);

  SMART(RequiredPlug);

  SMART(FrameBase);

  SMART(Topology);
  SMART(TopologyAtomInfo);

};


 template <>
struct gctools::GCInfo<chem::TopologyAtomInfo_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

 namespace chem {
   class TopologyAtomInfo_O : public core::CxxObject_O {
     LISP_CLASS(chem,ChemPkg,TopologyAtomInfo_O,"TopologyAtomInfo",core::CxxObject_O);
   public:
     core::T_sp _fftype;
     double     _AtomicCharge;
     ConstitutionAtomIndex0N _ConstitutionAtomIndex;
public:
   TopologyAtomInfo_O(core::T_sp fftype, double atomicCharge, ConstitutionAtomIndex0N index ) : _fftype(fftype), _AtomicCharge(atomicCharge), _ConstitutionAtomIndex(index) {};
   public:
     bool fieldsp() const { return true;};
     void fields(core::Record_sp node);
   public:
     CL_LISPIFY_NAME(make_topology_atom_info);
     CL_DEF_CLASS_METHOD static inline TopologyAtomInfo_sp make(core::T_sp fftype, double atomicCharge, chem::ConstitutionAtomIndex0N constitutionAtomIndex) {
     return gctools::GC<TopologyAtomInfo_O>::allocate(fftype,atomicCharge, constitutionAtomIndex);
   }
   public:
     inline CL_DEFMETHOD core::T_sp fftype() const { return this->_fftype; };
     inline CL_DEFMETHOD double atomicCharge() const { return this->_AtomicCharge; };
     inline CL_DEFMETHOD ConstitutionAtomIndex0N constitutionAtomIndex() const { return this->_ConstitutionAtomIndex; };
   };
 };




/*! @class A Topology describes one way that a Constitution can be connected to other Constitutions.
  It maintains a list of Plugs that describe how this Topology can plug into other Constitutions/Topologys.
  For example, a bis-amino acid can plug into two other bis-amino acids through DKPs on the leading end
  or the trailing end.  The same bis-amino acid could also have an incoming DKP and an outgoing peptide-amide
  and something on the outgoing ester.  The same bis-amino acid could also have an incoming peptide-amide
  an outgoing peptide-amide and something else connected to the trailing secondary amine and the leading 
  carbonyl group.  Each of these scenarios is described by a different Topology.
  A Topology also contains ExtractScaffold and ExtractFragment objects that know how to extract coordinates 
  from a Residue in the given Topology.

  A Topology also defines a map of atom names to atomic types and atomic charges. Atom types and charges will
  be topology dependent but the element is not and so the Constitution provides a map of atom names to atomic elements.

  I'm disabling the following for now - until I see a reason to turn it back on.
Each Topology is uniquely identified within a CandoDatabase by a TopologyIndex0N identifier that ranges
  between 0 and N-(the number of unique Topologys in the database).  This is to facilitate atom based lookups
  that describe the number of bonds between atoms in two interconnected Topologys.
*/
template <>
struct gctools::GCInfo<chem::Topology_O> {
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {
class Topology_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,Topology_O,"Topology",core::CxxObject_O); 
public:
public:
    void initialize();
private:
    typedef Plug_sp	plugType;
    typedef Plug_O	plugOType;
    typedef adapt::SymbolMap<Plug_O> Plugs;
private:
    core::Symbol_sp				_Name;
    Constitution_sp				_Constitution;
    bool					_SuppressTrainers;
    core::Vector_sp                             _AtomInfo;
    adapt::SymbolMap<StereoisomerAtoms_O>      	_StereoisomerAtomProperties;
    core::HashTableEq_sp			_Properties;
    core::Symbol_sp                             _DefaultOutPlugName;
    Plugs			_Plugs;
 public:
    bool fieldsp() const { return true;};
    void fields(core::Record_sp node);
public:
    CL_LISPIFY_NAME(make_topology);
    CL_DOCSTRING(R"doc(Create a topology and return it, after this the topology needs to be added to a constitution)doc");
    CL_DEF_CLASS_METHOD static Topology_sp make(core::Symbol_sp name, Constitution_sp constitution, core::List_sp plugs );
    
CL_LISPIFY_NAME(makeTopologyFromResidue);
CL_LAMBDA(residue topology_name &optional constitution);
 CL_DOCSTRING("Create a topology from a residue. The constitution may be NIL or a constitution to use to define the atoms");
 CL_DEF_CLASS_METHOD static Topology_mv makeTopologyFromResidue(chem::Residue_sp residue, core::Symbol_sp topologyName, core::T_sp constitution);
protected:
    void setFromMonomer(Monomer_sp mon);
public:
    void setConstitution(Constitution_sp c);
public:
    /* You can attach properties to the Topology later */
    core::HashTableEq_sp properties() const;

    string description() const;
    CL_DEFMETHOD Constitution_sp	getConstitution() const { return this->_Constitution; };
    MonomerContext_sp getMonomerContext(CandoDatabase_sp bdb);
    //! Return if we suppress trainers
CL_LISPIFY_NAME("suppressTrainers");
CL_DEFMETHOD     bool suppressTrainers() const { return this->_SuppressTrainers;};

 Residue_sp build_residue() const;

 string __repr__() const;
 
    //! Return all plugs as a Cons
    core::List_sp	plugsAsList();

    //! Return a Cons of Plugs that have Mates
    core::List_sp	plugsWithMatesAsList();

    //! Return a Cons of out plugs
    core::List_sp	outPlugsAsList();

    /// @brief Return all of the Incomplete Frames as Cons
    /// @return Cons of Incomplete Frames
    core::List_sp	incompleteFramesAsList();

    /// @brief Return the name of the Topology
CL_LISPIFY_NAME("getName");
CL_DEFMETHOD     core::Symbol_sp	getName() const {return (this->_Name);};
CL_LISPIFY_NAME("setName");
 CL_DEFMETHOD     void setName(core::Symbol_sp n) {this->_Name = n;};

    /// @brief Return true if this Topology has a plug named (name)
    bool	hasPlugNamed(core::Symbol_sp name);
    /// @brief Return the Plug with the name
    Plug_sp	plugNamed(core::Symbol_sp name);

    void setTemporaryObject(core::T_sp o);
    core::T_sp getTemporaryObject();

CL_LISPIFY_NAME("numberOfPlugs");
CL_DEFMETHOD     int	numberOfPlugs() { return this->_Plugs.size(); };

 CL_DEFMETHOD bool	hasInPlug();
 CL_DEFMETHOD plugType getInPlug();

 CL_DEFMETHOD core::Symbol_sp defaultOutPlugName() const { return this->_DefaultOutPlugName; };
 CL_DEFMETHOD void setf_default_out_plug_name(core::Symbol_sp outPlugName) {
   this->_DefaultOutPlugName = outPlugName;
 }

    /*! Return true if this Topology has all of the plugs in (plugSet) */
    bool	hasMatchingPlugsWithMates(adapt::SymbolSet_sp plugSet);

    bool	matchesMonomerEnvironment( Monomer_sp mon );
    RingClosingPlug_sp provideMissingRingClosingPlug( Monomer_sp mon );

    CL_DEFMETHOD core::T_sp atomInfo() const { return this->_AtomInfo; };
    CL_DEFMETHOD core::T_sp setf_atomInfo(core::Vector_sp ai) { this->_AtomInfo = ai; return ai; };

    void mapPlugs(std::function<void(Plug_sp)> );

    core::List_sp extractFragmentsAsList();

    void throwIfExtractFragmentsAreNotExclusive(ConstitutionAtoms_sp constitutionAtoms);

    CL_DEFMETHOD void addPlug(core::Symbol_sp nm, plugType op ) { this->_Plugs.set(nm, op);};
    bool		hasFlag(core::Symbol_sp f) const;
    bool		matchesTopology(Topology_sp cm);
    bool		matchesContext(MonomerContext_sp cm);

    /*! Return a StereoisomerAtoms_sp object for the stereoisomer name
      Either return the existing one or create a new one */
    StereoisomerAtoms_sp lookupOrCreateStereoisomerAtoms(core::Symbol_sp stereoisomerName);
    
 Topology_O() : _Name(_Nil<core::Symbol_O>()), _Constitution(_Unbound<chem::Constitution_O>()) {};
 Topology_O(core::Symbol_sp name, Constitution_sp constitution) : _Name(name), _Constitution(constitution) {};
};


 void connect_residues(Topology_sp prev_topology,
                               Residue_sp prev_residue,
                               core::Symbol_sp out_plug_name,
                               Topology_sp next_topology,
                               Residue_sp next_residue,
                       core::Symbol_sp in_plug_name);

};

#endif
