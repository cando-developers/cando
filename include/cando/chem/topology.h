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
//#include "kinematics/pointTemplate.fwd.h"
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
  SMART(ExtractCoreFragment);
  SMART(RequiredPlug);
  SMART(FrameBase);
  SMART(Topology);
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
  CL_DOCSTRING(R"(A Topology describes one way that a Constitution can be connected to other Constitutions.
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
)");
  LISP_CLASS(chem,ChemPkg,Topology_O,"Topology",core::CxxObject_O); 
public:
  void initialize();
public:
  typedef Plug_sp	plugType;
  typedef Plug_O	plugOType;
  typedef gc::SmallMap<core::Symbol_sp,Plug_sp> Plugs;
public:
  core::Symbol_sp			      _Name;
  Constitution_sp			      _Constitution;
  core::Symbol_sp                             _StereoisomerCoding; // kw::_sym_absolute or kw::_sym_coded
  gctools::Vec0<StereoisomerAtoms_sp>         _StereoisomerAtomProperties;
  core::List_sp 			      _Properties;
  core::Symbol_sp                             _DefaultOutPlugName;
  Plugs			                      _Plugs;
public:
  bool fieldsp() const { return true;};
  void fields(core::Record_sp node);
public:
  CL_LISPIFY_NAME(make_topology);
  CL_DOCSTRING(R"dx(Create a topology and return it, after this the topology needs to be added to a constitution)dx");
  CL_DEF_CLASS_METHOD static Topology_sp make(core::Symbol_sp name, Constitution_sp constitution, core::List_sp plugs );
    
  static Topology_mv makeTopologyFromResidue(chem::Residue_sp residue, core::Symbol_sp topologyName, core::T_sp constitution);
protected:
  void setFromMonomer(Monomer_sp mon);
public:
  void setConstitution(Constitution_sp c);
		/*! Remove the property from this Matters property list
		 */
  void	clearProperty(core::Symbol_sp propertySymbol);
		/*! Define/set the value of the property.
		 * If it already exists it is overwritten.
		 */
  void	setProperty(core::Symbol_sp propertySymbol, core::T_sp value);
		/*! Return the value of the property.
		 * Throw an exception if the property isn't defined.
		 */
  core::T_sp getProperty(core::Symbol_sp propertySymbol );
		/*! Return the value of the property or the
		 * default if it isn't defined.
		 */
  core::T_sp getPropertyOrDefault(core::Symbol_sp propertySymbol, core::T_sp defVal );
		/*! Return true if the property exists.
		 */
  bool	hasProperty(core::Symbol_sp propertySymbol );

public:

  void walkStereoisomerAtoms(core::Function_sp func);
  
  string description() const;
  CL_DEFMETHOD Constitution_sp	getConstitution() const { return this->_Constitution; };
  MonomerContext_sp getMonomerContext(CandoDatabase_sp bdb);

    Residue_sp buildResidueSingleName() const;
    Residue_sp buildResidueForIsomer(size_t isomer) const;
    Residue_sp buildResidueForMonomerName(core::Symbol_sp monomerName) const;
    Residue_sp buildResidueForIsoname(Isoname_sp isoname) const;

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
  CL_DEFMETHOD     core::Symbol_sp	getName() const {return (this->_Name);};
  CL_DEFMETHOD     void setName(core::Symbol_sp n) {this->_Name = n;};

    /// @brief Return true if this Topology has a plug named (name)
  bool	hasPlugNamed(core::Symbol_sp name);
    /// @brief Return the Plug with the name
  Plug_sp	plugNamed(core::Symbol_sp name);
  core::List_sp allOutPlugNamesThatMatchInPlugName(core::Symbol_sp inPlugName);
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
  core::T_sp provideMissingRingClosingPlug( Monomer_sp mon );

  void mapPlugs(std::function<void(Plug_sp)> );

  core::List_sp extractFragmentsAsList();

  void throwIfExtractFragmentsAreNotExclusive(ConstitutionAtoms_sp constitutionAtoms);

  CL_DEFMETHOD void addPlug(core::Symbol_sp nm, plugType op ) { this->_Plugs.set(nm, op);};
  bool		hasFlag(core::Symbol_sp f) const;
  bool		matchesPlugs(Topology_sp cm);
  bool		matchesContext(MonomerContext_sp cm);

  size_t numberOfStereoisomers() const { return this->_StereoisomerAtomProperties.size(); }
  void addStereoisomerAtoms(StereoisomerAtoms_sp sia);
  void setStereoisomerAtoms(core::Symbol_sp coding, core::List_sp stereoisomer_atoms);
  StereoisomerAtoms_sp getStereoisomerAtoms(core::Symbol_sp stereoisomerName) const;
  core::Symbol_sp getStereoisomerName(size_t index) const;
  core::T_mv lookupStereoisomerAtoms(Fixnum index);
    
  Topology_O() : _Name(nil<core::Symbol_O>()),
                 _Constitution(unbound<chem::Constitution_O>()),
//    _AtomInfo(unbound<core::Vector_O>()),
                 _StereoisomerCoding(kw::_sym_absolute),
                 _Properties(nil<core::T_O>()),
                 _DefaultOutPlugName(nil<core::T_O>())
  {};
  Topology_O(core::Symbol_sp name, Constitution_sp constitution) : _Name(name),
                                                                   _Constitution(constitution),
//    _AtomInfo(unbound<core::Vector_O>()),
                                                                   _StereoisomerCoding(kw::_sym_absolute),
                                                                   _Properties(nil<core::T_O>()),
                                                                   _DefaultOutPlugName(nil<core::T_O>())
  {};
};

void connect_residues(Topology_sp prev_topology,
                      Residue_sp prev_residue,
                      core::Symbol_sp out_plug_name,
                      Topology_sp next_topology,
                      Residue_sp next_residue,
                      core::Symbol_sp in_plug_name);

};

#endif
