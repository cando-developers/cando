       
       
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


/*! @class A Topology describes one way that a Constitution can be connected to other Constitutions.
  It maintains a list of Plugs that describe how this Topology can plug into other Constitutions/Topologys.
  For example, a bis-amino acid can plug into two other bis-amino acids through DKPs on the leading end
  or the trailing end.  The same bis-amino acid could also have an incoming DKP and an outgoing peptide-amide
  and something on the outgoing ester.  The same bis-amino acid could also have an incoming peptide-amide
  an outgoing peptide-amide and something else connected to the trailing secondary amine and the leading 
  carbonyl group.  Each of these scenarios is described by a different Topology.
  A Topology also contains ExtractScaffold and ExtractFragment objects that know how to extract coordinates 
  from a Residue in the given Topology.

  Each Topology is uniquely identified within a CandoDatabase by a TopologyIndex0N identifier that ranges
  between 0 and N-(the number of unique Topologys in the database).  This is to facilitate atom based lookups
  that describe the number of bonds between atoms in two interconnected Topologys.
*/

class Topology_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,Topology_O,"Topology",core::CxxObject_O); 
#if INIT_TO_FACTORIES
public:
    static Topology_sp make(core::Symbol_sp name, int netCharge, core::HashTableEq_sp properties, core::List_sp curPlugs, ConstitutionAtoms_sp residue ); // , kinematics::AtomTemplate_sp atomTreeTemplate, kinematics::ChiList_sp chiList);
#else
    DECLARE_INIT();
#endif
public:
    void initialize();
public:
//	void	archiveBase(core::ArchiveP node);
private:
    typedef Plug_sp	plugType;
    typedef Plug_O	plugOType;
public:
    typedef	gctools::SmallMap<core::Symbol_sp,Plug_sp>	Plugs;
private:
    Constitution_wp				_WeakConstitution;
    core::Symbol_sp				_Name;
    TopologyIndex0N				_TopologyIndex0N;
//	kinematics::AtomTemplate_sp		_AtomTreeTemplate;
//	kinematics::ChiList_sp			_ChiList;
    int					_ResidueNetCharge;
    gctools::SmallOrderedSet<core::Symbol_sp> 	_Flags;
    bool					_SuppressTrainers;
    adapt::SymbolMap<StereoisomerAtoms_O>       	_StereoisomerAtomProperties;
    core::HashTableEq_sp			_Properties;
private:	// Do not archive
    core::T_sp				_TemporaryObject;
 public:
    Plugs					_Plugs;
	
public:
    typedef gctools::Vec0<ExtractFragment_sp>::iterator	iterateExtractFragments;
//	typedef adapt::SymbolMap<ExtractScaffold_O>::iterator	iterateExtractScaffolds;
protected:
    void setFromMonomer(Monomer_sp mon);
public:
    void setConstitution(Constitution_sp c);
public:
//	Topology_O( Constitution_sp mold);


//	/ * ! Return the AtomTreeTemplate * /
//	kinematics::AtomTemplate_sp atomTreeTemplate() const {return this->_AtomTreeTemplate;};


    /* You can attach properties to the Topology later */
    core::HashTableEq_sp properties() const;

    string description() const;



    Constitution_sp	getConstitution();

    MonomerContext_sp getMonomerContext(CandoDatabase_sp bdb);

    //! Return if we suppress trainers
CL_NAME("suppressTrainers");
CL_DEFMETHOD     bool suppressTrainers() const { return this->_SuppressTrainers;};

    //! Return all plugs as a Cons
    core::List_sp	plugsAsCons();

    //! Return a Cons of Plugs that have Mates
    core::List_sp	plugsWithMatesAsCons();

    //! Return a Cons of out plugs
    core::List_sp	outPlugsAsCons();

    /// @brief Return all of the Incomplete Frames as Cons
    /// @return Cons of Incomplete Frames
    core::List_sp	incompleteFramesAsCons();

    /// @brief Return the name of the Topology
CL_NAME("getName");
CL_DEFMETHOD     core::Symbol_sp	getName() const {return (this->_Name);};

    /// @brief Return true if this Topology has a plug named (name)
    bool	hasPlugNamed(core::Symbol_sp name);
    /// @brief Return the Plug with the name
    Plug_sp	plugNamed(core::Symbol_sp name);

    void setTemporaryObject(core::T_sp o);
    core::T_sp getTemporaryObject();

CL_NAME("numberOfPlugs");
CL_DEFMETHOD     int	numberOfPlugs() { return this->_Plugs.size(); };

    bool	hasInPlug();
    plugType getInPlug();

    /*! Return true if this Topology has all of the plugs in (plugSet) */
    bool	hasMatchingPlugsWithMates(adapt::SymbolSet_sp plugSet);

    bool	matchesMonomerEnvironment( Monomer_sp mon );
    RingClosingPlug_sp provideMissingRingClosingPlug( Monomer_sp mon );

CL_NAME("getResidueNetCharge");
CL_DEFMETHOD     int	getResidueNetCharge() { return this->_ResidueNetCharge; };
    void	setResidueNetCharge(int nc) { this->_ResidueNetCharge = nc; };


    void mapPlugs(std::function<void(Plug_sp)> );

    core::List_sp extractFragmentsAsCons();

    void throwIfExtractFragmentsAreNotExclusive(ConstitutionAtoms_sp constitutionAtoms);

    void	addPlug(core::Symbol_sp nm, plugType op ) { this->_Plugs.set(nm, op);};

    bool		hasFlag(core::Symbol_sp f) const;

    bool		matchesTopology(Topology_sp cm);

    bool		matchesContext(MonomerContext_sp cm);

    /*! Return a StereoisomerAtoms_sp object for the stereoisomer name
      Either return the existing one or create a new one */
    StereoisomerAtoms_sp lookupOrCreateStereoisomerAtoms(core::Symbol_sp stereoisomerName);

    DEFAULT_CTOR_DTOR(Topology_O);
};


};

TRANSLATE(chem::Topology_O);
#endif
