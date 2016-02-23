/*
    File: constitution.h
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


#ifndef Constitution_H  //[
#define Constitution_H


/*
 *	constitution.h
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

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/adapt/stringList.fwd.h>
#include <cando/geom/vector2.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/quickDom.fwd.h>
#include <cando/chem/entity.h>
#include <cando/chem/constitution.fwd.h>
#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/chem/chemPackage.h>
#include <clasp/core/cons.fwd.h>// constitution.h wants Cons needs cons.fwd.h
#include <cando/chem/residue.fwd.h>



namespace chem
{

    SMART(RepresentedEntityNameSet);
    SMART(MonomerContext);
    SMART(CandoDatabase);
    SMART(Coupling);
    SMART(Oligomer);
    SMART(Constitution);
    SMART(StereoInformation);
    SMART(Monomer);
    SMART(Topology);
    SMART(Plug);
    SMART(RingClosingPlug);
    SMART(Stereoisomer);
    SMART(StereoInformation);
    SMART(Constitution);


    /*! @class A Constitution is like a Residue but its more general because it describes multiple
      stereoisomers. A Constitution stores a collection 
      of mutually bonded atoms (like a Residue) that connect to other collections of mutually bonded atoms.
      A Constitution also stores a collection of stereoisomers that have the same constitution.
      It also stores a collection of Topologys that describe different ways the Constitution can be connected
      to other Constitutions.  It also stores a collection of Plugs that describe how select atoms within
      the Constitution can be plugged into other Constitutions. Another thing that it keeps track of are a collection of 
      Frames that define coordinate systems within which internal Fragment atoms are built and that are exported
      to define coordinate systems within which atoms of other attaching Constitutions are built.  
      Finally, it stores a collection of Fragments that encompass all of the atoms in the Constitution and
      are built within coordinate systems defined by the Frames of each Topology.
      
      Each Constitution is assigned a unique ConstitutionIndex0N value when the CandoDatabase is created that uniquely 
      identifies it within all Constitutions of the CandoDatabase
    */

    class Constitution_O : public Entity_O
	{
	    LISP_CLASS(chem,ChemPkg,Constitution_O,"Constitution",Entity_O);
#if INIT_TO_FACTORIES
	public:
	    static Constitution_sp make(core::Symbol_sp name, const string& comment, core::Symbol_sp metaConstitution, ConstitutionAtoms_sp constitutionAtoms, StereoInformation_sp stereoInformation, core::List_sp plugs, core::List_sp topologies);
#else
	    DECLARE_INIT();
#endif
	public:
	    void initialize();
	public:
//	    void	archiveBase(core::ArchiveP node);

	    friend class Monomer_O;
	public:
	    typedef	gctools::SmallMap<core::Symbol_sp,Topology_sp> TopologyMap; //adapt::SymbolMap<Topology_O>	TopologyMap;
	    typedef    	gctools::SmallMap<core::Symbol_sp,Plug_sp> PlugMap; // adapt::SymbolMap<Plug_O>		PlugMap;
	private:
	    core::Symbol_sp			_Name;
	    ConstitutionIndex0N			_ConstitutionIndex0N;
	    core::Symbol_sp			_MetaConstitution;
	    string				_Comment;
	    /*! An object that that describes the atoms and bonding within this Constitution.
	      The index of the atom in the ConstitutionAtoms is used as the atom part of the AtomId 
	      when refered to using an AtomId object    */
	    ConstitutionAtoms_sp		_ConstitutionAtoms;
	    PlugMap				_PlugsByName;
	    TopologyMap				_Topologies;
	    StereoInformation_sp		_StereoInformation;
	private:
//	    Residue_sp	loadResidue(adapt::QDomNode_sp node);

	public:
	    typedef gctools::Vec0<Stereoisomer_sp>::iterator	stereoisomerIterator;
	    typedef gctools::Vec0<Stereoisomer_sp>::const_iterator	const_stereoisomerIterator;
	public:
	    virtual adapt::SymbolSet_sp	expandedNameSet() { return this->getMonomerNamesAsSymbolSet(); };
	public:



	    core::List_sp stereoisomersAsCons();
	    core::List_sp topologiesAsCons();
	    core::List_sp plugsAsCons();
	    core::List_sp plugsWithMatesAsCons();

	    stereoisomerIterator begin_Stereoisomers();
	    stereoisomerIterator end_Stereoisomers();

	    const_stereoisomerIterator begin_Stereoisomers() const;
	    const_stereoisomerIterator end_Stereoisomers() const;

	    TopologyMap::iterator	begin_Topologies() {
		return this->_Topologies.begin();};
	    TopologyMap::iterator	end_Topologies() {
		return this->_Topologies.end();};

	    PlugMap::iterator begin_Plugs() {
		return this->_PlugsByName.begin();};
	    PlugMap::iterator end_Plugs() {
		return this->_PlugsByName.end();};


	    bool	isSameAs(Constitution_sp& m) { return (this->getName()==m->getName());};

	    void	makeResidueConsistentWithStereoisomerNamed(Residue_sp res,
								   core::Symbol_sp stereoisomerName );

	    Topology_sp	getTopologyForContext(MonomerContext_sp c );

	    Topology_sp	getTopologyForMonomerEnvironment(Monomer_sp mon );

	    RingClosingPlug_sp getMissingRingClosingPlug(Monomer_sp mon, Monomer_sp mate);

	    RepresentativeList_sp expandedRepresentativeList() const;

	    bool hasStereoisomerWithName(core::Symbol_sp stereoisomerName);
	    Stereoisomer_sp getStereoisomerWithName(core::Symbol_sp stereoisomerName) const;

	    core::Symbol_sp	getName() const {return this->constitutionName();};

CL_NAME("constitutionName");
CL_DEFMETHOD 	    core::Symbol_sp	constitutionName() const	{return this->_Name;};
CL_NAME("setConstitutionName");
CL_DEFMETHOD 	    void	setConstitutionName(core::Symbol_sp tt) {this->_Name = tt;};

CL_NAME("getMetaConstitutionName");
CL_DEFMETHOD 	    core::Symbol_sp	getMetaConstitutionName()	{return this->_MetaConstitution;};
CL_NAME("setMetaConstitutionName");
CL_DEFMETHOD 	    void	setMetaConstitutionName(core::Symbol_sp tt) {this->_MetaConstitution = tt;};


	    /*! Add all of my Stereoisomers as Entities to the CandoDatabase */
	    void addStereoisomersToCandoDatabase(CandoDatabase_sp db);


	    RepresentedEntityNameSet_sp	asGroup();

CL_NAME("getConstitutionAtoms");
CL_DEFMETHOD 	    ConstitutionAtoms_sp getConstitutionAtoms() { return this->_ConstitutionAtoms;};

	    Topology_sp	simplestTopologyWithPlugNamed(core::Symbol_sp name);

	    Topology_sp topologyWithName(core::Symbol_sp name) const;


	    adapt::StringList_sp getMonomerNamesAsStringList();
	    adapt::SymbolSet_sp	getMonomerNamesAsSymbolSet();

	    adapt::SymbolSet_sp	getPlugNames();

	    adapt::StringList_sp getPdbNamesAsStringList();

	    //! Return the monomer name for the monomer with the nameOrPdb
	    core::Symbol_sp	nameFromNameOrPdb(core::Symbol_sp nm);
	    //! Return the pdb name for the monomer with the nameOrPdb
	    core::Symbol_sp	pdbFromNameOrPdb(core::Symbol_sp nm);


CL_NAME("setComment");
CL_DEFMETHOD 	    void	setComment( const string& sName ) { this->_Comment= sName; };
CL_NAME("getComment");
CL_DEFMETHOD 	    string	getComment( ) { return this->_Comment; };

	    Residue_sp	createResidueForStereoisomerName(core::Symbol_sp nameOrPdb);

CL_NAME("hasPlugNamed");
CL_DEFMETHOD 	    bool	hasPlugNamed( core::Symbol_sp name) { return this->_PlugsByName.contains(name); };
CL_NAME("getPlugNamed");
CL_DEFMETHOD 	    Plug_sp	getPlugNamed( core::Symbol_sp name) { return this->_PlugsByName[name]; };

//	void	testConsistency( std::ostream& sout );


//	set<string>	getValidInCouplingRuleNames();

	    bool	recognizesMonomerAndEnvironment(Monomer_sp mon);

	    bool	recognizesStereoisomerName(core::Symbol_sp nm );

	    string	description() const;

	    bool hasConstitution() { return true;};
	    Constitution_sp constitution() { return this->sharedThis<Constitution_O>();};
	    DEFAULT_CTOR_DTOR(Constitution_O);
	};


};


TRANSLATE(chem::Constitution_O);
#endif
