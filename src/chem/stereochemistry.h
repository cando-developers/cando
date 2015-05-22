       
       
//
// (C) 2004 Christian E. Schafmeister
//



/*
 *	stereochemistry.h
 *
 *
 */

#ifndef Stereochemistry_H
#define Stereochemistry_H
#include <clasp/core/common.h>

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/symbolSet.fwd.h>
#include <cando/candoBase/vector2.h>
#include <cando/candoBase/vector3.h>
#include <cando/chem/entity.h>


#include <cando/chem/chemPackage.h>


namespace chem
{


    SMART(RestraintList);
    SMART(Atom);
    SMART(Residue);
    SMART(Topology);
    SMART(ComplexRestraint);

    SMART(StereoConfiguration);
    class StereoConfiguration_O : public core::T_O
	{
	    LISP_BASE1(core::T_O);
	    LISP_CLASS(chem,ChemPkg,StereoConfiguration_O,"StereoConfiguration");
#if INIT_TO_FACTORIES
	public:
	    static StereoConfiguration_sp make(string atomName, string config);
#else
	    DECLARE_INIT();
#endif
	public:
	    static void lisp_initGlobals(core::Lisp_sp lisp);
	public:
//	    void	archive(core::ArchiveP node);

	private:
	    MatterName	_AtomName;
            core::Symbol_sp 	_Configuration;

	public:
	    static	core::Cons_sp stereochemicalPermutations(uint numberOfCenters, core::Lisp_sp e);
	    static	core::Cons_sp create_multiple(core::Cons_sp atomNames, core::Cons_sp configurations, core::Lisp_sp e);

	public:
	    MatterName	getAtomName()	{return this->_AtomName;};
	    void	setAtomName(MatterName n) {this->_AtomName = n;};
            core::Symbol_sp getConfiguration()	{return this->_Configuration;};
	    void	setConfiguration(core::Symbol_sp p) {this->_Configuration = p;};
	    int		getMoeConfiguration();

	    DEFAULT_CTOR_DTOR(StereoConfiguration_O);
	};








/*! Save a monomer name and pdb name pair
 */
    SMART(Stereoisomer);
    class Stereoisomer_O : public Entity_O
	{
	    LISP_BASE1(Entity_O);
	    LISP_CLASS(chem,ChemPkg,Stereoisomer_O,"Stereoisomer");
#if INIT_TO_FACTORIES
	public:
	    static Stereoisomer_sp make(core::Symbol_sp name, core::Symbol_sp pdb, core::Cons_sp configs);
#else
	    DECLARE_INIT();
#endif
	public:
	    void	initialize();
//	    void	archive(core::ArchiveP node);

	private:
	    Constitution_wp	_WeakConstitution;
	    core::Symbol_sp	_Name;
	    core::Symbol_sp	_Pdb;
	    core::Symbol_sp	_Enantiomer;
	    // Add configurations here
            gctools::Vec0<StereoConfiguration_sp>	_Configurations;
	public:
	    void lispFinalize();

	public:		// local implementations of CandoDatabaseDependent functions
	    bool	isTerminalName() { return true; };
	    bool	hasConstitution() { return true; };
	    Constitution_sp constitution() { return this->getConstitution(); };
	    core::SymbolSet_sp	expandedNameSet();
	    RepresentativeList_sp expandedRepresentativeList() const;

	    string getConfigurationForCenter(const string& centerName );
	public:
	    Constitution_sp	getConstitution();

	    core::Symbol_sp	getName() const {return this->_Name;};
	    void	setName(core::Symbol_sp n) {this->_Name = n;};
	    core::Symbol_sp	getPdb() {return this->_Pdb;};
	    void	setPdb(core::Symbol_sp p) {this->_Pdb = p;};
	    core::Symbol_sp	getEnantiomer()	{return this->_Enantiomer;};
	    void	setEnantiomer(core::Symbol_sp p) {this->_Enantiomer = p;};
	    void	addStereoConfiguration(StereoConfiguration_sp sc) {this->_Configurations.push_back(sc);};
	    string __repr__() const;
            gctools::Vec0<StereoConfiguration_sp>::iterator _Configurations_begin() { return this->_Configurations.begin();};
            gctools::Vec0<StereoConfiguration_sp>::iterator _Configurations_end() { return this->_Configurations.end();};



	    DEFAULT_CTOR_DTOR(Stereoisomer_O);
	};

    SMART(StereoInformation);
    class StereoInformation_O : public core::T_O
	{
	    LISP_BASE1(core::T_O);
	    LISP_CLASS(chem,ChemPkg,StereoInformation_O,"StereoInformation");
#if INIT_TO_FACTORIES
	public:
	    static StereoInformation_sp make(core::Cons_sp stereoisomers, core::Cons_sp restraints);
#else
	    DECLARE_INIT();
#endif
	public:
//	    void	archive(core::ArchiveP node);
	public:
//	    bool loadFinalize(core::ArchiveP node);
	private:
            gctools::Vec0<Stereoisomer_sp>		_Stereoisomers;
	    core::SymbolMap<Stereoisomer_O>	_NameOrPdbToStereoisomer;
            gctools::Vec0<ComplexRestraint_sp>	_ComplexRestraints;
	public:
	public:
            gctools::Vec0<Stereoisomer_sp>::iterator begin_Stereoisomers() { return this->_Stereoisomers.begin();};
            gctools::Vec0<Stereoisomer_sp>::iterator end_Stereoisomers() { return this->_Stereoisomers.end();};

            gctools::Vec0<Stereoisomer_sp>::const_iterator begin_Stereoisomers() const { return this->_Stereoisomers.begin();};
            gctools::Vec0<Stereoisomer_sp>::const_iterator end_Stereoisomers() const { return this->_Stereoisomers.end();};


#if 0
            gctools::Vec0<O_ProChiralCenter>::iterator begin_ProChiralCenters() { return this->_ProChiralCenters.begin();};
            gctools::Vec0<O_ProChiralCenter>::iterator end_ProChiralCenters() { return this->_ProChiralCenters.end();};
#endif
            gctools::Vec0<ComplexRestraint_sp>::iterator begin_ComplexRestraints() { return this->_ComplexRestraints.begin();};
            gctools::Vec0<ComplexRestraint_sp>::iterator end_ComplexRestraints() { return this->_ComplexRestraints.end();};

	    core::Cons_sp stereoisomersAsCons() { return core::Cons_O::createFromVec0(this->_Stereoisomers); };

	    void validate();

	    void addStereoisomer(Stereoisomer_sp s);
//    void addProChiralCenter(RPProChiralCenter s);
	    Stereoisomer_sp	getStereoisomer(core::Symbol_sp nameOrPdb);

	    //! Return true if this Constitution recognize the name or pdb name
	    bool	recognizesNameOrPdb(core::Symbol_sp nm)
	    {
		return ( this->_NameOrPdbToStereoisomer.contains(nm) != 0);
	    };
	    //! The name of the monomer with the nameOrPdb
	    core::Symbol_sp	nameFromNameOrPdb(core::Symbol_sp np)
	    {
		Stereoisomer_sp i;
                try {i = this->_NameOrPdbToStereoisomer[np];}
                catch (...) {CELL_ERROR(np);};
		ASSERTNOTNULL(i);
		return i->getName();
	    };
	    //! The pdb_name of the monomer with the nameOrPdb
	    core::Symbol_sp	pdbFromNameOrPdb(core::Symbol_sp np)
	    {
                Stereoisomer_sp i;
                try {i = this->_NameOrPdbToStereoisomer[np];}
                catch (...) {CELL_ERROR(np);};
		ASSERTNOTNULL(i);
		return i->getPdb();
	    };

	    //! Return the names of Monomers that can be created from this constitution
	    core::StringList_sp getMonomerNamesAsStringList();
	    core::SymbolSet_sp getMonomerNamesAsSymbolSet();
	    //! Return the pdb names of Monomers that can be created from this constitution
	    core::StringList_sp getPdbNamesAsStringList();
	    DEFAULT_CTOR_DTOR(StereoInformation_O);
	};


};

TRANSLATE(chem::StereoConfiguration_O);
TRANSLATE(chem::Stereoisomer_O);
TRANSLATE(chem::StereoInformation_O);
#endif
