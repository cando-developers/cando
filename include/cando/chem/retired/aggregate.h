       
       
//
// (C) 2004 Christian E. Schafmeister
//





#ifndef	AGGREGATE_H
#define	AGGREGATE_H

#include <clasp/core/common.h>
#include <string>
#include <cando/chem/matter.h>
#include <cando/chem/atomIdMap.fwd.h>
#include <cando/chem/molecule.h>



#include <cando/chem/chemPackage.h>

#include "core/lisp.fwd.h"// aggregate.h wants Lisp needs lisp.fwd.h

namespace chem {


SMART(Restraint);
SMART(RestraintList);
SMART(Lisp);

class  Aggregate_O : public  Matter_O 
{
    LISP_BASE1( Matter_O );
    LISP_CLASS(chem,ChemPkg,Aggregate_O,"Aggregate");
//    DECLARE_SERIALIZE();
public:

    void initialize();
    friend	class	Loop;
public:
    void	archiveBase(core::ArchiveP node );
    bool	loadFinalize(core::ArchiveP node);
    /*! Use AtomIds to quickly access Atoms in Aggregate */
//    AtomIdMap<Atom_sp>	_AtomIdMap;

public:
    static Aggregate_sp make(MatterName name);

public:
    static Aggregate_sp	open(const string& name);
public:

    /*! Build an AtomIdMap that maps moleculeId/residueId/atomId to atoms
     Call this after you make any changes to the structure of the Aggregate.
    This is automatically called whenever an Aggregate is loaded. */
    virtual AtomIdToAtomMap_sp buildAtomIdMap() const;


    virtual Atom_sp atomWithAtomId(AtomId_sp atomId) const;

    /*! Lookup an Atom using the atomId */
//    Atom_sp lookupAtom(const AtomId& atomId) const;


    virtual string __repr__() const;

    char	getMatterType() { return AGGREGATE_CLASS; };

    contentIterator begin_molecules() { return this->_contents.begin(); };
    contentIterator end_molecules() {return this->_contents.end(); };

    Atom_sp		atomWithId(int i);

    virtual bool equal(core::T_sp obj) const;
    virtual void	transferCoordinates(Matter_sp other);

    virtual void addMatter(Matter_sp matter);
    Aggregate_sp	addMolecule( Molecule_sp a ) {this->addMatter(a); return this->sharedThis<Aggregate_O>();};
    Aggregate_sp	addMoleculeRetainId( Molecule_sp a ) {this->addMatterRetainId(Matter_sp(a)); return this->sharedThis<Aggregate_O>();};
    Molecule_sp      firstMolecule();
    MatterName     firstMoleculeName();

//    Aggregate_sp	shallowCopy();
//    Aggregate_sp	deepCopy() const;

    /*! Return the atom according to the Chimera formated specification
     *  These have the format #[model]:[residueSpecification]@[atom name]
     *  eg: #0:12.B@CA
     */
    Atom_sp		atomWithChimeraSpecification(const string& name);
    bool		hasAtomWithChimeraSpecification(const string& name);


    Atom_sp		firstAtomWithName( MatterName name );

    void		perturbAtomPositions(double dist);

    void		removeMolecule( Molecule_sp mol );
//	void		dump();

    bool		testForVdwOverlap(Aggregate_sp other);
    bool		testAggregateConsistancy() {return this->testConsistancy(Matter_sp());};

    virtual uint	numberOfAtoms( );
    virtual bool isAggregate() { return true;};

    /*! Return a Cons with all molecules */
    core::Cons_sp molecules() { return this->contentsAsCons(); };

    /*! When a structure is loaded and all the atoms are dumped into the same
     * molecule and same residue then we may want to separate out the molecules
     * defined by their connectivity.
     *
     * This function identifies the separate molecules using spanning trees
     * and creates a separate molecule for each of them.
     *
     * Return the number of molecules found.
     */
    uint separateMolecules();

    /*! Sometimes we want to identify a molecule that contains an atom with
     * a specific name.
     */
    Molecule_sp firstMoleculeWithAtomNamed(MatterName name);


    void		writeToFile(const string& fileName);
    void		saveAs(const string& fileName) { _G(); this->writeToFile(fileName);};


    bool	canRender() { return true;}
#ifdef RENDER
    candoBase::Render_sp	rendered(core::Cons_sp opt);
#endif


    virtual string	subMatter() { return Molecule_O::static_className(); };
    virtual	string	description() const { stringstream ss; ss << "aggregate("<<this->getName()<<")@"<<std::hex<<this<<std::dec; return ss.str();};


    void		duplicate(const Aggregate_O * x );


public:
    virtual Matter_sp	copy();
    
protected:
    virtual Matter_sp copyDontRedirectAtoms();
	virtual void redirectAtoms();


    DEFAULT_CTOR_DTOR( Aggregate_O);
};


};
TRANSLATE(chem:: Aggregate_O);
#endif