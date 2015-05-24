       
       
//
// (C) 2004 Christian E. Schafmeister
//





#ifndef	MOLECULE_H
#define MOLECULE_H

#include <iostream>
#include <sstream>
#include <clasp/core/common.h>
#include <cando/chem/matter.h>
#include <cando/chem/atomIdMap.fwd.h>
#include <cando/chem/residue.h>

#include <cando/chem/chemPackage.h>


namespace chem {


class	Loop;

SMART(Molecule);
class Molecule_O : public Matter_O
{
    friend class Aggregate_O;
    LISP_BASE1(Matter_O);
    LISP_CLASS(chem,ChemPkg,Molecule_O,"Molecule");
//    DECLARE_SERIALIZE();
public:
	void initialize();
	void archiveBase(core::ArchiveP node);
friend	class	Loop;


public:
    static Molecule_sp make(MatterName name);

private:
	void	duplicate(const Molecule_O* m );
public:

	char	getMatterType() { return MOLECULE_CLASS; };
virtual	char	getClass() { return moleculeId;};

	Molecule_O(const Molecule_O& mol );




    virtual string __repr__() const;

	contentIterator begin_residues() { return this->_contents.begin(); };
	contentIterator end_residues() {return this->_contents.end(); };
	const_contentIterator begin_residues() const { return this->_contents.begin(); };
	const_contentIterator end_residues() const {return this->_contents.end(); };

	virtual bool equal(core::T_sp obj) const;
virtual void	transferCoordinates(Matter_sp other);

 Residue_sp	firstResidue() {return this->contentAt(0).as<Residue_O>();};
 Residue_sp	getResidue(int i) {return this->contentAt(i).as<Residue_O>();};
 Residue_sp	residueWithId( int lid ) { return this->contentWithId(lid).as<Residue_O>(); };
	bool		hasResidueWithId( int lid ) { return this->hasContentWithId(lid); };

	int		numberOfResiduesWithName(MatterName name);
	Residue_sp	getFirstResidueWithName(MatterName name);
    VectorResidue	getResiduesWithName(MatterName name);



	virtual bool isMolecule() { return true;};

	int		residueCount() {return (this->_contents.size());};

	void		moveAllAtomsIntoFirstResidue();

virtual string	subMatter() { return Residue_O::static_className(); };
 virtual	string	description() const { stringstream ss; ss << "molecule("<<_rep_(this->getName())<<")@"<<std::hex<<this<<std::dec; return ss.str();};

	bool		testMoleculeConsistancy() {IMPLEMENT_ME();};

	bool	canRender() { return true;}
#ifdef RENDER
	geom::Render_sp	rendered(core::Cons_sp options);
#endif

	virtual uint	numberOfAtoms( );

	/*! Build a map of AtomIds to Atoms */
	virtual AtomIdToAtomMap_sp buildAtomIdMap() const;

    virtual Atom_sp atomWithAtomId(AtomId_sp atomId) const;


public:
    virtual Matter_sp	copy();
    
protected:
    virtual Matter_sp copyDontRedirectAtoms();
	virtual void redirectAtoms();



	DEFAULT_CTOR_DTOR(Molecule_O);
};


//extern Molecule_sp	new_Molecule_sp();
//extern Molecule_sp	_copy_Molecule_sp(Molecule_O*);


};
TRANSLATE(chem::Molecule_O);
#endif

