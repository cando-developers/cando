/*
    File: atomIdMap.h
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
#ifndef	_ATOMIDMAP_H
#define _ATOMIDMAP_H

#include <clasp/core/common.h>
#include <cando/chem/atomId.h>

namespace chem
{



class AtomIdMapFunctor
{
public:
  virtual void operator()(const AtomId& atomId, Atom_sp atom) = 0;
};



template <class T>
class AtomIdMap
{
public:
  typedef	vector<T>		AtomMap;
  typedef	vector<AtomMap>		ResidueMap;
  typedef	vector<ResidueMap>	MoleculeMap;

public:	// Functions here

	/*! Number of molecules */
  int numberOfMolecules() const { return this->_MoleculeMap.size(); };

	/*! Number of residues */
  int numberOfResidues(int mid) const { return this->_MoleculeMap[mid].size(); };

	/*! Number of atoms */
  int numberOfAtoms(int mid,int rid) const { return this->_MoleculeMap[mid][rid].size(); };

	/*! Iterate over every atom in the AtomIdMap and call the Functor with the AtomId and the atom */
  void iterate(AtomIdMapFunctor& functor) const
  {
    for ( int mid=0; mid<(int)this->_MoleculeMap.size(); mid++ )
    {
      for ( int rid=0; rid<(int)this->_MoleculeMap[mid].size(); rid++ )
      {
        for ( int aid=0; aid<(int)this->_MoleculeMap[mid][rid].size(); aid++ )
        {
          AtomId atomId(mid,rid,aid);
          functor(atomId,this->_MoleculeMap[mid][rid][aid]);
        }
      }
    }
  }


	/*! Resize the number of molecules */
  void resize(uint nmols)
  {
    HARD_ASSERTF(nmols>0,BF("Illegal number of molecules[%d]") % nmols );
    this->_MoleculeMap.resize(nmols);
  };
	/*! Resize the number of residues in a molecule */
  void resize(uint mol, uint nres)
  {
    HARD_ASSERTF(mol<this->_MoleculeMap.size(),BF("Illegal molecule[%d]") % mol );
    HARD_ASSERTF(nres>=0,BF("Illegal number of residues[%d]") % nres );
    this->_MoleculeMap[mol].resize(nres);
  }

	/*! Resize the number of atoms in a residue of a molecule */
  void resize(uint mol, uint res, uint natoms)
  {
    HARD_ASSERTF(mol<this->_MoleculeMap.size(),BF("Illegal molecule[%d]") % mol );
    HARD_ASSERTF(res<this->_MoleculeMap[mol].size(),BF("Illegal residue[%d]") % res );
    HARD_ASSERTF(natoms>=0,BF("Illegal number of atoms[%d]") % natoms );
    this->_MoleculeMap[mol][res].resize(natoms);
  }

  void throwIfInvalidMoleculeId(const AtomId& atomId) const
  {_G();
    if ( atomId.moleculeId() < (int)this->_MoleculeMap.size() ) return;
    SIMPLE_ERROR(BF("moleculeId[%d] is out of range - it must be less than %d") % atomId.moleculeId() % this->_MoleculeMap.size() );
  }

  void throwIfInvalidResidueId(const AtomId& atomId) const
  {_G();
    ASSERT(atomId.moleculeId() < (int)this->_MoleculeMap.size() );
    int numResidues = this->_MoleculeMap[atomId.moleculeId()].size();
    if ( atomId.residueId() < numResidues ) return;
    SIMPLE_ERROR(BF("residueId[%d] is out of range - it must be less than %d") % atomId.residueId() % numResidues );
  }

  void throwIfInvalidAtomId(const AtomId& atomId) const
  {_G();
    ASSERT(atomId.moleculeId() < (int)this->_MoleculeMap.size() );
    ASSERT(atomId.residueId() < (int)this->_MoleculeMap[atomId.moleculeId()].size());
    int numAtoms = this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()].size();
    if ( atomId.atomId() < numAtoms ) return;
    SIMPLE_ERROR(BF("atomId[%d] is out of range - it must be less than %d") % atomId.atomId() % numAtoms );
  }



public: // Access
  T const& operator[](const AtomId& atomId) const
  {
    HARD_ASSERTF((atomId.moleculeId()!=AtomId::UndefId) && (atomId.residueId()!=AtomId::UndefId) && (atomId.atomId()!=AtomId::UndefId),
                 BF("AtomId is undefined"));
    HARD_ASSERTF(atomId.moleculeId()<(int)this->_MoleculeMap.size(),
                 BF("Molecule index[%d] is out of range - number of molecules[%d]")
                 % atomId.moleculeId() % this->_MoleculeMap.size() );
    HARD_ASSERTF(atomId.residueId()<(int)this->_MoleculeMap[atomId.moleculeId()].size(),
                 BF("Residue index[%d] is out of range - number of residues[%d] in molecule[%d]")
                 % atomId.residueId()
                 % this->_MoleculeMap[atomId.moleculeId()].size()
                 % atomId.moleculeId() );
    HARD_ASSERTF(atomId.atomId()<(int)this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()].size(),
                 BF("Atom index[%d] is out of range - number of atoms[%d] in residue[%d] in molecule[%d]")
                 % atomId.atomId()
                 % this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()].size()
                 % atomId.residueId()
                 % atomId.moleculeId() );
    return this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()][atomId.atomId()];
  }
    
  T& operator[](const AtomId& atomId)
  {
    HARD_ASSERTF((atomId.moleculeId()!=AtomId::UndefId) && (atomId.residueId()!=AtomId::UndefId) && (atomId.atomId()!=AtomId::UndefId), BF("AtomId is undefined"));
    HARD_ASSERTF(atomId.moleculeId()<(int)this->_MoleculeMap.size(),
                 BF("Molecule index[%d] is out of range - number of molecules[%d]")
                 % atomId.moleculeId() % this->_MoleculeMap.size() );
    HARD_ASSERTF(atomId.residueId()<(int)this->_MoleculeMap[atomId.moleculeId()].size(),
                 BF("Residue index[%d] is out of range - number of residues[%d] in molecule[%d]")
                 % atomId.residueId()
                 % this->_MoleculeMap[atomId.moleculeId()].size()
                 % atomId.moleculeId() );
    HARD_ASSERTF(atomId.atomId()<(int)this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()].size(),
                 BF("Atom index[%d] is out of range - number of atoms[%d] in residue[%d] in molecule[%d]")
                 % atomId.atomId()
                 % this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()].size()
                 % atomId.residueId()
                 % atomId.moleculeId() );
    return this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()][atomId.atomId()];
  }

  ResidueMap const& operator[](int moleculeId) const
  {
    HARD_ASSERTF(moleculeId!=AtomId::UndefId,BF("Illegal moleculeId[%d]")% moleculeId );
    HARD_ASSERTF(moleculeId<(int)this->_MoleculeMap.size(),
                 BF("Molecule index[%d] is out of range - number of molecules[%d]")
                 % moleculeId % this->_MoleculeMap.size() );
    return this->_MoleculeMap[moleculeId];
  }

  ResidueMap& operator[](int moleculeId)
  {
    HARD_ASSERTF(moleculeId!=AtomId::UndefId,BF("Illegal moleculeId[%d]")% moleculeId );
    HARD_ASSERTF(moleculeId<(int)this->_MoleculeMap.size(),
                 BF("Molecule index[%d] is out of range - number of molecules[%d]")
                 % moleculeId % this->_MoleculeMap.size() );
    return this->_MoleculeMap[moleculeId];
  }

  T safeLookup(const AtomId& atomId) const
  {
    if (atomId.moleculeId()<this->_MoleculeMap.size()) {
      if (atomId.residueId()<this->_MoleculeMap[atomId.moleculeId()].size()) {
        if (atomId.atomId()<this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()].size()) {
          return this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()][atomId.atomId()];
        }
        SIMPLE_ERROR(BF("atom index %lu is out of range %lu") % atomId.atomId() % this->_MoleculeMap[atomId.moleculeId()][atomId.residueId()].size() );
      }
      SIMPLE_ERROR(BF("residue index %lu is out of range %lu") % atomId.residueId() % this->_MoleculeMap[atomId.moleculeId()].size() );
    }
    SIMPLE_ERROR(BF("molecule index %lu is out of range %lu") % atomId.moleculeId() % this->_MoleculeMap.size());
  }


public:
	// instance variables here
  T		_DefaultValue;
  MoleculeMap	_MoleculeMap;
};




FORWARD(AtomIdToAtomMap);
class AtomIdToAtomMap_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,AtomIdToAtomMap_O,"AtomIdToAtomMap",core::CxxObject_O);
private:
  AtomIdMap<Atom_sp>	_AtomIdMap;
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
  void initialize();
//	void archiveBase(core::ArchiveP node);

public:
  void resize(int numMolecules);
  void resize(int mol, int numRes);
  void resize(int mol, int res, int numAtoms);
  void set(AtomId const& atomId, Atom_sp atom);
  Atom_sp lookupAtom(const AtomId& atomId) const ;
  void walk(core::Function_sp callback);
  DEFAULT_CTOR_DTOR(AtomIdToAtomMap_O);
};



}; /* chem */

TRANSLATE(chem::AtomIdToAtomMap_O);

#endif /* _ATOMIDMAP_H */
