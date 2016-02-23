/*
    File: mol2.cc
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
       
#define DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//





/*
 *	mol2.cc
 *
 *	Functions that read and write TRIPOS Mol2 files.
 *
 */

#include <cando/chem/mol2.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <stdio.h>
#include <map>
#include <clasp/core/foundation.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/pathname.h>
#include <cando/adapt/indexedObjectBag.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/loop.h>
#include <cando/chem/molecule.h>
#include <clasp/core/str.h>
#include <cando/chem/chemPackage.fwd.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/chemInfo.h>
#include <cando/chem/typeAssignmentRules.h>

namespace chem {


static chem::TypeAssignmentRules_sp	sybylRules;


class	Mol2File 
{
private:
  core::T_sp	fIn;
  bool _eof;
  std::stringstream	mLine;
public:
  bool eof() { return this->_eof;};
  void	advanceLine() {
    bool ignoreNextCr;
    char c;
    this->mLine.str("");
    ignoreNextCr = false;
    while (1) {
      c = core::clasp_read_char(this->fIn);
      if ( c == EOF ) {
        this->_eof = true;
        break;
      }
      if ( c == 10 ) {
        if ( ignoreNextCr ) {
          ignoreNextCr = false;
        } else break;
      }
      if ( c == '\\' ) {
        ignoreNextCr = true;
      } else {
        if ( c >= ' ' ) {
          this->mLine.put(c);
        }
      }
    }
  };
  stringstream&	line() { return this->mLine; };
  std::queue<string>	splitLine() {
    std::queue<string>	qWords;
    core::queueSplitString(this->mLine.str(), qWords, " \t" );
    return qWords;
  };
  bool	hasDataLine() {
    if ( this->eof() ) return false;
    if ( this->mLine.str().size() == 0 ) return false;
    if ( this->mLine.str()[0] != '@' ) return true;
    return false;
  };
  void	openFileName(core::T_sp pn)
  {_OF();
    this->fIn = core::clasp_openRead(pn);
    if ( this->fIn.nilp() ) {
      SIMPLE_ERROR(BF("Mol2 file doesn't exist: %s" ) % _rep_(pn) );
    }
    this->advanceLine();
  }
  Mol2File() : _eof(false) {};
  ~Mol2File() {
    cl__close(fIn);
  }


};

struct TriposMolecule {
	string	mMoleculeName;
	uint	mAtoms;
	uint	mBonds;
	uint	mResidues;
	uint	m_Unknown1;
	uint	m_Unknown2;
	string	m_Comment1_Unknown;
	string	m_Comment2_Unknown;
};

struct TriposAtom {
	uint	mIndex;
	string	atom_name;
	double	mX;
	double	mY;
	double	mZ;
	string	mType;
	uint	mSubstId;
	string	mSubstName;
	double	mCharge;
	Atom_sp	mAtom;
};

struct TriposBond {
	uint	mId;
	uint	mAtom1Id;
	uint	mAtom2Id;
	string	mOrder;
};


struct TriposSubstructure {
	uint	mId;
	string	subst_name;
	uint	root_atom;
	string	subst_type;
	uint	dict_type;
	string	chain;
	string	sub_type;
	uint	mInterBonds;
	string	mStatus;
	string	mComment;
	uint	file_sequence_number;
};



void	splitStringIntoNameAndNumber(const string& nameNum, string& name, uint& seqNum )
{
size_t pos = nameNum.find_last_not_of("0123456789");
size_t seqPos, nameLen;
    if ( pos == string::npos )
    {
        seqNum = -1;
	name = nameNum;
    }
    else
    {
    	seqPos = pos + 1;
        if ( nameNum[pos] == '_' )
	{
	    nameLen = pos;
	} else
	{
	    nameLen = pos + 1;
	}
	name = nameNum.substr(0,nameLen);
	seqNum = atoi(nameNum.substr(seqPos,9999).c_str());
    }
}





void	_calculateElementAndHybridization(Atom_sp a)
{
    MatterName	mnm;
    Element el;
    Hybridization geom;
    uint	neighbors;
    mnm = a->getName();
    string nm = mnm->symbolName()->get();
    
    // Set the element from the name HACK HACK HACK
    if ( nm.substr(0,2) == "CU" || nm.substr(0,2) == "Cu" ) {
        el = element_Cu;
    } else if ( nm.substr(0,2) == "CL" || nm.substr(0,2) == "Cl" ) {
        el = element_Cl;
    } else if ( nm.substr(0,2) == "BR" || nm.substr(0,2) == "Br" ) {
        el = element_Br;
    } else if ( nm.substr(0,2) == "RU" || nm.substr(0,2) == "Ru" ) {
        el = element_Ru;
    } else if ( nm.substr(0,1) == "C" ) {
        el = element_C;
    } else if ( nm.substr(0,1) == "N" ) {
        el = element_N;
    } else if ( nm.substr(0,1) == "O" ) {
        el = element_O;
    } else if ( nm.substr(0,1) == "H" ) {
        el = element_H;
    } else if ( nm.substr(0,1) == "S" ) {
        el = element_S;
    } else if ( nm.substr(0,1) == "P" ) {
        el = element_P;
    } else if ( nm.substr(0,1) == "B" ) {
        el = element_B;
    } else if ( nm.substr(0,1) == "F" ) {
        el = element_F;
    } else if ( nm.substr(0,1) == "I" ) {
        el = element_I;
    } else {
        el = element_Undefined;
    }
    geom = hybridization_sp3;
    neighbors = a->numberOfBonds();
    if ( el == element_C ) {
	switch (neighbors) {
        case 2:
            geom = hybridization_sp;
            break;
        case 3:
            geom = hybridization_sp2;
            break;
	}
    } else if ( el == element_N ) {
	if ( neighbors == 2 ) {
	    geom = hybridization_sp2;
	}
    } else if ( el == element_O || el == element_S ) {
	if ( neighbors == 1 ) {
	    geom = hybridization_sp2;
	}
    }
    a->setElement(el);
    a->setHybridization(geom);
}


void	mol2ReadAggregateFromFileName( Aggregate_sp aggregate, core::T_sp sFileName )
{
    string				line;
    TriposMolecule			mol;
    vector<TriposMolecule>		molecules;
    map<uint,TriposAtom>		atoms;
    vector<TriposBond>		bonds;
    vector<TriposSubstructure>	substructures;
    gctools::SmallMap<MatterName,Molecule_sp>		chains;
    gctools::SmallMap<uint,Residue_sp>		residues;
    std::stringstream			sin;
    std::queue<string>			words;
    string				sAtoms, oneWord;
    TriposAtom			oneAtom;
    TriposBond			oneBond;
    TriposSubstructure		oneSubst;
    Molecule_sp 			m;
    Residue_sp 				res;
    Atom_sp 				a;
    uint				lastmId;
    bool				sawSubStructure = false;
    Mol2File fIn;
    fIn.openFileName(sFileName);
    while ( !fIn.eof() ) {
	line = fIn.line().str();
	if ( fIn.hasDataLine() ) {
	    SIMPLE_ERROR(BF("Mol2 file encountered data line when expecting command line in top level: %s") % fIn.line().str() );
	}
	fIn.advanceLine();
	if ( line.size() == 0 ) continue;
//	LOG(BF("line: %s") % line.c_str()  );
	if ( line[0] == '#' ) continue;
	if ( line == "@<TRIPOS>MOLECULE" ) {
	    LOG(BF("Got @TRIPOS.MOLECULE") );
	    mol.mMoleculeName = fIn.line().str();
	    fIn.advanceLine(); // num_atoms
	    words = fIn.splitLine();
	    if ( words.empty() ) continue;
	    mol.mAtoms = atoi(words.front().c_str()); words.pop();
	    LOG(BF("Got num_atoms=%d") % mol.mAtoms  );
	    if ( words.empty() ) continue;
	    mol.mBonds = atoi(words.front().c_str()); words.pop();
	    LOG(BF("Got num_bonds=%d") % mol.mBonds  );
	    if ( words.empty() ) continue;
	    mol.mBonds = atoi(words.front().c_str()); words.pop();
	    fIn.advanceLine(); // mol_type
	    if ( !fIn.hasDataLine() ) {
	        SIMPLE_ERROR(BF("Mol2 file encountered cmd line when expecting data line: %s") % fIn.line().str() );
	    }
	    fIn.advanceLine(); // charge_type
	    while ( fIn.hasDataLine() ) {
		fIn.advanceLine();
	    }
	    molecules.push_back(mol);
	} else if ( line == "@<TRIPOS>ATOM" ) {
	    _BLOCK_TRACE("Reading @<TRIPOS>ATOM");
	    while ( fIn.hasDataLine() ) {
		line = fIn.line().str();
		if ( line[0] == '#' ) {
		    fIn.advanceLine();
		    continue;
		}
		words = fIn.splitLine();
		oneAtom.mIndex = atoi(words.front().c_str()); words.pop();
		oneAtom.atom_name = words.front(); words.pop();
		oneAtom.mX = atof(words.front().c_str()); words.pop();
		oneAtom.mY = atof(words.front().c_str()); words.pop();
		oneAtom.mZ = atof(words.front().c_str()); words.pop();
		oneAtom.mType = words.front(); words.pop();
		oneAtom.mSubstId = 0;
		oneAtom.mCharge = 0.0;
		oneAtom.mSubstName = "";
		if ( !words.empty() ) {
		    oneAtom.mSubstId = atoi(words.front().c_str());
		    LOG(BF("Read mSubstId=%d") % oneAtom.mSubstId  );
		    words.pop();
		    if ( !words.empty() ) {
			oneAtom.mSubstName = words.front();
		        LOG(BF("Read mSubstName=%s") % oneAtom.mSubstName.c_str() );
			words.pop();
			if ( !words.empty() ) {
			    oneWord = words.front();
			    LOG(BF("Charge word = |%s|") % oneWord.c_str()  );
			    oneAtom.mCharge = atof(oneWord.c_str());
		            LOG(BF("Read charge =%lf") % oneAtom.mCharge  );
			}
		    }
		}
		atoms[oneAtom.mIndex] = oneAtom;
		LOG(BF("Added atom with id: %d name: %s") % oneAtom.mIndex % oneAtom.atom_name.c_str()  );
		fIn.advanceLine();
	    }
	} else if ( line == "@<TRIPOS>BOND" ) {
	    _BLOCK_TRACE("Reading @<TRIPOS>BOND");
	    while ( fIn.hasDataLine() ) {
		line = fIn.line().str();
		LOG(BF("line len(%d) line: %s") % line.size() % line.c_str()  );
		if ( line[0] == '#' ) {
		    fIn.advanceLine();
		    continue;
		}
		words = fIn.splitLine();
		oneBond.mId = atoi(words.front().c_str()); words.pop();
		oneBond.mAtom1Id = atoi(words.front().c_str()); words.pop();
		oneBond.mAtom2Id = atoi(words.front().c_str()); words.pop();
		oneBond.mOrder= words.front(); words.pop();
		bonds.push_back(oneBond);
		fIn.advanceLine();
	    }
	} else if ( line == "@<TRIPOS>SUBSTRUCTURE" ) {
	    _BLOCK_TRACE("Reading @<TRIPOS>SUBSTRUCTURE");
	    sawSubStructure = true;
	    lastmId = -1;
	    string	nameNum;
	    while ( fIn.hasDataLine() ) {
		line = fIn.line().str();
		LOG(BF("line = (%s)") % line.c_str()  );
		if ( line[0] == '#' ) {
		    fIn.advanceLine();
		    continue;
		}
		words = fIn.splitLine();
		oneSubst.mId = atoi(words.front().c_str()); words.pop();
		if ( oneSubst.mId == lastmId ) {
		    SIMPLE_ERROR(BF("First entry of SUBSTRUCTURE records are not valid"));
		}
		lastmId = oneSubst.mId;
		nameNum = words.front(); words.pop();
		splitStringIntoNameAndNumber( nameNum, oneSubst.subst_name,oneSubst.file_sequence_number);
//		oneSubst.subst_name = words.front(); words.pop();
		LOG(BF("oneSubst.root_atom str(%s)") % words.front().c_str()  );
		oneSubst.root_atom = atoi(words.front().c_str()); words.pop();
		oneSubst.dict_type = 0;
		oneSubst.chain = "";
		oneSubst.sub_type = oneSubst.subst_name;
		oneSubst.mInterBonds = 0;
		oneSubst.mStatus = "";
		oneSubst.mComment = "";
		if ( !words.empty() ) {
		    oneSubst.subst_type = words.front(); words.pop();
		    if ( !words.empty() ) {
			oneSubst.dict_type = atoi(words.front().c_str());
			words.pop();
			if ( !words.empty() ) {
			    oneSubst.chain = words.front(); words.pop();
			    if ( !words.empty() ) {
				oneSubst.sub_type = words.front();
				LOG(BF("Got sub_type = %s") % oneSubst.sub_type.c_str()  );
				words.pop();
				if ( !words.empty() ) {
				    oneSubst.mInterBonds = atoi(words.front().c_str());
				    words.pop();
				    if ( !words.empty() ) {
					oneSubst.mStatus = words.front();
					words.pop();
					while ( !words.empty() ) {
					    oneSubst.mComment = oneSubst.mComment +
						words.front() + " ";
					    words.pop();
					}
				    }
				}
			    }
			}
		    }
		}
		fIn.advanceLine();
		LOG(BF("Read substituent  mId=%d subst_type(%s)") % oneSubst.mId % oneSubst.subst_type.c_str() );
		substructures.push_back(oneSubst);
	    }
	} else
	{
	    _lisp->print(BF("Ignoring mol2 command: %s") % line.c_str() );
	    while ( fIn.hasDataLine() ) {
		fIn.advanceLine();
	    }
	}
    }
    //
    // First create all of the Aggregate
    //	(Tripos calls this a molecule and it can have multiple chains.)
    //
    vector<TriposMolecule>::iterator	mi;
    if ( molecules.size() != 1 ) {
	SIMPLE_ERROR(BF("There can be only one @<TRIPOS>MOLECULE entry"));
    }
    mi=molecules.begin();
    aggregate->setName(chemkw_intern(mi->mMoleculeName));

    // Now create all of the atoms and the residues and put the
    // atoms in the residues
    map<uint,TriposAtom>::iterator	ai;
    Vector3				pos;
    string				nm, el;
    for ( ai=atoms.begin(); ai!= atoms.end(); ai++ ) {
	GC_ALLOCATE(Atom_O, tempa ); // new_Atom_sp();
	a = tempa;
	a->setName(chemkw_intern(ai->second.atom_name));

	LOG(BF("Creating atom with id: %d name: %s element: %s charge: %8.2f") % ai->second.mIndex % ai->second.atom_name.c_str() % el.c_str() % ai->second.mCharge  );
	pos.set(ai->second.mX,ai->second.mY,ai->second.mZ);
	a->setPosition(pos);
	a->setType(chemkw_intern(ai->second.mType));
	a->setCharge(ai->second.mCharge);
	LOG(BF(" atom info: %s") % a->description().c_str()  );
	ai->second.mAtom = a;
	//
	// If we haven't seen this residue before
	// we need to create it
	//
	Residue_sp res;
	if ( residues.count(ai->second.mSubstId)==0 ) {
	    GC_ALLOCATE(Residue_O, tempRes ); // new_Residue_sp();
	    res = tempRes;
	    // April 06/2007 decided that its best to set the 
	    // residue name from the substructure records below
	    //
	    // But I'll set them here incase there isn't 
	    // a substructure entry
	    res->setName(chemkw_intern(ai->second.mSubstName));
	    res->setPdbName(chemkw_intern(ai->second.mSubstName));
//            printf("%s:%d setFileSequenceNumber --> %d\n", __FILE__, __LINE__, ai->second.mSubstId);
            res->setFileSequenceNumber(ai->second.mSubstId);
	    LOG(BF("Creating residue with id(%d) and name(%s)") % ai->second.mSubstId % ai->second.mSubstName.c_str()  );
	    residues[ai->second.mSubstId] = res;
	    LOG(BF("New residue %s") % res->description().c_str()  );
	} else {
	    res = residues[ai->second.mSubstId];
	    LOG(BF("Reusing residue with mSubstId=%d %s") % ai->second.mSubstId % res->description().c_str()  );
	}
	LOG(BF("Adding atom(%s) to residue(%s) SubstId(%d)") % a->getName().c_str() % res->getPdbName().c_str() % ai->second.mSubstId  );
	res->addAtom(a);
    }

    //
    // Now form all the bonds
    //
    vector<TriposBond>::iterator	bi;
    Atom_sp				a1,a2;
    BondOrder				bo;
    for ( bi=bonds.begin(); bi!=bonds.end(); bi++ ) {
	LOG(BF("Creating bond between atom id: %d - %d") % bi->mAtom1Id % bi->mAtom2Id  );
	a1 = atoms[bi->mAtom1Id].mAtom;
	a2 = atoms[bi->mAtom2Id].mAtom;
	LOG(BF("Creating bond between %s - %s") % a1->description().c_str() % a2->description().c_str()  );
	if ( bi->mOrder=="1" ) bo = singleBond;
	else if ( bi->mOrder=="2" ) bo = doubleBond;
	else if ( bi->mOrder=="3" ) bo = tripleBond;
	else if ( bi->mOrder=="ar" ) bo = aromaticBond;
	else bo = singleBond;
	a1->bondTo(a2,bo);
    }

    //
    // Now create all of the molecules
    //
    if ( sawSubStructure )
    {
	vector<TriposSubstructure>::iterator	si;
	{_BLOCK_TRACE("Creating molecules");
	    for ( si=substructures.begin(); si!=substructures.end(); si++ ) {
                core::Symbol_sp schain = chemkw_intern(si->chain);
		if ( chains.count(schain) == 0 ) {
		    LOG(BF("Creating molecule: %s") % si->chain.c_str()  );
		    GC_ALLOCATE(Molecule_O, tempm );
		    m = tempm;
		    LOG(BF("Setting name") );
		    m->setName(schain);
		    LOG(BF("Adding molecule to aggregate") );
		    aggregate->addMolecule(m);
		    LOG(BF("setting chains index=%s") % si->chain.c_str() );
		    chains[schain] = m;
		    LOG(BF("Done with molecule") );
		}
		m = chains[schain];
		if ( m.nilp() )
		{
		    SIMPLE_ERROR(BF("Could not find molecule with chain=%s") % _rep_(schain));
		}
		LOG(BF("Looking up residue with mId=%d") % si->mId );
		res = residues[si->mId];
		if ( res.nilp() )
		{
		    stringstream serr;
		    serr << "Could not find residue with si->mId=";
		    serr << si->mId;
		    serr << std::endl;
		    SIMPLE_ERROR(BF(serr.str()));
		}
		m->addMatter(residues[si->mId]);
		LOG(BF("Setting residue name at index: %d") % si->mId  );
		LOG(BF("Setting residue name to: %s") % si->subst_name.c_str()  );
		residues[si->mId]->setName(chemkw_intern(si->subst_name));
		LOG(BF("Setting residue PDB name to: %s") % si->sub_type.c_str()  );
		residues[si->mId]->setPdbName(chemkw_intern(si->sub_type));
//		printf("%s:%d Setting residue file_sequence_number(%d)", __FILE__, __LINE__, si->file_sequence_number );
//		residues[si->mId]->setFileSequenceNumber(si->file_sequence_number);
		LOG(BF("Adding residue number: %d name(%s) pdbName(%s) to molecule: %s") % si->mId % residues[si->mId]->getName().c_str() % residues[si->mId]->getPdbName().c_str() % si->chain.c_str()  );
	    }
	}
    } else
    {
	// No SUBSTRUCTURE entry was found, toss all of
	// the residues into one "molecule"
	GC_ALLOCATE(Molecule_O, mol );
        gctools::SmallMap<uint,Residue_sp>::iterator rit;
	for ( rit=residues.begin(); rit!=residues.end(); rit++ )
	{
	    mol->addMatter(rit->second);
	}
	aggregate->addMolecule(mol);
    }

    //
    // Determine the elements and hybridization of each atom
    //
    for ( ai=atoms.begin(); ai!= atoms.end(); ai++ ) {
	_calculateElementAndHybridization(ai->second.mAtom);
    }
}





void	mol2WriteAggregateStream( Aggregate_sp agg, std::ostream &out, bool useSybylTypes )
{
Loop		loop, lRes;
uint		atomCount, bondCount, residueCount;
Atom_sp		a1, a2, a;
Residue_sp	r;
Molecule_sp	m;
gctools::Vec0<AtomInfo> atomList;
gctools::Vec0<ResidueOut> residueList;
ResidueOut	oneResOut;
AtomInfo	one;

 	//
	// Count the atoms
	//
    atomCount = agg->numberOfAtoms();
 	//
	// Count the bonds
	//
    bondCount = 0;
    loop.loopTopGoal(agg,BONDS);
    while ( loop.advanceLoopAndProcess() ) bondCount++;

 	//
	// Count the residues
	//
    residueCount = 0;
    loop.loopTopGoal(agg,RESIDUES);
    while ( loop.advanceLoopAndProcess() ) residueCount++;

    out << "@<TRIPOS>MOLECULE" << std::endl;
    out << agg->getName()->symbolNameAsString() << std::endl;
    out << atomCount << " "
	<< bondCount << " "
	<< residueCount << " 0 0" << std::endl;
    out << "POLYMER" << std::endl;
    out << "CHARGES" << std::endl;
    out << "" << std::endl;
    out << "" << std::endl;

	//
	// Assign an ID to every atom
	//
    LOG(BF("Assigning ID to every atom") );
    uint atomId = 1;
    uint resId = 1;
    lRes.loopTopGoal(agg,RESIDUES);
    while ( lRes.advanceLoopAndProcess() ) {
	r = lRes.getResidue();
	r->setTempFileId(resId);
	oneResOut.res = r;
	oneResOut.firstAtom = atomId;
	residueList.push_back(oneResOut);
	resId++;
	loop.loopTopGoal(r,ATOMS);
	while ( loop.advanceLoopAndProcess() ) {
	    a = loop.getAtom();
	    a->setTempFileId(atomId);
	    one._Atom = a;
            if ( useSybylTypes ) {
              chem::TypeAssignmentRules_sp sybylRules = gc::As<TypeAssignmentRules_sp>(chem::_sym_STARsybyl_type_assignment_rulesSTAR->symbolValue());
              one._Type = sybylRules->calculateType(a);
            } else {
              if ( a->getType() ) {
                one._Type = a->getType();
              } else {
                one._Type = _Nil<core::T_O>();
              }
            }
	    atomList.push_back(one);
	    atomId++;
	}
    }

	//
	// Set the name of every molecule to the chainID.
	//
    LOG(BF("Setting molecule chainID") );
    loop.loopTopGoal(agg,MOLECULES);
    char chainId = 'A';
    string nameId;
    while ( loop.advanceLoopAndProcess() ) {
	m = loop.getMolecule();
	nameId = chainId;
	m->setName(chemkw_intern(nameId));
	chainId++;
        if ( chainId > 'Z' ) chainId = 'A'; // recycle chain ids
    }

	//
	// Write out the atoms
	//
    LOG(BF("Writing ATOMs") );
    out << "@<TRIPOS>ATOM" << std::endl;
    Vector3 pos;
    gctools::Vec0<AtomInfo>::iterator	ai;
    for ( ai=atomList.begin(); ai!= atomList.end(); ai++ ) {
	a = ai->_Atom;
	out << a->getTempFileId() << " ";
	out << a->getName()->symbolNameAsString() << " ";
	pos = a->getPosition();
	out << pos.getX() << " ";
	out << pos.getY() << " ";
	out << pos.getZ() << " ";
	LOG(BF("Writing mol2 atom(%s) pos(%s)") % a->description() % a->getPosition().asString() );
	out << ai->_Type->symbolNameAsString() << " ";
	out << a->containedBy()->getTempFileId() << " ";
	out << a->containedBy()->getName()->symbolNameAsString() << "_"
		<< a->containedBy()->getTempFileId() << " ";
	out << std::setiosflags(std::ios::fixed)
            << std::setw(8)
            << std::setprecision(3);
	out << a->getCharge() << std::endl;
    }
	//
	// Write out the bonds
	//
    LOG(BF("Writing BOND") );
    out << "@<TRIPOS>BOND" << std::endl;
    loop.loopTopGoal(agg,BONDS);
    uint id = 1;
    BondOrder		o;
    string		orderName;
    while ( loop.advanceLoopAndProcess() ) {
	a1 = loop.getBondA1();
	a2 = loop.getBondA2();
	o = loop.getBondOrder();
	orderName = "1";
	if ( o == singleBond ) orderName = "1";
	else if ( o == doubleBond ) orderName = "2";
	else if ( o == tripleBond ) orderName = "3";
	else if ( o == aromaticBond ) orderName = "ar";
	out << id << " ";
	out << a1->getTempFileId() << " ";
	out << a2->getTempFileId() << " ";
	out << orderName << " ";
	out << std::endl;
	id++;
    }
	//
	// Write out the substructures
	//
    LOG(BF("Writing SUBSTRUCTURE") );
    out << "@<TRIPOS>SUBSTRUCTURE" << std::endl;
    gctools::Vec0<ResidueOut>::iterator	ri;
    id = 1;
    for ( ri=residueList.begin(); ri!=residueList.end(); ri++ ) {
	r = ri->res;
//	out << id << " ";
	out << r->getTempFileId() << " ";
	out << r->getName()->symbolNameAsString() << "_" << r->getTempFileId() << " ";
	out << ri->firstAtom << " ";
	out << "RESIDUE 1 ";
	out << downcast<Molecule_O>(r->containedBy())->getName()->symbolNameAsString() << " ";
	out << r->getName()->symbolNameAsString() << " 1";
	out << std::endl;
	id++;
    }
    LOG(BF("Returning") );
}



void	mol2WriteAggregateToFileName( Aggregate_sp agg, core::T_sp fname, bool useSybylTypes )
{
  core::Str_sp sname = gc::As<core::Str_sp>(core::cl__namestring(fname));
  std::ofstream	fout;
  fout.open(sname->get().c_str(),std::ios::out);
  mol2WriteAggregateStream( agg, fout, useSybylTypes );
  fout.close();
}



void	mol2WriteMatterToFileName(Matter_sp matter, core::T_sp fileName, bool useSybylTypes )
{
  if ( Aggregate_sp agg = matter.asOrNull<Aggregate_O>() ) {
    mol2WriteAggregateToFileName(agg,fileName,useSybylTypes);
    return;
  }
  if ( Molecule_sp mol = matter.asOrNull<Molecule_O>() ) {
    Aggregate_sp agg = Aggregate_O::create();
    agg->addMolecule(mol);
    mol2WriteAggregateToFileName(agg,fileName,useSybylTypes);
    return;
  }
  SIMPLE_ERROR(BF("You must pass a Molecule or Aggregate"));
}



void	initialize_mol2_type_rules()
{
  LOG(BF("Initializing sybyl types") );
  GC_ALLOCATE(chem::WildElementDict_O,wilds);
	// WILDATOM XX C N O S P
  SYMBOL_EXPORT_SC_(ChemKwPkg,XX);
  SYMBOL_EXPORT_SC_(ChemKwPkg,XA);
  SYMBOL_EXPORT_SC_(ChemKwPkg,XB);
  SYMBOL_EXPORT_SC_(ChemKwPkg,XC);
  SYMBOL_EXPORT_SC_(ChemKwPkg,XD);

  wilds->addWildName(chemkw::_sym_XX);
  wilds->addWildNameMap(chemkw::_sym_XX,chemkw::_sym_C);
  wilds->addWildNameMap(chemkw::_sym_XX,chemkw::_sym_N);
  wilds->addWildNameMap(chemkw::_sym_XX,chemkw::_sym_O);
  wilds->addWildNameMap(chemkw::_sym_XX,chemkw::_sym_S);
  wilds->addWildNameMap(chemkw::_sym_XX,chemkw::_sym_P);
	// WILDATOM XA O S
  wilds->addWildName(chemkw::_sym_XA);
  wilds->addWildNameMap(chemkw::_sym_XA,chemkw::_sym_O);
  wilds->addWildNameMap(chemkw::_sym_XA,chemkw::_sym_S);
	// WILDATOM XB N P
  wilds->addWildName(chemkw::_sym_XB);
  wilds->addWildNameMap(chemkw::_sym_XB,chemkw::_sym_N);
  wilds->addWildNameMap(chemkw::_sym_XB,chemkw::_sym_P);
	// WILDATOM XC F Cl Br I
  wilds->addWildName(chemkw::_sym_XC);
  wilds->addWildNameMap(chemkw::_sym_XC,chemkw::_sym_F);
  wilds->addWildNameMap(chemkw::_sym_XC,chemkw::_sym_Cl);
  wilds->addWildNameMap(chemkw::_sym_XC,chemkw::_sym_Br);
  wilds->addWildNameMap(chemkw::_sym_XC,chemkw::_sym_I);
	// WILDATOM XD S P
  wilds->addWildName(chemkw::_sym_XD);
  wilds->addWildNameMap(chemkw::_sym_XD,chemkw::_sym_S);
  wilds->addWildNameMap(chemkw::_sym_XD,chemkw::_sym_P);
  chem::TypeAssignmentRules_sp sybylRules = chem::TypeAssignmentRules_O::create();
  sybylRules->setWildElementDict(wilds);
  const char*	sybylTypeRules[] = {
      "ATD  C.3   *   6   4   &",
      "ATD  C.cat *   6   3   *   *   *  		(N3,N3,N3)  &",
      "ATD  C.ar  *   6   3   *   *   [AR1]    &",
      "ATD  C.2   *   6   3   & ",
      "ATD  C.1   *   6   2   & ",
      "ATD  C.1   *   6   1   & ",
      "ATD  H     *   1   &  ",
      "ATD  F     *   9   &",
      "ATD  Cl    *   17  &",
      "ATD  Br    *   35  &      ",
      "ATD  I     *   53  & ",
      "ATD  P.3   *   15  & ",
      "ATD  N.4   *   7   4   & ",
      "ATD  N.am  *   7   3   *   *   *                (C3(XA1))       &",
      "ATD  N.pl3 *   7   3   *   *   *       		(O1,O1)		& ",
      "ATD  N.pl3 *   7   3   *   *   [AR1.AR2.AR3]    &",
      "ATD  N.3   *   7   3   & ",
      "ATD  N.ar  *   7   2   *   *   [AR1]    &",
      "ATD  N.2   *   7   2   & ",
      "ATD  N.1   *   7   1   & ",
      "ATD  O.co2 *   8   1   *   *   *       		(C3(O1))	& ",
      "ATD  O.co2 *   8   1   *   *   *       		(P(O1))		& ",
      "ATD  O.2   *   8   1   & ",
      "ATD  O.3   *   8   2   &",
      "ATD  S.2   *   16  1   & ",
      "ATD  S.3   *   16  2   & ",
      "ATD  S.o   *   16  3   *   *   *       		(O1[DB'])  	&     ",
      "ATD  S.o2  *   16  4   *   *   *       		(O1[DB'],O1[DB'])  	&     ",
      "ATD  s.3   *   16  &  ",
      "ATD  Li    *   3   & ",
      "ATD  Be    *   4   & ",
      "ATD  B     *   5   & ",
      "ATD  Na    *   11  & ",
      "ATD  Mg    *   12  & ",
      "ATD  Al    *   13  & ",
      "ATD  Si    *   14  & ",
      "ATD  K     *   19  & ",
      "ATD  Ca    *   20  & ",
      "ATD  Sr    *   38  & ",
      "ATD  Ba    *   56  & ",
      "ATD  Sc    *   21  & ",
      "ATD  Ti    *   22  & ",
      "ATD  V     *   23  & ",
      "ATD  Cr    *   24  & ",
      "ATD  Mn    *   25  & ",
      "ATD  Fe    *   26  & ",
      "ATD  Co    *   27  & ",
      "ATD  Ni    *   28  & ",
      "ATD  Cu    *   29  & ",
      "ATD  Zn    *   30  & ",
      "ATD  Ga    *   31  & ",
      "ATD  Ge    *   32  & ",
      "ATD  As    *   33  & ",
      "ATD  Se    *   34  & ",
      "ATD  Ru    *   44  & ",
      "ATD  Rh    *   45  & ",
      "ATD  Pd    *   46  & ",
      "ATD  Ag    *   47  & ",
      "ATD  Cd    *   48  & ",
      "ATD  Pt    *   78  & ",
      "ATD  Au    *   79  & ",
      "ATD  Hg    *   80  & ",
      "ATD  Tl    *   81  & ",
      "ATD  Pb    *   82  & ",
      "ATD  lp    *   0   1   &",
      "ATD  ANY   & ",
      ""
  };

  _sym_STARparserNodeHolderSTAR->defparameter(adapt::IndexedObjectBag_O::create());
  for ( uint i=0; strlen(sybylTypeRules[i])>0; i++ )
  {
    const char* typeRuleString = sybylTypeRules[i];
    chem::OneTypeRule_sp rule = chem::OneTypeRule_O::create();
    chem::ChemInfo_sp match = chem::ChemInfo_O::create();
    match->compileAntechamber(typeRuleString,wilds);
    rule->setAssignTypeName(match->getAssignType());
    rule->setCode(match);
    sybylRules->appendRule(rule);
  }
  SYMBOL_EXPORT_SC_(ChemPkg,STARsybyl_type_assignment_rulesSTAR);
  chem::_sym_STARsybyl_type_assignment_rulesSTAR->defparameter(sybylRules);

}




};
