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
       
#define DEBUG_LEVEL_NONE


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
#include <clasp/core/ql.h>
#include <clasp/core/pathname.h>
#include <cando/adapt/indexedObjectBag.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/bformat.h>
#include <clasp/core/hashTableEql.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/loop.h>
#include <cando/chem/molecule.h>
#include <clasp/core/array.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/chemPackage.fwd.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/chemInfo.h>

namespace chem {


std::string residue_molecule_name(const AtomInfo& ai) {
  stringstream ss;
  ss << ai._Residue->getName()->symbolNameAsString();
  ss << "_" << ai._MoleculeIndex << "/" << ai._ResidueIndex;
  return ss.str();
}


void	Mol2File::advanceLine() {
  bool ignoreNextCr;
  char c;
  this->mLine.str("");
  ignoreNextCr = false;
  while (1) {
    c = core::stream_read_char(this->fIn);
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
  if (chem__verbose(2)) {
    core::clasp_write_string(fmt::format("advanceLine mLine = {}\n" , this->mLine.str()));
  }
};
std::queue<string>	Mol2File::splitLine() {
  std::queue<string>	qWords;
  if (chem__verbose(1)) {
    core::clasp_write_string(fmt::format("Splitting read line: {}\n" , this->mLine.str()));
  }
  core::queueSplitString(this->mLine.str(), qWords, " \t" );
  return qWords;
};
bool	Mol2File::hasDataLine() {
  if ( this->eof() ) return false;
  if ( this->mLine.str().size() == 0 ) return false;
  if ( this->mLine.str()[0] != '@' ) return true;
  return false;
};
void	Mol2File::openFileName(core::T_sp pn)
{
  this->fIn = cl__open(pn);
  if ( this->fIn.nilp() ) {
    SIMPLE_ERROR("Mol2 file doesn't exist: {}" , _rep_(pn) );
  }
  this->advanceLine();
}
Mol2File::~Mol2File() {
  cl__close(fIn);
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
  uint  mArrayIndex;
  string	atom_name;
  double	mX;
  double	mY;
  double	mZ;
  string	mType;
  uint	mSubstId;
  string	mSubstName;
  double	mCharge;
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
  size_t pos = nameNum.find("_");
  size_t seqPos, nameLen;
  if ( pos == string::npos ) {
    seqNum = -1;
    name = nameNum;
  } else {
    seqPos = pos + 1;
    nameLen = pos;
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
  string nm = mnm->symbolName()->get_std_string();
    
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


/*! If the mol2 file contains multiple @<TRIPOS>MOLECULE entries then return a list of aggregates and a hash-table with atom-types.
Otherwise return (values aggregate atom-types) */
core::T_mv mol2Read(Mol2File& fIn)
{
  core::HashTable_sp                    atom_types = core::HashTable_O::createEq();
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
  std::string moleculeName;
  uint				lastmId;
  bool				sawSubStructure = false;
  bool  firstMolecule = true;
  core::List_sp boundingBox = nil<core::T_O>();
  core::T_sp solventName = nil<core::T_O>();
  size_t atom_counter = 0;
  while ( !fIn.eof() ) {
    line = fIn.line().str();
    if ( fIn.hasDataLine() ) {
      SIMPLE_ERROR("Mol2 file encountered data line when expecting command line in top level: {}" , fIn.line().str() );
    }
    if (chem__verbose(1)) {
      core::clasp_write_string(fmt::format("Read mol2 line: {}\n" , line));
    }
    if ( line == "@<TRIPOS>MOLECULE" ) {
      if (!firstMolecule) break; // Read only one @<TRIPOS>MOLECULE
      firstMolecule = false;
    }
    fIn.advanceLine();
    if ( line.size() == 0 ) continue;
//	LOG("line: {}" , line.c_str()  );
    if ( line[0] == '#' ) continue;
    if ( line == "@<TRIPOS>MOLECULE" ) {
      LOG("Got @TRIPOS.MOLECULE" );
      moleculeName = fIn.line().str();
      firstMolecule = false;
      mol.mMoleculeName = fIn.line().str();
      fIn.advanceLine(); // num_atoms
      words = fIn.splitLine();
      if ( words.empty() ) continue;
      mol.mAtoms = atoi(words.front().c_str()); words.pop();
      LOG("Got num_atoms={}" , mol.mAtoms  );
      if ( words.empty() ) continue;
      mol.mBonds = atoi(words.front().c_str()); words.pop();
      LOG("Got num_bonds={}" , mol.mBonds  );
      if ( words.empty() ) continue;
      mol.mBonds = atoi(words.front().c_str()); words.pop();
      fIn.advanceLine(); // mol_type
      if ( !fIn.hasDataLine() ) {
        SIMPLE_ERROR("Mol2 file encountered cmd line when expecting data line: {}" , fIn.line().str() );
      }
      fIn.advanceLine(); // charge_type
      while ( fIn.hasDataLine() ) {
        fIn.advanceLine();
      }
      molecules.push_back(mol);
    } else if ( line == "@<TRIPOS>ATOM" ) {
      
      while ( fIn.hasDataLine() ) {
        line = fIn.line().str();
        if ( line[0] == '#' ) {
          fIn.advanceLine();
          continue;
        }
        words = fIn.splitLine();
        oneAtom.mIndex = atoi(words.front().c_str()); words.pop();
        oneAtom.mArrayIndex = atom_counter++;
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
          LOG("Read mSubstId={}" , oneAtom.mSubstId  );
          words.pop();
          if ( !words.empty() ) {
            oneAtom.mSubstName = words.front();
            if (chem__verbose(1)) {
              core::clasp_write_string(fmt::format("Read atom mSubstName: {}\n" , oneAtom.mSubstName));
            }
            LOG("Read mSubstName={}" , oneAtom.mSubstName.c_str() );
            words.pop();
            if ( !words.empty() ) {
              oneWord = words.front();
              LOG("Charge word = |{}|" , oneWord.c_str()  );
              oneAtom.mCharge = atof(oneWord.c_str());
              LOG("Read charge ={}" , oneAtom.mCharge  );
            }
          }
        }
        atoms[oneAtom.mIndex] = oneAtom;
        //        printf("%s:%d Adding atom[%lu] -> %s\n", __FILE__, __LINE__, oneAtom.mArrayIndex, oneAtom.atom_name.c_str());
        LOG("Added atom with id: {} name: {}" , oneAtom.mArrayIndex , oneAtom.atom_name.c_str()  );
        fIn.advanceLine();
      }
    } else if ( line == "@<TRIPOS>BOND" ) {
      
      while ( fIn.hasDataLine() ) {
        line = fIn.line().str();
        LOG("line len({}) line: {}" , line.size() , line.c_str()  );
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
    } else if ( line == "@<TRIPOS>FF_PBC" ) {
      line = fIn.line().str();
      if (chem__verbose(1)) {
        core::clasp_write_string("@<TRIPOSE>FF_PBC\n");
        core::clasp_write_string(fmt::format("{}\n" , line));
      }
      words = fIn.splitLine();
      fIn.advanceLine();
      if (words.empty()) {
        SIMPLE_ERROR("No data followed @<TRIPOSE>FF_PBC");
      }
      words.pop(); // get rid of version
      int pbc_type = atoi(words.front().c_str());
      words.pop(); // pop pbc type
      if (pbc_type == 0) {
        // Skip the pbc - there is none
        for ( size_t xx = 0; xx<6; ++xx) words.pop();
        boundingBox = nil<core::T_O>();
      } else if (pbc_type ==1) {
        // Rectanglar pbc
        double values[6];
        for ( size_t xx = 0; xx<6; ++xx) {
          values[xx] = atof(words.front().c_str());
          words.pop();
        }
        boundingBox = core::Cons_O::createList(core::clasp_make_double_float(values[3]-values[0]),
                                               core::clasp_make_double_float(values[4]-values[1]),
                                               core::clasp_make_double_float(values[5]-values[2]));
      } else {
        SIMPLE_ERROR("Unknown solvent box type in FF_PBC: {}" , line);
      }
      // Now get the solvent name
      if (words.size()>0) {
        solventName = _lisp->intern(words.front(),KeywordPkg);
        words.pop();
      }
    } else if ( line == "@<TRIPOS>SUBSTRUCTURE" ) {
      
      sawSubStructure = true;
      lastmId = -1;
      string	nameNum;
      while ( fIn.hasDataLine() ) {
        line = fIn.line().str();
        LOG("line = ({})" , line.c_str()  );
        if ( line[0] == '#' ) {
          fIn.advanceLine();
          continue;
        }
        words = fIn.splitLine();
        oneSubst.mId = atoi(words.front().c_str()); words.pop();
        if ( oneSubst.mId == lastmId ) {
          SIMPLE_ERROR("First entry of SUBSTRUCTURE records are not valid");
        }
        lastmId = oneSubst.mId;
        nameNum = words.front(); words.pop();
        splitStringIntoNameAndNumber( nameNum, oneSubst.subst_name,oneSubst.file_sequence_number);
        LOG("oneSubst.root_atom str({})" , words.front().c_str()  );
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
                LOG("Got sub_type = {}" , oneSubst.sub_type.c_str()  );
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
        LOG("Read substituent  mId={} subst_type({})" , oneSubst.mId , oneSubst.subst_type.c_str() );
        substructures.push_back(oneSubst);
      }
    } else
    {
      core::clasp_write_string(fmt::format("Ignoring mol2 command: {}\n" , line ));
      while ( fIn.hasDataLine() ) {
        fIn.advanceLine();
      }
    }
  }
  if (molecules.size()==0) return Values( nil<core::T_O>(), nil<core::T_O>());
    //
    // First create all of the Aggregate
    //	(Tripos calls this a molecule and it can have multiple chains.)
    //
  Aggregate_sp aggregate = Aggregate_O::create();
  if (boundingBox.notnilp()) {
    BoundingBox_sp bb = BoundingBox_O::make(boundingBox,nil<core::T_O>(),nil<core::T_O>());
    aggregate->setBoundingBox(bb);
  }
  aggregate->setName(chemkw_intern(moleculeName));
    
    // Now create all of the atoms and the residues and put the
    // atoms in the residues
  map<uint,TriposAtom>::iterator	ai;
  Vector3				pos;
  string				nm, el;
  core::SimpleVector_sp atoms_sv = core::SimpleVector_O::make(atoms.size());
  for ( ai=atoms.begin(); ai!= atoms.end(); ai++ ) {
    auto  tempa  = gctools::GC<Atom_O>::allocate_with_default_constructor();
    a = tempa;
    a->setName(chemkw_intern(ai->second.atom_name));
    LOG("Creating atom with id: {} name: {} element: {} charge: %8.2f" , ai->second.mIndex , ai->second.atom_name.c_str() , el.c_str() , ai->second.mCharge  );
    pos.set(ai->second.mX,ai->second.mY,ai->second.mZ);
    a->setPosition(pos);
    atom_types->setf_gethash(a,chemkw_intern(ai->second.mType));
    a->setCharge(ai->second.mCharge);
    LOG(" atom info: {}" , a->description().c_str()  );
    if (ai->second.mArrayIndex<atoms.size()) {
      (*atoms_sv)[ai->second.mArrayIndex] = a;
    } else {
      SIMPLE_ERROR("Bad index {} for atom atoms.size() -> {}\n" , ai->second.mArrayIndex , atoms.size());
    }
	//
	// If we haven't seen this residue before
	// we need to create it
	//
    Residue_sp res;
    if ( residues.count(ai->second.mSubstId)==0 ) {
      auto  tempRes  = gctools::GC<Residue_O>::allocate_with_default_constructor();
      res = tempRes;
	    // April 06/2007 decided that its best to set the 
	    // residue name from the substructure records below
	    //
	    // But I'll set them here incase there isn't 
	    // a substructure entry
      if (chem__verbose(1)) {
        core::clasp_write_string(fmt::format("Setting residue name: {} with id {}\n" , ai->second.mSubstName , ai->second.mSubstId ));
      }
      res->setName(chemkw_intern(ai->second.mSubstName));
      res->setPdbName(chemkw_intern(ai->second.mSubstName));
//            printf("%s:%d setFileSequenceNumber --> %d\n", __FILE__, __LINE__, ai->second.mSubstId);
      res->setFileSequenceNumber(ai->second.mSubstId);
      LOG("Creating residue with id({}) and name({})" , ai->second.mSubstId , ai->second.mSubstName.c_str()  );
      residues[ai->second.mSubstId] = res;
      LOG("New residue {}" , res->description().c_str()  );
    } else {
      res = residues[ai->second.mSubstId];
      LOG("Reusing residue with mSubstId={} {}" , ai->second.mSubstId , res->description().c_str()  );
    }
    LOG("Adding atom({}) to residue({}) SubstId({})" , a->getName().c_str() , res->getPdbName().c_str() , ai->second.mSubstId  );
    res->addAtom(a);
  }

    //
    // Now form all the bonds
    //
  vector<TriposBond>::iterator	bi;
  Atom_sp				a1,a2;
  BondOrder				bo;
  for ( bi=bonds.begin(); bi!=bonds.end(); bi++ ) {
    LOG("Creating bond between atom id: {} - {}" , bi->mAtom1Id , bi->mAtom2Id  );
    TriposAtom& at1 = atoms[bi->mAtom1Id];
    TriposAtom& at2 = atoms[bi->mAtom2Id];
    a1 = gc::As<Atom_sp>((*atoms_sv)[at1.mArrayIndex]);
    a2 = gc::As<Atom_sp>((*atoms_sv)[at2.mArrayIndex]);
    // printf("%s:%d Bonding atom[%lu] to atom[%lu]\n", __FILE__, __LINE__, bi->mAtom1Id, bi->mAtom2Id);
    LOG("Creating bond between {} - {}" , a1->description().c_str() , a2->description().c_str()  );
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
    {
      for ( si=substructures.begin(); si!=substructures.end(); si++ ) {
        core::Symbol_sp schain = chemkw_intern(si->chain);
        if ( chains.count(schain) == 0 ) {
          LOG("Creating molecule: {}" , si->chain.c_str()  );
          auto  tempm  = gctools::GC<Molecule_O>::allocate_with_default_constructor();
          m = tempm;
          LOG("Setting name" );
          m->setName(schain);
          LOG("Adding molecule to aggregate" );
          aggregate->addMolecule(m);
          LOG("setting chains index={}" , si->chain.c_str() );
          chains[schain] = m;
          LOG("Done with molecule" );
        }
        m = chains[schain];
        if ( m.nilp() )
        {
          SIMPLE_ERROR("Could not find molecule with chain={}" , _rep_(schain));
        }
        LOG("Looking up residue with mId={}" , si->mId );
        if (!residues.contains(si->mId)) {
          if (residues.size()==1) {
            // If there is only one residue - then use that.
            gctools::SmallMap<uint,Residue_sp>::iterator it = residues.begin();
            res = it->second;
          } else {
            stringstream ss;
            for (auto it = residues.begin(); it<residues.end(); ++it ) {
              if (it!=residues.begin()) {
                ss << ", ";
              }
              ss << it->first;
            }
            SIMPLE_ERROR("Could not find residue with subst_id {} - there are {} total residues with subst_ids {}" , si->mId , residues.size() , ss.str());
          }
        } else {
          res = residues[si->mId];
        }
        if (chem__verbose(1)) {
          core::clasp_write_string(fmt::format("About to add residue id {} to molecule - number of residues {}\n" , si->mId , residues.size()));
          core::clasp_write_string(fmt::format("About to add residue {} to molecule\n" , _rep_(res)));
        }
        m->addMatter(res);
        LOG("Setting residue name at index: {}" , si->mId  );
        LOG("Setting residue name to: {}" , si->subst_name.c_str()  );
        core::T_sp sub_name = chemkw_intern(si->subst_name);
        if (chem__verbose(1)) {
          core::clasp_write_string(fmt::format("Setting substructure residue name: {}\n" , _rep_(sub_name)));
        }
        res->setName(sub_name);
        LOG("Setting residue PDB name to: {}" , si->sub_type.c_str()  );
        core::T_sp sub_type = chemkw_intern(si->sub_type);
        res->setPdbName(sub_type);
        if (solventName.notnilp() && (sub_name==solventName || sub_type==solventName))
        {
          //          printf("%s:%d Setting molecule %s to solvent\n", __FILE__, __LINE__, _rep_(m).c_str());
          m->setf_molecule_type(kw::_sym_solvent);
        }
        LOG("Adding residue number: {} name({}) pdbName({}) to molecule: {}" , si->mId , res->getName().c_str() , res->getPdbName().c_str() , si->chain.c_str()  );
      }
    }
  } else
  {
	// No SUBSTRUCTURE entry was found, toss all of
	// the residues into one "molecule"
    auto  mol  = gctools::GC<Molecule_O>::allocate_with_default_constructor();
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
    Atom_sp aa = gc::As<Atom_sp>((*atoms_sv)[ai->second.mArrayIndex]);
    _calculateElementAndHybridization(aa);
  }
  return Values(aggregate,atom_types);
}

SYMBOL_EXPORT_SC_(ChemPkg,initialize_sybyl_type_rules);
SYMBOL_EXPORT_SC_(ChemPkg,STARsybyl_type_assignment_rulesSTAR);


core::HashTable_sp assignSybylTypes( Matter_sp matter) {
  if (chem::_sym_STARsybyl_type_assignment_rulesSTAR->symbolValue().nilp()) {
    SIMPLE_ERROR("The chem:*sybyl-type-assignment-rules* dynamic variable needs to be bound to a chem:fftypes-db database for sybyl types");
  }
  FFTypesDb_sp sybylRules = gc::As<FFTypesDb_sp>(chem::_sym_STARsybyl_type_assignment_rulesSTAR->symbolValue());
  core::HashTable_sp atom_types = core::HashTable_O::createEq();
  Loop lAtoms;
  lAtoms.loopTopGoal(matter,ATOMS);
  while ( lAtoms.advanceLoopAndProcess() ) {
    Atom_sp a = lAtoms.getAtom();
    core::Symbol_sp type = sybylRules->assignType(a);
    atom_types->setf_gethash(a,type);
  }
  return atom_types;
}

core::HashTable_sp	mol2WriteAggregateStream( Aggregate_sp agg, std::ostream &out, core::HashTable_sp atom_types )
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
  core::HashTable_sp ht = core::HashTable_O::createEql();
 
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
  LOG("Assigning ID to every atom" );
  uint atomId = 1;
  uint resId = 1;
  size_t moleculeIndex = 0;
  chem::FFTypesDb_sp sybylRules;
  Loop lMol;
  lMol.loopTopGoal(agg,MOLECULES);
  while (lMol.advanceLoopAndProcess()) {
    Molecule_sp mol = lMol.getMolecule();
    oneResOut.mol = mol;
    size_t residueIndex = 0;
    lRes.loopTopGoal(mol,RESIDUES);
    while ( lRes.advanceLoopAndProcess() ) {
      r = lRes.getResidue();
      ht->setf_gethash(r,core::clasp_make_fixnum(resId));
      oneResOut.res = r;
      oneResOut.firstAtom = atomId;
      residueList.push_back(oneResOut);
      resId++;
      loop.loopTopGoal(r,ATOMS);
      while ( loop.advanceLoopAndProcess() ) {
        a = loop.getAtom();
        ht->setf_gethash(a,core::clasp_make_fixnum(atomId));
        one._Atom = a;
        one._MoleculeIndex = moleculeIndex;
        one._ResidueIndex = residueIndex;
        if ( a->getType(atom_types).notnilp() ) {
          core::Symbol_sp type = gc::As<core::Symbol_sp>(a->getType(atom_types));
          one._Type = type;
        } else {
          one._Type = nil<core::Symbol_O>();
        }
        one._Residue = r;
        atomList.push_back(one);
        atomId++;
      }
      residueIndex++;
    }
    moleculeIndex++;
  }


        // Figure out if there is a solvent residue and a bounding box

  core::Symbol_sp solventName = nil<core::T_O>();
  core::T_sp boundingBox = nil<core::T_O>();
  {
    Loop lMol;
    loop.loopTopGoal(agg,MOLECULES);
    while ( loop.advanceLoopAndProcess() ) {
      m = loop.getMolecule();
      if (m->molecule_type() == kw::_sym_solvent) {
        size_t numResidues = m->contentSize();
        if (numResidues==1) {
          solventName = m->contentAt(0)->getName();
          printf("%s:%d  Read solvent name: %s\n", __FILE__, __LINE__, _rep_(solventName).c_str());
          goto GOT_SOLVENT;
        } else {
          printf("%s:%d There is more than one residue in a solvent molecule %s\n", __FILE__, __LINE__, _rep_(m).c_str());
        }
      }
    }
  GOT_SOLVENT:
    if (agg->boundingBoxBoundP()) {
      boundingBox = agg->boundingBox();
    }
  }

  if (boundingBox.notnilp()) {
    out << "@<TRIPOS>FF_PBC" << std::endl;
    out << "0.1 ";
    if (boundingBox.nilp()) {
      out << "0 0.0 0.0 0.0 0.0 0.0 0.0 ";
    } else {
      BoundingBox_sp bb = gc::As<BoundingBox_sp>(boundingBox);
      out << "1 0.0 0.0 0.0 ";
      if (bb->cuboidp()) {
        out << bb->get_x_width() << " ";
        out << bb->get_y_width() << " ";
        out << bb->get_z_width() << " ";
      } else {
        SIMPLE_WARN("Mol2 only supports cuboid bounding-box - you have {}", _rep_(bb));
      }
    }
    if (solventName.notnilp()) {
      out << solventName->symbolNameAsString() << " ";
    }
    out << std::endl;
    out << std::endl;
  }

#if 0
    // Don't change the names of molecules
    
	//
	// Set the name of every molecule to the chainID.
	//
  LOG("Setting molecule chainID" );
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
#endif
	//
	// Write out the atoms
	//
  LOG("Writing ATOMs" );
  out << "@<TRIPOS>ATOM" << std::endl;
  Vector3 pos;
  gctools::Vec0<AtomInfo>::iterator	ai;
  for ( ai=atomList.begin(); ai!= atomList.end(); ai++ ) {
    a = ai->_Atom;
    int atom_temp_int = ht->gethash(a).unsafe_fixnum();
    out << atom_temp_int << " ";
    out << a->getName()->symbolNameAsString() << " ";
    pos = a->getPosition();
    out << fmt::format("{:5.3f} {:5.3f} {:5.3f} " , pos.getX() , pos.getY() , pos.getZ() );
    LOG("Writing mol2 atom({}) pos({})" , a->description() , a->getPosition().asString() );
    out << ai->_Type->symbolNameAsString() << " ";
    int residue_temp_int = ht->gethash(ai->_Residue).unsafe_fixnum();
    out << residue_temp_int << " ";
    out << residue_molecule_name(*ai);
    out << " ";
    out << fmt::format("{:5.4}" , a->getCharge()) << std::endl;
  }
	//
	// Write out the bonds
	//
  LOG("Writing BOND" );
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
    int a1_ti = ht->gethash(a1).unsafe_fixnum();
    int a2_ti = ht->gethash(a2).unsafe_fixnum();
    out << a1_ti << " ";
    out << a2_ti << " ";
    out << orderName << " ";
    out << std::endl;
    id++;
  }
	//
	// Write out the substructures
	//
  LOG("Writing SUBSTRUCTURE" );
  out << "@<TRIPOS>SUBSTRUCTURE" << std::endl;
  gctools::Vec0<ResidueOut>::iterator	ri;
  id = 1;
  for ( ri=residueList.begin(); ri!=residueList.end(); ri++ ) {
    r = ri->res;
//	out << id << " ";
    int r_ti = ht->gethash(r).unsafe_fixnum();
    out << r_ti << " ";
    out << r->getName()->symbolNameAsString() << " "; // << "_" << r_ti << " ";
    out << ri->firstAtom << " ";
    out << "RESIDUE 1 ";
    out << (ri->mol->getName()->symbolNameAsString()) << " ";
    out << r->getName()->symbolNameAsString() << " 1";
    out << std::endl;
    id++;
  }
  LOG("Returning" );
  return ht;
}



void	mol2WriteAggregateToFileName( Aggregate_sp agg, core::T_sp fname, core::HashTable_sp atom_types )
{
  std::string sname = gc::As<core::String_sp>(core::cl__namestring(fname))->get_std_string();
  std::ofstream	fout;
  fout.open(sname.c_str(),std::ios::out);
  mol2WriteAggregateStream( agg, fout, atom_types );
  fout.close();
}



void	mol2WriteMatterToFileName(Matter_sp matter, core::T_sp fileName, core::HashTable_sp atom_types )
{
  if ( Aggregate_sp agg = matter.asOrNull<Aggregate_O>() ) {
    mol2WriteAggregateToFileName(agg,fileName,atom_types);
    return;
  }
  if ( Molecule_sp mol = matter.asOrNull<Molecule_O>() ) {
    Aggregate_sp agg = Aggregate_O::create();
    agg->addMolecule(mol);
    mol2WriteAggregateToFileName(agg,fileName,atom_types);
    return;
  }
  SIMPLE_ERROR("You must pass a Molecule or Aggregate");
}





};
