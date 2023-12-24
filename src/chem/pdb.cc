/*
    File: pdb.cc
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
#define	DEBUG_LEVEL_NONE
#include <iomanip>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/chem/pdb.h>
#include <clasp/core/pathname.h>
#include <cando/chem/chemPackage.fwd.h>
#include <clasp/core/lispStream.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/atom.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/mbbCoreTools.h>
#include <cando/chem/residue.h>
#include <cando/chem/bond.h>
#include <cando/chem/loop.h>
#include <clasp/core/wrappers.h>


namespace chem
{


struct ConnectPdbRec
{
  int	_atom1;
  int	_bonded[4];

  void write(core::T_sp fout);
  void parse(const string& line);
};







struct	EntirePdbRec
{
  gctools::Vec0<AtomPdbRec>	_atoms;
  gctools::Vec0<Molecule_sp>	_molecules;
  gctools::Vec0<Residue_sp>	_residues;
  std::vector<ConnectPdbRec>	_connects;
  std::map<int,int>		_atomRecIndexFromSerial;

  void startMolecule();
  void addAtomPdbRec(AtomPdbRec& atom);
  void addConnect(const ConnectPdbRec& connect);
  Aggregate_sp createAggregate();
};



string pdb_substr(const string& str, int firstChar1, int lastChar1)
{
  int strLen = str.size();
  stringstream chrs;
  int	firstChar = firstChar1 - 1;
  int lastChar = lastChar1 - 1;
  chrs.str("");
  for ( int x=firstChar; x<=lastChar && x<strLen; x++ )
  {
    if ( str[x] > ' ' )
    {  
      chrs << str[x];
    }
  }
  return chrs.str();
}

double pdb_substr_double(const string& str, int firstChar1, int lastChar1, double dflt )
{
  string part;
  part = pdb_substr(str,firstChar1,lastChar1);
  if ( part == "" ) return dflt;
  return atof(part.c_str());
}

int pdb_substr_int(const string& str, int firstChar1, int lastChar1, int dflt )
{
  string part;
  part = pdb_substr(str,firstChar1,lastChar1);
  if ( part == "" ) return dflt;
  return atoi(part.c_str());
}



void AtomPdbRec::write(core::T_sp fout)
{
  string name;
  if ( this->_name->symbolName()->length() > 3 )
  {
    string sname = this->_name->symbolName()->get_std_string();
    name = sname; // sname.substr(3,1)+sname.substr(0,3);
  } else
  {
    name = " " + this->_name->symbolName()->get_std_string();
  }
  core::clasp_write_string(fmt::format("ATOM  {:5d} {:<4} {:3} {:1}{:4}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}           {:2}\n",
                                       this->_serial, name, /* altLoc, */ this->_resName->symbolNameAsString(), this->_chainId,
                                       this->_resSeq, /* iCode, */ this->_x, this->_y, this->_z, this->_occupancy,
                                       this->_tempFactor, this->_element /* charge */),
                           fout);
}

void	EntirePdbRec::addAtomPdbRec(AtomPdbRec& atom)
{
  this->_atoms.push_back(atom);
  this->_atomRecIndexFromSerial[atom._serial] = this->_atoms.size()-1;
}


void	EntirePdbRec::addConnect(const ConnectPdbRec& connect)
{
  this->_connects.push_back(connect);
}





void	AtomPdbRec::parse(const string& line)
{
  this->_line = line;
  this->_recordName = pdb_substr(line,1,6);
  this->_serial = pdb_substr_int(line,7,11,-1);
  string name3 = pdb_substr(line,14,16);
  string name1 = pdb_substr(line,13,13);
  this->_name = chemkw_intern(name3+name1);
  this->_altLoc = pdb_substr(line,17,17);
  this->_resName = chemkw_intern(pdb_substr(line,18,20));
  this->_chainId = pdb_substr(line,22,22);
  this->_resSeq = pdb_substr_int(line,23,26,0);
  this->_iCode = pdb_substr(line,27,27);
  this->_x = pdb_substr_double(line,31,38,0.0);
  this->_y = pdb_substr_double(line,39,46,0.0);
  this->_z = pdb_substr_double(line,47,54,0.0);
  this->_occupancy = pdb_substr_double(line,55,60,0.0);
  this->_tempFactor = pdb_substr_double(line,61,66,0.0);
  this->_element = pdb_substr(line,77,78);
  this->_charge = pdb_substr(line,79,80);
}


Atom_sp AtomPdbRec::createAtom()
{
  Atom_sp atom = Atom_O::create();
  atom->setName(this->_name);
  Vector3 pos;
  pos.set(this->_x,this->_y,this->_z);
  atom->setPosition(pos);
  if ( this->_element != "" )
  {
    atom->setElementFromString(this->_element);
  } else
  {
    atom->setElement(elementFromAtomNameStringCaseInsensitive(this->_name->symbolName()->get_std_string()));
  }
  return atom;
}





void	ConnectPdbRec::parse(const string& line)
{
  this->_atom1   = pdb_substr_int(line,  7, 11, -1 );
  this->_bonded[0] = pdb_substr_int(line, 12, 16, -1 );
  this->_bonded[1] = pdb_substr_int(line, 17, 21, -1 );
  this->_bonded[2] = pdb_substr_int(line, 22, 26, -1 );
  this->_bonded[3] = pdb_substr_int(line, 27, 31, -1 );
}


void	ConnectPdbRec::write(core::T_sp stream)
{
  stringstream fout;
  fout << "CONECT";
  fout << std::setw(5);
  fout << this->_atom1;
  for ( uint i=0; i< 4; i++ )
  {
    if ( this->_bonded[i]>0 )
    {
      fout << std::setw(5) << this->_bonded[i];
    }
  }
  fout << std::endl;
  core::clasp_write_string(fout.str(),stream);
}


Aggregate_sp	EntirePdbRec::createAggregate()
{
  Aggregate_sp agg = Aggregate_O::create();
  {
    int moleculeIdx = -1;
    int residueSeq = -1;
    int currentResidueIdx = -1;
    gctools::Vec0<AtomPdbRec>::iterator	ai;
    for ( ai=this->_atoms.begin();ai!=this->_atoms.end(); ai++ )
    {
		// If the moleculeIdx has not been seen then
		// create a new molecule entry for this index
      if ( moleculeIdx != ai->_moleculeIdx )
      { 
        moleculeIdx = ai->_moleculeIdx;
        ASSERT((unsigned)(moleculeIdx)==this->_molecules.size());
        this->_molecules.resize(moleculeIdx+1,nil<Molecule_O>());
        Molecule_sp mol = Molecule_O::create();
        mol->setName(chemkw_intern(ai->_chainId));
        this->_molecules[moleculeIdx] = mol;
		    // whenever a new molecule is created
		    // the residueSeq should be reset
        residueSeq = -1;	
        agg->addMatter(mol);
      }
		// If the residueSeq has not been seen then
		// create a new residue entry and assign its residueIdx
		// and put it in its molecule
      LOG("residueSeq({})   ai->_resSeq({})" , residueSeq , ai->_resSeq  );
      if ( residueSeq != ai->_resSeq)
      { 
        residueSeq = ai->_resSeq;
        Residue_sp res = Residue_O::create();
        res->setName(ai->_resName);
        res->setPdbName(ai->_resName);
//        res->setId(ai->_resSeq);
        res->setFileSequenceNumber(ai->_resSeq);
        Molecule_sp mol = this->_molecules[ai->_moleculeIdx];
        mol->addMatter(res);
        ai->_residueIdx = this->_residues.size();
        currentResidueIdx = ai->_residueIdx;
        this->_residues.push_back(res);
      }
	    	// Now create the atom and put it in its residue
      Atom_sp atom = ai->createAtom();
      ai->_atom = atom;
      LOG("currentResidueIdx = {}" , currentResidueIdx  );
      ASSERTF(currentResidueIdx>=0,("residue index must be positive - line--> {}") , ai->_line );
      ASSERT_lessThan(currentResidueIdx,(int)(this->_residues.size()));
      ai->_residueIdx = currentResidueIdx;
      Residue_sp res = this->_residues[ai->_residueIdx];
      res->addMatter(atom);
    }
  }
  {
    LOG("There are {} bonds" , this->_connects.size()  );
    for ( unsigned ci=0; ci<this->_connects.size(); ci++ )
    { 
      LOG("CONECT indices: {} {} {} {} {}" , this->_connects[ci]._atom1 , this->_connects[ci]._bonded[0] , this->_connects[ci]._bonded[1] , this->_connects[ci]._bonded[2] , this->_connects[ci]._bonded[3]  );
      AtomPdbRec*	atomPtr;
      AtomPdbRec*	bondedPtr;
      int bondedIdx, bondedAtomRecIdx;
      LOG("Looking atom with index: {}" , this->_connects[ci]._atom1  );
      int atomRecIndex = this->_atomRecIndexFromSerial[this->_connects[ci]._atom1];
      atomPtr = &(this->_atoms[atomRecIndex]);
      LOG("Found atomPdbRec pointer at {}" , atomPtr  );
      LOG("Found atom entry with serial#{}" , atomPtr->_serial  );
      LOG("Atom.get = {}" , atomPtr->_atom.get()  );
      ANN(atomPtr->_atom);
      ASSERTP(atomPtr->_atom.notnilp(), "Atom is nil");
      bondedIdx = this->_connects[ci]._bonded[0];
      if ( bondedIdx >= 0 )
      {
        LOG("Creating bond2 between atoms with indices {}-{}" , atomRecIndex , bondedIdx );
        bondedAtomRecIdx = this->_atomRecIndexFromSerial[bondedIdx];
        bondedPtr = &(this->_atoms[bondedAtomRecIdx]);
        ANN(bondedPtr->_atom);
        ASSERTP(bondedPtr->_atom.notnilp(), "Atom is nil");
        atomPtr->_atom->bondTo(bondedPtr->_atom,singleBond);
      }
      bondedIdx = this->_connects[ci]._bonded[1];
      if ( bondedIdx >= 0 )
      {
        LOG("Creating bond3 between atoms with indices {}-{}" , atomRecIndex , bondedIdx );
        bondedAtomRecIdx = this->_atomRecIndexFromSerial[bondedIdx];
        bondedPtr = &(this->_atoms[bondedAtomRecIdx]);
        ANN(bondedPtr->_atom);
        ASSERTP(bondedPtr->_atom.notnilp(), "Atom is nil");
        atomPtr->_atom->bondTo(bondedPtr->_atom,singleBond);
      }
      bondedIdx = this->_connects[ci]._bonded[2];
      if ( bondedIdx >= 0 )
      {
        LOG("Creating bond4 between atoms with indices {}-{}" , atomRecIndex , bondedIdx );
        bondedAtomRecIdx = this->_atomRecIndexFromSerial[bondedIdx];
        bondedPtr = &(this->_atoms[bondedAtomRecIdx]);
        ANN(bondedPtr->_atom);
        ASSERTP(bondedPtr->_atom.notnilp(), "Atom is nil");
        atomPtr->_atom->bondTo(bondedPtr->_atom,singleBond);
      }
      bondedIdx = this->_connects[ci]._bonded[3];
      if ( bondedIdx >= 0 )
      {
        LOG("Creating bond5 between atoms with indices {}-{}" , atomRecIndex , bondedIdx );
        bondedAtomRecIdx = this->_atomRecIndexFromSerial[bondedIdx];
        bondedPtr = &(this->_atoms[bondedAtomRecIdx]);
        ANN(bondedPtr->_atom);
        ASSERTP(bondedPtr->_atom.notnilp(), "Atom is nil");
        atomPtr->_atom->bondTo(bondedPtr->_atom,singleBond);
      }
    }
  }
  return agg;
}



Aggregate_sp PdbReader_O::loadPdb(core::T_sp fileName)
{
  PdbReader_sp pdb = PdbReader_O::create();
  core::T_sp stream = cl__open(fileName,
                               core::StreamDirection::input,
                               cl::_sym_character,
                               core::StreamIfExists::nil, false,
                               core::StreamIfDoesNotExist::nil, false,
                               kw::_sym_default,
                               nil<T_O>());
  Aggregate_sp agg = pdb->parse(stream);
  return agg;
}

Aggregate_sp PdbReader_O::loadPdbFromStreamConnectAtoms(core::T_sp stream)
{
  PdbReader_sp pdb = PdbReader_O::create();
  Aggregate_sp agg = pdb->parse(stream);
  chem__connectAtomsInMatterInCovalentContact(agg);
  return agg;
}

Aggregate_sp PdbReader_O::loadPdbConnectAtoms(core::T_sp fileName)
{
  PdbReader_sp pdb = PdbReader_O::create();
  core::T_sp stream = cl__open(fileName,
                               core::StreamDirection::input,
                               cl::_sym_character,
                               core::StreamIfExists::nil, false,
                               core::StreamIfDoesNotExist::nil, false,
                               kw::_sym_default,
                               nil<T_O>());
  Aggregate_sp agg = pdb->parse(stream);
  chem__connectAtomsInMatterInCovalentContact(agg);
  return agg;
}



    
    
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__simple_load_pdb(core::String_sp fileName)
{
  Aggregate_sp agg = PdbReader_O::loadPdbConnectAtoms(fileName);
  return agg;
}

DOCGROUP(cando);
CL_DEFUN core::T_sp chem__load_pdb_from_stream(core::T_sp stream)
{
  Aggregate_sp agg = PdbReader_O::loadPdbFromStreamConnectAtoms(stream);
  return agg;
}

void	PdbReader_O::initialize()
{
  this->Base::initialize();
}

#ifdef XML_ARCHIVE
void	PdbReader_O::archive(core::ArchiveP node)
{
  IMPLEMENT_ME();
}
#endif



Aggregate_sp PdbReader_O::parse(core::T_sp stream)
{
  EntirePdbRec pdbRec;
  int moleculeIdx = 0;
  {
    while ( true ) {
      core::T_mv line = cl__read_line(stream,nil<core::T_O>(),nil<core::T_O>(),nil<core::T_O>());
      if (line.nilp()) break;
      string oneLine = gc::As<core::String_sp>(line)->get_std_string();
      LOG("PDB line |{}|" , oneLine.c_str()  );
      string recordName = pdb_substr(oneLine,1,6);
      if ( recordName == "ATOM" ) {
        AtomPdbRec atomRec;
        atomRec.parse(oneLine);
        atomRec._moleculeIdx = moleculeIdx;
        pdbRec.addAtomPdbRec(atomRec);
      } else if ( recordName == "HETATM" ) {
        AtomPdbRec atomRec;
        atomRec.parse(oneLine);
        atomRec._moleculeIdx = moleculeIdx;
        pdbRec.addAtomPdbRec(atomRec);
      } else if ( recordName == "CONECT" )
      {
        ConnectPdbRec connectRec;
        connectRec.parse(oneLine);
        pdbRec._connects.push_back(connectRec);
      } else if ( recordName == "TER" )
      {
        moleculeIdx++;
      } else if ( recordName == "END" )
      {
        moleculeIdx++;
      }
    }
  }
  return pdbRec.createAggregate();
}

#if INIT_TO_FACTORIES

#define ARGS_PdbWriter_O_make "(file_name)"
#define DECL_PdbWriter_O_make ""
#define DOCS_PdbWriter_O_make "make PdbWriter"
PdbWriter_sp PdbWriter_O::make(core::T_sp fileName)
{
  auto  me  = gctools::GC<PdbWriter_O>::allocate_with_default_constructor();
  me->open(fileName);
  return me;
};

#else

core::T_sp PdbWriter_O::__init__(core::Function_sp exec, core::Cons_sp args,
                                 core::Environment_sp env)
{
  string fn = env->lookup(_lisp->internWithPackageName(ChemPkg,"fileName")).as<core::Str_O>()->get();
  this->open(fn);
  return _Nil<core::T_O>();
}

#endif


void	PdbWriter_O::initialize()
{
  this->Base::initialize();
  this->_Out = nil<core::T_O>();
}



size_t _setupAtomAndConnectRecordsForOneMolecule(Molecule_sp mol, 
                                               gctools::Vec0<AtomPdbRec>& pdbAtoms, 
                                               vector<ConnectPdbRec>& pdbConnects, 
                                               char chainId, size_t serialStart )
{
  pdbAtoms.clear();
  pdbConnects.clear();
  size_t serial = serialStart;
  uint resSeq = 1;
  Loop lResidues;
  lResidues.loopTopGoal(mol,RESIDUES);
  while ( lResidues.advance() )
  {
    Residue_sp res = lResidues.getResidue();
    Loop lAtoms;
    lAtoms.loopTopGoal(res,ATOMS);
    while ( lAtoms.advance() )
    {
      Atom_sp a = lAtoms.getAtom();
      AtomPdbRec atom;
      atom._atom = a;
      atom._name = a->getName();
      atom._element = a->getElementAsString();
      atom._recordName = "ATOM";
      atom._serial = serial++;
      atom._resName = res->getPdbName().notnilp() ? res->getPdbName() : res->getName();
      atom._chainId = chainId;
      atom._resSeq = resSeq;
      atom._x = a->getPosition().getX();
      atom._y = a->getPosition().getY();
      atom._z = a->getPosition().getZ();
      atom._occupancy = 1.0;
      atom._tempFactor = 1.0;
      a->setTempInt(pdbAtoms.size());
      pdbAtoms.push_back(atom);
    }
    resSeq++;
  }
  for ( gctools::Vec0<AtomPdbRec>::iterator ai=pdbAtoms.begin(); ai!=pdbAtoms.end(); ai++ )
  {
    gctools::Vec0<Atom_sp>	bonded = ai->_atom->getBondedAtoms();
    gctools::Vec0<Atom_sp>	uniqueBonded;
    for ( gctools::Vec0<Atom_sp>::iterator bi=bonded.begin(); bi!=bonded.end(); bi++ ) 
    {
      if ( ai->_atom->getTempInt()<(*bi)->getTempInt() )
      {
        uniqueBonded.push_back(*bi);
      }
    }
    if ( uniqueBonded.size() > 0 )
    {
      ConnectPdbRec connect;
      connect._atom1 = ai->_atom->getTempInt()+1;
      for ( uint i=0; i<uniqueBonded.size(); i++ )
      {
        connect._bonded[i] = uniqueBonded[i]->getTempInt()+1;
      }
      for ( uint j=uniqueBonded.size(); j<4; j++ )
      {
        connect._bonded[j] = 0;
      }
      pdbConnects.push_back(connect);
    }
  }
  return serial;
}

void _writeAtomAndConnectRecords( core::T_sp fout, gctools::Vec0<AtomPdbRec>& pdbAtoms, 
                                  vector<ConnectPdbRec>& pdbConnects )
{
  gctools::Vec0<AtomPdbRec>::iterator ai;
  for ( ai=pdbAtoms.begin(); ai!=pdbAtoms.end(); ai++ )
  {
    ai->write(fout);
  }
#if 0
  vector<ConnectPdbRec>::iterator ci;
  for ( ci=pdbConnects.begin(); ci!=pdbConnects.end(); ci++ )
  {
    ci->write(fout);
  }
#endif
  core::clasp_write_string("TER\n", fout );
}

CL_LISPIFY_NAME("pdb-open");
CL_DEFMETHOD void PdbWriter_O::open(core::T_sp pathDesignator)
{
  this->_Out = cl__open(pathDesignator,
                        core::StreamDirection::output,
                        cl::_sym_character,
                        core::StreamIfExists::supersede, true,
                        core::StreamIfDoesNotExist::nil, false,
                        kw::_sym_default,
                        nil<T_O>());
  if ( this->_Out.nilp() ) {
    SIMPLE_ERROR("Could not open file: {}" , _rep_(pathDesignator));
  }
}

void	PdbWriter_O::write(Matter_sp matter)
{
  size_t serialStart = 1;
  if ( matter.isA<Aggregate_O>() )
  {
    Loop lMolecules;
    lMolecules.loopTopGoal(matter,MOLECULES);
    char chainId = 'A';
    while ( lMolecules.advance() )
    {
      gctools::Vec0<AtomPdbRec> pdbAtoms;
      vector<ConnectPdbRec> pdbConnects;
      Molecule_sp mol = lMolecules.getMolecule();
      serialStart = _setupAtomAndConnectRecordsForOneMolecule(mol,pdbAtoms,pdbConnects,chainId,serialStart);
      _writeAtomAndConnectRecords(this->_Out,pdbAtoms,pdbConnects);
      chainId++;
      if ( chainId > 'Z' ) chainId = 'A'; // recycle chain-ids if there are too many
    }
  } else
  {
    gctools::Vec0<AtomPdbRec> pdbAtoms;
    vector<ConnectPdbRec> pdbConnects;
    serialStart = _setupAtomAndConnectRecordsForOneMolecule(matter.as<Molecule_O>(),pdbAtoms,pdbConnects,'A',serialStart);
    _writeAtomAndConnectRecords(this->_Out,pdbAtoms,pdbConnects);
  }
}


CL_LISPIFY_NAME("writeModel");
CL_DEFMETHOD     void PdbWriter_O::writeModel(Matter_sp matter, int model)
{
  core::clasp_write_string(fmt::format("MODEL     {}\n" , model ), this->_Out);
  this->write(matter);
  core::clasp_write_string("ENDMDL\n",this->_Out);
}

CL_LISPIFY_NAME("pdb-close");
CL_DEFMETHOD     void PdbWriter_O::close()
{
  if ( this->_Out.notnilp() ) {
    core::clasp_write_string("END\n", this->_Out);
    core::cl__close(this->_Out);
  }
}

CL_LISPIFY_NAME(save-pdb);
DOCGROUP(cando);
CL_DEFUN void	PdbWriter_O::savePdb(Matter_sp matter, core::T_sp fileName )
{
  PdbWriter_sp writer = PdbWriter_O::create();
  writer->open(fileName);
  writer->write(matter);
  writer->close();
}











};
