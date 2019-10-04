/*
    File: ffTypesDb.cc
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
//
// (C) 2004 Christian E. Schafmeister
//

/*
 *	ffTypes.cc
 *
 *	Maintains databases and structures to store types and type assignement
 *	rules.
 */
#include <clasp/core/foundation.h>
#include <clasp/core/hashTable.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/record.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/topology.h>
#include <cando/chem/constitution.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/loop.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>
#include <clasp/core/bformat.h>





namespace chem {


void FFTypeRule_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,test), this->_Test);
  node->field(INTERN_(kw,type), this->_Type);
}




void	FFTypesDb_O::initialize()
{
    this->Base::initialize();
    this->_TypeAssignmentRules.clear();
}


SYMBOL_EXPORT_SC_(ChemPkg,assignType);

CL_LISPIFY_NAME("assignType");
CL_LAMBDA((types-db !) atom)
CL_DEFMETHOD core::Symbol_sp FFTypesDb_O::assignType(chem::Atom_sp atom) {
  LOG(BF("Got atom") );
  LOG(BF("atom name: %s") % atom->getName().c_str() );
  LOG(BF("Assigning type for atom: %s") % atom->description().c_str()  );
  if (chem__verbose(0)) {
    core::write_bf_stream(BF("Assigning type for atom: %s\n") % _rep_(atom));
  }

  {_BLOCK_TRACE("Testing every type rule");
    for ( auto it=this->_TypeAssignmentRules.begin();
          it!=this->_TypeAssignmentRules.end(); it++ ) 
    {_BLOCK_TRACEF(BF("Testing rule code(%s)") % (*it)->getCode().c_str() );
      Root_sp root = (*it)->_Test;
      core::T_mv matches_mv = chem::chem__chem_info_match(root,atom);
      ChemInfoMatch_sp match = gc::As<ChemInfoMatch_sp>(matches_mv.second());
      if ( matches_mv.notnilp() ) {
        LOG(BF("Rule MATCH!!!") );
        if (chem__verbose(2)) core::write_bf_stream(BF("Matched %s type-> %s\n") % _rep_(root) % _rep_((*it)->_Type));
        return (*it)->_Type;
      } else {
        if (chem__verbose(2)) core::write_bf_stream(BF("Did not match %s type-> %s\n") % _rep_(root) % _rep_((*it)->_Type));
      }
      LOG(BF("Rule does not match, keep going") );
    }
  }
  return _Nil<core::Symbol_O>();
}


CL_LISPIFY_NAME("assignTypes");
CL_DEFMETHOD void    FFTypesDb_O::assignTypes(chem::Matter_sp matter)
{ 
  chem::Loop    				lAtoms;
  chem::Atom_sp  				atom;
  chem::Matter_sp				c;
  chem::Molecule_sp                     	mol;
  chem::Residue_sp				res;
  string 		                        name;
  size_t missing_types = 0;
  size_t total_atoms = 0;
  // first clear out old atom types
  lAtoms.loopTopGoal(matter,ATOMS);
  while (lAtoms.advanceLoopAndProcess()) {
    atom = lAtoms.getAtom();
    atom->setType(_Nil<core::T_O>());
    total_atoms++;
  }
  // Now assign the types first checking for residues and then
  // applying atom type rules for missing types.
  if (chem__verbose(2)) {
    core::write_bf_stream(BF("Assigning atom types to %s\n") % _rep_(matter));
  }
  Loop lRes;
  lRes.loopTopGoal(matter,RESIDUES);
  while (lRes.advanceLoopAndProcess()) {
    Residue_sp res = lRes.getResidue();
    core::T_sp name = res->getName();
    if (chem__verbose(2)) {
      core::write_bf_stream(BF("Looking for residue: %s\n") % _rep_(name));
    }
    if (name.notnilp()) {
      core::T_mv result_mv = chem__findTopology(name,false);
      if (chem__verbose(2)) {
        core::write_bf_stream(BF("chem__findTopology -> %s %s\n") % _rep_(result_mv) % _rep_(result_mv.second()));
      }
      if (result_mv.second().notnilp()) {
        if (chem__verbose(2)) {
          core::write_bf_stream(BF("Found topology for residue name: %s\n") % _rep_(name));
        }
        Topology_sp topology = gc::As<Topology_sp>(result_mv);
        Constitution_sp constitution = topology->getConstitution();
        ConstitutionAtoms_sp constitution_atoms = constitution->getConstitutionAtoms();
        // Use the Topology to assign atom types
        lAtoms.loopTopGoal(res,ATOMS);
        while (lAtoms.advanceLoopAndProcess()) {
          atom = lAtoms.getAtom();
          core::T_sp atom_name = atom->getName();
          if (chem__verbose(2)) {
            core::write_bf_stream(BF("Looking for atom with name: %s\n") % _rep_(atom_name));
          }
          core::T_mv constitution_atom_mv = constitution_atoms->atomWithName(atom_name,false);
          if (constitution_atom_mv.second().notnilp()) {
            core::T_sp single = constitution_atom_mv;
            ConstitutionAtom_sp constitution_atom = gc::As_unsafe<ConstitutionAtom_sp>(single);
            core::T_sp type = constitution_atom->_AtomType;
            atom->setType(type);
            if (chem__verbose(2)) {
              core::write_bf_stream(BF("Assigned atom type %s using topology %s\n") % _rep_(type) % _rep_(name));
            }
          } else {
            if (chem__verbose(2)) {
              core::write_bf_stream(BF("   Not found\n"));
            }
          }
        }
      }
    }
    lAtoms.loopTopGoal(res,ATOMS);
    while (lAtoms.advanceLoopAndProcess()) {
      atom = lAtoms.getAtom();
      if (atom->getType().nilp()) {
        if (this->_TypeAssignmentRules.size()!=0) {
          core::T_sp type = this->assignType(atom);
          if (chem__verbose(2)) {
            core::write_bf_stream(BF("Assigned atom type %s using type rules\n") % _rep_(type));
          }
          atom->setType(type);
        }
      }
      if (atom->getType().nilp()) missing_types++;
    }
  }
  if (missing_types>0) {
    SIMPLE_WARN(BF("There were %lu missing types of %lu atoms") % missing_types % total_atoms);
  }
}


string FFTypeRule_O::__repr__() const {
  stringstream ss;
  ss << "#<FFTypeRule :type " << _rep_(this->_Type) << " :test " << _rep_(this->_Test);
  return ss.str();
}

CL_DOCSTRING("Return the atom types as a hash-table with the atoms as keys and types as values");
CL_DEFMETHOD core::HashTableEq_sp    FFTypesDb_O::atomTypes(chem::Matter_sp matter)
{ 
  chem::Loop    				lAtoms;
  chem::Atom_sp  				atom;
  core::HashTableEq_sp atomTypes = core::HashTableEq_O::create_default();
  if (this->_TypeAssignmentRules.size()==0)  return atomTypes;
  lAtoms.loopTopGoal(matter,ATOMS);
  LOG(BF("defined loop") );
  while ( lAtoms.advanceLoopAndProcess() ) {
    LOG(BF("Getting container") );
    atom = lAtoms.getAtom();
    core::Symbol_sp type = this->assignType(atom);
    atomTypes->setf_gethash(atom,type);
  }
  return atomTypes;
}


void FFTypesDb_O::fields(core::Record_sp node)
{	
  node->field(INTERN_(kw,type_rules), this->_TypeAssignmentRules);
  this->FFBaseDb_O::fields(node);
}



CL_LISPIFY_NAME("FFTypes-getRule");
CL_DEFMETHOD chem::FFTypeRule_sp FFTypesDb_O::getRule(uint index)
{_OF();
  if ( index < this->_TypeAssignmentRules.size() ) {
    return this->_TypeAssignmentRules[index];
  }
  SIMPLE_ERROR(BF("Rule index is out of bounds"));
}


CL_DEFUN FFTypesDb_sp chem__make_fftypes_db()
{
  return FFTypesDb_O::create();
}


};
