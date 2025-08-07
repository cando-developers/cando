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
#include <cctype>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/character.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/pathname.h>
#include <clasp/core/hashTable.h>
#include <clasp/core/hashTableEqual.h>
#include <clasp/core/string.h>
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
#include <cando/chem/pdb.h>
#include <cando/chem/write_mmcif.h>
#include <clasp/core/wrappers.h>

// ========================================================================================
// Implementation
// ========================================================================================
#include <fmt/format.h>

namespace chem
{

SYMBOL_EXPORT_SC_(ChemPkg,formula);
SYMBOL_EXPORT_SC_(ChemPkg,formula_weight);
SYMBOL_EXPORT_SC_(ChemPkg,formula_string);

/* -- open ---------------------------------------------------------------------- */
void MmcifWriter_O::open(core::T_sp pathDesignator)
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

/* -- close --------------------------------------------------------------------- */
void MmcifWriter_O::close()
{
  if ( this->_Out.notnilp() ) {
    core::clasp_write_string("#\n", this->_Out);
    core::cl__close(this->_Out);
    this->_Out = nil<core::T_O>();
  }
}

/* -- collectAtoms -------------------------------------------------------------- */
size_t MmcifWriter_O::collectAtoms(Matter_sp matter,
                                   std::vector<AtomCifRec>& atoms,
                                   std::map<std::string,std::string>& resNameMap,
                                   char chainIdStart)
{
  // Build unique residue name set first – exactly like PdbWriter does
  std::set<std::string> longNames;
  Loop resLoop; resLoop.loopTopGoal(matter, RESIDUES);
  while (resLoop.advance()) {
    Residue_sp res = resLoop.getResidue();
    longNames.insert(res->getName()->symbolNameAsString());
  }
  // Use PdbWriter helper to map -> three‑letter codes
  std::vector<std::string> longNamesVec(longNames.begin(), longNames.end());
  std::vector<std::string> pdbNames   = make_pdb_tags(longNamesVec);
  for (size_t i=0;i<longNamesVec.size();++i) resNameMap[longNamesVec[i]] = pdbNames[i];

  // Now collect atoms
  char chainId = chainIdStart;
  Loop molLoop; molLoop.loopTopGoal(matter, MOLECULES);
  while ( molLoop.advance() )
    {
      Molecule_sp mol = molLoop.getMolecule();
      Loop resLoop; resLoop.loopTopGoal(mol, RESIDUES);
      int resSeq = 1;
      while (resLoop.advance())
        {
          Residue_sp res = resLoop.getResidue();
          Loop atLoop; atLoop.loopTopGoal(res, ATOMS);
          while (atLoop.advance())
            {
              Atom_sp a = atLoop.getAtom();
              AtomCifRec rec;
              rec._atom      = a;
              rec._name      = a->getName();
              rec._element   = a->getElementAsString();
              rec._resName   = res->getName();
              rec._pdbResName = resNameMap[ res->getName()->symbolNameAsString() ];
              rec._chainId   = chainId;
              rec._resSeq    = resSeq;
              rec._x         = a->getPosition().getX();
              rec._y         = a->getPosition().getY();
              rec._z         = a->getPosition().getZ();
              rec._occupancy = 1.0;
              rec._biso      = 1.0;
              atoms.push_back(rec);
            }
          ++resSeq;
        }
      ++chainId;
    }
  // Assign tempInt serial numbers for bonding reference
  for (size_t i=0;i<atoms.size();++i) atoms[i]._atom->setTempInt(i+1);
  return atoms.size();
}

/* -- collectResidues -------------------------------------------------------------- */
core::HashTableEq_sp MmcifWriter_O::collectResidues(Matter_sp matter )
{
  // Build unique residue name set first – exactly like PdbWriter does
  core::HashTableEq_sp residues = core::HashTableEq_O::create_default();
  Loop resLoop;
  resLoop.loopTopGoal(matter, RESIDUES);
  while (resLoop.advance()) {
    Residue_sp res = resLoop.getResidue();
    if (residues->gethash(res->getName()).nilp()) {
      residues->setf_gethash(res->getName(), res );
    }
  }
  return residues;
}

/* -- writeAllResidues ----  */
void MmcifWriter_O::writeAllResidues(core::HashTableEq_sp residues) {
  core::clasp_write_string(fmt::format("loop_\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp.id\n"), this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp.name\n"), this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp.type\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_type\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.formula\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.formula_weight\n"), this->_Out );
  residues->maphash( [this] (core::T_sp name, core::T_sp tresidue ) {
    Residue_sp residue = gc::As<Residue_sp>(tresidue);
    core::T_sp tformula = core::eval::funcall( chem::_sym_formula, residue);
    core::T_sp tformulaString = core::eval::funcall( chem::_sym_formula_string, tformula );
    core::T_sp tformulaWeight = core::eval::funcall( chem::_sym_formula_weight, tformula );
    double formulaWeight = core::clasp_to_double(tformulaWeight);
    core::clasp_write_string(fmt::format("{} \"{}\" {} {} {} {:.3f}\n",
                                         residue->getName()->symbolNameAsString(),
                                         residue->getName()->symbolNameAsString(),
                                         "other",
                                         "HETAIN",
                                         _rep_(tformulaString),
                                         formulaWeight),
                             this->_Out);
  } );
  core::clasp_write_string(fmt::format("#\n"), this->_Out );
  core::clasp_write_string(fmt::format("loop_\n"),this->_Out);
  
  core::clasp_write_string(fmt::format("_chem_comp_atom.comp_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.atom_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.alt_atom_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.type_symbol\n"),this->_Out);
  // core::clasp_write_string(fmt::format("_chem_comp_atom.charge\n"),this->_Out);
  // core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_align\n"),this->_Out);
  //  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_aromatic_flag\n"),this->_Out);
  //  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_leaving_atom_flag\n"),this->_Out);
  //  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_stereo_config\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.model_Cartn_x\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.model_Cartn_y\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.model_Cartn_z\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_model_Cartn_x_ideal\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_model_Cartn_y_ideal\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_model_Cartn_z_ideal\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_ordinal\n"),this->_Out);
    
  residues->maphash( [this] (core::T_sp name, core::T_sp tresidue ) {
    Residue_sp residue = gc::As<Residue_sp>(tresidue);
    Loop lAtoms;
    lAtoms.loopTopGoal(residue, ATOMS);
    size_t ordinal = 0;
    while (lAtoms.advance()) {
      Atom_sp atm = lAtoms.getAtom();
      core::clasp_write_string(fmt::format("{} {} {} {} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:d}\n",
                                           residue->getName()->symbolNameAsString(),
                                           atm->getName()->symbolNameAsString(),
                                           atm->getName()->symbolNameAsString(),
                                           atm->getElementAsString(),
                                           atm->getPosition().getX(),
                                           atm->getPosition().getY(),
                                           atm->getPosition().getZ(),
                                           atm->getPosition().getX(),
                                           atm->getPosition().getY(),
                                           atm->getPosition().getZ(),
                                           ++ordinal),
                               this->_Out);
    }
  } );
  core::clasp_write_string(fmt::format("#\n"),this->_Out);
  core::clasp_write_string(fmt::format("loop_\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.comp_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.atom_id_1\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.atom_id_2\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.value_order\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.pdbx_ordinal\n"),this->_Out);
  size_t bordinal = 0;
  residues->maphash( [&bordinal,this] (core::T_sp name, core::T_sp tresidue ) {
    Residue_sp residue = gc::As<Residue_sp>(tresidue);
    core::HashTable_sp atomToResidue = residue->atomToResidueMap();
    Loop lBonds;
    lBonds.loopTopGoal(residue, BONDS);
    while (lBonds.advance()) {
      Atom_sp atm1 = lBonds.getAtom1();
      Atom_sp atm2 = lBonds.getAtom2();
      if (atomToResidue->gethash(atm1)==residue && atomToResidue->gethash(atm2)==residue) {
        BondOrder bo = lBonds.getBondOrder();
        string order;
        switch (bo) {
        case singleBond: order = "SING";
          break;
        case doubleBond: order = "DOUB";
          break;
        case tripleBond: order = "TRIP";
          break;
        default: order = "?";
        }
        core::clasp_write_string(fmt::format("{} {} {} {} {:d}\n",
                                             residue->getName()->symbolNameAsString(),
                                             atm1->getName()->symbolNameAsString(),
                                             atm2->getName()->symbolNameAsString(),
                                             order,
                                             ++bordinal),
                                 this->_Out);
      }
    }
  } );
  core::clasp_write_string(fmt::format("#\n"),this->_Out);
}

/* -- writeResidue ----  */
void MmcifWriter_O::writeResidue(Residue_sp residue) {
  //  core::clasp_write_string(fmt::format("data_{}\n#\n",residue->getName()->symbolNameAsString()), this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp.id {}\n",residue->getName()->symbolNameAsString()), this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp.name \"{}\"\n",residue->getName()->symbolNameAsString()), this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp.type other\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_type HETAIN\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.formula ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.mon_nstd_parent_comp_id ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_synonyms ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_formal_charge ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_initial_date ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_modified_date ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_ambiguous_flag ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_release_status ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_replaced_by ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_replaces ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.formula_weight ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.one_letter_code ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.three_letter_code ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_model_coordinates_details ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_model_coordinates_missing_flag N\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_ideal_coordinates_details ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_ideal_coordinates_missing_flag N\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_model_coordinates_db_code ?\n"), this->_Out );
  core::clasp_write_string(fmt::format("_chem_comp.pdbx_processing_site RCSB\n"), this->_Out );

  core::clasp_write_string(fmt::format("#\n"),this->_Out);
  core::clasp_write_string(fmt::format("loop_\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.comp_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.atom_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.alt_atom_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.type_symbol\n"),this->_Out);
  // core::clasp_write_string(fmt::format("_chem_comp_atom.charge\n"),this->_Out);
  // core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_align\n"),this->_Out);
  //  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_aromatic_flag\n"),this->_Out);
  //  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_leaving_atom_flag\n"),this->_Out);
  //  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_stereo_config\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.model_Cartn_x\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.model_Cartn_y\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.model_Cartn_z\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_model_Cartn_x_ideal\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_model_Cartn_y_ideal\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_model_Cartn_z_ideal\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_atom.pdbx_ordinal\n"),this->_Out);
  Loop lAtoms;
  lAtoms.loopTopGoal(residue, ATOMS);
  size_t ordinal = 0;
  while (lAtoms.advance()) {
    Atom_sp atm = lAtoms.getAtom();
    core::clasp_write_string(fmt::format("{} {} {} {} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:d}\n",
                                         residue->getName()->symbolNameAsString(),
                                         atm->getName()->symbolNameAsString(),
                                         atm->getName()->symbolNameAsString(),
                                         atm->getElementAsString(),
                                         atm->getPosition().getX(),
                                         atm->getPosition().getY(),
                                         atm->getPosition().getZ(),
                                         atm->getPosition().getX(),
                                         atm->getPosition().getY(),
                                         atm->getPosition().getZ(),
                                         ++ordinal),
                             this->_Out);
  }
  core::clasp_write_string(fmt::format("#\n"),this->_Out);
  core::clasp_write_string(fmt::format("loop_\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.comp_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.atom_id_1\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.atom_id_2\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.value_order\n"),this->_Out);
  core::clasp_write_string(fmt::format("_chem_comp_bond.pdbx_ordinal\n"),this->_Out);
  size_t bordinal = 0;
  core::HashTable_sp atomToResidue = residue->atomToResidueMap();
  Loop lBonds;
  lBonds.loopTopGoal(residue, BONDS);
  while (lBonds.advance()) {
    Atom_sp atm1 = lBonds.getAtom1();
    Atom_sp atm2 = lBonds.getAtom2();
    if (atomToResidue->gethash(atm1)==residue && atomToResidue->gethash(atm2)==residue) {
      BondOrder bo = lBonds.getBondOrder();
      string order;
      switch (bo) {
      case singleBond: order = "SING";
        break;
      case doubleBond: order = "DOUB";
        break;
      case tripleBond: order = "TRIP";
        break;
      default: order = "?";
      }
      core::clasp_write_string(fmt::format("{} {} {} {} {:d}\n",
                                           residue->getName()->symbolNameAsString(),
                                           atm1->getName()->symbolNameAsString(),
                                           atm2->getName()->symbolNameAsString(),
                                           order,
                                           ++bordinal),
                               this->_Out);
    }
  }
  core::clasp_write_string(fmt::format("#\n"),this->_Out);
}

/* -- write --------------------------------------------------------------------- */
core::HashTableEqual_sp
MmcifWriter_O::write(Matter_sp matter, core::List_sp metadata) {
  // ---------------------------------------------------------------- cell/symmetry
  std::string dataName("model1");
  core::clasp_write_string(fmt::format("data_{}\n#\n",dataName), this->_Out);

  // ---------------------------------------------------------------- collect residues
  core::HashTableEq_sp residues = collectResidues(matter);
  this->writeAllResidues(residues);
  core::HashTableEq_sp atomToResidue = core::HashTableEq_O::create_default();
  core::HashTableEq_sp atomToIndex = core::HashTableEq_O::create_default();
  core::HashTableEq_sp residueToIndex = core::HashTableEq_O::create_default();
  core::HashTableEq_sp residueToMolecule = core::HashTableEq_O::create_default();
  core::HashTableEq_sp moleculeToChainId = core::HashTableEq_O::create_default();
  Loop lMolecules;
  lMolecules.loopTopGoal(matter,MOLECULES);
  char chainId = 'A';
  while ( lMolecules.advance() ) {
    Molecule_sp mol = lMolecules.getMolecule();
    moleculeToChainId->setf_gethash(mol,core::clasp_make_standard_character(chainId++));
    size_t atomIndex = 0;
    size_t residueIndex = 0;
    Loop lResidues;
    lResidues.loopTopGoal(mol,RESIDUES);
    while (lResidues.advance()) {
      Residue_sp res = lResidues.getResidue();
      residueToMolecule->setf_gethash(res,mol);
      residueToIndex->setf_gethash(res,core::make_fixnum(++residueIndex));
      Loop lAtoms;
      lAtoms.loopTopGoal(res,ATOMS);
      while (lAtoms.advance()) {
        Atom_sp atm = lAtoms.getAtom();
        atomToResidue->setf_gethash(atm,res);
        atomToIndex->setf_gethash(atm,core::make_fixnum(atomIndex++));
      }
    }
  }
  if (metadata.consp()) {
    core::T_sp val = gc::As_unsafe<core::Cons_sp>(metadata)->getf(kw::_sym_cryst1,nil<core::T_O>());
    if (val.consp()) {
      core::Cons_sp cval = gc::As_unsafe<core::Cons_sp>(val);
      double a = core::clasp_to_double(cval->elt(0));
      double b = core::clasp_to_double(cval->elt(1));
      double c = core::clasp_to_double(cval->elt(2));
      double alpha = core::clasp_to_double(cval->elt(3));
      double beta  = core::clasp_to_double(cval->elt(4));
      double gamma = core::clasp_to_double(cval->elt(5));
      std::string spgr = gc::As<core::String_sp>(cval->elt(6))->get_std_string();
      size_t z_pdb;
      if (cval->elt(7).fixnump()) {
        z_pdb = cval->elt(7).unsafe_fixnum();
      } else {
        SIMPLE_ERROR("You must provide the number of asymmetric units as a fixnum");
      }
      core::clasp_write_string(fmt::format("_cell.length_a    {:.3f}\n",a), this->_Out);
      core::clasp_write_string(fmt::format("_cell.length_b    {:.3f}\n",b), this->_Out);
      core::clasp_write_string(fmt::format("_cell.length_c    {:.3f}\n",c), this->_Out);
      core::clasp_write_string(fmt::format("_cell.angle_alpha {:.2f}\n",alpha), this->_Out);
      core::clasp_write_string(fmt::format("_cell.angle_beta  {:.2f}\n",beta), this->_Out);
      core::clasp_write_string(fmt::format("_cell.angle_gamma {:.2f}\n",gamma), this->_Out);
      core::clasp_write_string(fmt::format("_symmetry.space_group_name_H-M   '{} '\n",spgr), this->_Out);
      core::clasp_write_string(fmt::format("_cell.Z_PDB {}\n#\n", z_pdb), this->_Out);
    }
  }

  // ---------------------------------------------------------------- collect atoms
  std::vector<AtomCifRec> atoms;
  std::map<std::string,std::string> resNameMap;
  collectAtoms(matter, atoms, resNameMap);

  // ---------------------------------------------------------------- atom_site loop header
  core::clasp_write_string("loop_\n", this->_Out);
  core::clasp_write_string("_atom_site.group_PDB\n", this->_Out);
  core::clasp_write_string("_atom_site.id\n", this->_Out);
  core::clasp_write_string("_atom_site.type_symbol\n", this->_Out);
  core::clasp_write_string("_atom_site.label_atom_id\n", this->_Out);
  core::clasp_write_string("_atom_site.label_comp_id\n", this->_Out);
  core::clasp_write_string("_atom_site.label_asym_id\n", this->_Out);
  core::clasp_write_string("_atom_site.label_seq_id\n", this->_Out);
  core::clasp_write_string("_atom_site.Cartn_x\n", this->_Out);
  core::clasp_write_string("_atom_site.Cartn_y\n", this->_Out);
  core::clasp_write_string("_atom_site.Cartn_z\n", this->_Out);
  core::clasp_write_string("_atom_site.occupancy\n", this->_Out);
  core::clasp_write_string("_atom_site.B_iso_or_equiv\n", this->_Out);
  core::clasp_write_string("#\n", this->_Out);

  // ---------------------------------------------------------------- atom_site rows
  size_t serial = 1;
  for (const auto& rec : atoms) {
    core::clasp_write_string(
        fmt::format("{} {} {} {} {} {} {} {:.3f} {:.3f} {:.3f} {:.2f} {:.2f}\n",
                    "ATOM",
                    serial,
                    rec._element,
                    rec._name->symbolNameAsString(),
                    rec._resName->symbolNameAsString(),
                    rec._chainId,
                    rec._resSeq,
                    rec._x, rec._y, rec._z,
                    rec._occupancy,
                    rec._biso),
        this->_Out);
    ++serial;
  }
  core::clasp_write_string("#\n", this->_Out);


  core::clasp_write_string(fmt::format("loop_\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.conn_type_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.ptnr1_label_asym_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.ptnr1_label_seq_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.ptnr1_label_comp_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.ptnr1_label_atom_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.ptnr2_label_asym_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.ptnr2_label_seq_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.ptnr2_label_comp_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.ptnr2_label_atom_id\n"),this->_Out);
  core::clasp_write_string(fmt::format("_struct_conn.pdbx_value_order\n"),this->_Out);
  Loop lBonds;
  lBonds.loopTopGoal(matter,BONDS);
  size_t lnk = 0;
  while (lBonds.advance()) {
    Atom_sp a1 = lBonds.getAtom1();
    Atom_sp a2 = lBonds.getAtom2();
    if (atomToResidue->gethash(a1)!=atomToResidue->gethash(a2)) {
      BondOrder bo = lBonds.getBondOrder();
      string order;
      switch (bo) {
      case singleBond: order = "SING";
        break;
      case doubleBond: order = "DOUB";
        break;
      case tripleBond: order = "TRIP";
        break;
      default: order = "?";
      }
      Residue_sp res1 = gc::As<Residue_sp>(atomToResidue->gethash(a1));
      Residue_sp res2 = gc::As<Residue_sp>(atomToResidue->gethash(a2));
      Molecule_sp mol1 = gc::As<Molecule_sp>(residueToMolecule->gethash(res1));
      Molecule_sp mol2 = gc::As<Molecule_sp>(residueToMolecule->gethash(res2));
      core::Character_sp cid1 = (core::Character_sp)moleculeToChainId->gethash(mol1);
      char cchainId1 = unbox_character(cid1);
      core::Character_sp cid2 = (core::Character_sp)moleculeToChainId->gethash(mol2);
      char cchainId2 = unbox_character(cid2);
      size_t seqId1 = residueToIndex->gethash(res1).unsafe_fixnum();
      size_t seqId2 = residueToIndex->gethash(res2).unsafe_fixnum();
      core::clasp_write_string(fmt::format("LNK{:d} {} {} {} {} {} {} {} {} {} {}\n",
                                           ++lnk,
                                           "covale",
                                           cchainId1,
                                           seqId1,
                                           res1->getName()->symbolNameAsString(),
                                           a1->getName()->symbolNameAsString(),
                                           cchainId2,
                                           seqId2,
                                           res2->getName()->symbolNameAsString(),
                                           a2->getName()->symbolNameAsString(),
                                           order),
                               this->_Out);
      //       LNK1 covale A 1 CGLY N   A 2 PRO4SS CAC SING
    }
  }

  // Build residue map identical to PdbWriter
  core::HashTableEqual_sp mapping = core::HashTableEqual_O::create_default();
  for (auto& kv : resNameMap) {
    mapping->setf_gethash(core::SimpleBaseString_O::make(kv.first),
                          core::SimpleBaseString_O::make(kv.second));
  }
  return mapping;
}

/* -- writeModel ---------------------------------------------------------------- */
void MmcifWriter_O::writeModel(Matter_sp matter,int model,core::List_sp metadata)
{
  core::clasp_write_string(fmt::format("data_model{}\n",model), this->_Out);
  this->write(matter,metadata);
}

} // namespace chem

namespace chem
{

CL_LISPIFY_NAME(save-mmcif);
DOCGROUP(cando);
CL_LAMBDA(matter filename &key metadata);
CL_DEFUN core::HashTableEqual_sp MmcifWriter_O::saveMmcif(Matter_sp matter, core::T_sp fileName, core::List_sp metadata=nil<core::T_O>() )
{
  MmcifWriter_sp writer = MmcifWriter_O::create();
  writer->open(fileName);
  core::HashTableEqual_sp pdb_residue_name_map = writer->write(matter,metadata);
  writer->close();
  return pdb_residue_name_map;
}


};
