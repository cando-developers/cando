
/*
    File: candoScript.cc
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

#include <stdlib.h>
#include <boost/filesystem.hpp>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/ql.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/symbolTable.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/fileSystem.h>
//#include "core/xmlLoadArchive.h"
//#include "core/xmlSaveArchive.h"
#include <cando/adapt/stringList.h>
#include <clasp/core/designators.h>
#include <clasp/core/bundle.h>
#include <clasp/core/evaluator.h>
#include <cando/geom/coordinateArray.h>
#include <cando/chem/candoScript.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/monomer.h>
#include <cando/chem/oligomer.h>
#include <cando/chem/mol2.h>
#include <cando/chem/loop.h>
//#include "monomerCoordinates.h"
//#include "groupPart.h"
//#include "fragmentCoordinates.h"
#include <cando/chem/calculatePosition.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitution.h>
#include <cando/chem/chemPackage.h>
#include <clasp/core/wrappers.h>

#ifdef	READLINE
extern "C" char *readline( const char* prompt);
extern "C" void add_history(char* line);
#endif

namespace chem {


/*! Everything that depends on SMARTS parsing is initialized here */

DOCGROUP(cando);
CL_DEFUN void chem__initialize_smarts_users() {
//  printf("%s:%d:%s\n", __FILE__, __LINE__, __FUNCTION__);
  energyFunction_initializeSmarts();
}




/*!
 * Replace the alchemist database of this Lisp environment.
 * If there already was one then tell it to give all of its
 * dependent objects to the new CandoDatabase
 */
void setCandoDatabase(CandoDatabase_sp bdb, core::LispPtr lisp)
{
    LOG("Setting *DATABASE* to database@{}" , bdb.get() );
    _sym_candoDatabase->defparameter(bdb);
}

/*
__BEGIN_DOC( candoScript.general.loadMol2, subsection, loadMol2 )
\scriptCmdRet{loadMol2}{Text::fileName}{Matter::result}

Loads the Tripos ``mol2'' file and returns the Aggregate within it.
__END_DOC
 */

 DOCGROUP(cando);
CL_DEFUN geom::SimpleVectorCoordinate_sp chem__make_simple_vector_coordinate_from_atom_sequence(core::T_sp atoms)
{
  geom::SimpleVectorCoordinate_sp coords;
  size_t num = core::cl__length(atoms);
  if (gc::IsA<core::List_sp>(atoms)) {
    core::List_sp latoms = gc::As_unsafe<core::List_sp>(atoms);
    coords = geom::SimpleVectorCoordinate_O::make(num);
    size_t idx(0);
    for ( auto cur : latoms ) {
      Atom_sp a = gc::As<chem::Atom_sp>(oCar(cur));
      Vector3 v = a->getPosition();
      (*coords)[idx] = v;
      ++idx;
    }
  } else if (gc::IsA<core::Array_sp>(atoms)) {
    core::Array_sp vatoms = gc::As_unsafe<core::Array_sp>(atoms);
    coords = geom::SimpleVectorCoordinate_O::make(num);
    size_t idx(0);
    for ( size_t ii = 0; ii<num; ii++ ) {
      Atom_sp a = gc::As<chem::Atom_sp>(vatoms->rowMajorAref(ii));
      Vector3 v = a->getPosition();
      (*coords)[idx] = v;
      ++idx;
    }
  } else {
    TYPE_ERROR(atoms,core::Cons_O::createList(cl::_sym_list,cl::_sym_vector));
  }
  return coords;
}

CL_DOCSTRING(R"dx(Load the only or the first aggregate from the mol2 file.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_mv chem__read_mol2(core::T_sp stream)
{
  Mol2File fin(stream);
  return mol2Read(fin);
}

CL_DOCSTRING(R"dx(Read the only or the first aggregate from the mol2 file.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_mv chem__load_mol2(core::T_sp fileName)
{
  Mol2File fin;
  fin.openFileName(fileName);
  return mol2Read(fin);
}


core::T_sp chem__read_mol2_list_common(Mol2File& fin, core::T_sp number_to_load)
{
  core::T_sp make_progress;
  core::T_sp progress_advance;
  core::T_sp progress_done;
  if (chem__verbose(0)) {
    make_progress = core::eval::funcall(cl::_sym_findSymbol,core::SimpleBaseString_O::make("MAKE-PROGRESS-BAR"),
                                        core::SimpleBaseString_O::make("CANDO"));
    progress_advance = core::eval::funcall(cl::_sym_findSymbol,core::SimpleBaseString_O::make("PROGRESS-ADVANCE"),
                                           core::SimpleBaseString_O::make("CANDO"));
    progress_done = core::eval::funcall(cl::_sym_findSymbol,core::SimpleBaseString_O::make("PROGRESS-DONE"),
                                        core::SimpleBaseString_O::make("CANDO"));
    if (make_progress.nilp()||progress_advance.nilp()||progress_done.nilp()) {
      SIMPLE_ERROR("Could not get progress bar functions make-progress, progress-advance, progress-done");
    }
  }
  core::T_sp progress_bar = nil<core::T_O>();
  if (chem__verbose(0)) {
    if (number_to_load.nilp()) {
      progress_bar = core::eval::funcall(make_progress,
                                         INTERN_(kw,total), core::stream_length(fin.fIn));
    } else if (gc::IsA<core::Real_sp>(number_to_load)){
      progress_bar = core::eval::funcall(make_progress,
                                         INTERN_(kw,total), number_to_load);
    }
  }
  ql::list result;
  size_t count = ~0;
  if (number_to_load.fixnump()) {
    count = number_to_load.unsafe_fixnum();
  }
  size_t num = 0;
  while (!fin.eof()) {
    if (count) {
      core::T_mv res = mol2Read(fin);
      core::T_sp tagg = res;
      core::MultipleValues &values = core::lisp_multipleValues();
      values.second(res.number_of_values());
      core::T_sp atom_types = values.second(res.number_of_values());
      if (gc::IsA<Aggregate_sp>(tagg)) {
        result << core::Cons_O::createList(tagg,atom_types);
        if (chem__verbose(0)) {
          if (number_to_load.nilp()) {
            core::eval::funcall(progress_advance, progress_bar, core::stream_position(fin.fIn));
          } else {
            core::eval::funcall(progress_advance, progress_bar, core::make_fixnum((size_t)num));
          }
        }
      }
    } else {
      break;
    }
    ++num;
    --count;
  }
  if (chem__verbose(0)) {
    core::eval::funcall(progress_done,progress_bar);
  }
  return result.cons();
}

CL_DOCSTRING(R"dx(Read all or a number of aggregates from the mol2 stream.)dx");
CL_LAMBDA(stream &optional number-to-load);
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__read_mol2_list(core::T_sp stream, core::T_sp number_to_load)
{
  Mol2File fin(stream);
  return chem__read_mol2_list_common(fin, number_to_load);
}

CL_DOCSTRING(R"dx(Load all or a number of aggregates from the mol2 file.)dx");
CL_LAMBDA(file-name &optional number-to-load);
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__load_mol2_list(core::T_sp fileName, core::T_sp number_to_load)
{
  Mol2File fin;
  fin.openFileName(fileName);
  return chem__read_mol2_list_common(fin, number_to_load);
}

/*
__BEGIN_DOC( candoScript.general.saveMol2, subsection, saveMol2)
\scriptCmd{saveMol2}{ Matter::structure Text::fileName}

Saves \sa{structure} to the Tripos "mol2" file with the name: \sa{fileName}.
__END_DOC
 */


CL_LAMBDA(matter dest-desig &optional use-sybyl-types atom-types);
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__save_mol2(Matter_sp matter, core::T_sp destDesig, bool useSybylTypes, core::T_sp atom_types)
{
  core::HashTable_sp ht;
  if (useSybylTypes) {
    ht = assignSybylTypes(matter);
  } else {
    ht = gc::As<core::HashTable_sp>(atom_types);
  }
  mol2WriteMatterToFileName(matter,destDesig,ht);
  return nil<core::T_O>();
}


CL_LAMBDA(matter &optional use-sybyl-types atom-types);
DOCGROUP(cando);
CL_DEFUN std::string chem__aggregate_as_mol2_string(Aggregate_sp matter, bool useSybylTypes, core::T_sp atom_types)
{
  core::HashTable_sp ht;
  if (useSybylTypes) {
    ht = assignSybylTypes(matter);
  } else {
    ht = gc::As<core::HashTable_sp>(atom_types);
  }
  stringstream ss;
  mol2WriteAggregateStream(matter,ss,ht);
  return ss.str();
}

CL_LAMBDA(matter &optional use-sybyl-types atom-types);
DOCGROUP(cando);
CL_DEFUN std::string chem__matter_as_mol2_string(Matter_sp matter, bool useSybylTypes, core::T_sp atom_types)
{
  core::HashTable_sp ht;
  if (useSybylTypes) {
    ht = assignSybylTypes(matter);
  } else {
    ht = gc::As<core::HashTable_sp>(atom_types);
  }
  stringstream ss;
  if (gc::IsA<Aggregate_sp>(matter)) {
    Aggregate_sp agg = gc::As_unsafe<Aggregate_sp>(matter);
    mol2WriteAggregateStream(agg,ss,ht);
    return ss.str();
  } else if (gc::IsA<Molecule_sp>(matter)) {
    Aggregate_sp agg = Aggregate_O::create();
    agg->addMatter(gc::As_unsafe<Molecule_sp>(matter));
    mol2WriteAggregateStream(agg,ss,ht);
    return ss.str();
  }
  SIMPLE_ERROR("Matter must be an aggregate or molecule - it is {}" , _rep_(matter));
}

/*!
 * Look for a residue using an identifier,
 * this looks for residues with a sequence number if the identifier is an int
 * or a residue name or a residue alias
 */

Residue_sp findResidue(Matter_sp matter, core::T_sp residueIdentifier )
{
  Molecule_sp molecule;
  core::Fixnum_sp	residueSequenceNumber;

	    // downcast the identifer to a Symbol and an Int object
	    // one of them will be nil and the other will have a value
  core::Symbol_sp resIdName = residueIdentifier.asOrNull<core::Symbol_O>();
  core::Fixnum_sp resIdSeqNum = residueIdentifier.asOrNull<core::Fixnum_O>();
  Loop l;
  l.loopTopGoal(matter,RESIDUES);
  while ( l.advanceLoopAndProcess() )
  {
    Residue_sp res = l.getResidue();
    LOG("Looking at residue with sequence number: {}" , res->getFileSequenceNumber()  );
    bool foundResidue = false;
    if ( resIdName )
    {
      LOG("Checking if residue has name({}) that matches({})" , res->getName().c_str() , _rep_(resIdName) );
      if ( resIdName==res->getName() )
      {
        LOG("Found residue with sequence name: {}" , _rep_(resIdName) );
        foundResidue = true;
      }
#ifdef USE_TOPOLOGY
      else if ( res->recognizesMonomerAlias(resIdName) )
      {
        LOG("Found residue with Monomer alias: {}" , _rep_(resIdName) );
        foundResidue = true;
      }
#endif
    } else if ( resIdSeqNum )
    {
      LOG("Checking if residue has fileSequenceNumber({}) that matches({})" , res->getFileSequenceNumber() , resIdSeqNum->get()  );
      if ( (int)(res->getFileSequenceNumber()) == resIdSeqNum.unsafe_fixnum())
      {
        LOG("Found residue with sequence number: {}" , resIdSeqNum->get()  );
        foundResidue = true;
      }
    }
    if ( foundResidue )
    {
      return res;
    }
  }
  return nil<Residue_O>();
}

/*
__BEGIN_DOC( candoScript.general.findResidue, subsection, findResidue)
\scriptCmdRet{findResidue}{ structure:Molecule residueIdentifier}{Residue}\par
\scriptCmdRet{findResidue}{ structure:Aggregate chain:text residueIdentifier}{Residue}

Return the Residue in the Aggregate or Molecule: \scriptArg{structure}. The residueIdentifier can either be a sequence number of a residue name.
__END_DOC
 */




#define ARGS_chem__find_residue "(&rest args)"
#define DECL_chem__find_residue ""
#define DOCS_chem__find_residue "findResidue"
CL_LAMBDA(&rest args);
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__find_residue(core::List_sp args)
{
  Molecule_sp molecule;
  core::Fixnum_sp	residueSequenceNumber;
  core::T_sp residueIdentifier;
  if ( core::cl__length(args)==3 ) 
  {
    Aggregate_sp agg = args.asCons()->onth(0).as<Aggregate_O>();
    core::String_sp chain = args.asCons()->onth(1).as<core::String_O>();
    residueIdentifier = args.asCons()->onth(2).as<core::T_O>();
    molecule = gc::As<Molecule_sp>(agg->contentWithName(chemkw_intern(chain->get_std_string())));
  } else if ( core::cl__length(args)==2 ) 
  {
    molecule = args.asCons()->onth(0).as<Molecule_O>();
    residueIdentifier = args.asCons()->onth(1).as<core::T_O>();
  } else
  {
    SIMPLE_ERROR("You must provide a molecule, residueId");
  }
  Residue_sp res = findResidue(molecule,residueIdentifier);
  if ( res.nilp() )
  {
    SIMPLE_ERROR("Molecule does not contain residue with identifier: {}", residueIdentifier);
  }
  return res;
}






/*
__BEGIN_DOC( candoScript.general.atomPos, subsection, atomPos)
\scriptCmdRet{atomPos}{ structure:Molecule residueNumber::int atomName:text }{vector}\par
\scriptCmdRet{atomPos}{ structure:Aggregate chain:text residueNumber::int atomName:text}{vector}

Return the position of the atom in the Aggregate or Molecule: \scriptArg{structure}. Identify the atom with \scriptArg{residueNumber} and \sa{atomName} and optionally the \sa{molecule/chain}. This function returns the position as a three element list of numbers.
__END_DOC
 */



#define ARGS_chem__atom_pos "(&rest args)"
#define DECL_chem__atom_pos ""
#define DOCS_chem__atom_pos "atomPos"
CL_LAMBDA(&rest args);
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__atom_pos(core::List_sp args)
{
  Molecule_sp molecule;
  core::Fixnum_sp	residueSequenceNumber;
  core::Symbol_sp	atomName;
  core::T_sp residueIdentifier;
  if ( core::cl__length(args)==4 ) 
  {
    Aggregate_sp agg = args.asCons()->onth(0).as<Aggregate_O>();
    core::String_sp chain = args.asCons()->onth(1).as<core::String_O>();
    residueIdentifier = args.asCons()->onth(2).as<core::T_O>();
    atomName = args.asCons()->onth(3).as<core::Symbol_O>();
    molecule = gc::As<Molecule_sp>(agg->contentWithName(chemkw_intern(chain->get_std_string())));
  } else if ( core::cl__length(args)==3 ) 
  {
    molecule = args.asCons()->onth(0).as<Molecule_O>();
    residueIdentifier = args.asCons()->onth(1).as<core::T_O>();
    atomName = args.asCons()->onth(2).as<core::Symbol_O>();
  } else
  {
    SIMPLE_ERROR("You must provide a molecule, residueId and atom name");
  }
  Residue_sp foundResidue = findResidue(molecule,residueIdentifier);
  if ( foundResidue.notnilp() )
  {
    if ( foundResidue->hasAtomWithName(atomName)) {
      LOG("Found atom with name: {}" , _rep_(atomName)  );
      Vector3 pos = gc::As_unsafe<Atom_sp>(foundResidue->atomWithName(atomName))->getPosition();
      geom::OVector3_sp v = geom::OVector3_O::create();
      v->setAll3(pos.getX(),pos.getY(),pos.getZ());
      return v;
    }
    SIMPLE_ERROR("Residue does not contain atom named: {}" , _rep_(atomName));
  }
  SIMPLE_ERROR("Molecule does not contain residue with identifier: {}" , _rep_(residueIdentifier));
}



DOCGROUP(cando);
CL_DEFUN bool chem__overlap_solvent(NVector_sp solute_xvec,
                                    NVector_sp solute_yvec,
                                    NVector_sp solute_zvec,
                                    Molecule_sp solvent,
                                    double atom_max_diameter_squared)
{
  for ( size_t resi=0, endResi(solvent->contentSize()) ; resi<endResi; resi++ ) {
    Residue_sp solvres = gc::As_unsafe<Residue_sp>(solvent->contentAt(resi));
    for ( size_t ati=0, endAti(solvres->contentSize()); ati<endAti; ati++) {
      Atom_sp solvatom = gc::As_unsafe<Atom_sp>(solvres->contentAt(ati));
      Vector3 solvent_pos = solvatom->getPosition();
      double x_solvent = solvent_pos.getX();
      double y_solvent = solvent_pos.getY();
      double z_solvent = solvent_pos.getZ();
      for ( size_t idx=0, endIdx(solute_xvec->length()); idx<endIdx; ++idx ) {
        double x_solute = (*solute_xvec)[idx];
        double x_delta =  x_solute - x_solvent;
        double x_delta_squared = x_delta*x_delta;
        if (x_delta_squared < atom_max_diameter_squared) {
          double y_solute = (*solute_yvec)[idx];
          double y_delta =  y_solute - y_solvent;
          double y_delta_squared = y_delta*y_delta;
          if (y_delta_squared < atom_max_diameter_squared) {
            double z_solute = (*solute_zvec)[idx];
            double z_delta =  z_solute - z_solvent;
            double z_delta_squared = z_delta*z_delta;
            if (z_delta_squared < atom_max_diameter_squared) {
              double delta_squared = x_delta_squared+y_delta_squared+z_delta_squared;
              if (delta_squared < atom_max_diameter_squared) {
                return false;
              }
            }
          }
        }
      }
    }
  }
  return true;
}

SYMBOL_EXPORT_SC_(KeywordPkg,macrocycle);
CL_LAMBDA(&optional (name :macrocycle))
CL_DEFUN core::T_sp chem__sample_aggregate(core::T_sp name) {
  if (name==kw::_sym_macrocycle) {
    core::T_sp fin = core::StringInputStream_O::make(R"big(@<TRIPOS>MOLECULE
NIL
240 258 1 0 0
POLYMER
CHARGES


@<TRIPOS>FF_PBC
0.1 1 0.0 0.0 0.0 29.9341 25.152 27.7661 

@<TRIPOS>ATOM
1 O 3.458 1.716 -9.712 O.2 1 NIL_1 0.0000
2 C_2091621551 2.656 1.853 -8.803 C.2 1 NIL_1 0.0000
3 C_2091621549 1.143 1.927 -9.014 C.3 1 NIL_1 0.0000
4 NH 3.004 1.866 -7.465 N.am 1 NIL_1 0.0000
5 N 0.670 2.562 -10.344 N.3 1 NIL_1 0.0000
6 C_2091621550 0.469 0.466 -8.949 C.3 1 NIL_1 0.0000
7 C_2091621746 4.373 2.003 -6.946 C.3 1 NIL_1 0.0000
8 C_2091621547 -0.586 1.821 -10.858 C.3 1 NIL_1 0.0000
9 C_2091621577 0.574 4.033 -10.387 C.3 1 NIL_1 0.0000
10 C_2091621542 -0.844 0.554 -9.878 C.3 1 NIL_1 0.0000
11 C_2091621744 4.472 1.115 -5.684 C.2 1 NIL_1 0.0000
12 C_2091621747 4.678 3.506 -6.685 C.3 1 NIL_1 0.0000
13 C_2091621578 1.896 4.736 -10.106 C.2 1 NIL_1 0.0000
14 C_2091621752 -2.055 0.828 -8.943 C.2 1 NIL_1 0.0000
15 NH_0 -1.037 -0.682 -10.676 N.am 1 NIL_1 0.0000
16 O_0 3.689 0.201 -5.471 O.2 1 NIL_1 0.0000
17 NH_1 5.463 1.440 -4.772 N.am 1 NIL_1 0.0000
18 C_2091621749 6.177 3.879 -6.763 C.3 1 NIL_1 0.0000
19 C_2091621579 1.979 5.576 -9.045 C.2 1 NIL_1 0.0000
20 C_2091621583 2.977 4.618 -10.932 C.2 1 NIL_1 0.0000
21 N_1 -2.796 -0.285 -8.569 N.am 1 NIL_1 0.0000
22 O_1 -2.294 1.935 -8.491 O.2 1 NIL_1 0.0000
23 C_2091621544 -1.587 -1.851 -10.158 C.2 1 NIL_1 0.0000
24 C_2091621742 5.779 0.600 -3.611 C.3 1 NIL_1 0.0000
25 NH2 6.720 3.670 -8.109 N.3 1 NIL_1 0.0000
26 C_2091621580 3.106 6.279 -8.804 C.2 1 NIL_1 0.0000
27 C_2091621582 4.165 5.177 -10.580 C.2 1 NIL_1 0.0000
28 C_2091621605 2.948 3.989 -12.128 C.2 1 NIL_1 0.0000
29 C_2091621751 -4.224 -0.263 -8.007 C.3 1 NIL_1 0.0000
30 C_2091621539 -2.304 -1.725 -8.800 C.3 1 NIL_1 0.0000
31 O_2 -1.461 -2.920 -10.728 O.2 1 NIL_1 0.0000
32 C_2091621741 6.589 1.367 -2.543 C.3 1 NIL_1 0.0000
33 C_2091621581 4.201 6.072 -9.566 C.2 1 NIL_1 0.0000
34 C_2091621602 5.315 4.740 -11.160 C.2 1 NIL_1 0.0000
35 C_2091621604 4.077 3.775 -12.848 C.2 1 NIL_1 0.0000
36 C_2091621750 -4.652 -1.805 -7.747 C.3 1 NIL_1 0.0000
37 C_2091621540 -3.588 -2.683 -8.593 C.3 1 NIL_1 0.0000
38 C_2091621736 6.628 0.620 -1.185 C.3 1 NIL_1 0.0000
39 N_2 5.237 4.017 -12.265 N.2 1 NIL_1 0.0000
40 NH_2 -6.038 -2.064 -8.187 N.am 1 NIL_1 0.0000
41 C_2091621557 -4.483 -2.114 -6.224 C.2 1 NIL_1 0.0000
42 C_2091621735 7.581 1.298 -0.169 C.2 1 NIL_1 0.0000
43 NH_3 5.267 0.497 -0.630 N.am 1 NIL_1 0.0000
44 C_2091621554 -7.031 -1.622 -7.332 C.2 1 NIL_1 0.0000
45 O_3 -3.382 -2.261 -5.721 O.2 1 NIL_1 0.0000
46 NH_4 -5.633 -2.266 -5.444 N.am 1 NIL_1 0.0000
47 O_4 8.653 1.773 -0.503 O.2 1 NIL_1 0.0000
48 NH_5 7.194 1.302 1.167 N.am 1 NIL_1 0.0000
49 C_2091621738 5.020 0.146 0.689 C.2 1 NIL_1 0.0000
50 C_2091621555 -7.000 -2.392 -5.999 C.3 1 NIL_1 0.0000
51 O_5 -7.655 -0.593 -7.527 O.2 1 NIL_1 0.0000
52 C_2091621718 6.153 0.406 1.715 C.3 1 NIL_1 0.0000
53 O_6 3.966 -0.379 1.007 O.2 1 NIL_1 0.0000
54 C_2091621560 -8.048 -1.929 -4.962 C.3 1 NIL_1 0.0000
55 C_2091621719 5.508 1.057 3.053 C.3 1 NIL_1 0.0000
56 C_2091621722 6.769 -1.011 2.190 C.3 1 NIL_1 0.0000
57 C_2091621561 -9.506 -2.170 -5.405 C.3 1 NIL_1 0.0000
58 N_6 5.399 -0.066 4.092 N.am 1 NIL_1 0.0000
59 C_2091621721 6.044 -1.368 3.588 C.3 1 NIL_1 0.0000
60 NH_6 -10.432 -1.800 -4.331 N.am 1 NIL_1 0.0000
61 C_2091621723 4.947 0.088 5.395 C.2 1 NIL_1 0.0000
62 C_2091621726 6.972 -1.972 4.659 C.2 1 NIL_1 0.0000
63 C_2091621563 -11.047 -0.561 -4.263 C.2 1 NIL_1 0.0000
64 O_7 4.113 0.931 5.682 O.2 1 NIL_1 0.0000
65 C_2091621724 5.516 -0.907 6.444 C.3 1 NIL_1 0.0000
66 NH_7 6.759 -1.567 5.973 N.am 1 NIL_1 0.0000
67 O_8 7.825 -2.798 4.383 O.2 1 NIL_1 0.0000
68 O_9 -11.505 -0.003 -5.247 O.2 1 NIL_1 0.0000
69 C_2091621565 -10.980 0.084 -2.869 C.3 1 NIL_1 0.0000
70 C_2091621514 5.845 -0.137 7.831 C.3 1 NIL_1 0.0000
71 C_2091621731 4.397 -2.012 6.807 C.3 1 NIL_1 0.0000
72 C_2091621566 -11.767 1.412 -2.748 C.3 1 NIL_1 0.0000
73 NH_9 -9.560 0.278 -2.532 N.am 1 NIL_1 0.0000
74 N_9 5.042 -0.830 8.949 N.3 1 NIL_1 0.0000
75 C_2091621730 3.976 -1.763 8.343 C.3 1 NIL_1 0.0000
76 C_2091621568 -13.295 1.209 -2.839 C.3 1 NIL_1 0.0000
77 C_2091621701 -8.979 -0.157 -1.354 C.2 1 NIL_1 0.0000
78 C_2091621570 4.605 0.016 10.071 C.3 1 NIL_1 0.0000
79 C_2091621732 3.833 -3.079 9.112 C.2 1 NIL_1 0.0000
80 NH2_9 -14.008 2.479 -2.680 N.3 1 NIL_1 0.0000
81 C_2091621699 -7.769 0.669 -0.898 C.3 1 NIL_1 0.0000
82 O_10 -9.394 -1.116 -0.724 O.2 1 NIL_1 0.0000
83 C_2091621571 5.777 0.637 10.817 C.2 1 NIL_1 0.0000
84 O_11 4.794 -3.779 9.387 O.2 1 NIL_1 0.0000
85 NH_11 2.523 -3.445 9.372 N.am 1 NIL_1 0.0000
86 N_11 -8.089 1.445 0.398 N.3 1 NIL_1 0.0000
87 C_2091621700 -6.531 -0.305 -0.540 C.3 1 NIL_1 0.0000
88 C_2091621572 6.147 1.908 10.516 C.2 1 NIL_1 0.0000
89 C_2091621576 6.515 -0.060 11.729 C.2 1 NIL_1 0.0000
90 C_2091621715 2.140 -4.687 10.069 C.3 1 NIL_1 0.0000
91 C_2091621697 -7.557 0.670 1.618 C.3 1 NIL_1 0.0000
92 C_2091621584 -9.393 2.125 0.508 C.3 1 NIL_1 0.0000
93 C_2091621692 -6.415 -0.344 1.073 C.3 1 NIL_1 0.0000
94 C_2091621573 7.221 2.474 11.106 C.2 1 NIL_1 0.0000
95 C_2091621575 7.480 0.571 12.449 C.2 1 NIL_1 0.0000
96 C_2091621601 6.378 -1.386 11.957 C.2 1 NIL_1 0.0000
97 C_2091621713 1.024 -5.360 9.247 C.2 1 NIL_1 0.0000
98 C_2091621716 1.721 -4.456 11.547 C.3 1 NIL_1 0.0000
99 C_2091621585 -9.443 3.081 1.693 C.2 1 NIL_1 0.0000
100 C_2091621691 -5.017 0.154 1.531 C.2 1 NIL_1 0.0000
101 NH_12 -6.719 -1.704 1.577 N.am 1 NIL_1 0.0000
102 C_2091621574 7.882 1.809 12.077 C.2 1 NIL_1 0.0000
103 C_2091621598 7.935 0.015 13.604 C.2 1 NIL_1 0.0000
104 C_2091621600 7.042 -2.010 12.960 C.2 1 NIL_1 0.0000
105 O_12 -0.044 -5.706 9.726 O.2 1 NIL_1 0.0000
106 NH_13 1.301 -5.490 7.894 N.am 1 NIL_1 0.0000
107 C_2091621569 2.868 -3.974 12.465 C.3 1 NIL_1 0.0000
108 C_2091621586 -10.117 2.807 2.849 C.2 1 NIL_1 0.0000
109 C_2091621590 -8.729 4.232 1.613 C.2 1 NIL_1 0.0000
110 N_13 -4.144 -0.809 2.021 N.am 1 NIL_1 0.0000
111 O_13 -4.692 1.327 1.457 O.2 1 NIL_1 0.0000
112 C_2091621694 -5.776 -2.721 1.680 C.2 1 NIL_1 0.0000
113 N_14 7.675 -1.261 13.843 N.2 1 NIL_1 0.0000
114 C_2091621711 0.293 -5.290 6.848 C.3 1 NIL_1 0.0000
115 NH2_14 3.124 -2.538 12.335 N.3 1 NIL_1 0.0000
116 C_2091621587 -10.211 3.755 3.819 C.2 1 NIL_1 0.0000
117 C_2091621594 -10.685 1.613 3.135 C.2 1 NIL_1 0.0000
118 C_2091621589 -8.680 5.091 2.653 C.2 1 NIL_1 0.0000
119 C_2091621687 -3.133 -0.573 3.148 C.3 1 NIL_1 0.0000
120 C_2091621689 -4.291 -2.301 1.687 C.3 1 NIL_1 0.0000
121 O_14 -6.115 -3.891 1.722 O.2 1 NIL_1 0.0000
122 C_2091621710 0.502 -3.950 6.106 C.3 1 NIL_1 0.0000
123 C_2091621588 -9.426 4.856 3.752 C.2 1 NIL_1 0.0000
124 C_2091621597 -11.177 3.657 4.772 C.2 1 NIL_1 0.0000
125 C_2091621595 -11.446 1.431 4.242 C.2 1 NIL_1 0.0000
126 C_2091621686 -2.523 -2.015 3.566 C.3 1 NIL_1 0.0000
127 C_2091621690 -3.391 -3.110 2.755 C.3 1 NIL_1 0.0000
128 C_2091621705 -0.626 -3.647 5.088 C.3 1 NIL_1 0.0000
129 N_15 -11.794 2.498 4.936 N.2 1 NIL_1 0.0000
130 NH_15 -2.680 -2.204 5.023 N.am 1 NIL_1 0.0000
131 C_2091621707 -1.030 -2.108 3.149 C.2 1 NIL_1 0.0000
132 C_2091621704 -1.921 -3.109 5.756 C.2 1 NIL_1 0.0000
133 NH_16 -0.158 -2.705 4.053 N.am 1 NIL_1 0.0000
134 O_16 -0.634 -1.716 2.065 O.2 1 NIL_1 0.0000
135 O_17 -2.296 -3.474 6.857 O.2 1 NIL_1 0.0000
136 H_2091621549 0.713 2.514 -8.192 H 1 NIL_1 0.0000
137 HH 2.283 1.615 -6.800 H 1 NIL_1 0.0000
138 H_20916215501 1.166 -0.271 -9.371 H 1 NIL_1 0.0000
139 H_20916215502 0.243 0.163 -7.918 H 1 NIL_1 0.0000
140 H_2091621746 5.074 1.606 -7.689 H 1 NIL_1 0.0000
141 H_20916215471 -0.410 1.482 -11.886 H 1 NIL_1 0.0000
142 H_20916215472 -1.461 2.484 -10.863 H 1 NIL_1 0.0000
143 H_20916215771 -0.182 4.356 -9.660 H 1 NIL_1 0.0000
144 H_20916215772 0.211 4.342 -11.375 H 1 NIL_1 0.0000
145 H_20916217471 4.150 4.112 -7.430 H 1 NIL_1 0.0000
146 H_20916217472 4.269 3.806 -5.711 H 1 NIL_1 0.0000
147 HH_17 -0.923 -0.618 -11.681 H 1 NIL_1 0.0000
148 HH_18 5.841 2.379 -4.795 H 1 NIL_1 0.0000
149 H_20916217491 6.289 4.937 -6.496 H 1 NIL_1 0.0000
150 H_20916217492 6.756 3.292 -6.042 H 1 NIL_1 0.0000
151 H_2091621579 1.118 5.713 -8.393 H 1 NIL_1 0.0000
152 H_20916217421 6.344 -0.279 -3.944 H 1 NIL_1 0.0000
153 H_20916217422 4.834 0.253 -3.175 H 1 NIL_1 0.0000
154 HH21 6.216 4.258 -8.772 H 1 NIL_1 0.0000
155 HH22 7.693 3.973 -8.132 H 1 NIL_1 0.0000
156 H_2091621580 3.165 6.959 -7.957 H 1 NIL_1 0.0000
157 H_2091621605 2.004 3.700 -12.583 H 1 NIL_1 0.0000
158 H_20916217511 -4.885 0.196 -8.754 H 1 NIL_1 0.0000
159 H_20916217512 -4.265 0.324 -7.081 H 1 NIL_1 0.0000
160 H_2091621539 -1.568 -1.943 -8.016 H 1 NIL_1 0.0000
161 H_20916217411 6.150 2.362 -2.389 H 1 NIL_1 0.0000
162 H_20916217412 7.608 1.523 -2.920 H 1 NIL_1 0.0000
163 H_2091621581 5.145 6.546 -9.303 H 1 NIL_1 0.0000
164 H_2091621602 6.290 4.874 -10.699 H 1 NIL_1 0.0000
165 H_2091621604 4.063 3.576 -13.916 H 1 NIL_1 0.0000
166 H_20916215401 -4.017 -2.924 -9.574 H 1 NIL_1 0.0000
167 H_20916215402 -3.326 -3.628 -8.100 H 1 NIL_1 0.0000
168 H_2091621736 7.012 -0.395 -1.351 H 1 NIL_1 0.0000
169 HH_19 -6.138 -3.034 -8.459 H 1 NIL_1 0.0000
170 HH_20 4.480 0.746 -1.217 H 1 NIL_1 0.0000
171 HH_21 -5.550 -2.078 -4.452 H 1 NIL_1 0.0000
172 HH_22 7.544 2.043 1.761 H 1 NIL_1 0.0000
173 H_2091621555 -7.170 -3.453 -6.224 H 1 NIL_1 0.0000
174 H_20916215601 -7.874 -2.468 -4.021 H 1 NIL_1 0.0000
175 H_20916215602 -7.901 -0.863 -4.744 H 1 NIL_1 0.0000
176 H_20916217191 6.169 1.836 3.455 H 1 NIL_1 0.0000
177 H_20916217192 4.522 1.495 2.854 H 1 NIL_1 0.0000
178 H_20916217221 7.848 -0.899 2.360 H 1 NIL_1 0.0000
179 H_20916217222 6.627 -1.802 1.442 H 1 NIL_1 0.0000
180 H_20916215611 -9.734 -1.589 -6.307 H 1 NIL_1 0.0000
181 H_20916215612 -9.649 -3.233 -5.636 H 1 NIL_1 0.0000
182 H_2091621721 5.233 -2.084 3.402 H 1 NIL_1 0.0000
183 HH_23 -10.425 -2.366 -3.491 H 1 NIL_1 0.0000
184 HH_24 7.545 -1.600 6.610 H 1 NIL_1 0.0000
185 H_2091621565 -11.409 -0.636 -2.160 H 1 NIL_1 0.0000
186 H_20916215141 6.918 -0.161 8.060 H 1 NIL_1 0.0000
187 H_20916215142 5.538 0.913 7.754 H 1 NIL_1 0.0000
188 H_20916217311 4.812 -3.022 6.690 H 1 NIL_1 0.0000
189 H_20916217312 3.522 -1.942 6.149 H 1 NIL_1 0.0000
190 H_20916215661 -11.537 1.873 -1.778 H 1 NIL_1 0.0000
191 H_20916215662 -11.434 2.118 -3.520 H 1 NIL_1 0.0000
192 HH_25 -8.982 0.796 -3.182 H 1 NIL_1 0.0000
193 H_2091621730 3.010 -1.240 8.342 H 1 NIL_1 0.0000
194 H_20916215681 -13.556 0.769 -3.809 H 1 NIL_1 0.0000
195 H_20916215682 -13.616 0.513 -2.054 H 1 NIL_1 0.0000
196 H_20916215701 3.942 0.804 9.692 H 1 NIL_1 0.0000
197 H_20916215702 4.011 -0.583 10.771 H 1 NIL_1 0.0000
198 HH21_25 -15.010 2.318 -2.772 H 1 NIL_1 0.0000
199 HH22_25 -13.744 3.111 -3.435 H 1 NIL_1 0.0000
200 H_2091621699 -7.490 1.362 -1.702 H 1 NIL_1 0.0000
201 HH_26 1.803 -2.991 8.823 H 1 NIL_1 0.0000
202 H_20916217001 -6.662 -1.307 -0.970 H 1 NIL_1 0.0000
203 H_20916217002 -5.612 0.112 -0.972 H 1 NIL_1 0.0000
204 H_2091621572 5.592 2.469 9.767 H 1 NIL_1 0.0000
205 H_2091621715 2.991 -5.379 10.053 H 1 NIL_1 0.0000
206 H_20916216971 -8.374 0.081 2.052 H 1 NIL_1 0.0000
207 H_20916216972 -7.173 1.350 2.389 H 1 NIL_1 0.0000
208 H_20916215841 -10.202 1.388 0.555 H 1 NIL_1 0.0000
209 H_20916215842 -9.560 2.705 -0.409 H 1 NIL_1 0.0000
210 H_2091621573 7.519 3.489 10.852 H 1 NIL_1 0.0000
211 H_2091621601 5.790 -2.008 11.289 H 1 NIL_1 0.0000
212 H_20916217161 0.877 -3.756 11.601 H 1 NIL_1 0.0000
213 H_20916217162 1.361 -5.410 11.954 H 1 NIL_1 0.0000
214 HH_27 -7.637 -1.869 1.970 H 1 NIL_1 0.0000
215 H_2091621574 8.677 2.305 12.631 H 1 NIL_1 0.0000
216 H_2091621598 8.416 0.599 14.385 H 1 NIL_1 0.0000
217 H_2091621600 7.211 -3.084 12.973 H 1 NIL_1 0.0000
218 HH_28 2.272 -5.557 7.612 H 1 NIL_1 0.0000
219 H_20916215691 2.593 -4.182 13.507 H 1 NIL_1 0.0000
220 H_20916215692 3.782 -4.540 12.243 H 1 NIL_1 0.0000
221 H_2091621590 -8.159 4.451 0.712 H 1 NIL_1 0.0000
222 H_20916217111 0.365 -6.118 6.131 H 1 NIL_1 0.0000
223 H_20916217112 -0.706 -5.326 7.298 H 1 NIL_1 0.0000
224 HH21_28 3.858 -2.264 12.987 H 1 NIL_1 0.0000
225 HH22_28 2.288 -2.020 12.601 H 1 NIL_1 0.0000
226 H_2091621594 -10.481 0.735 2.528 H 1 NIL_1 0.0000
227 H_2091621589 -8.094 6.005 2.587 H 1 NIL_1 0.0000
228 H_20916216871 -3.656 -0.116 3.998 H 1 NIL_1 0.0000
229 H_20916216872 -2.347 0.112 2.808 H 1 NIL_1 0.0000
230 H_2091621689 -3.890 -2.459 0.679 H 1 NIL_1 0.0000
231 H_20916217101 1.465 -3.999 5.581 H 1 NIL_1 0.0000
232 H_20916217102 0.578 -3.125 6.827 H 1 NIL_1 0.0000
233 H_2091621588 -9.477 5.600 4.545 H 1 NIL_1 0.0000
234 H_2091621597 -11.532 4.515 5.336 H 1 NIL_1 0.0000
235 H_2091621595 -11.639 0.449 4.666 H 1 NIL_1 0.0000
236 H_20916216901 -4.041 -3.659 3.450 H 1 NIL_1 0.0000
237 H_20916216902 -2.752 -3.847 2.252 H 1 NIL_1 0.0000
238 H_2091621705 -0.892 -4.576 4.568 H 1 NIL_1 0.0000
239 HH_29 -3.257 -1.544 5.531 H 1 NIL_1 0.0000
240 HH_30 0.816 -2.427 4.048 H 1 NIL_1 0.0000
@<TRIPOS>BOND
1 2 3 1 
2 2 1 2 
3 2 4 1 
4 3 5 1 
5 3 136 1 
6 4 7 1 
7 4 137 1 
8 5 8 1 
9 5 9 1 
10 6 3 1 
11 6 10 1 
12 6 138 1 
13 6 139 1 
14 7 140 1 
15 8 10 1 
16 8 141 1 
17 8 142 1 
18 9 143 1 
19 9 144 1 
20 10 14 1 
21 10 15 1 
22 11 16 2 
23 11 7 1 
24 11 17 1 
25 12 7 1 
26 13 9 1 
27 14 21 1 
28 17 24 1 
29 18 12 1 
30 18 25 1 
31 19 13 2 
32 20 27 2 
33 20 13 1 
34 20 28 1 
35 21 29 1 
36 22 14 2 
37 23 15 1 
38 23 30 1 
39 24 32 1 
40 26 19 1 
41 27 33 1 
42 27 34 1 
43 28 35 2 
44 29 36 1 
45 30 21 1 
46 31 23 2 
47 32 38 1 
48 33 26 2 
49 35 39 1 
50 36 40 1 
51 36 41 1 
52 37 30 1 
53 37 36 1 
54 38 42 1 
55 39 34 2 
56 41 46 1 
57 42 48 1 
58 43 38 1 
59 44 40 1 
60 45 41 2 
61 46 50 1 
62 47 42 2 
63 49 43 1 
64 50 44 1 
65 51 44 2 
66 52 48 1 
67 52 49 1 
68 53 49 2 
69 54 50 1 
70 55 52 1 
71 56 59 1 
72 56 52 1 
73 57 54 1 
74 58 55 1 
75 59 58 1 
76 59 62 1 
77 60 57 1 
78 61 58 1 
79 61 64 2 
80 61 65 1 
81 63 60 1 
82 66 62 1 
83 66 65 1 
84 67 62 2 
85 68 63 2 
86 69 63 1 
87 70 65 1 
88 71 75 1 
89 71 65 1 
90 72 69 1 
91 73 69 1 
92 74 70 1 
93 74 78 1 
94 75 74 1 
95 76 72 1 
96 76 80 1 
97 77 81 1 
98 77 73 1 
99 79 75 1 
100 81 86 1 
101 81 200 1 
102 82 77 2 
103 83 78 1 
104 83 89 1 
105 84 79 2 
106 85 90 1 
107 85 79 1 
108 86 91 1 
109 86 92 1 
110 87 81 1 
111 87 93 1 
112 87 202 1 
113 87 203 1 
114 88 83 2 
115 88 94 1 
116 89 95 2 
117 89 96 1 
118 90 97 1 
119 91 93 1 
120 91 206 1 
121 91 207 1 
122 93 100 1 
123 93 101 1 
124 95 102 1 
125 95 103 1 
126 96 104 2 
127 97 106 1 
128 98 90 1 
129 98 107 1 
130 99 92 1 
131 100 110 1 
132 101 214 1 
133 102 94 2 
134 104 113 1 
135 105 97 2 
136 106 114 1 
137 107 115 1 
138 108 99 2 
139 108 117 1 
140 109 118 2 
141 109 99 1 
142 110 119 1 
143 111 100 2 
144 112 101 1 
145 112 120 1 
146 113 103 2 
147 114 122 1 
148 116 108 1 
149 116 124 1 
150 118 123 1 
151 119 126 1 
152 119 228 1 
153 119 229 1 
154 120 110 1 
155 120 230 1 
156 121 112 2 
157 122 128 1 
158 123 116 2 
159 124 129 2 
160 125 117 2 
161 126 131 1 
162 127 120 1 
163 127 126 1 
164 127 236 1 
165 127 237 1 
166 128 132 1 
167 129 125 1 
168 130 126 1 
169 130 132 1 
170 130 239 1 
171 131 133 1 
172 133 128 1 
173 134 131 2 
174 135 132 2 
175 145 12 1 
176 146 12 1 
177 147 15 1 
178 148 17 1 
179 149 18 1 
180 150 18 1 
181 151 19 1 
182 152 24 1 
183 153 24 1 
184 154 25 1 
185 155 25 1 
186 156 26 1 
187 157 28 1 
188 158 29 1 
189 159 29 1 
190 160 30 1 
191 161 32 1 
192 162 32 1 
193 163 33 1 
194 164 34 1 
195 165 35 1 
196 166 37 1 
197 167 37 1 
198 168 38 1 
199 169 40 1 
200 170 43 1 
201 171 46 1 
202 172 48 1 
203 173 50 1 
204 174 54 1 
205 175 54 1 
206 176 55 1 
207 177 55 1 
208 178 56 1 
209 179 56 1 
210 180 57 1 
211 181 57 1 
212 182 59 1 
213 183 60 1 
214 184 66 1 
215 185 69 1 
216 186 70 1 
217 187 70 1 
218 188 71 1 
219 189 71 1 
220 190 72 1 
221 191 72 1 
222 192 73 1 
223 193 75 1 
224 194 76 1 
225 195 76 1 
226 196 78 1 
227 197 78 1 
228 198 80 1 
229 199 80 1 
230 201 85 1 
231 204 88 1 
232 205 90 1 
233 208 92 1 
234 209 92 1 
235 210 94 1 
236 211 96 1 
237 212 98 1 
238 213 98 1 
239 215 102 1 
240 216 103 1 
241 217 104 1 
242 218 106 1 
243 219 107 1 
244 220 107 1 
245 221 109 1 
246 222 114 1 
247 223 114 1 
248 224 115 1 
249 225 115 1 
250 226 117 1 
251 227 118 1 
252 231 122 1 
253 232 122 1 
254 233 123 1 
255 234 124 1 
256 235 125 1 
257 238 128 1 
258 240 133 1 
@<TRIPOS>SUBSTRUCTURE
1 NIL_1 1 RESIDUE 1 NIL NIL 1
)big");
    Mol2File mfin(fin);
    return mol2Read(mfin);
  } else {
    SIMPLE_ERROR("Illegal sample-aggregate option {}", _rep_(name) );
  }
}


};
