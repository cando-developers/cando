
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
CL_DEFUN geom::SimpleVectorCoordinate_sp chem__make_simple_vector_coordinate_from_atom_list(core::List_sp atoms)
{
  size_t num = core::cl__length(atoms);
  geom::SimpleVectorCoordinate_sp coords = geom::SimpleVectorCoordinate_O::make(num);
  int idx(0);
  for ( auto cur : atoms ) {
    Atom_sp a = gc::As<chem::Atom_sp>(oCar(cur));
    Vector3 v = a->getPosition();
    (*coords)[idx] = v;
    ++idx;
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
                                         INTERN_(kw,total), core::cl__file_length(fin.fIn));
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
            core::eval::funcall(progress_advance, progress_bar, core::cl__file_position(fin.fIn,nil<core::T_O>()));
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
CL_DEFUN bool chem__overlap_solvent(core::SimpleVector_double_sp solute_xvec,
                                    core::SimpleVector_double_sp solute_yvec,
                                    core::SimpleVector_double_sp solute_zvec,
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








};
