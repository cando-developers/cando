/*
    File: readAmberParameters.cc
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

#include <clasp/core/foundation.h>
#include <cando/chem/readAmberParameters.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/numerics.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/bformat.h>
#include <clasp/core/readtable.h>
#include <clasp/core/primitives.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/array.h>
#include <cando/chem/chemInfo.h>
#include <cando/chem/units.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/forceField.h>
#include <clasp/core/wrappers.h>


namespace chem {


CL_DEFUN ReadAmberParameters_sp chem__make_read_amber_parameters()
{
  return ReadAmberParameters_O::create();
}



CL_LISPIFY_NAME("readTypes");
CL_DEFMETHOD void	ReadAmberParameters_O::readTypes(core::T_sp fin)
{
  this->_Types = this->parseTypeRules(fin);
}


CL_LISPIFY_NAME("readParameters");
CL_DEFMETHOD void	ReadAmberParameters_O::readParameters(core::T_sp fin)
{
  this->_ForceField = this->parseAmberFormattedForceFieldOrFrcMod(fin);
}



CL_LISPIFY_NAME("getForceField");
CL_DEFMETHOD ForceField_sp ReadAmberParameters_O::getForceField()
{
  if ( this->_ForceField.unboundp() ) {
    GC_ALLOCATE(ForceField_O, temp );
    this->_ForceField = temp;
  }
  if (!this->_Types.unboundp()) this->_ForceField->setFFTypeDb(this->_Types);
  return this->_ForceField;
}


/*! Check if the stream contains a FrcMod file or a regular Amber Parameter file.
    If it's a FrcMod file then frcmod == T.
    If it contains a MASS entry then mass == T.
    if it contains a NONBOND entry then NONBOND == T.
    Return (values frcmod mass nonbond).
*/

core::T_mv ReadAmberParameters_O::determineParmSetFrcModType(core::T_sp fin)
{
  bool bNew = false;
  bool bMass = false;
  bool bNonBond = false;
  while ( true ) {
    core::T_mv mv_line = core::cl__read_line(fin,_Nil<T_O>());
    if ( mv_line.nilp() ) break;
    string line = mv_line.as<core::Str_O>()->get();
    LOG(BF("Read line(%s)") % line  );
    if ( line.size()>=4 && line.substr(0,4) == "MASS") {
      bNew = true;
      bMass = true;
    } else if ( line.size()>=4 && line.substr(0,4) == "BOND") {
      bNew = true;
    } else if ( line.size()>=4 && line.substr(0,4) == "ANGL") {
      bNew = true;
    } else if ( line.size()>=4 && line.substr(0,4) == "DIHE") {
      bNew = true;
    } else if ( line.size()>=4 && line.substr(0,4) == "IMPR") {
      bNew = true;
    } else if ( line.size()>=4 && line.substr(0,4) == "HBON") {
      bNew = true;
    } else if ( line.size()>=4 && line.substr(0,4) == "NONB") {
      bNew = true;
      bNonBond = true;
    }
  }
  clasp_file_position_set(fin,core::clasp_make_fixnum(0));
  return Values( _lisp->_boolean(bNew),_lisp->_boolean(bMass),_lisp->_boolean(bNonBond));
}

FFTypesDb_sp ReadAmberParameters_O::parseTypeRules(core::T_sp fin)
{
  GC_ALLOCATE(WildElementDict_O, wildCardElementDictionary );
  GC_ALLOCATE(FFTypesDb_O, ffTypesDb );
  uint lineno = 1;
  vector< pair< uint, string> > entries;
  while ( true ) {
    core::T_mv mv_line = core::cl__read_line(fin,_Nil<T_O>());
    if ( mv_line.nilp() ) break;
    string line = mv_line.as<core::Str_O>()->get();
    LOG(BF("Read line(%s)") % line  );
    if ( line.substr(0,8) == "WILDATOM" ) {
      vector<string> names = core::split(line," \t");
      string wildName = names[1];
      wildCardElementDictionary->addWildName(chemkw_intern(wildName));
      for ( vector<string>::iterator it=names.begin()+1; it!=names.end(); it++ )
      {
        core::Symbol_sp wildSymbol = chemkw_intern(wildName);
        core::Symbol_sp wildOtherSymbol = chemkw_intern(*it);
        wildCardElementDictionary->addWildNameMap(wildSymbol,wildOtherSymbol);
      }
    } else if ( line.substr(0,3) == "ATD" )
    {
      pair<uint,string> entry;
      entry.first = lineno;
      entry.second = line;
      entries.push_back( entry );
    }
  }
  vector< pair<uint,string> >::iterator ei;
  for ( ei=entries.begin(); ei!=entries.end(); ei++ )
  {
    ChemInfo_sp typeRule = ChemInfo_O::create();
    typeRule->compileAntechamber(ei->second,wildCardElementDictionary);
    if ( typeRule->compileSucceeded() ) {
      ffTypesDb->add(typeRule);
    } else {
      SIMPLE_ERROR(BF("Antechamber compile failed on: "+ei->second+"\n"+typeRule->compilerMessage() ));
    }
  }
  return ffTypesDb;
}




FFNonbondDb_sp ReadAmberParameters_O::parseMasses(core::T_sp fin, FFNonbondDb_sp ffNonbondDb )
{
  core::DynamicScopeManager scope(cl::_sym_STARpackageSTAR, _lisp->keywordPackage());
  core::ReadTable_sp readtable = gctools::As<core::ReadTable_sp>(cl::_sym_STARreadtableSTAR->symbolValue());
  core::ReadTable_sp copy_readtable = readtable->copyReadTable(_Nil<core::T_O>());
  copy_readtable->setf_readtable_case(kw::_sym_preserve);
  scope.bind(cl::_sym_STARreadtableSTAR,copy_readtable);
  bool done = false;
  while ( not done ) {
    core::T_mv ol = core::cl__read_line(fin);
    string line = gc::As<core::String_sp>(ol)->get_std_string();
    if ( line.size() == 0 ) {
      done = true;
    } else {
      core::T_sp linestream = core::cl__make_string_input_stream(ol,core::make_fixnum(0),_Nil<core::T_O>());
      LOG(BF("Parsing line|%s|") % line.c_str()  );
      printf("%s:%d:%s parseMasses line: %s\n", __FILE__, __LINE__, __FUNCTION__, line.c_str());
      core::Symbol_sp typeSymbol = gc::As<core::Symbol_sp>(core::cl__read(linestream,_Nil<core::T_O>()));
      core::T_sp omass = core::cl__read(linestream,_Nil<core::T_O>());
      core::T_sp maybePolarizability = core::cl__read(linestream,_Nil<core::T_O>());
      double mass = core::clasp_to_double(omass);
      double polarizability = 0.0;
      if ( maybePolarizability.fixnump() || gc::IsA<core::Float_sp>(maybePolarizability)) {
        polarizability = core::clasp_to_double(gc::As_unsafe<core::Float_sp>(maybePolarizability));
      }
      FFNonbond_sp ffNonbond;
      core::T_sp nonbond = ffNonbondDb->findType(typeSymbol);
      bool newp = false;
      if (nonbond.nilp()) {
        ffNonbond = gctools::GC<FFNonbond_O>::allocate_with_default_constructor();
        newp = true;
      } else if (gc::IsA<FFNonbond_sp>(nonbond)) {
        ffNonbond = gc::As_unsafe<FFNonbond_sp>(nonbond);
      } else {
        SIMPLE_ERROR(BF("Illegal atom type %s") % _rep_(nonbond));
      }
      ffNonbond->setType(typeSymbol);
      ffNonbond->setMass(mass);
      ffNonbond->setPolarizability(polarizability);
      printf("%s:%d Adding nonbond type with name: %s\n", __FILE__, __LINE__, _rep_(typeSymbol).c_str());
      if (newp) ffNonbondDb->add(ffNonbond);
    }
  }
  return ffNonbondDb;
}
 
FFStretchDb_sp ReadAmberParameters_O::parseStretchDb(core::T_sp fin)
{
    // read stretch terms 
  GC_ALLOCATE(FFStretchDb_O, ffStretchDb );
  bool done = false;
  while ( not done )
  {
    string line = core::cl__read_line(fin).as<core::Str_O>()->get();
    if ( line.size() == 0 )
    {
      done = true;
    } else {
      LOG(BF("Parsing line|%s|") % line.c_str()  );
      string types = line.substr(0,5);
      vector<string>typeParts = core::split(types,"-");
      string type1Name = core::trimWhiteSpace(typeParts[0]);
      string type2Name = core::trimWhiteSpace(typeParts[1]);
      GC_ALLOCATE(FFStretch_O, ffStretch );
      ffStretch->_Type1 = chemkw_intern(type1Name);
      ffStretch->_Type2 = chemkw_intern(type2Name);
      string parms = line.substr(6);
      vector<string> parmsParts = core::split(parms);
      double kb_kJPerNanometerSquared = kCalPerAngstromSquared_to_kJPerNanometerSquared(atof(parmsParts[0].c_str()));
      double r0_Nanometer = angstrom_to_nanometer(atof(parmsParts[1].c_str()));
      ffStretch->setKb_kJPerNanometerSquared(kb_kJPerNanometerSquared);
      ffStretch->setR0_Nanometer(r0_Nanometer);
//            print "(%s)-(%s) k=%lf len=%lf"%(ffStretch._Type1,ffStretch._Type2,ffStretch._kb, ffStretch._r0)
      ffStretchDb->add(ffStretch);
    }
  }
  return ffStretchDb;
}


FFAngleDb_sp ReadAmberParameters_O::parseAngleDb(core::T_sp fin)
{
    //
    // read angle terms 
    //
  GC_ALLOCATE(FFAngleDb_O, ffAngleDb );
  bool done = false;
  while ( not done )
  {
    string line = core::cl__read_line(fin).as<core::Str_O>()->get();
    if (line.size()== 0)
    {
      done = true;
    } else
    {
      LOG(BF("Parsing line|%s|") % line.c_str()  );
      string types = line.substr(0,8);
      vector<string>typeParts = core::split(types,"-");
      string t1 = core::trimWhiteSpace(typeParts[0]);
      string t2 = core::trimWhiteSpace(typeParts[1]);
      string t3 = core::trimWhiteSpace(typeParts[2]);
      GC_ALLOCATE(FFAngle_O, ffAngle );
      ffAngle->_Type1 = chemkw_intern(t1);
      ffAngle->_Type2 = chemkw_intern(t2);
      ffAngle->_Type3 = chemkw_intern(t3);
//      printf("%s:%d  Parsing line: %s\n", __FILE__, __LINE__, line.c_str());
      fflush(stdout);
      string parms = line.substr(9);
      vector<string> parmsParts = core::split(parms);
      ffAngle->_K2__kJPerRadianSquared = kCalPerRadianSquared_to_kJPerRadianSquared(atof(parmsParts[0].c_str()));
      ffAngle->_AngRad = core::radians(atof(parmsParts[1].c_str()));
//            print "(%s)-(%s)-(%s) k=%lf ang(deg)=%lf"%(ffAngle._Type1,ffAngle._Type2, ffAngle._Type3,ffAngle._K2, ffAngle._AngRad/0.0174533)
      ffAngleDb->add(ffAngle);
    }
  }
  return ffAngleDb;
};



FFPtorDb_sp ReadAmberParameters_O::parsePtorDb(core::T_sp fin)
{
    //
    // read Ptor terms 
    //
  GC_ALLOCATE(FFPtorDb_O, ffPtorDb );
  bool done = false;
  while (not done )
  {
    string line = core::cl__read_line(fin).as<core::Str_O>()->get();
    if (line.size()== 0)
    {
      done = true;
    } else
    {
      LOG(BF("Parsing line|%s|") % line.c_str()  );
      string types = line.substr(0,12);
      vector<string>typeParts = core::split(types,"-");
      string t1 = core::trimWhiteSpace(typeParts[0]);
      string t2 = core::trimWhiteSpace(typeParts[1]);
      string t3 = core::trimWhiteSpace(typeParts[2]);
      string t4 = core::trimWhiteSpace(typeParts[3]);
      GC_ALLOCATE(FFPtor_O, ffPtor );
      core::Symbol_sp st1;
      if ( t1 != "X" ) {
        st1 = chemkw_intern(t1);
      } else {
        st1 = _Nil<core::T_O>();
      }
      core::Symbol_sp st2 = chemkw_intern(t2);
      core::Symbol_sp st3 = chemkw_intern(t3);
      core::Symbol_sp st4;
      if ( t4 != "X" ) {
        st4 = chemkw_intern(t4);
      } else {
        st4 = _Nil<core::T_O>();
      }
      ffPtor->setTypes(st1,st2,st3,st4);
      string parms = line.substr(13);
      vector<string> parmsParts = core::split(parms);
      double idivf =  atof(parmsParts[0].c_str());
      double pk = atof(parmsParts[1].c_str());
      double phaseRad = atof(parmsParts[2].c_str())*0.0174533;
      int pn = abs(int(atof(parmsParts[3].c_str())));
      ASSERTF(pn>=1 && pn<=6,BF("pn[%d] must range between [1,6]") % pn);
      ffPtor->setV_kCal(pn,pk/idivf);
      ffPtor->setPhaseRad(pn,phaseRad);
      ffPtorDb->add(ffPtor);
    }
  }
  return ffPtorDb;
}





FFItorDb_sp ReadAmberParameters_O::parseItorDb(core::T_sp fin)
{
    //
    // read Itor terms 
    //
  GC_ALLOCATE(FFItorDb_O, ffItorDb );
  bool done = false;
  while (not done )
  {
    string line = core::cl__read_line(fin).as<core::Str_O>()->get();
//      printf("%s:%d line: %s\n", __FILE__, __LINE__, line.c_str());
    if (line.size()== 0)
    {
      done = true;
    } else
    {
      LOG(BF("Parsing line|%s|") % line.c_str()  );
      string types = line.substr(0,12);
      vector<string>typeParts = core::split(types,"-");
      string t1 = core::trimWhiteSpace(typeParts[0]);
      string t2 = core::trimWhiteSpace(typeParts[1]);
      string t3 = core::trimWhiteSpace(typeParts[2]);
      string t4 = core::trimWhiteSpace(typeParts[3]);
      GC_ALLOCATE(FFItor_O, ffItor );
      core::Symbol_sp st1;
      if ( t1 != "X" ) {
        st1 = chemkw_intern(t1);
      } else {
        st1 = _Nil<core::T_O>();
      }
      core::Symbol_sp st2;
      if ( t2 != "X" ) {
        st2 = chemkw_intern(t2);
      } else {
        st2 = _Nil<core::T_O>();
      }
      core::Symbol_sp st3 = chemkw_intern(t3);
      core::Symbol_sp st4;
      if ( t4 != "X" ) {
        st4 = chemkw_intern(t4);
      } else {
        st4 = _Nil<core::T_O>();
      }
      ffItor->setTypes(st1,st2,st3,st4);
      string parms = line.substr(13);
      vector<string> parmsParts = core::split(parms);
      double pk = atof(parmsParts[0].c_str());
      double phaseRad = atof(parmsParts[1].c_str())*0.0174533;
      int pn = abs(int(atof(parmsParts[2].c_str())));
      ASSERTF(pn>=1 && pn<=FFItor_O::IMaxPeriodicity,BF("Illegal pn[%d] must be in [1,%d]")
              % pn % FFItor_O::IMaxPeriodicity );
      ffItor->setV_kCal(pn,pk);
      ffItor->setPhaseRad(pn,phaseRad);
      ffItorDb->add(ffItor);
    }
  }
  return ffItorDb;
}



string ReadAmberParameters_O::parseNonbondLabelKindNB(core::T_sp fin)
{
  string line = core::cl__read_line(fin).as<core::Str_O>()->get();
  vector<string>parts = core::split(line);
  if (parts.size() == 0) {
    printf("%s:%d Could not split line: %s\n", __FILE__, __LINE__, line.c_str());
  }
  string label = parts[0];
  string kindnb = parts[1];
  if ( kindnb != "RE" )
  {
    SIMPLE_ERROR(BF("Nonbond parameters must be of kindnb=RE this file has kindnb="+kindnb));
  }
  return kindnb;
}

/*!
 * This modifies entries in the ffNonbondDb database
 */
void ReadAmberParameters_O::parseNonbondDb(core::T_sp fin, FFNonbondDb_sp ffNonbondDb)
{
    bool done = false;
    while ( 1 ) {
      string line = core::cl__read_line(fin).as<core::Str_O>()->get();
//      printf("%s:%d:%s line: %s\n", __FILE__, __LINE__, __FUNCTION__, line.c_str());
      line = core::trimWhiteSpace(line);
      if ( line.size() == 0 ) break;
      else {
        LOG(BF("Parsing line|%s|") % line.c_str()  );
        FFNonbond_sp ffNonbond;
        string type = core::trimWhiteSpace(line.substr(0,4));
        core::Symbol_sp stype = chemkw_intern(type);
        if ( ffNonbondDb->hasType(stype) ) {
          ffNonbond = ffNonbondDb->findType(stype);
        } else {
          SIMPLE_ERROR(BF("Could not find type: %s") % type);
        }
        string parms = line.substr(5);
        vector<string>parmsParts = core::split(parms);
        double radius  = atof(parmsParts[0].c_str());
        double edep = atof(parmsParts[1].c_str());
//	    print "parseNonbondDb::",
//	    print ffNonbond.asXml().asString()
        ffNonbond->setRadius_Angstroms(radius);
        ffNonbond->setEpsilon_kCal(edep);
//	    print "parseNonbondDb::",
//	    print ffNonbond.asXml().asString()
//            print "parseNonbondDb::(%s) radius=%lf edep=%lf"%(ffNonbond._Type, ffNonbond._Radius,ffNonbond._Well)
      }
    }
}


void ReadAmberParameters_O::parseAtomEquivalences(core::T_sp fin)
{
  SIMPLE_WARN(BF("Warning!  Skipping force field atom equivalences\n"));
  while (1) {
    string line = core::cl__read_line(fin).as<core::Str_O>()->get();
    if (line=="") return;
  }
}



ForceField_sp ReadAmberParameters_O::parseAmberFormattedForceField(core::T_sp fin)
{_OF();
    string line = core::cl__read_line(fin).as<core::Str_O>()->get();
    GC_ALLOCATE(FFNonbondDb_O, ffNonbondsDb );
    ffNonbondsDb = this->parseMasses(fin,ffNonbondsDb);
    FFStretchDb_sp ffStretchesDb = this->parseStretchDb(fin);
    FFAngleDb_sp ffAnglesDb = this->parseAngleDb(fin);
    FFPtorDb_sp ffPtorsDb = this->parsePtorDb(fin);
    FFItorDb_sp ffItorsDb = this->parseItorDb(fin);
    SIMPLE_WARN(BF("Warning!  Skipping 10-12 hbond\n"));
    core::cl__read_line(fin); // skp 10-12 hbond
    core::cl__read_line(fin); // blank
    this->parseAtomEquivalences(fin);
    string kindNb = this->parseNonbondLabelKindNB(fin);
    this->parseNonbondDb(fin,ffNonbondsDb);
    LOG(BF("All parameters read fine") );
    GC_ALLOCATE(ForceField_O, ff );
    ff->setTitle(line);
    ff->setFFStretchDb(ffStretchesDb);
//    ff->setFFTypeDb( _Nil<FFTypesDb_O>());
    ff->setFFAngleDb(ffAnglesDb);
    ff->setFFPtorDb(ffPtorsDb);
    ff->setFFItorDb(ffItorsDb);
    ff->setFFNonbondDb(ffNonbondsDb);
    //
    // parameters from  Antechamber paper
    //
    ff->_Angles->addZConstant(chemkw_intern("H"), 0.784);
    ff->_Angles->addZConstant(chemkw_intern("C"), 1.183);
    ff->_Angles->addZConstant(chemkw_intern("N"), 1.212);
    ff->_Angles->addZConstant(chemkw_intern("O"), 1.219);
    ff->_Angles->addZConstant(chemkw_intern("F"), 1.166);
    ff->_Angles->addZConstant(chemkw_intern("Cl"), 1.272);
    ff->_Angles->addZConstant(chemkw_intern("Br"), 1.378);
    ff->_Angles->addZConstant(chemkw_intern("I"), 1.398);
    ff->_Angles->addZConstant(chemkw_intern("P"), 1.620);
    ff->_Angles->addZConstant(chemkw_intern("S"), 1.280);
    ff->_Angles->addCConstant(chemkw_intern("C"), 1.339);
    ff->_Angles->addCConstant(chemkw_intern("N"), 1.300);
    ff->_Angles->addCConstant(chemkw_intern("O"), 1.249);
    ff->_Angles->addCConstant(chemkw_intern("P"), 0.906);
    ff->_Angles->addCConstant(chemkw_intern("S"), 1.448);

    //
    // Parameters from Antechamber for estimating stretching force constants
    //
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("H"), 0.738, 4.661 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("C"), 1.526, 7.643 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("N"), chemkw_intern("N"), 1.441, 7.634 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("O"), chemkw_intern("O"), 1.460, 7.561 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("F"), chemkw_intern("F"), 1.406, 7.358 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("Cl"), chemkw_intern("Cl"), 2.031, 8.648 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("Br"), chemkw_intern("Br"), 2.337, 9.012 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("I"), chemkw_intern("I"), 2.836, 9.511 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("P"), chemkw_intern("P"), 2.324, 8.805 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("S"), chemkw_intern("S"), 2.038, 8.316 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("C"), 1.090, 6.217 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("N"), 1.010, 6.057 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("O"), 0.960, 5.794 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("F"), 0.920, 5.600 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("Cl"), 1.280, 6.937 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("Br"), 1.410, 7.301 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("I"), 1.600, 7.802 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("P"), 1.410, 7.257 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("H"), chemkw_intern("S"), 1.340, 7.018 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("N"), 1.470, 7.504 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("O"), 1.440, 7.347 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("F"), 1.370, 7.227 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("Cl"), 1.800, 8.241 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("Br"), 1.940, 8.478 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("I"), 2.160, 8.859 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("P"), 1.830, 8.237 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("C"), chemkw_intern("S"), 1.820, 8.117 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("N"), chemkw_intern("O"), 1.420, 7.526 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("N"), chemkw_intern("F"), 1.420, 7.475 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("N"), chemkw_intern("Cl"), 1.750, 8.266 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("N"), chemkw_intern("Br"), 1.930, 8.593 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("N"), chemkw_intern("I"), 2.120, 8.963 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("N"), chemkw_intern("P"), 1.720, 8.212 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("N"), chemkw_intern("S"), 1.690, 8.073 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("O"), chemkw_intern("F"), 1.410, 7.375 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("O"), chemkw_intern("Cl"), 1.700, 8.097 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("O"), chemkw_intern("Br"), 1.790, 8.276 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("O"), chemkw_intern("I"), 2.110, 8.854 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("O"), chemkw_intern("P"), 1.640, 7.957 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("O"), chemkw_intern("S"), 1.650, 7.922 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("F"), chemkw_intern("Cl"), 1.648, 7.947 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("Cl"), chemkw_intern("I"), 2.550, 9.309 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("Br"), chemkw_intern("I"), 2.671, 9.380 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("F"), chemkw_intern("P"), 1.500, 7.592 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("F"), chemkw_intern("S"), 1.580, 7.733 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("Cl"), chemkw_intern("P"), 2.040, 8.656 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("Cl"), chemkw_intern("S"), 2.030, 8.619 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("Br"), chemkw_intern("P"), 2.240, 8.729 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("Br"), chemkw_intern("S"), 2.210, 8.728 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("I"), chemkw_intern("P"), 2.490, 9.058 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("I"), chemkw_intern("S"), 2.560, 9.161 );
    ff->_Stretches->addEstimateStretch( chemkw_intern("P"), chemkw_intern("S"), 2.120, 8.465 );
    return ff;
}

ForceField_sp ReadAmberParameters_O::parseFrcModFile(core::T_sp fin) {
  GC_ALLOCATE(ForceField_O, ff);
  GC_ALLOCATE(FFNonbondDb_O, nonbondDb);
  int lastSegment = 0;
  bool readMasses = false;
  bool readNonbond = false;
  core::T_mv mv_title = core::cl__read_line(fin,_Nil<core::T_O>());
  if (mv_title.notnilp()) {
    ff->setTitle(mv_title.as<core::Str_O>()->get());
  }
  while ( true ) {
    core::T_mv mv_line = core::cl__read_line(fin,_Nil<T_O>());
    if ( mv_line.nilp() ) break;
    string line = mv_line.as<core::Str_O>()->get();
    LOG(BF("Read line(%s)") % line  );
    if ( line.size()>=4 && line.substr(0,4) == "MASS") {
      nonbondDb = this->parseMasses(fin,nonbondDb);
      lastSegment = 1;
      readMasses = true;
    } else if ( line.size()>=4 && line.substr(0,4) == "BOND") {
      ff->setFFStretchDb(this->parseStretchDb(fin));
      lastSegment = 2;
    } else if ( line.size()>=4 && line.substr(0,4) == "ANGL") {
      ff->setFFAngleDb(this->parseAngleDb(fin));
      lastSegment = 3;
    } else if ( line.size()>=4 && line.substr(0,4) == "DIHE") {
      ff->setFFPtorDb(this->parsePtorDb(fin));
      lastSegment = 4;
    } else if ( line.size()>=4 && line.substr(0,4) == "IMPR") {
      ff->setFFItorDb(this->parseItorDb(fin));
      lastSegment = 5;
    } else if ( line.size()>=4 && line.substr(0,4) == "HBON") {
      SIMPLE_WARN(BF("Skipping HBON FrcMod term"));
    } else if ( line.size()>=4 && line.substr(0,4) == "NONB") {
      this->parseNonbondDb(fin,nonbondDb);
      readNonbond = true;
      lastSegment = 7;
    } else {
      if (line.size()>0) {
        SIMPLE_WARN(BF("Ignoring %s in FrcMod") % line );
      }
    }
  }
  if (readMasses&&readNonbond) {
    ff->setFFNonbondDb(nonbondDb);
  }
  return ff;
}


ForceField_sp ReadAmberParameters_O::parseAmberFormattedForceFieldOrFrcMod(core::T_sp strm) {
  core::T_mv test_mv = this->determineParmSetFrcModType(strm);
  core::T_sp masses = test_mv.second();
  core::T_sp nonbond = test_mv.third();
  if (test_mv.nilp()) {
    return this->parseAmberFormattedForceField(strm);
  } else {
    if (masses.nilp() != nonbond.nilp()) {
      SIMPLE_ERROR(BF("Modified force field files must contain both a MASS and a NONB entry, or neither"));
    }
    return this->parseFrcModFile(strm);
  }
}







};
