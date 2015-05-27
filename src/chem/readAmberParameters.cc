#define	DEBUG_LEVEL_FULL

#include <clasp/core/foundation.h>
#include <cando/chem/readAmberParameters.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/numerics.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/str.h>
#include <cando/chem/chemInfo.h>
#include <cando/chem/units.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/forceField.h>
#include <clasp/core/wrappers.h>


namespace chem {





void	ReadAmberParameters_O::initialize()
{
    this->Base::initialize();
    this->_Types = _Nil<FFTypesDb_O>();
    this->_ForceField = _Nil<ForceField_O>();
}


void	ReadAmberParameters_O::readTypes(core::T_sp fin)
{_G();
    this->_Types = this->parseTypeRules(fin);
}


void	ReadAmberParameters_O::readParameters(core::T_sp fin)
{
    this->_ForceField = this->parseAmberFormattedForceField(fin);
}



ForceField_sp ReadAmberParameters_O::getForceField()
{_G();
    if ( this->_Types.nilp() ) {
	SIMPLE_ERROR(BF("In getForceField this->_Types cannot be nil"));
    }
    if ( this->_ForceField.notnilp() )
    {
	this->_ForceField->setFFTypeDb(this->_Types);
    } else
    if ( this->_ForceField.nilp() )
    {
	GC_ALLOCATE(ForceField_O, temp );
	this->_ForceField = temp;
    }
    this->_ForceField->setFFTypeDb(this->_Types);
    return this->_ForceField;
}




FFTypesDb_sp ReadAmberParameters_O::parseTypeRules(core::T_sp fin)
{
  GC_ALLOCATE(WildElementDict_O, wildCardElementDictionary );
  GC_ALLOCATE(FFTypesDb_O, ffTypesDb );
  uint lineno = 1;
  vector< pair< uint, string> > entries;
  while ( true ) {
    core::T_mv mv_line = core::cl_read_line(fin,_Nil<T_O>());
    if ( mv_line.nilp() ) break;
    string line = mv_line.as<core::Str_O>()->get();
    LOG(BF("Read line(%s)") % line  );
    if ( line.substr(0,8) == "WILDATOM" ) {
      vector<string> names = core::split(line," \t");
      string wildName = names[1];
      wildCardElementDictionary->addWildName(_lisp->intern(wildName,ChemKwPkg));
      for ( vector<string>::iterator it=names.begin()+1; it!=names.end(); it++ )
      {
        core::Symbol_sp wildSymbol = _lisp->intern(wildName,ChemKwPkg);
        core::Symbol_sp wildOtherSymbol = _lisp->intern(*it,ChemKwPkg);
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
    if ( typeRule->compileSucceeded() )
    {
      ffTypesDb->add(typeRule);
    } else
    {
      SIMPLE_ERROR(BF("Antechamber compile failed on: "+ei->second+"\n"+typeRule->compilerMessage() ));
    }
  }
  return ffTypesDb;
}




FFNonbondDb_sp ReadAmberParameters_O::parseFirstNonbondDb(core::T_sp fin)
{_G();
  GC_ALLOCATE(FFNonbondDb_O, ffNonbondDb );
  bool done = false;
  while ( not done )
  {
    core::T_mv ol = core::cl_read_line(fin);
    string line = gc::As<core::Str_sp>(ol)->get();
    if ( line.size() == 0 )
    {
      done = true;
    } else {
      LOG(BF("Parsing line|%s|") % line.c_str()  );
      vector<string> parts = core::split(line," \t");
      string typeName = parts[0];
      double mass = atof(parts[1].c_str());
      double polarizability = atof(parts[2].c_str());
      GC_ALLOCATE(FFNonbond_O, ffNonbond );
      ffNonbond->setType(_lisp->intern(typeName,ChemKwPkg));
      ffNonbond->setMass(mass);
      ffNonbond->setPolarizability(polarizability);
      ffNonbondDb->add(ffNonbond);
    }
  }
  core::cl_read_line(fin); // blank line
  return ffNonbondDb;
}
 
FFStretchDb_sp ReadAmberParameters_O::parseStretchDb(core::T_sp fin)
{_G();
    // read stretch terms 
    GC_ALLOCATE(FFStretchDb_O, ffStretchDb );
    bool done = false;
    while ( not done )
    {
      string line = core::cl_read_line(fin).as<core::Str_O>()->get();
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
{_G();
    //
    // read angle terms 
    //
    GC_ALLOCATE(FFAngleDb_O, ffAngleDb );
    bool done = false;
    while ( not done )
    {
      string line = core::cl_read_line(fin).as<core::Str_O>()->get();
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
            ffAngle->_Type1 = _lisp->intern(t1,ChemKwPkg);
            ffAngle->_Type2 = _lisp->intern(t2,ChemKwPkg);
            ffAngle->_Type3 = _lisp->intern(t3,ChemKwPkg);
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
{_G();
    //
    // read Ptor terms 
    //
    GC_ALLOCATE(FFPtorDb_O, ffPtorDb );
    bool done = false;
    while (not done )
    {
      string line = core::cl_read_line(fin).as<core::Str_O>()->get();
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
            if ( t1 == "X" ) t1 = "";
            if ( t4 == "X" ) t4 = "";
            GC_ALLOCATE(FFPtor_O, ffPtor );
            core::Symbol_sp st1 = _lisp->intern(t1,ChemKwPkg);
            core::Symbol_sp st2 = _lisp->intern(t2,ChemKwPkg);
            core::Symbol_sp st3 = _lisp->intern(t3,ChemKwPkg);
            core::Symbol_sp st4 = _lisp->intern(t4,ChemKwPkg);
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
{_G();
    //
    // read Itor terms 
    //
    GC_ALLOCATE(FFItorDb_O, ffItorDb );
    bool done = false;
    while (not done )
    {
      string line = core::cl_read_line(fin).as<core::Str_O>()->get();
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
            if ( t1 == "X" ) t1 = "";
            if ( t2 == "X" ) t2 = "";
            if ( t4 == "X" ) t4 = "";
            GC_ALLOCATE(FFItor_O, ffItor );
            core::Symbol_sp st1 = _lisp->intern(t1,ChemKwPkg);
            core::Symbol_sp st2 = _lisp->intern(t2,ChemKwPkg);
            core::Symbol_sp st3 = _lisp->intern(t3,ChemKwPkg);
            core::Symbol_sp st4 = _lisp->intern(t4,ChemKwPkg);
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



/*!
 * This modifies entries in the ffNonbondDb database
 */
void ReadAmberParameters_O::parseNonbondDb(core::T_sp fin, FFNonbondDb_sp ffNonbondDb)
{_G();
    bool done = false;
    string line = core::cl_read_line(fin).as<core::Str_O>()->get();
    vector<string>parts = core::split(line);
    string label = parts[0];
    string kindnb = parts[1];
    if ( kindnb != "RE" )
    {
	SIMPLE_ERROR(BF("Nonbond parameters must be of kindnb=RE this file has kindnb="+kindnb));
    }
    while ( not done )
    {
      string line = core::cl_read_line(fin).as<core::Str_O>()->get();
      if ( line.size() == 0 )
      {
        done = true;
      } else
      {
        LOG(BF("Parsing line|%s|") % line.c_str()  );
        FFNonbond_sp ffNonbond;
        string type = core::trimWhiteSpace(line.substr(0,4));
        core::Symbol_sp stype = _lisp->intern(type,ChemKwPkg);
        if ( ffNonbondDb->hasType(stype) ) {
          ffNonbond = ffNonbondDb->findType(stype);
        } else
        {
          SIMPLE_ERROR(BF("Could not find type: "+type));
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




ForceField_sp ReadAmberParameters_O::parseAmberFormattedForceField(core::T_sp fin)
{_OF();
    string line = core::cl_read_line(fin).as<core::Str_O>()->get();
    FFNonbondDb_sp ffNonbondsDb = this->parseFirstNonbondDb(fin);
    FFStretchDb_sp ffStretchesDb = this->parseStretchDb(fin);
    FFAngleDb_sp ffAnglesDb = this->parseAngleDb(fin);
    FFPtorDb_sp ffPtorsDb = this->parsePtorDb(fin);
    FFItorDb_sp ffItorsDb = this->parseItorDb(fin);
    core::cl_read_line(fin); // skp 10-12 hbond
    core::cl_read_line(fin); // blank
    core::cl_read_line(fin); // skip equivalence
    core::cl_read_line(fin); // blank
    this->parseNonbondDb(fin,ffNonbondsDb);
    LOG(BF("All parameters read fine") );
    GC_ALLOCATE(ForceField_O, ff );
    ff->setTitle(line);
    ff->setFFStretchDb(ffStretchesDb);
    ff->setFFTypeDb( _Nil<FFTypesDb_O>());
    ff->setFFAngleDb(ffAnglesDb);
    ff->setFFPtorDb(ffPtorsDb);
    ff->setFFItorDb(ffItorsDb);
    ff->setFFNonbondDb(ffNonbondsDb);
    //
    // parameters from  Antechamberpaper
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



    void ReadAmberParameters_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<ReadAmberParameters_O>("makeReadAmberParameters")
	.def("getTypeRules",&ReadAmberParameters_O::getTypeRules)
	.def("getForceField",&ReadAmberParameters_O::getForceField)
	.def("readTypes",&ReadAmberParameters_O::readTypes)
	.def("readParameters",&ReadAmberParameters_O::readParameters)
	;
}

    void ReadAmberParameters_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ReadAmberParameters,"","",_lisp)
	.def("getTypeRules",&ReadAmberParameters_O::getTypeRules)
	.def("getForceField",&ReadAmberParameters_O::getForceField)
	.def("readTypes",&ReadAmberParameters_O::readTypes)
	.def("readParameters",&ReadAmberParameters_O::readParameters)
    ;
#endif
#endif //]
}



    EXPOSE_CLASS(chem,ReadAmberParameters_O);
};




