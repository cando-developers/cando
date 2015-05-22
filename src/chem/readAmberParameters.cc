#define	DEBUG_LEVEL_FULL

#include "readAmberParameters.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "core/numerics.h"
#include "chemInfo.h"
#include "units.h"
#include "ffTypesDb.h"
#include "forceField.h"
#include "core/wrappers.h"


namespace chem {





void	ReadAmberParameters_O::initialize()
{
    this->Base::initialize();
    this->_Types = _Nil<FFTypesDb_O>();
    this->_ForceField = _Nil<ForceField_O>();
}


void	ReadAmberParameters_O::readTypes(const string& fileName)
{_G();
    this->_Types = this->parseTypeRules(fileName);
}


void	ReadAmberParameters_O::readParameters(const string& fileName)
{
    this->_ForceField = this->parseAmberFormattedForceField(fileName);
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




FFTypesDb_sp ReadAmberParameters_O::parseTypeRules(const string& fileName)
{_G();
    GC_ALLOCATE(WildElementDict_O, wildCardElementDictionary );
    GC_ALLOCATE(FFTypesDb_O, ffTypesDb );
    uint lineno = 1;
    vector< pair< uint, string> > entries;
    ifstream fin;
    fin.open(fileName.c_str(),ios::in);
    if ( !fin.is_open() )
    {
	SIMPLE_ERROR(BF("Could not open file: "+fileName ));
    }
    char	buffer[1024];
    while ( fin.good() )
    {
	fin.getline(buffer,sizeof(buffer)-1);
	LOG(BF("Read line(%s)") % buffer  );
	string line = buffer;
        if ( line.substr(0,8) == "WILDATOM" )
	{
            vector<string> names = core::split(line," \t");
            string wildName = names[1];
            wildCardElementDictionary->addWildName(wildName);
            for ( vector<string>::iterator it=names.begin()+1; it!=names.end(); it++ )
	    {
                wildCardElementDictionary->addWildNameMap(wildName,*it);
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
    fin.close();
    return ffTypesDb;
}




FFNonbondDb_sp ReadAmberParameters_O::parseFirstNonbondDb(ifstream& fin)
{_G();
    GC_ALLOCATE(FFNonbondDb_O, ffNonbondDb );
    bool done = false;
    char buffer[1024];
    while ( not done )
    {
	fin.getline(buffer,sizeof(buffer)-1);
	string line = buffer;
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
            ffNonbond->setType(typeName);
            ffNonbond->setMass(mass);
            ffNonbond->setPolarizability(polarizability);
            ffNonbondDb->add(ffNonbond);
	}
    }
    fin.getline(buffer,sizeof(buffer)-1);
    return ffNonbondDb;
}
 
FFStretchDb_sp ReadAmberParameters_O::parseStretchDb(ifstream& fin)
{_G();
    // read stretch terms 
    GC_ALLOCATE(FFStretchDb_O, ffStretchDb );
    bool done = false;
    char buffer[1024];
    while ( not done )
    {
	fin.getline(buffer,sizeof(buffer)-1);
	string line = buffer;
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
            ffStretch->_Type1 = type1Name;
            ffStretch->_Type2 = type2Name;
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


FFAngleDb_sp ReadAmberParameters_O::parseAngleDb(ifstream& fin)
{_G();
    //
    // read angle terms 
    //
    GC_ALLOCATE(FFAngleDb_O, ffAngleDb );
    bool done = false;
    char buffer[1024];
    while ( not done )
    {
	fin.getline(buffer,sizeof(buffer)-1);
	string line = buffer;
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
            ffAngle->_Type1 = t1;
            ffAngle->_Type2 = t2;
            ffAngle->_Type3 = t3;
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



FFPtorDb_sp ReadAmberParameters_O::parsePtorDb(ifstream& fin)
{_G();
    //
    // read Ptor terms 
    //
    GC_ALLOCATE(FFPtorDb_O, ffPtorDb );
    bool done = false;
    char buffer[1024];
    while (not done )
    {
	fin.getline(buffer,sizeof(buffer)-1);
	string line = buffer;
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
            ffPtor->setTypes(t1,t2,t3,t4);
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





FFItorDb_sp ReadAmberParameters_O::parseItorDb(ifstream& fin)
{_G();
    //
    // read Itor terms 
    //
    GC_ALLOCATE(FFItorDb_O, ffItorDb );
    bool done = false;
    char buffer[1024];
    while (not done )
    {
	fin.getline(buffer,sizeof(buffer)-1);
	string line = buffer;
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
            ffItor->setTypes(t1,t2,t3,t4);
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
void ReadAmberParameters_O::parseNonbondDb(ifstream& fin, FFNonbondDb_sp ffNonbondDb)
{_G();
    char buffer[1024];
    bool done = false;
    fin.getline(buffer,sizeof(buffer)-1);
    string kind = buffer;
    vector<string>parts = core::split(kind);
    string label = parts[0];
    string kindnb = parts[1];
    if ( kindnb != "RE" )
    {
	SIMPLE_ERROR(BF("Nonbond parameters must be of kindnb=RE this file has kindnb="+kindnb));
    }
    while ( not done )
    {
        fin.getline(buffer,sizeof(buffer)-1);
        string line = buffer;
        if ( line.size() == 0 )
	{
            done = true;
	} else
	{
	    LOG(BF("Parsing line|%s|") % line.c_str()  );
	    FFNonbond_sp ffNonbond;
            string type = core::trimWhiteSpace(line.substr(0,4));
	    if ( ffNonbondDb->hasType(type) )
	    {
		ffNonbond = ffNonbondDb->findType(type);
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




ForceField_sp ReadAmberParameters_O::parseAmberFormattedForceField(const string& parametersFileName)
{_OF();
    char buffer[1024];
    ifstream fin;
    fin.open(parametersFileName.c_str(),ios::in);
    if ( !fin.is_open() )
    {
	SIMPLE_ERROR(BF("Could not open file: "+parametersFileName ));
    }
    fin.getline(buffer,sizeof(buffer)-1);
    string title = buffer;
    FFNonbondDb_sp ffNonbondsDb = this->parseFirstNonbondDb(fin);
    FFStretchDb_sp ffStretchesDb = this->parseStretchDb(fin);
    FFAngleDb_sp ffAnglesDb = this->parseAngleDb(fin);
    FFPtorDb_sp ffPtorsDb = this->parsePtorDb(fin);
    FFItorDb_sp ffItorsDb = this->parseItorDb(fin);
    fin.getline(buffer,sizeof(buffer)-1); // skp 10-12 hbond
    fin.getline(buffer,sizeof(buffer)-1); // blank
    fin.getline(buffer,sizeof(buffer)-1); // skip equivalence
    fin.getline(buffer,sizeof(buffer)-1); // blank
    this->parseNonbondDb(fin,ffNonbondsDb);
    LOG(BF("All parameters read fine") );
    GC_ALLOCATE(ForceField_O, ff );
    ff->setTitle(title);
    ff->setFFStretchDb(ffStretchesDb);
    ff->setFFTypeDb( _Nil<FFTypesDb_O>());
    ff->setFFAngleDb(ffAnglesDb);
    ff->setFFPtorDb(ffPtorsDb);
    ff->setFFItorDb(ffItorsDb);
    ff->setFFNonbondDb(ffNonbondsDb);
    //
    // parameters from  Antechamberpaper
    //
    ff->_Angles->addZConstant("H", 0.784);
    ff->_Angles->addZConstant("C", 1.183);
    ff->_Angles->addZConstant("N", 1.212);
    ff->_Angles->addZConstant("O", 1.219);
    ff->_Angles->addZConstant("F", 1.166);
    ff->_Angles->addZConstant("Cl", 1.272);
    ff->_Angles->addZConstant("Br", 1.378);
    ff->_Angles->addZConstant("I", 1.398);
    ff->_Angles->addZConstant("P", 1.620);
    ff->_Angles->addZConstant("S", 1.280);
    ff->_Angles->addCConstant("C", 1.339);
    ff->_Angles->addCConstant("N", 1.300);
    ff->_Angles->addCConstant("O", 1.249);
    ff->_Angles->addCConstant("P", 0.906);
    ff->_Angles->addCConstant("S", 1.448);

    //
    // Parameters from Antechamber for estimating stretching force constants
    //
    ff->_Stretches->addEstimateStretch( "H", "H", 0.738, 4.661 );
    ff->_Stretches->addEstimateStretch( "C", "C", 1.526, 7.643 );
    ff->_Stretches->addEstimateStretch( "N", "N", 1.441, 7.634 );
    ff->_Stretches->addEstimateStretch( "O", "O", 1.460, 7.561 );
    ff->_Stretches->addEstimateStretch( "F", "F", 1.406, 7.358 );
    ff->_Stretches->addEstimateStretch( "Cl", "Cl", 2.031, 8.648 );
    ff->_Stretches->addEstimateStretch( "Br", "Br", 2.337, 9.012 );
    ff->_Stretches->addEstimateStretch( "I", "I", 2.836, 9.511 );
    ff->_Stretches->addEstimateStretch( "P", "P", 2.324, 8.805 );
    ff->_Stretches->addEstimateStretch( "S", "S", 2.038, 8.316 );
    ff->_Stretches->addEstimateStretch( "H", "C", 1.090, 6.217 );
    ff->_Stretches->addEstimateStretch( "H", "N", 1.010, 6.057 );
    ff->_Stretches->addEstimateStretch( "H", "O", 0.960, 5.794 );
    ff->_Stretches->addEstimateStretch( "H", "F", 0.920, 5.600 );
    ff->_Stretches->addEstimateStretch( "H", "Cl", 1.280, 6.937 );
    ff->_Stretches->addEstimateStretch( "H", "Br", 1.410, 7.301 );
    ff->_Stretches->addEstimateStretch( "H", "I", 1.600, 7.802 );
    ff->_Stretches->addEstimateStretch( "H", "P", 1.410, 7.257 );
    ff->_Stretches->addEstimateStretch( "H", "S", 1.340, 7.018 );
    ff->_Stretches->addEstimateStretch( "C", "N", 1.470, 7.504 );
    ff->_Stretches->addEstimateStretch( "C", "O", 1.440, 7.347 );
    ff->_Stretches->addEstimateStretch( "C", "F", 1.370, 7.227 );
    ff->_Stretches->addEstimateStretch( "C", "Cl", 1.800, 8.241 );
    ff->_Stretches->addEstimateStretch( "C", "Br", 1.940, 8.478 );
    ff->_Stretches->addEstimateStretch( "C", "I", 2.160, 8.859 );
    ff->_Stretches->addEstimateStretch( "C", "P", 1.830, 8.237 );
    ff->_Stretches->addEstimateStretch( "C", "S", 1.820, 8.117 );
    ff->_Stretches->addEstimateStretch( "N", "O", 1.420, 7.526 );
    ff->_Stretches->addEstimateStretch( "N", "F", 1.420, 7.475 );
    ff->_Stretches->addEstimateStretch( "N", "Cl", 1.750, 8.266 );
    ff->_Stretches->addEstimateStretch( "N", "Br", 1.930, 8.593 );
    ff->_Stretches->addEstimateStretch( "N", "I", 2.120, 8.963 );
    ff->_Stretches->addEstimateStretch( "N", "P", 1.720, 8.212 );
    ff->_Stretches->addEstimateStretch( "N", "S", 1.690, 8.073 );
    ff->_Stretches->addEstimateStretch( "O", "F", 1.410, 7.375 );
    ff->_Stretches->addEstimateStretch( "O", "Cl", 1.700, 8.097 );
    ff->_Stretches->addEstimateStretch( "O", "Br", 1.790, 8.276 );
    ff->_Stretches->addEstimateStretch( "O", "I", 2.110, 8.854 );
    ff->_Stretches->addEstimateStretch( "O", "P", 1.640, 7.957 );
    ff->_Stretches->addEstimateStretch( "O", "S", 1.650, 7.922 );
    ff->_Stretches->addEstimateStretch( "F", "Cl", 1.648, 7.947 );
    ff->_Stretches->addEstimateStretch( "Cl", "I", 2.550, 9.309 );
    ff->_Stretches->addEstimateStretch( "Br", "I", 2.671, 9.380 );
    ff->_Stretches->addEstimateStretch( "F", "P", 1.500, 7.592 );
    ff->_Stretches->addEstimateStretch( "F", "S", 1.580, 7.733 );
    ff->_Stretches->addEstimateStretch( "Cl", "P", 2.040, 8.656 );
    ff->_Stretches->addEstimateStretch( "Cl", "S", 2.030, 8.619 );
    ff->_Stretches->addEstimateStretch( "Br", "P", 2.240, 8.729 );
    ff->_Stretches->addEstimateStretch( "Br", "S", 2.210, 8.728 );
    ff->_Stretches->addEstimateStretch( "I", "P", 2.490, 9.058 );
    ff->_Stretches->addEstimateStretch( "I", "S", 2.560, 9.161 );
    ff->_Stretches->addEstimateStretch( "P", "S", 2.120, 8.465 );
    fin.close();
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




