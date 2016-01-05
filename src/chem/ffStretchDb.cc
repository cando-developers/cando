       
#define	 DEBUG_LEVEL_NONE	// _FULL/_SILENT/_OFF

//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffStretch.cc
 *
 *	Maintains a database of stretch types
 */
#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/units.h>
#include <cando/chem/symbolTable.h>
#include <clasp/core/wrappers.h>



namespace chem
{


    EXPOSE_CLASS(chem, FFStretch_O );
    EXPOSE_CLASS(chem, FFStretchDb_O );




#ifdef XML_ARCHIVE
adapt::QDomNode_sp	EstimateStretch::asXml()
{
    adapt::QDomNode_sp	node;
    node = adapt::QDomNode_O::create("EstimateStretch");
    node->addAttributeSymbol("ti",this->_ti);
    node->addAttributeSymbol("tj",this->_tj);
    node->addAttributeDoubleScientific("rij",this->_rij);
    node->addAttributeDoubleScientific("lnKij",this->_lnKij);
    return node;
}


void	EstimateStretch::parseFromXml(adapt::QDomNode_sp node)
{
    this->_ti = node->getAttributeString("ti");
    this->_tj = node->getAttributeString("tj");
    this->_rij = node->getAttributeDouble("rij");
    this->_lnKij = node->getAttributeDouble("lnKij");
}
#endif



    void FFStretchDb_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<FFStretchDb_O>()
	    .def("findTerm",&FFStretchDb_O::findTerm)
	    ;
    }

    void FFStretchDb_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,FFStretchDb,"","",_lisp)
	    .def("findTerm",&FFStretchDb_O::findTerm)
	    ;
#endif
    }


void	FFStretchDb_O::initialize()
{
    this->Base::initialize();
    this->clearEstimateStretch();
    this->_Lookup.clear();
    this->_Terms.clear();
}





    core::Symbol_sp stretchKey(core::Symbol_sp t1
                               , core::Symbol_sp t2 )
    {
        return chemkw_intern(t1->symbolName()->get()+"-"
                           + t2->symbolName()->get());
    }


void    FFStretchDb_O::add(FFStretch_sp term)
{
    core::Symbol_sp key;
    this->_Terms.push_back(term);
    key = stretchKey(term->_Type1,term->_Type2);        // forwards
    this->_Lookup.set(key,term);
    key = stretchKey(term->_Type2,term->_Type1);        // backwards
    this->_Lookup.set(key,term);
}


CL_LISPIFY_NAME("findTerm");
CL_DEFMETHOD FFStretch_sp	FFStretchDb_O::findTerm(chem::Atom_sp a1, chem::Atom_sp a2)
{_G();
FFStretch_sp			match;
 core::Symbol_sp                          key;
core::Symbol_sp				t1, t2;

    t1 = a1->getType();
    t2 = a2->getType();
    LOG(BF("Looking for stretch between types (%s)-(%s)") % t1.c_str() % t2.c_str() );
    key = stretchKey(t1,t2);
    if ( this->_Lookup.count(key) != 0 ) {
        return this->_Lookup.get(key);
    }
    key = stretchKey(t2,t1);
    if ( this->_Lookup.count(key) != 0 ) {
        return this->_Lookup.get(key);
    }
#if 1
    stringstream ss;
    ss << "PARAMETER_PROBLEM Can't find stretch term for types(" << t1 << ")-("<<t2<<")";
    ss << " atoms("<<a1->description()<<") - (" << a2->description()<<")";
    SIMPLE_ERROR(BF("%s")%ss.str());
#endif
    match = _Nil<FFStretch_O>();
    return match;
}


#if 0 //[
bool    FFStretchDb_O::hasTerm(string t1, string t2 )
{
FFStretch_sp			match;
string                          key;
    key = stretchKey(t1,t2);
    if ( this->_Lookup.count(key) != 0 ) {
        return true;
    }
    key = stretchKey(t2,t1);
    if ( this->_Lookup.count(key) != 0 ) {
        return true;
    }
    return false;
}


void    FFStretchDb_O::cantFind(string t1, string t2 )
{
    stringstream ss;
    ss << "Can't find stretch term for types(" << t1 << ")-("<<t2<<")";
    SIMPLE_ERROR(BF("%s")%ss.str());
}

#endif //]


void	FFStretchDb_O::clearEstimateStretch()
{
    this->_EstimateStretch.clear();
}

void	FFStretchDb_O::_addEstimateStretch(const EstimateStretch& esi)
{
  core::Symbol_sp		tt;
EstimateStretch	es;
    es = esi;
    if ( core::SymbolComparer::order(es._ti,es._tj) == 1 ) {
      es._ti.swap(es._tj);
#if 0
      tt = es._ti;
	es._ti = es._tj;
	es._tj = tt;
#endif
    }
    tt = chemkw_intern(es._ti->symbolNameAsString()+"-"+es._tj->symbolNameAsString());
    this->_EstimateStretch[tt] = es;
}




void	FFStretchDb_O::addEstimateStretch(core::Symbol_sp ti, core::Symbol_sp tj, double rij, double lnKij ) {
EstimateStretch	es;
    es._ti = ti;
    es._tj = tj;
    es._rij = rij;
    es._lnKij = lnKij;
    this->_addEstimateStretch(es);
}



void	FFStretch_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
  node->field( INTERN_(kw,type1), this->_Type1 );
  node->field( INTERN_(kw,type2), this->_Type2 );
  node->/*pod_*/field( INTERN_(kw,r0), this->_R0_Nanometer );
  node->/*pod_*/field( INTERN_(kw,kb), this->_Kb_kJPerNanometerSquared );
}

void	FFStretchDb_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
  node->field( INTERN_(kw,stretches),this->_Terms );
  node->/*pod_*/field( INTERN_(kw,map),this->_Lookup );
}

    void FFStretch_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<FFStretch_O>()
	    ;
    }

    void FFStretch_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,FFStretch,"","",_lisp)
	    ;
#endif
    }





void FFStretch_O::initialize()
{
    this->Base::initialize();
    this->_Type1 = _Nil<core::T_O>();
    this->_Type2 = _Nil<core::T_O>();
    this->_R0_Nanometer = 0.0;
    this->_Kb_kJPerNanometerSquared = 0.0;
//    this->_K3 = 0.0;
//    this->_K4 = 0.0;
//    this->_Bci = 0.0;
}



double FFStretch_O::getR0_Nanometer() const
{
    return this->_R0_Nanometer;
}

double FFStretch_O::getR0_Angstrom() const
{
    return this->_R0_Nanometer*10.0;
}

double FFStretch_O::getKb_kJPerNanometerSquared() const
{
    return this->_Kb_kJPerNanometerSquared;
}


void FFStretch_O::setKb_kJPerNanometerSquared(double val)
{
    this->_Kb_kJPerNanometerSquared = val;
}


void FFStretch_O::setR0_Nanometer(double val)
{
    this->_R0_Nanometer = val;
}


double FFStretch_O::getKb_kCalPerAngstromSquared() const
{
    return kJPerNanometerSquared_to_kCalPerAngstromSquared(this->_Kb_kJPerNanometerSquared);
}





string	FFStretch_O::levelDescription()
{
stringstream	desc;
    desc << this->FFParameter_O::levelDescription();
    desc << " stretch between types "<<this->_Type1<< " "
		<< this->_Type2;
    return desc.str();
}


};

