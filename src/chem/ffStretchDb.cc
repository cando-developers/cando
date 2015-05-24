       
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
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>



namespace chem
{


    EXPOSE_CLASS(chem, FFStretch_O );
    EXPOSE_CLASS(chem, FFStretchDb_O );





geom::QDomNode_sp	EstimateStretch::asXml(core::Lisp_sp env)
{
    geom::QDomNode_sp	node;
    node = geom::QDomNode_O::create(env,"EstimateStretch");
    node->addAttributeString("ti",this->_ti);
    node->addAttributeString("tj",this->_tj);
    node->addAttributeDoubleScientific("rij",this->_rij);
    node->addAttributeDoubleScientific("lnKij",this->_lnKij);
    return node;
}


void	EstimateStretch::parseFromXml(geom::QDomNode_sp node)
{
    this->_ti = node->getAttributeString("ti");
    this->_tj = node->getAttributeString("tj");
    this->_rij = node->getAttributeDouble("rij");
    this->_lnKij = node->getAttributeDouble("lnKij");
}




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
        return chem_intern(t1->symbolName()->get()+"-"
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


FFStretch_sp	FFStretchDb_O::findTerm(chem::Atom_sp a1, chem::Atom_sp a2)
{_G();
FFStretch_sp			match;
string                          key;
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
string		tt;
EstimateStretch	es;
    es = esi;
    if ( es._ti>es._tj ) {
	tt = es._ti;
	es._ti = es._tj;
	es._tj = tt;
    }
    tt = es._ti+"-"+es._tj;
    this->_EstimateStretch[tt] = es;
}




void	FFStretchDb_O::addEstimateStretch(string ti, string tj, double rij, double lnKij ) {
string		tt;
EstimateStretch	es;
    es._ti = ti;
    es._tj = tj;
    es._rij = rij;
    es._lnKij = lnKij;
    tt = ti+"-"+tj;
    this->_addEstimateStretch(es);
}



    void	FFStretch_O::archiveBase(core::ArchiveP node)
    {
	node->attribute( "type1", this->_Type1 );
	node->attribute( "type2", this->_Type2 );
	node->attribute( "r0", this->_R0_Nanometer );
	node->attribute( "kb", this->_Kb_kJPerNanometerSquared );
    }


    void	FFStretchDb_O::archiveBase(core::ArchiveP node)
    {
	this->FFBaseDb_O::archiveBase(node);
	node->attributeVector0<FFStretch_O>(KW("stretches"),this->_Terms );
	node->attributeStringMap<FFStretch_O>(KW("map"),this->_Lookup );
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
    this->_Type1 = "";
    this->_Type2 = "";
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
