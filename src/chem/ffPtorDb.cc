       
#define	TURN_DEBUG_OFF

//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffPtor.cc
 *
 *	Maintains a database of ptor terms
 */
#include <cando/chem/ffPtorDb.h>
//#include "core/archiveNode.h"
#include <clasp/core/str.h>
#include <clasp/core/wrappers.h>
#include <cando/chem/units.h>

namespace chem
{
    const int FFPtor_O::MaxPeriodicity;

    REGISTER_CLASS(chem, FFPtor_O );
    REGISTER_CLASS(chem, FFPtorDb_O );


    void	FFPtor_O::archiveBase(core::ArchiveP node)
{
    node->attribute( "type1", this->_T1 );
    node->attribute( "type2", this->_T2 );
    node->attribute( "type3", this->_T3 );
    node->attribute( "type4", this->_T4 );
    for ( int i=0; i<MaxPeriodicity; i++ )
    {
	node->attributeIfDefined( KW((BF("v%d")%i).str()), this->_hasPeriodicity[i], this->_Vs_kJ[i] );
	node->attributeIfDefined( KW((BF("ph%d")%i).str()), this->_hasPeriodicity[i],  this->_PhaseRads[i] );
    }
}

#ifdef CONSPACK
    void	FFPtorDb_O::archiveBase(core::ArchiveP node)
    {
	this->FFBaseDb_O::archiveBase(node);
	node->attributeVector0<FFPtor_O>(KW("ptors"),this->_Terms );
	node->attributeStringMap<FFPtor_O>(KW("map"),this->_Lookup );
    }
#endif



void	FFPtorDb_O::initialize()
{
    this->Base::initialize();
    this->_Lookup.clear();
    this->_Terms.clear();
}






void    FFPtorDb_O::cantFind(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{_OF();
    stringstream ss;
    ss << "Can't find ptor term for ("<<t1<<")-("<<t2<<")-("<<t3<<")-("<<t4<<")";
    SIMPLE_ERROR(BF("%s")%ss.str());
}


void FFPtor_O::initialize()
{
    this->Base::initialize();
    int i;
    this->_T1 = _Nil<core::Symbol_O>();
    this->_T2 = _Nil<core::Symbol_O>();
    this->_T3 = _Nil<core::Symbol_O>();
    this->_T4 = _Nil<core::Symbol_O>();
    for ( i=0; i<MaxPeriodicity; i++ )
    {
	this->_hasPeriodicity[i] = false;
	this->_Vs_kJ[i] = 0.0;
	this->_PhaseRads[i] = 0.0;
    }
}



void    FFPtor_O::setTypes(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
    this->_T1 = t1;
    this->_T2 = t2;
    this->_T3 = t3;
    this->_T4 = t4;
}

bool	FFPtor_O::hasPeriodicity(int idx ) const
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal periodicity for hasPeriodicity"));
    }
    return this->_hasPeriodicity[idx-1];
}




double  FFPtor_O::getV_kJ(int idx) const
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal periodicity for getV"));
    }
    return this->_Vs_kJ[idx-1];
}

void    FFPtor_O::setV_kJ(int idx, double val)
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal periodicity for setV"));
    }
    this->_hasPeriodicity[idx-1] = true;
    this->_Vs_kJ[idx-1] = val;
}

double FFPtor_O::getV_kCal(int idx) const
{_OF();
    return kJ_to_kCal(this->getV_kJ(idx));
}

void FFPtor_O::setV_kCal(int idx, double val)
{_OF();
    this->setV_kJ(idx,kCal_to_kJ(val));
}


double  FFPtor_O::getPhaseRad(int idx) const
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal periodicity[%d] for getPhaseRad - must be in [1,%d]") % idx % MaxPeriodicity );
    }
    return this->_PhaseRads[idx-1];
}

void    FFPtor_O::setPhaseRad(int idx, double val)
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal periodicity for setPhaseRad"));
    }
    this->_PhaseRads[idx-1] = val;
}


void    FFPtor_O::mergeWith(FFPtor_sp ptor)
{
int     i;
    for ( i=0; i<MaxPeriodicity; i++ ) {
	if ( ptor->_hasPeriodicity[i] ) {
            this->_hasPeriodicity[i] = ptor->_hasPeriodicity[i];
            this->_Vs_kJ[i] = ptor->_Vs_kJ[i];
            this->_PhaseRads[i] = ptor->_PhaseRads[i];
        }
    }
}

#define keyString(t1,t2,t3,t4) chemkw_intern(t1->symbolName()->get()+"-"+t2->symbolName()->get()+"-"+t3->symbolName()->get()+"-"+t4->symbolName()->get())

void	FFPtorDb_O::add( FFPtor_sp ptor )
{
core::Symbol_sp key;
FFPtor_sp        ptorOld;
    if (this->hasExactTerm(ptor->_T1,ptor->_T2,ptor->_T3,ptor->_T4) ){
        ptorOld = this->findExactTerm(ptor->_T1,ptor->_T2,ptor->_T3,ptor->_T4);
        ptorOld->mergeWith(ptor);
    } else {
        this->_Terms.push_back(ptor);
        key = keyString(ptor->_T1,ptor->_T2,ptor->_T3,ptor->_T4);
        this->_Lookup.set(key,ptor);
        key = keyString(ptor->_T4,ptor->_T3,ptor->_T2,ptor->_T1);
        this->_Lookup.set(key,ptor);
    }
}

bool FFPtorDb_O::hasExactTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
core::Symbol_sp key;
FFPtor_sp        ptor;
    key = keyString(t1,t2,t3,t4);
    if ( this->_Lookup.count(key)!=0 ) return true;
    key = keyString(t4,t3,t2,t1);
    if ( this->_Lookup.count(key)!=0 ) return true;
    return false;
}

FFPtor_sp FFPtorDb_O::findExactTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
core::Symbol_sp key;
FFPtor_sp        ptor;
    key = keyString(t1,t2,t3,t4);
    if ( this->_Lookup.count(key)!=0 ) {
        return this->_Lookup.get(key);
    }
    key = keyString(t4,t3,t2,t1);
    if ( this->_Lookup.count(key)!=0 ) {
        return this->_Lookup.get(key);
    }
    ptor = _Nil<FFPtor_O>();
    return ptor;
}

FFPtor_sp FFPtorDb_O::findBestTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
core::Symbol_sp key,bl;
FFPtor_sp        ptor;
    key = keyString(t1,t2,t3,t4);
    if ( this->_Lookup.count(key)!=0 ) {
        return this->_Lookup.get(key);
    }
    key = keyString(t4,t3,t2,t1);
    if ( this->_Lookup.count(key)!=0 ) {
        return this->_Lookup.get(key);
    }
    bl = _Nil<core::Symbol_O>();
    key = keyString(bl,t2,t3,bl);
    if ( this->_Lookup.count(key)!=0 ) {
        return this->_Lookup.get(key);
    }
    key = keyString(bl,t3,t2,bl);
    if ( this->_Lookup.count(key)!=0 ) {
        return this->_Lookup.get(key);
    }
    ptor = _Nil<FFPtor_O>();
    return ptor;
}


bool    FFPtorDb_O::hasBestTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
string          key,bl;
FFPtor_sp        ptor;
    ptor = this->findBestTerm(t1,t2,t3,t4);
    return ptor.notnilp();
}



string	FFPtor_O::levelDescription()
{
stringstream	desc;
    desc << this->FFParameter_O::levelDescription();
    desc << " ptor between types "<<this->_T1<< " "
		<< this->_T2 << " "
		<< this->_T3 << " "
		<< this->_T4;
    return desc.str();
}

};

