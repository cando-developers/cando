       
//
// (C) 2004 Christian E. Schafmeister
//

#define	TURN_DEBUG_OFF

/*
 *	ffItor.cc
 *
 *	Maintains a database of itor terms
 */
#include "core/common.h"
#include "core/str.h"
#include "ffItorDb.h"
//#include "core/archiveNode.h"
#include "core/wrappers.h"
#include "units.h"


namespace chem
{

    const int FFItor_O::IMaxPeriodicity;


    REGISTER_CLASS(chem, FFItor_O );
    REGISTER_CLASS(chem, FFItorDb_O );

void	FFItor_O::initialize()
{
    this->Base::initialize();
    int i;
    this->_T1 = _Nil<core::Symbol_O>();
    this->_T2 = _Nil<core::Symbol_O>();
    this->_T3 = _Nil<core::Symbol_O>();
    this->_T4 = _Nil<core::Symbol_O>();
    for ( i=0; i<FFItor_O::IMaxPeriodicity; i++ )
    {
	this->_hasPeriodicity[i] = false;
	this->_Vs_kJ[i] = 0.0;
	this->_PhaseRads[i] = 0.0;
    }
}

    void	FFItor_O::archiveBase(core::ArchiveP node)
    {
	node->attribute( "type1", this->_T1 );
	node->attribute( "type2", this->_T2 );
	node->attribute( "type3", this->_T3 );
	node->attribute( "type4", this->_T4 );
	for ( int i=0; i<FFItor_O::IMaxPeriodicity; i++ )
	{
	    node->attributeIfDefined( KW((BF("v%d")%i).str()), this->_hasPeriodicity[i], this->_Vs_kJ[i]);
	    node->attributeIfDefined( KW((BF("ph%d")%i).str()), this->_hasPeriodicity[i], this->_PhaseRads[i]);
	}
    }


#ifdef CONSPACK
    void	FFItorDb_O::archiveBase(core::ArchiveP node)
{
    this->FFBaseDb_O::archiveBase(node);
    node->attributeVector0<FFItor_O>(KW("ptors"),this->_Terms );
    node->attributeStringMap<FFItor_O>(KW("map"),this->_Lookup );
}
#endif




void    FFItor_O::setTypes(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
    this->_T1 = t1;
    this->_T2 = t2;
    this->_T3 = t3;
    this->_T4 = t4;
}

double  FFItor_O::getV_kJ(int period) const
{_OF();
    if ( period < 1 || period > IMaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal index for getV"));
    }
    return this->_Vs_kJ[period-1];
}

void    FFItor_O::setV_kJ(int period, double val)
{_OF();
    if ( period < 1 || period > IMaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal index for setV"));
    }
    this->_Vs_kJ[period-1] = val;
    this->_hasPeriodicity[period-1] = true;
}


bool FFItor_O::hasPeriodicity(int period) const
{_OF();
    ASSERTF(period>=1 && period<=FFItor_O::IMaxPeriodicity,
	    BF("hasPeriodicityperiod[%d] must be in [1,%d]") % period % IMaxPeriodicity );
    return this->_hasPeriodicity[period-1];
}


double  FFItor_O::getV_kCal(int period) const
{_OF();
    return kJ_to_kCal(this->getV_kJ(period));
}

void    FFItor_O::setV_kCal(int period, double val)
{_OF();
    this->setV_kJ(period,kCal_to_kJ(val));
}


double  FFItor_O::getPhaseRad(int period) const
{_OF();
    if ( period < 1 || period > IMaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal index for getPhaseRad"));
    }
    return this->_PhaseRads[period-1];
}

void    FFItor_O::setPhaseRad(int period, double val)
{_OF();
    if ( period < 1 || period > IMaxPeriodicity ) {
        SIMPLE_ERROR(BF("Illegal index for setPhaseRad"));
    }
    this->_PhaseRads[period-1] = val;
    this->_hasPeriodicity[period-1] = true;
}


void    FFItor_O::mergeWith(FFItor_sp itor)
{
int     i;
    for ( i=1; i<=IMaxPeriodicity; i++ ) {
        if ( itor->hasPeriodicity(i))
	{
	    this->_hasPeriodicity[i-1] = itor->_hasPeriodicity[i-1];
            this->_Vs_kJ[i-1] = itor->_Vs_kJ[i-1];
            this->_PhaseRads[i-1] = itor->_PhaseRads[i-1];
        }
    }
}

#define keyString(t1,t2,t3,t4) chemkw_intern(t1->symbolName()->get()+"-"+t2->symbolName()->get()+"-"+t3->symbolName()->get()+"-"+t4->symbolName()->get())









void	FFItorDb_O::add( FFItor_sp itor )
{
core::Symbol_sp          key;
FFItor_sp        itorOld;
    if (this->hasExactTerm(itor->_T1,itor->_T2,itor->_T3,itor->_T4) ){
        itorOld = this->findExactTerm(itor->_T1,itor->_T2,itor->_T3,itor->_T4);
        itorOld->mergeWith(itor);
    } else {
        this->_Terms.push_back(itor);
        key = keyString(itor->_T1,itor->_T2,itor->_T3,itor->_T4);
#ifdef DEBUG_ON
	if ( itor->_T3=="c" && itor->_T4=="o" ) {
	    LOG(BF("FFItorDb::add adding term with key: %s") % key.c_str()  );
	}
#endif
        this->_Lookup.set(key,itor);
    }
}

bool FFItorDb_O::hasExactTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
core::Symbol_sp         key;
FFItor_sp        itor;
    key = keyString(t1,t2,t3,t4);
    if ( this->_Lookup.count(key)!=0 ) return true;
    return false;
}

FFItor_sp FFItorDb_O::findExactTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
core::Symbol_sp         key;
FFItor_sp        itor;
    key = keyString(t1,t2,t3,t4);
    if ( this->_Lookup.count(key)!=0 ) {
        return this->_Lookup.get(key);
    }
    itor = _Nil<FFItor_O>();
    return itor;
}

FFItor_sp FFItorDb_O::findBestTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{_G();
core::Symbol_sp        key,key1,key2,key3,key4,bl;
FFItor_sp        itor;
    key4 = keyString(t1,t2,t3,t4);
    if ( this->_Lookup.count(key4)!=0 ) {
        itor = this->_Lookup.get(key4);
	key = key4;
	goto DONE;
    }
    bl = _Nil<core::Symbol_O>();
    key3 = keyString(bl,t2,t3,t4);
    if ( this->_Lookup.count(key3)!=0 ) {
        itor = this->_Lookup.get(key3);
	key = key3;
	goto DONE;
    }
    key2 = keyString(bl,bl,t3,t4);
    if ( this->_Lookup.count(key2)!=0 ) {
        itor = this->_Lookup.get(key2);
	key = key2;
	goto DONE;
    }
    key1 = keyString(bl,bl,t3,bl);
    if ( this->_Lookup.count(key1)!=0 ) {
        itor = this->_Lookup.get(key1);
	key = key1;
	goto DONE;
    }
    itor = _Nil<FFItor_O>();
DONE:
#ifdef	DEBUG_ON
    if ( t3=="c" && t4=="o" ) {
	LOG(BF("FFITorDb::findBestTerm for types %s-%s-%s-%s") % t1.c_str() % t2.c_str() % t3.c_str() % t4.c_str()  );
	if ( itor.notnilp() )
	{
	    LOG(BF("Found term with key: %s") % key.c_str()  );
	} else {
	    LOG(BF("Could not find term with key: %s") % key4.c_str()  );
	    LOG(BF("Could not find term with key: %s") % key3.c_str()  );
	    LOG(BF("Could not find term with key: %s") % key2.c_str()  );
	    LOG(BF("Could not find term with key: %s") % key1.c_str()  );
	    {_BLOCK_TRACE("itor dictionary");
	        for ( gctools::SmallMap<string,FFItor_sp>::iterator ii=this->_Lookup.begin();
				ii!=this->_Lookup.end();ii++ ) {
		    LOG(BF("Entry key=%s") % ii->first.c_str()  );
		}
	    }
	}
    }
#endif
    return itor;
}


bool    FFItorDb_O::hasBestTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
string          key,bl;
FFItor_sp        itor;
    itor = this->findBestTerm(t1,t2,t3,t4);
    return itor.notnilp();
}


void    FFItorDb_O::cantFind(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{_OF();
    stringstream ss;
    ss << "Can't find itor term for ("<<t1<<")-("<<t2<<")-("<<t3<<")-("<<t4<<")";
    SIMPLE_ERROR(BF("%s")%ss.str());
}



string	FFItor_O::levelDescription()
{
stringstream	desc;
    desc << this->FFParameter_O::levelDescription();
    desc << " itor between types "<<this->_T1<< " "
		<< this->_T2 << " "
		<< this->_T3 << " "
		<< this->_T4;
    return desc.str();
}





};

