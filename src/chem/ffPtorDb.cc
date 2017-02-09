/*
    File: ffPtorDb.cc
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
#include <clasp/core/array.h>
#include <clasp/core/wrappers.h>
#include <cando/chem/units.h>

namespace chem
{
    const int FFPtor_O::MaxPeriodicity;



void	FFPtor_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,type1), this->_T1 );
  node->field( INTERN_(kw,type2), this->_T2 );
  node->field( INTERN_(kw,type3), this->_T3 );
  node->field( INTERN_(kw,type4), this->_T4 );
  ASSERT(MaxPeriodicity==6);
  node->field_if_defined(INTERN_(kw,v0) ,this->_hasPeriodicity[0],this->_Vs_kJ[0]);
  node->field_if_defined(INTERN_(kw,ph0),this->_hasPeriodicity[0],this->_PhaseRads[0]);
  node->field_if_defined(INTERN_(kw,v1) ,this->_hasPeriodicity[1],this->_Vs_kJ[1]);
  node->field_if_defined(INTERN_(kw,ph1),this->_hasPeriodicity[1],this->_PhaseRads[1]);
  node->field_if_defined(INTERN_(kw,v2) ,this->_hasPeriodicity[2],this->_Vs_kJ[2]);
  node->field_if_defined(INTERN_(kw,ph2),this->_hasPeriodicity[2],this->_PhaseRads[2]);
  node->field_if_defined(INTERN_(kw,v3) ,this->_hasPeriodicity[3],this->_Vs_kJ[3]);
  node->field_if_defined(INTERN_(kw,ph3),this->_hasPeriodicity[3],this->_PhaseRads[3]);
  node->field_if_defined(INTERN_(kw,v4) ,this->_hasPeriodicity[4],this->_Vs_kJ[4]);
  node->field_if_defined(INTERN_(kw,ph4),this->_hasPeriodicity[4],this->_PhaseRads[4]);
  node->field_if_defined(INTERN_(kw,v5) ,this->_hasPeriodicity[5],this->_Vs_kJ[5]);
  node->field_if_defined(INTERN_(kw,ph5),this->_hasPeriodicity[5],this->_PhaseRads[5]);
  this->Base::fields(node);
}

void FFPtorDb_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
}


void	FFPtorDb_O::initialize()
{
    this->Base::initialize();
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
#if 0
 printf("%s:%d:%s  Adding ptor %s-%s-%s-%s\n", __FILE__, __LINE__, __FUNCTION__,
        _rep_(ptor->_T1).c_str(),
        _rep_(ptor->_T2).c_str(),
        _rep_(ptor->_T3).c_str(),
        _rep_(ptor->_T4).c_str() );
#endif
 core::T_sp tptor = this->findExactTerm(ptor->_T1,ptor->_T2,ptor->_T3,ptor->_T4);
 if (tptor.notnilp()) {
   ptorOld = gc::As<FFPtor_sp>(tptor);
   ptorOld->mergeWith(ptor);
 } else {
   key = keyString(ptor->_T1,ptor->_T2,ptor->_T3,ptor->_T4);
   this->_Parameters->setf_gethash(key,ptor);
 }
}


core::T_sp FFPtorDb_O::findExactTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
core::Symbol_sp key;
    core::T_sp parm;
    key = keyString(t1,t2,t3,t4);
    parm = this->_Parameters->gethash(key);
    if (parm.notnilp()) return parm;
    key = keyString(t4,t3,t2,t1);
    parm = this->_Parameters->gethash(key);
    if (parm.notnilp()) return parm;
    parm = _Nil<core::T_O>();
    return parm;
}

core::T_sp FFPtorDb_O::findBestTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
  core::Symbol_sp key,bl;
  core::T_sp ptor;
  key = keyString(t1,t2,t3,t4);
  ptor = this->_Parameters->gethash(key);
  if (ptor.notnilp()) return ptor;
  key = keyString(t4,t3,t2,t1);
  ptor = this->_Parameters->gethash(key);
  if (ptor.notnilp()) return ptor;
  bl = _Nil<core::Symbol_O>();
  key = keyString(bl,t2,t3,bl);
  ptor = this->_Parameters->gethash(key);
  if (ptor.notnilp()) return ptor;
  key = keyString(bl,t3,t2,bl);
  ptor = this->_Parameters->gethash(key);
  if (ptor.notnilp()) return ptor;
  ptor = _Nil<core::T_O>();
  return ptor;
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
