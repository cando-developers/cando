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
#include <clasp/core/evaluator.h>
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
  node->field_if_defined(INTERN_(kw,v0) ,this->_HasPeriodicity[0],this->_Vs_kJ[0]);
  node->field_if_defined(INTERN_(kw,ph0),this->_HasPeriodicity[0],this->_PhaseRads[0]);
  node->field_if_defined(INTERN_(kw,v1) ,this->_HasPeriodicity[1],this->_Vs_kJ[1]);
  node->field_if_defined(INTERN_(kw,ph1),this->_HasPeriodicity[1],this->_PhaseRads[1]);
  node->field_if_defined(INTERN_(kw,v2) ,this->_HasPeriodicity[2],this->_Vs_kJ[2]);
  node->field_if_defined(INTERN_(kw,ph2),this->_HasPeriodicity[2],this->_PhaseRads[2]);
  node->field_if_defined(INTERN_(kw,v3) ,this->_HasPeriodicity[3],this->_Vs_kJ[3]);
  node->field_if_defined(INTERN_(kw,ph3),this->_HasPeriodicity[3],this->_PhaseRads[3]);
  node->field_if_defined(INTERN_(kw,v4) ,this->_HasPeriodicity[4],this->_Vs_kJ[4]);
  node->field_if_defined(INTERN_(kw,ph4),this->_HasPeriodicity[4],this->_PhaseRads[4]);
  node->field_if_defined(INTERN_(kw,v5) ,this->_HasPeriodicity[5],this->_Vs_kJ[5]);
  node->field_if_defined(INTERN_(kw,ph5),this->_HasPeriodicity[5],this->_PhaseRads[5]);
  node->field_if_defined(INTERN_(kw,scee),this->_HasScee,this->_Scee);
  node->field_if_defined(INTERN_(kw,scnb),this->_HasScnb,this->_Scnb);
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
{
    stringstream ss;
    ss << "Can't find ptor term for ("<<t1<<")-("<<t2<<")-("<<t3<<")-("<<t4<<")";
    SIMPLE_ERROR(("%s") , ss.str());
}


void FFPtor_O::initialize()
{
    this->Base::initialize();
    int i;
    this->_T1 = nil<core::Symbol_O>();
    this->_T2 = nil<core::Symbol_O>();
    this->_T3 = nil<core::Symbol_O>();
    this->_T4 = nil<core::Symbol_O>();
    for ( i=0; i<MaxPeriodicity; i++ )
    {
	this->_HasPeriodicity[i] = false;
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
{
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for hasPeriodicity") , idx);
    }
    return this->_HasPeriodicity[idx-1];
}




double  FFPtor_O::getV_kJ(int idx) const
{
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for getV") , idx);
    }
    return this->_Vs_kJ[idx-1];
}

CL_DEFMETHOD void    FFPtor_O::setV_kJ(int idx, double val)
{
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for setV") , idx);
    }
    this->_HasPeriodicity[idx-1] = true;
    this->_Vs_kJ[idx-1] = val;
}

CL_DEFMETHOD double FFPtor_O::getV_kCal(int idx) const
{
    return kJ_to_kCal(this->getV_kJ(idx));
}

CL_DEFMETHOD void FFPtor_O::setV_kCal(int idx, double val)
{
    this->setV_kJ(idx,kCal_to_kJ(val));
}


CL_DEFMETHOD double  FFPtor_O::getPhaseRad(int idx) const
{
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity[%d] for getPhaseRad - must be in [1,%d]") , idx , MaxPeriodicity );
    }
    return this->_PhaseRads[idx-1];
}

CL_DEFMETHOD void    FFPtor_O::setPhaseRad(int idx, double val)
{
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for setPhaseRad") , idx);
    }
    this->_PhaseRads[idx-1] = val;
}


void    FFPtor_O::mergeWith(FFPtor_sp ptor)
{
int     i;
    for ( i=0; i<MaxPeriodicity; i++ ) {
	if ( ptor->_HasPeriodicity[i] ) {
            this->_HasPeriodicity[i] = ptor->_HasPeriodicity[i];
            this->_Vs_kJ[i] = ptor->_Vs_kJ[i];
            this->_PhaseRads[i] = ptor->_PhaseRads[i];
        }
    }
}

#define keyString(t1,t2,t3,t4) chemkw_intern(t1->symbolName()->get_std_string()+"-"+t2->symbolName()->get_std_string()+"-"+t3->symbolName()->get_std_string()+"-"+t4->symbolName()->get_std_string())

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

CL_DEFMETHOD void FFPtorDb_O::addFFPtor(FFPtor_sp obj) {
  this->add(obj);
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
    parm = nil<core::T_O>();
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
  bl = nil<core::Symbol_O>();
  key = keyString(bl,t2,t3,bl);
  ptor = this->_Parameters->gethash(key);
  if (ptor.notnilp()) return ptor;
  key = keyString(bl,t3,t2,bl);
  ptor = this->_Parameters->gethash(key);
  if (ptor.notnilp()) return ptor;
  ptor = nil<core::T_O>();
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


void FFPtorDb_O::forceFieldMerge(FFBaseDb_sp bother)
{
  FFPtorDb_sp other = gc::As<FFPtorDb_sp>(bother);
  other->_Parameters->maphash([this] (core::T_sp key, core::T_sp value) {
      core::Symbol_sp skey = gc::As<core::Symbol_sp>(key);
      this->_Parameters->hash_table_setf_gethash(skey,value);
    } );
}

};
