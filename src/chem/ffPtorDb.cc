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

CL_LISPIFY_NAME(FFPtorDb/make);
CL_DEF_CLASS_METHOD FFPtorDb_sp FFPtorDb_O::make()
{
  return FFPtorDb_O::create();
}


void	FFPtor_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,type1), this->_T1 );
  node->field( INTERN_(kw,type2), this->_T2 );
  node->field( INTERN_(kw,type3), this->_T3 );
  node->field( INTERN_(kw,type4), this->_T4 );
  ASSERT(MaxPeriodicity==6);
  node->field_if_defined(INTERN_(kw,id0) ,this->_HasPeriodicity[0],this->_idivf[0]);
  node->field_if_defined(INTERN_(kw,v0) ,this->_HasPeriodicity[0],this->_Vs_kj[0]);
  node->field_if_defined(INTERN_(kw,ph0),this->_HasPeriodicity[0],this->_PhaseDegrees[0]);
  node->field_if_defined(INTERN_(kw,id1) ,this->_HasPeriodicity[1],this->_idivf[1]);
  node->field_if_defined(INTERN_(kw,v1) ,this->_HasPeriodicity[1],this->_Vs_kj[1]);
  node->field_if_defined(INTERN_(kw,ph1),this->_HasPeriodicity[1],this->_PhaseDegrees[1]);
  node->field_if_defined(INTERN_(kw,id2) ,this->_HasPeriodicity[2],this->_idivf[2]);
  node->field_if_defined(INTERN_(kw,v2) ,this->_HasPeriodicity[2],this->_Vs_kj[2]);
  node->field_if_defined(INTERN_(kw,ph2),this->_HasPeriodicity[2],this->_PhaseDegrees[2]);
  node->field_if_defined(INTERN_(kw,id3) ,this->_HasPeriodicity[3],this->_idivf[3]);
  node->field_if_defined(INTERN_(kw,v3) ,this->_HasPeriodicity[3],this->_Vs_kj[3]);
  node->field_if_defined(INTERN_(kw,ph3),this->_HasPeriodicity[3],this->_PhaseDegrees[3]);
  node->field_if_defined(INTERN_(kw,id4) ,this->_HasPeriodicity[4],this->_idivf[4]);
  node->field_if_defined(INTERN_(kw,v4) ,this->_HasPeriodicity[4],this->_Vs_kj[4]);
  node->field_if_defined(INTERN_(kw,ph4),this->_HasPeriodicity[4],this->_PhaseDegrees[4]);
  node->field_if_defined(INTERN_(kw,id5) ,this->_HasPeriodicity[5],this->_idivf[5]);
  node->field_if_defined(INTERN_(kw,v5) ,this->_HasPeriodicity[5],this->_Vs_kj[5]);
  node->field_if_defined(INTERN_(kw,ph5),this->_HasPeriodicity[5],this->_PhaseDegrees[5]);
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
{_OF();
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
	this->_Vs_kj[i] = 0.0;
	this->_PhaseDegrees[i] = 0.0;
    }
}



CL_LISPIFY_NAME(FFPtor/getTypes);
CL_DEFMETHOD
core::T_mv FFPtor_O::getTypes() const
{
  return Values(this->_T1,this->_T2,this->_T3,this->_T4);
}

CL_LISPIFY_NAME(FFPtor/setTypes);
CL_DEFMETHOD
void    FFPtor_O::setTypes(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
    this->_T1 = t1;
    this->_T2 = t2;
    this->_T3 = t3;
    this->_T4 = t4;
}

CL_LISPIFY_NAME(FFPtor/hasPeriodicity);
CL_DEFMETHOD
bool	FFPtor_O::hasPeriodicity(int idx ) const
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for hasPeriodicity") , idx);
    }
    return this->_HasPeriodicity[idx-1];
}

CL_LISPIFY_NAME(FFPtor/maxPeriodicity);
CL_DEFMETHOD
size_t FFPtor_O::maxPeriodicity() const
{
  for (size_t periodicity=MaxPeriodicity; periodicity>0; periodicity-- ) {
    size_t idx = periodicity-1;
    if (this->_HasPeriodicity[idx]) return periodicity;
  }
  return 0;
}


CL_LISPIFY_NAME(FFPtor/getIdivf);
CL_DEFMETHOD
int FFPtor_O::getIdivf(int idx) const
{
  if ( idx < 1 || idx > MaxPeriodicity ) {
    SIMPLE_ERROR(("Illegal periodicity (%d) for getIdivf") , idx);
  }
  return this->_idivf[idx-1];
}

CL_LISPIFY_NAME(FFPtor/setIdivf);
CL_DEFMETHOD
void FFPtor_O::setIdivf(int idx, int idivf)
{
  if ( idx < 1 || idx > MaxPeriodicity ) {
    SIMPLE_ERROR(("Illegal periodicity (%d) for getIdivf") , idx);
  }
  this->_idivf[idx-1] = idivf;
}



CL_LISPIFY_NAME(FFPtor/getV_kj);
CL_DEFMETHOD
double  FFPtor_O::getV_kj(int idx) const
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for getV") , idx);
    }
    return this->_Vs_kj[idx-1];
}

CL_LISPIFY_NAME(FFPtor/setV_kj);
CL_DEFMETHOD
void    FFPtor_O::setV_kj(int idx, double val)
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for setV") , idx);
    }
    this->_HasPeriodicity[idx-1] = true;
    this->_Vs_kj[idx-1] = val;
}

CL_LISPIFY_NAME(FFPtor/getV_kcal);
CL_DEFMETHOD
double FFPtor_O::getV_kcal(int idx) const
{_OF();
    return kj_to_kcal(this->getV_kj(idx));
}

CL_LISPIFY_NAME(FFPtor/setV_kcal);
CL_DEFMETHOD
void FFPtor_O::setV_kcal(int idx, double val)
{_OF();
    this->setV_kj(idx,kcal_to_kj(val));
}


CL_LISPIFY_NAME(FFPtor/getPhaseRadians);
CL_DEFMETHOD
double  FFPtor_O::getPhaseRadians(int idx) const
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity[%d] for getPhaseRadians - must be in [1,%d]") , idx , MaxPeriodicity );
    }
    return this->_PhaseDegrees[idx-1]*0.0174533;
}

CL_LISPIFY_NAME(FFPtor/setPhaseRadians);
CL_DEFMETHOD
void    FFPtor_O::setPhaseRadians(int idx, double val)
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for setPhaseRadians") , idx);
    }
    this->_PhaseDegrees[idx-1] = val/0.0174533;
}


CL_LISPIFY_NAME(FFPtor/getPhaseDegrees);
CL_DEFMETHOD
double  FFPtor_O::getPhaseDegrees(int idx) const
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity[%d] for getPhaseDegrees - must be in [1,%d]") , idx , MaxPeriodicity );
    }
    return this->_PhaseDegrees[idx-1];
}

CL_LISPIFY_NAME(FFPtor/setPhaseDegrees);
CL_DEFMETHOD
void    FFPtor_O::setPhaseDegrees(int idx, double val)
{_OF();
    if ( idx < 1 || idx > MaxPeriodicity ) {
        SIMPLE_ERROR(("Illegal periodicity (%d) for setPhaseDegrees") , idx);
    }
    this->_PhaseDegrees[idx-1] = val;
}


void    FFPtor_O::mergeWith(FFPtor_sp ptor)
{
int     i;
    for ( i=0; i<MaxPeriodicity; i++ ) {
	if ( ptor->_HasPeriodicity[i] ) {
            this->_HasPeriodicity[i] = ptor->_HasPeriodicity[i];
            this->_idivf[i] = ptor->_idivf[i];
            this->_Vs_kj[i] = ptor->_Vs_kj[i];
            this->_PhaseDegrees[i] = ptor->_PhaseDegrees[i];
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
