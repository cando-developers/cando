/*
    File: ffItorDb.cc
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
       
//
// (C) 2004 Christian E. Schafmeister
//

#define	TURN_DEBUG_OFF

/*
 *	ffItor.cc
 *
 *	Maintains a database of itor terms
 */
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/ffItorDb.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>
#include <cando/chem/units.h>


namespace chem
{

    const int FFItor_O::IMaxPeriodicity;



void	FFItor_O::initialize()
{
    this->Base::initialize();
    int i;
    this->_T1 = nil<core::Symbol_O>();
    this->_T2 = nil<core::Symbol_O>();
    this->_T3 = nil<core::Symbol_O>();
    this->_T4 = nil<core::Symbol_O>();
    for ( i=0; i<FFItor_O::IMaxPeriodicity; i++ )
    {
	this->_hasPeriodicity[i] = false;
	this->_Vs_kJ[i] = 0.0;
	this->_PhaseRads[i] = 0.0;
    }
}

void	FFItor_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,type1), this->_T1 );
  node->field( INTERN_(kw,type2), this->_T2 );
  node->field( INTERN_(kw,type3), this->_T3 );
  node->field( INTERN_(kw,type4), this->_T4 );
  ASSERT(IMaxPeriodicity==6);
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

#define keyString(t1,t2,t3,t4) chemkw_intern(t1->symbolName()->get_std_string()+"-"+t2->symbolName()->get_std_string()+"-"+t3->symbolName()->get_std_string()+"-"+t4->symbolName()->get_std_string())





void	FFItorDb_O::initialize()
{
    this->Base::initialize();
}

void	FFItorDb_O::add( FFItor_sp itor )
{
  core::Symbol_sp          key;
  FFItor_sp        itorOld;
  core::T_sp titorOld = this->findExactTerm(itor->_T1,itor->_T2,itor->_T3,itor->_T4);
  if (titorOld.notnilp()) {
    itorOld = gc::As<FFItor_sp>(titorOld);
    itorOld->mergeWith(itor);
  } else {
    key = keyString(itor->_T1,itor->_T2,itor->_T3,itor->_T4);
    this->_Parameters->setf_gethash(key,itor);
#ifdef DEBUG_ON
    if ( itor->_T3=="c" && itor->_T4=="o" ) {
      LOG(BF("FFItorDb::add adding term with key: %s") % key.c_str()  );
    }
#endif
  }
}

core::T_sp FFItorDb_O::findExactTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
core::Symbol_sp         key;
 core::T_sp itor;
 key = keyString(t1,t2,t3,t4);
 itor = this->_Parameters->gethash(key);
 return itor;
}

// Order the two atoms according to their atom type and atom name
// If the atom types are the same then order on atom name (alphabetical)
// If the atom types are different then order on atom type (alphabetical)
// Atom names and atom types are stored as interned names (Symbols)
bool ordered(Atom_sp a1, Atom_sp a2) {
  if (a1->getType() == a2->getType()) {
    return a1->getName()->symbolNameAsString() < a2->getName()->symbolNameAsString();
  }
  return gc::As<core::Symbol_sp>(a1->getType())->symbolNameAsString() < gc::As<core::Symbol_sp>(a2->getType())->symbolNameAsString();
}

// Sort the three atoms of an improper a1-a2-a3-a4 where a3 is the central atom.
void FFItorDb_O::improperAtomSort(Atom_sp& a1, Atom_sp& a2, Atom_sp& a4)
{
  Atom_sp ta;
  if (!ordered(a1,a2)) {ta = a1; a1 = a2; a2 = ta;}
  if (!ordered(a2,a4)) {ta = a2; a2 = a4; a4 = ta;}
  if (!ordered(a1,a2)) {ta = a1; a1 = a2; a2 = ta;}
}

CL_DEFUN core::T_mv chem__improperAtomSort(Atom_sp a1, Atom_sp a2, Atom_sp a3)
{
  FFItorDb_O::improperAtomSort(a1,a2,a3);
  return Values(a1,a2,a3);
}

core::T_sp FFItorDb_O::findBestTerm( core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{
  core::Symbol_sp        key,key1,key2,key3,key4,bl;
  core::T_sp itor;
  key4 = keyString(t1,t2,t3,t4);
  itor = this->_Parameters->gethash(key4);
  if (itor.notnilp()) {
    key = key4;
    goto DONE;
  }
  bl = nil<core::Symbol_O>();
  key3 = keyString(bl,t2,t3,t4);
  itor = this->_Parameters->gethash(key3);
  if (itor.notnilp()) {
    key = key3;
    goto DONE;
  }
  key2 = keyString(bl,bl,t3,t4);
  itor = this->_Parameters->gethash(key2);
  if (itor.notnilp()) {
    key = key2;
    goto DONE;
  }
  key1 = keyString(bl,bl,t3,bl);
  itor = this->_Parameters->gethash(key1);
  if (itor.notnilp()) {
    key = key1;
    goto DONE;
  }
  itor = nil<core::T_O>();
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



void    FFItorDb_O::cantFind(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 )
{_OF();
    stringstream ss;
    ss << "Can't find itor term for ("<<t1<<")-("<<t2<<")-("<<t3<<")-("<<t4<<")";
    SIMPLE_ERROR(BF("%s")%ss.str());
}



void FFItorDb_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
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


void FFItorDb_O::forceFieldMerge(FFBaseDb_sp bother)
{
  //SIMPLE_WARN(BF("Merging FFItorDb terms - but terms with different type orders will create duplicates!"));
  FFItorDb_sp other = gc::As<FFItorDb_sp>(bother);
  other->_Parameters->maphash([this] (core::T_sp key, core::T_sp value) {
      core::Symbol_sp skey = gc::As<core::Symbol_sp>(key);
      this->_Parameters->hash_table_setf_gethash(skey,value);
    } );
}




};
