/*
    File: ffAngleDb.cc
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
       
//#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffAngle.cc
 *
 *	Maintains a database of stretch types
 */
#include <clasp/core/foundation.h>
#include <clasp/core/array.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/units.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>


namespace chem
{



void FFAngle_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,type1),this->_Type1);
  node->field(INTERN_(kw,type2),this->_Type2);
  node->field(INTERN_(kw,type3),this->_Type3);
  node->field(INTERN_(kw,angRad), this->_AngRad);
  node->field(INTERN_(kw,k2),this->_K2__kJPerRadianSquared);
  this->Base::fields(node);
}

void	FFAngleDb_O::fields(core::Record_sp node)
{
    node->field(INTERN_(kw,angleFunction),this->_AngleFunction);
    node->field(INTERN_(kw,zConstants),this->_ZConstants);
    node->field(INTERN_(kw,cConstants),this->_CConstants);
    this->Base::fields(node);
}


void	FFAngleDb_O::initialize()
{
    this->Base::initialize();
    this->_ZConstants.clear();
    this->_CConstants.clear();
}




void FFAngle_O::setK2_kCalPerRadianSquared(double val)
{
    this->_K2__kJPerRadianSquared = kCalPerRadianSquared_to_kJPerRadianSquared(val);
}


void FFAngle_O::setK2_kJPerRadianSquared(double val)
{
    this->_K2__kJPerRadianSquared = val;
}


double FFAngle_O::getK2_kCalPerRadianSquared() const
{
    return kJPerRadianSquared_to_kCalPerRadianSquared(this->_K2__kJPerRadianSquared);
}


double FFAngle_O::getK2_kJPerRadianSquared() const
{
    return this->_K2__kJPerRadianSquared;
}






    core::Symbol_sp angleKey(core::Symbol_sp t1
                             , core::Symbol_sp t2
                             , core::Symbol_sp t3)
    {
        return chemkw_intern(t1->symbolName()->get_std_string()+"-"
                           + t2->symbolName()->get_std_string()+"-"
                           + t3->symbolName()->get_std_string());
    }



void    FFAngleDb_O::add(FFAngle_sp ang)
{
    core::Symbol_sp  key;
    key = angleKey(ang->_Type1,ang->_Type2,ang->_Type3);
    this->_Parameters->setf_gethash(key,ang);
}


FFAngle_sp FFAngleDb_O::findTerm(FFStretchDb_sp ffstretch, chem::Atom_sp a1, chem::Atom_sp a2, chem::Atom_sp a3 )
{ 
FFAngle_sp       match;
core::Symbol_sp          key;
core::Symbol_sp		t1, t2, t3;
    t1 = a1->getType();
    t2 = a2->getType();
    t3 = a3->getType();
    key = angleKey(t1,t2,t3);
    core::T_sp parm;
    parm = this->_Parameters->gethash(key);
    if (parm.notnilp()) return gc::As_unsafe<FFAngle_sp>(parm);
    key = angleKey(t3,t2,t1);
    parm = this->_Parameters->gethash(key);
    if (parm.notnilp()) return gc::As_unsafe<FFAngle_sp>(parm);
    match = this->estimateTerm(ffstretch,a1,a2,a3);
    return match;
}



void	FFAngle_O::initialize()
{
    this->Base::initialize();
    this->_Type1 = nil<core::Symbol_O>();
    this->_Type2 = nil<core::Symbol_O>();
    this->_Type3 = nil<core::Symbol_O>();
    this->_AngRad = 0.0;
    this->_K2__kJPerRadianSquared = 0.0;
//    this->_K3 = 0.0;
//    this->_K4 = 0.0;
//    this->_Ub_k = 0.0;
//    this->_Ub_len = 0.0;
}

SYMBOL_EXPORT_SC_(ChemPkg,warn_estimated_angle_term);

/*! Estimate the angle term according to Wang et al. J. Comput. Chem 25, 1157-1174 (2004) */
FFAngle_sp	FFAngleDb_O::estimateTerm(FFStretchDb_sp ffstretch, chem::Atom_sp a1, chem::Atom_sp a2, chem::Atom_sp a3 )
{
FFStretch_sp	r12,r32;
FFAngle_sp	ff_121;
FFAngle_sp	ff_323;
double		d, z1, c2, z3, angRad, k;
core::Symbol_sp		element1, element2, element3;
core::Symbol_sp		t1, t2, t3;
    LOG(BF("status") );
    t1 = a1->getType();
    t2 = a2->getType();
    t3 = a3->getType();
    LOG(BF("status") );
    // To avoid an std::endless loop of estimating terms just go straight to guess if t1==t3
    if ( t1==t3 ) goto GUESS;
    LOG(BF("status") );
    ff_121 = this->findTerm(ffstretch,a1,a2,a1);
    ff_323 = this->findTerm(ffstretch,a3,a2,a3);
    LOG(BF("status") );
    if ( ff_121.nilp() ) goto GUESS;
    LOG(BF("status") );
    if ( ff_323.nilp() ) goto GUESS;
    LOG(BF("status") );
    angRad = (ff_121->_AngRad + ff_323->_AngRad)/2.0;
    element1 = a1->getElementAsSymbol();
    element2 = a2->getElementAsSymbol();
    element3 = a3->getElementAsSymbol();
    LOG(BF("status") );
    if ( this->_ZConstants.count(element1) == 0 ) goto GUESS;
    LOG(BF("status") );
    if ( this->_CConstants.count(element2) == 0 ) goto GUESS;
    LOG(BF("status") );
    if ( this->_ZConstants.count(element3) == 0 ) goto GUESS;
    LOG(BF("status") );
    z1 = this->_ZConstants.get(element1)->get();
    c2 = this->_CConstants.get(element2)->get();
    z3 = this->_ZConstants.get(element3)->get();
    r12 = gc::As_unsafe<FFStretch_sp>(ffstretch->findTerm(a1,a2));
    r32 = gc::As_unsafe<FFStretch_sp>(ffstretch->findTerm(a3,a2));
    LOG(BF("status") );
    if ( r12.nilp() ) goto GUESS;
    LOG(BF("status") );
    if ( r32.nilp() ) goto GUESS;
    LOG(BF("status") );
    d = (r12->getR0_Angstrom()-r32->getR0_Angstrom())/(r12->getR0_Angstrom()+r32->getR0_Angstrom());
    d = d*d;
    k = 143.9*z1*c2*z3*(1.0/((r12->getR0_Angstrom()+r32->getR0_Angstrom())*(angRad*angRad)))*exp(-2.0*d);
    {
      FFAngle_sp angle = FFAngle_O::create();
      angle->_Level = estimated;
      angle->_Type1 = t1;
      angle->_Type2 = t2;
      angle->_Type3 = t3;
      angle->_AngRad = angRad;
      angle->setK2_kCalPerRadianSquared(k);
      core::eval::funcall(_sym_warn_estimated_angle_term,t1,t2,t3,
                          core::clasp_make_double_float(angle->getAngle_Radian()),
                          core::clasp_make_double_float(angle->getK2_kCalPerRadianSquared()));
      return angle;
    }
GUESS:
    {
      FFAngle_sp angle = FFAngle_O::create();
      angle->_Level = rough;
      angle->_Type1 = t1;
      angle->_Type2 = t2;
      angle->_Type3 = t3;
      if ( a2->getHybridization() == hybridization_sp3 ) {
	angle->_AngRad = 109.5*0.0174533;
      } else if ( a2->getHybridization() == hybridization_sp2 ) {
	angle->_AngRad = 120.0*0.0174533;
      } else if ( a2->getHybridization() == hybridization_sp ) {
	angle->_AngRad = 180.0*0.0174533;
      } else {
	angle->_AngRad = 109.5*0.0174533;
      }
		// See Gaff paper for where I get these force constants
      if ( a1->getElement() == element_H && a3->getElement() == element_H ) {
	angle->_K2__kJPerRadianSquared = kCalPerRadianSquared_to_kJPerRadianSquared(32.5);
      } else if ( a1->getElement() == element_H || a3->getElement() == element_H ) {
	angle->_K2__kJPerRadianSquared = kCalPerRadianSquared_to_kJPerRadianSquared(50.0);
      } else {
	angle->_K2__kJPerRadianSquared = kCalPerRadianSquared_to_kJPerRadianSquared(70.0);
      }
      core::eval::funcall(_sym_warn_estimated_angle_term,t1,t2,t3,core::clasp_make_double_float(angle->getAngle_Radian()),core::clasp_make_double_float(angle->getK2_kCalPerRadianSquared()));
      return angle;
    }
}


string	FFAngle_O::levelDescription()
{
stringstream	desc;
    desc << this->FFParameter_O::levelDescription();
    desc << " angle between types "<<this->_Type1<< " "
		<< this->_Type2 << " "
		<< this->_Type3;
    return desc.str();
}

void FFAngleDb_O::forceFieldMerge(FFBaseDb_sp bother)
{
  //SIMPLE_WARN(BF("Merging angle terms - but terms with different type orders will create duplicates!"));
  FFAngleDb_sp other = gc::As<FFAngleDb_sp>(bother);
  other->_Parameters->maphash([this] (core::T_sp key, core::T_sp value) {
      core::Symbol_sp skey = gc::As<core::Symbol_sp>(key);
      this->_Parameters->hash_table_setf_gethash(skey,value);
    } );
}



};
