/*
    File: ffStretchDb.cc
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
#include <clasp/core/evaluator.h>
#include <clasp/core/array.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/units.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>

SYMBOL_EXPORT_SC_(ChemPkg,warn_estimated_stretch_term);


namespace chem
{

CL_LISPIFY_NAME(FFStretchDb/make);
CL_DEF_CLASS_METHOD FFStretchDb_sp FFStretchDb_O::make()
{
  return FFStretchDb_O::create();
}





core::Symbol_sp stretchKey(core::Symbol_sp t1
                           , core::Symbol_sp t2 )
{
  return chemkw_intern(t1->symbolName()->get_std_string()+"-"
                       + t2->symbolName()->get_std_string());
}

CL_NAME(CHEM:FFSTRETCH-DB/ADD);
CL_DEFMETHOD void    FFStretchDb_O::add(FFStretch_sp term)
{
  core::Symbol_sp key;
  key = stretchKey(term->_Type1,term->_Type2);        // forwards
  this->_Parameters->setf_gethash(key,term);
}

//  (chem:canonical-stretch-key (chem:get-type atom1) (chem:get-type atom2))
DOCGROUP(cando);
CL_DEFUN core::Symbol_sp chem__canonicalStretchKey(core::Symbol_sp type1, core::Symbol_sp type2) {
  std::string type1str = type1->symbolNameAsString();
  std::string type2str = type2->symbolNameAsString();
  if (type1str<=type2str) {
    return chemkw_intern(type1str+"-"+type2str);
  }
  return chemkw_intern(type2str+"-"+type1str);
}

CL_LISPIFY_NAME("findTerm");
CL_DOCSTRING(R"dx(Return the stretch term or NIL if none was found)dx")
CL_DEFMETHOD core::T_sp	FFStretchDb_O::findTermForTypes(core::Symbol_sp t1, core::Symbol_sp t2)
{
  FFStretch_sp	match;
  core::Symbol_sp key;
  LOG("Looking for stretch between types ({})-({})" , t1.c_str() , t2.c_str() );
  key = stretchKey(t1,t2); // forwards
  core::T_sp parm = this->_Parameters->gethash(key);
  if (parm.notnilp()) return parm;
  key = stretchKey(t2,t1); // backwards
  parm = this->_Parameters->gethash(key);
  if (parm.notnilp()) return parm;
  FFStretch_sp missing = FFStretch_O::create_missing(t1,t2);
  return missing;
}



void	FFStretchDb_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
}




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
  node->field( INTERN_(kw,type1), this->_Type1 );
  node->field( INTERN_(kw,type2), this->_Type2 );
  node->/*pod_*/field( INTERN_(kw,r0), this->_R0_Nanometer );
  node->/*pod_*/field( INTERN_(kw,kb), this->_Kb_kjPerNanometerSquared );
  this->Base::fields(node);
}

void FFStretch_O::initialize()
{
  this->Base::initialize();
  this->_Type1 = nil<core::T_O>();
  this->_Type2 = nil<core::T_O>();
  this->_R0_Nanometer = 0.0;
  this->_Kb_kjPerNanometerSquared = 0.0;
  //    this->_K3 = 0.0;
  //    this->_K4 = 0.0;
  //    this->_Bci = 0.0;
}


CL_LISPIFY_NAME(FFStretch/getTypes);
CL_DEFMETHOD
core::T_mv FFStretch_O::getTypes() const
{
  return Values(this->_Type1,this->_Type2);
}

CL_LISPIFY_NAME(FFStretch/getR0_Nanometer);
CL_DEFMETHOD
double FFStretch_O::getR0_Nanometer() const
{
  return this->_R0_Nanometer;
}

CL_LISPIFY_NAME(FFStretch/getR0_Angstrom);
CL_DEFMETHOD
double FFStretch_O::getR0_Angstrom() const
{
  return this->_R0_Nanometer*10.0;
}

CL_LISPIFY_NAME(FFStretch/getKb_kjPerNanometerSquared);
CL_DEFMETHOD
double FFStretch_O::getKb_kjPerNanometerSquared() const
{
  return this->_Kb_kjPerNanometerSquared;
}


CL_LISPIFY_NAME(FFStretch/setKb_kjPerNanometerSquared);
CL_DEFMETHOD
void FFStretch_O::setKb_kjPerNanometerSquared(double val)
{
  this->_Kb_kjPerNanometerSquared = val;
}


CL_LISPIFY_NAME(FFStretch/setR0_Nanometer);
CL_DEFMETHOD
void FFStretch_O::setR0_Nanometer(double val)
{
  this->_R0_Nanometer = val;
}


CL_LISPIFY_NAME(FFStretch/getKb_kcalPerAngstromSquared);
CL_DEFMETHOD
double FFStretch_O::getKb_kcalPerAngstromSquared() const
{
  return kjPerNanometerSquared_to_kcalPerAngstromSquared(this->_Kb_kjPerNanometerSquared);
}





string	FFStretch_O::levelDescription()
{
  stringstream	desc;
  desc << this->FFParameter_O::levelDescription();
  desc << " stretch between types "<<this->_Type1<< " "
       << this->_Type2;
  return desc.str();
}


void FFStretchDb_O::forceFieldMerge(FFBaseDb_sp bother)
{
  FFStretchDb_sp other = gc::As<FFStretchDb_sp>(bother);
  other->_Parameters->maphash([this] (core::T_sp key, core::T_sp value) {
    core::Symbol_sp skey = gc::As<core::Symbol_sp>(key);
    this->_Parameters->hash_table_setf_gethash(skey,value);
  } );
}

};
