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
#include <clasp/core/array.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/units.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>



namespace chem
{






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

    core::Symbol_sp stretchKey(core::Symbol_sp t1
                               , core::Symbol_sp t2 )
    {
        return chemkw_intern(t1->symbolName()->get()+"-"
                           + t2->symbolName()->get());
    }

CL_NAME(CHEM:FFSTRETCH-DB-ADD);
CL_DEFMETHOD void    FFStretchDb_O::add(FFStretch_sp term)
{
    core::Symbol_sp key;
    key = stretchKey(term->_Type1,term->_Type2);        // forwards
    this->_Parameters->setf_gethash(key,term);
}


CL_LISPIFY_NAME("findTerm");
CL_DOCSTRING("Return the stretch term or NIL if none was found");
CL_DEFMETHOD core::T_sp	FFStretchDb_O::findTerm(chem::Atom_sp a1, chem::Atom_sp a2)
{
  FFStretch_sp	match;
  core::Symbol_sp key;
  core::Symbol_sp t1, t2;
  t1 = a1->getType();
  t2 = a2->getType();
  LOG(BF("Looking for stretch between types (%s)-(%s)") % t1.c_str() % t2.c_str() );
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
  node->/*pod_*/field( INTERN_(kw,kb), this->_Kb_kJPerNanometerSquared );
  this->Base::fields(node);
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
