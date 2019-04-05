/*
    File: ffNonbondDb.cc
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

#define TURN_DEBUG_OFF

/*
 *	ffNonbonds.cc
 *
 *	Maintains databases and structures to store types and type assignement
 *	rules.
 */
#include <clasp/core/foundation.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/array.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/loop.h>
#include <cando/chem/units.h>
#include <clasp/core/wrappers.h>

namespace chem {
void FFNonbondCrossTermTable_O::fillUsingFFNonbondDb( FFNonbondDb_sp db )
{
  uint numberOfTypes = db->numberOfTypes();
  this->_NumberOfTypes = numberOfTypes;
  for ( uint it1 = 0; it1 < numberOfTypes; it1++ )
  {
    for ( uint it2 = 0; it2 < numberOfTypes; it2++ )
    {
      FFNonbondCrossTerm term;
      term._Type1Index = it1;
      term._Type2Index = it2;
      FFNonbond_sp ffNonbond1 = db->getFFNonbondUsingTypeIndex(it1);
      FFNonbond_sp ffNonbond2 = db->getFFNonbondUsingTypeIndex(it2);
      term._RStar = ffNonbond1->getRadius_Angstroms() + ffNonbond2->getRadius_Angstroms();
      double epsilonij = sqrt(ffNonbond1->getEpsilon_kCal()*ffNonbond2->getEpsilon_kCal());
      term._A = epsilonij*pow(term._RStar,12.0);
      term._C = 2.0*epsilonij*pow(term._RStar,6.0);
      this->_CrossTerms.push_back(term);
    }
  }
}

uint FFNonbondCrossTermTable_O::typeMajorIndex(uint typeIndex)
{
  ASSERT_lt(typeIndex,this->_NumberOfTypes);
  return typeIndex*this->_NumberOfTypes;
}

FFNonbondCrossTerm& FFNonbondCrossTermTable_O::nonbondCrossTermFromTypes(uint type1, uint type2 )
{
  ASSERT_lt(type1,this->_NumberOfTypes);
  ASSERT_lt(type2,this->_NumberOfTypes);
  uint idx = type1*this->_NumberOfTypes + type2;
  ASSERT_lt(idx,this->_CrossTerms.size());
  return this->_CrossTerms[idx];
}

FFNonbondCrossTerm& FFNonbondCrossTermTable_O::nonbondCrossTerm(uint idx)
{_OF();
  ASSERT_lt(idx,this->_CrossTerms.size());
  return this->_CrossTerms[idx];
}

uint FFNonbondCrossTermTable_O::nonbondCrossTermIndexFromTypes(uint type1, uint type2 )
{_OF();
  ASSERT_lt(type1,this->_NumberOfTypes);
  ASSERT_lt(type2,this->_NumberOfTypes);
  uint idx = type1*this->_NumberOfTypes + type2;
  return idx;
}

core::NullTerminatedEnumAssociation daEnum[] = {
    { "daDonor", daDonor },
    { "daAcceptor", daAcceptor },
    { "daNeither", daNeither },
    { "", -1 }
};

core::NullTerminatedEnumAssociation dielectricEnum[] = {
    { "edConstant", edConstant},
    { "edDistance", edDistance},
    { "", -1 }
};

core::NullTerminatedEnumAssociation vmrEnum[] = {
    { "vmrAverage", vmrAverage},
    { "vmrGeometric", vmrGeometric},
    { "vmrMmff94", vmrMmff94},
    { "", -1 }
};

core::NullTerminatedEnumAssociation vmwEnum[] = {
    { "vmwAverage", vmwAverage},
    { "vmwMmff94", vmwMmff94},
    { "", -1 }
};

SYMBOL_EXPORT_SC_(ChemPkg,nonbond_force_field_component_merge);

SYMBOL_EXPORT_SC_(ChemPkg,STARDonorAcceptorEnumConverterSTAR);
SYMBOL_EXPORT_SC_(ChemPkg,STAREleDielectricEnumConverterSTAR);
SYMBOL_EXPORT_SC_(ChemPkg,STARVdwMixRadiusEnumConverterSTAR);
SYMBOL_EXPORT_SC_(ChemPkg,STARVdwMixWellEnumConverterSTAR);



SYMBOL_EXPORT_SC_(KeywordPkg,daDonor);
SYMBOL_EXPORT_SC_(KeywordPkg,daAcceptor);
SYMBOL_EXPORT_SC_(KeywordPkg,daNeither);
SYMBOL_EXPORT_SC_(KeywordPkg,edConstant);
SYMBOL_EXPORT_SC_(KeywordPkg,edDistance);
SYMBOL_EXPORT_SC_(KeywordPkg,vmrAverage);
SYMBOL_EXPORT_SC_(KeywordPkg,vmrGeometric);
SYMBOL_EXPORT_SC_(KeywordPkg,vmrMmff94);
SYMBOL_EXPORT_SC_(KeywordPkg,vmwAverage);
SYMBOL_EXPORT_SC_(KeywordPkg,vmwMmff94);

CL_BEGIN_ENUM(DonorAcceptorEnum,_sym_STARDonorAcceptorEnumConverterSTAR,"DonorAcceptorEnum");
CL_VALUE_ENUM( kw::_sym_daDonor, daDonor );
CL_VALUE_ENUM( kw::_sym_daAcceptor, daAcceptor );
CL_VALUE_ENUM( kw::_sym_daNeither, daNeither );
CL_END_ENUM(_sym_STARDonorAcceptorEnumConverterSTAR);

CL_BEGIN_ENUM(EleDielectricEnum,_sym_STAREleDielectricEnumConverterSTAR,"EleDielectricEnum");
CL_VALUE_ENUM( kw::_sym_edConstant, edConstant);
CL_VALUE_ENUM( kw::_sym_edDistance, edDistance);
CL_END_ENUM(_sym_STAREleDielectricEnumConverterSTAR);

CL_BEGIN_ENUM(VdwMixRadiusEnum,_sym_STARVdwMixRadiusEnumConverterSTAR,"VdwMixRadiusEnum");
CL_VALUE_ENUM( kw::_sym_vmrAverage, vmrAverage);
CL_VALUE_ENUM( kw::_sym_vmrGeometric, vmrGeometric);
CL_VALUE_ENUM( kw::_sym_vmrMmff94, vmrMmff94);
CL_END_ENUM(_sym_STARVdwMixRadiusEnumConverterSTAR);

CL_BEGIN_ENUM(VdwMixWellEnum,_sym_STARVdwMixWellEnumConverterSTAR,"VdwMixWellEnum");
CL_VALUE_ENUM( kw::_sym_vmwAverage, vmwAverage);
CL_VALUE_ENUM( kw::_sym_vmwMmff94, vmwMmff94);
CL_END_ENUM(_sym_STARVdwMixWellEnumConverterSTAR);


CL_LAMBDA(type &key (radius_nanometers 0.0) (epsilon_kj 0.0) (apol 0.0) (neff 0.0) (mass 0.0) (polarizability 0.0) (initial_charge 0.0) (fcadj 0.0) (pbci 0.0) (donor_acceptor :da-neither));
CL_DEF_CLASS_METHOD
FFNonbond_sp FFNonbond_O::make_FFNonbond(core::Symbol_sp type,
                                         double radius_nanometers,
                                         double epsilon_kj,
                                         double apol,
                                         double neff,
                                         double mass,
                                         double polarizability,
                                         double initial_charge,
                                         double fcadj,
                                         double pbci,
                                         DonorAcceptorEnum donor_acceptor) {
  GC_ALLOCATE(FFNonbond_O, res);
  res->_Type = type;
  res->_Radius_Nanometers = radius_nanometers;
  res->_Epsilon_kJ = epsilon_kj;
  res->_Apol = apol;
  res->_Neff = neff;
  res->_Mass = mass;
  res->_Polarizability = polarizability;
  res->_InitialCharge = initial_charge;
  res->_Fcadj = fcadj;
  res->_Pbci = pbci;
  res->_DonorAcceptor = donor_acceptor;
  return res;
}

void	FFNonbond_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,type),this->_Type );
  node->field(INTERN_(kw,radius),this->_Radius_Nanometers);
  node->field(INTERN_(kw,well),this->_Epsilon_kJ);
  node->field_if_not_default(INTERN_(kw,apol),this->_Apol,0.0);
  node->field_if_not_default(INTERN_(kw,neff),this->_Neff,0.0);
  node->field(INTERN_(kw,mass),this->_Mass);
  node->field_if_not_default(INTERN_(kw,pol),this->_Polarizability,0.0);
  node->field_if_not_default(INTERN_(kw,initCharge),this->_InitialCharge,0.0);
  node->field_if_not_default(INTERN_(kw,fcadj),this->_Fcadj,0.0);
  node->field_if_not_default(INTERN_(kw,pbci),this->_Pbci,0.0);
  node->field(INTERN_(kw,da),this->_DonorAcceptor );
  node->field(INTERN_(kw,sameParms),this->_SameParms);
  this->Base::fields(node);
}

string FFNonbond_O::__repr__() const {
  stringstream ss;
  ss << "#<" << cl__class_of(this->asSmartPtr())->_classNameAsString() << " :type " << _rep_(this->_Type) << ">";
  return ss.str();
}

void FFNonbondDb_O::forceFieldMergeGlobalParameters(FFNonbondDb_sp other)
{
  if (other->EleDielectricValueDefined) {
    this->set_EleDielectricValue(other->EleDielectricValue);
  }
  if (other->EleBufferDefined) {
    this->set_EleBuffer(other->EleBuffer);
  }
  if (other->EleScale14Defined) {
    this->set_EleScale14(other->EleScale14);
  }
  if (other->VdwScale14Defined) {
    this->set_VdwScale14(other->VdwScale14);
  }
  if (other->VdwScaleBufferADefined) {
    this->set_VdwScaleBufferA(other->VdwScaleBufferA);
  }
  if (other->VdwScaleBufferBDefined) {
    this->set_VdwScaleBufferB(other->VdwScaleBufferB);
  }
  if (other->EleDielectricCodeDefined) {
    this->set_EleDielectricCode(other->EleDielectricCode);
  }
  if (other->VdwMixRadiusDefined) {
    this->set_VdwMixRadius(other->VdwMixRadius);
  }
  if (other->VdwMixWellDefined) {
    this->set_VdwMixWell(other->VdwMixWell);
  }
}

void FFNonbondDb_O::forceFieldMerge(FFBaseDb_sp bother)
{
  FFNonbondDb_sp other = gc::As<FFNonbondDb_sp>(bother);
  this->forceFieldMergeGlobalParameters(other);
  // Merge the terms and overwrite old ones with the same name
  for ( size_t i(0), iEnd(other->_Terms.size()); i<iEnd; ++i ) {
    bool new_type = true;
    for (size_t j(0), jEnd(this->_Terms.size()); j<jEnd; ++j ) {
      if (this->_Terms[j]->getType() == other->_Terms[i]->getType() ) {
        this->_Terms[j] = other->_Terms[i];
        new_type = false;
        break;
      }
    }
    if (new_type) this->_Terms.push_back(other->_Terms[i]);
  }
  this->_Parameters->clrhash();
  for ( size_t i(0), iEnd(this->_Terms.size()); i<iEnd; ++i ) {
//    printf("%s:%d this->_Parameters->_HashTableCount @%p-> %lu\n", __FILE__, __LINE__, (void*)&this->_Parameters->_HashTableCount, this->_Parameters->_HashTableCount );
    this->_Parameters->hash_table_setf_gethash(this->_Terms[i]->getType(),core::make_fixnum(i));
  }
  this->Base::Base::forceFieldMerge(other);
}


CL_DEFUN void chem__FFNonbondDb_merge(FFNonbondDb_sp dest, FFNonbondDb_sp source)
{
  dest->forceFieldMerge(source);
}


#if 0
CL_LISPIFY_NAME(FFNonbondDb-concatenate);
CL_DOCSTRING("Join the other FFNonbondDb_sp to the current one - signal an error if there are any type conflicts");
CL_DEFMETHOD void FFNonbondDb_O::FFNonbondDb_concatenate(FFNonbondDb_sp other, core::Symbol_sp other_name) {
  // What do we do with the nonbond parameters?
  // There needs to be something to control how nonbond interactions work across the whole system
  core::HashTable_sp terms_ht = core::HashTableEq_O::create_default();
  for ( size_t myi=0; myi<this->_Terms.size(); ++myi ) {
    terms_ht->setf_gethash(this->_Terms[myi]->_Type,_lisp->_true());
  }
  for ( size_t otheri=0; otheri<other->_Terms.size(); ++otheri ) {
    FFNonbond_sp term = this->_Terms[otheri];
    core::T_sp type = term->_Type;
    core::T_sp found = terms_ht->gethash(type,_Nil<core::T_O>());
    if (found.notnilp()) {
      SIMPLE_ERROR(BF("The type %s exists in the force-field %s and another force-field - there can be only one") % _rep_(type) % _rep_(other_name));
    }
    this->_Terms.push_back(term);
  }
}
#endif

CL_DOCSTRING("Pass an alist of (nonbond-force-field . force-field-name) and a nonbond-force-field will be returned that combines them all. The global_nonbond defines the global parameters only.");
CL_DEFUN FFNonbondDb_sp chem__combine_nonbond_force_fields(FFNonbondDb_sp global_nonbond, core::List_sp nonbond_force_fields) {
  // What do we do with the nonbond parameters?
  // There needs to be something to control how nonbond interactions work across the whole system
  FFNonbondDb_sp conc = FFNonbondDb_O::create();
  conc->forceFieldMergeGlobalParameters(global_nonbond);
  core::HashTable_sp terms_ht = core::HashTableEq_O::create_default();
  for ( auto cur : nonbond_force_fields ) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    FFNonbondDb_sp nb = gc::As<FFNonbondDb_sp>(CONS_CAR(pair));
    core::T_sp nb_name = CONS_CDR(pair);
    for ( size_t otheri=0; otheri<nb->_Terms.size(); ++otheri ) {
      FFNonbond_sp term = nb->_Terms[otheri];
      core::T_sp type = term->_Type;
      core::T_sp found = terms_ht->gethash(type,_Nil<core::T_O>());
      if (found.notnilp()) {
        SIMPLE_ERROR(BF("The type %s exists in the force-field %s and %s - there can be only one") % _rep_(type) % _rep_(nb_name) % _rep_(found));
      }
      terms_ht->setf_gethash(type,nb_name);
      conc->_Terms.push_back(term);
    }
  }
  return conc;
}
  


void	FFNonbondDb_O::fields(core::Record_sp node)
{
  node->field_if_not_nil(INTERN_(kw,eleChargeFcn), this->EleChargeFcn);
  node->field_if_defined(INTERN_(kw,eleDielectricValue), this->EleDielectricValueDefined, this->EleDielectricValue);
  node->field_if_defined(INTERN_(kw,eleBuffer), this->EleBufferDefined, this->EleBuffer);
  node->field_if_defined(INTERN_(kw,eleScale14), this->EleScale14Defined, this->EleScale14);
  node->field_if_defined(INTERN_(kw,vdwScale14), this->VdwScale14Defined, this->VdwScale14);
  node->field_if_defined(INTERN_(kw,vdwScaleBufferA), this->VdwScaleBufferADefined, this->VdwScaleBufferA);
  node->field_if_defined(INTERN_(kw,vdwScaleBufferB), this->VdwScaleBufferBDefined, this->VdwScaleBufferB);
  node->field_if_defined(INTERN_(kw,eleDielectricCode), this->EleDielectricCodeDefined, this->EleDielectricCode);
  node->field_if_defined(INTERN_(kw,vdwMixRadius), this->VdwMixRadiusDefined, this->VdwMixRadius);
  node->field_if_defined(INTERN_(kw,vdwMixWell), this->VdwMixWellDefined, this->VdwMixWell);
  node->field(INTERN_(kw,terms),this->_Terms);
  this->Base::fields(node);
}

CL_LISPIFY_NAME(FFNonbondDb_add);
CL_DEFMETHOD void    FFNonbondDb_O::add(FFNonbond_sp nb)
{
  uint index = this->_Terms.size();
  this->_Terms.push_back(nb);
  this->_Parameters->setf_gethash(nb->getType(),core::clasp_make_fixnum(index));
}

CL_DEFMETHOD bool    FFNonbondDb_O::hasType(core::Symbol_sp type)
{
  return this->_Parameters->gethash(type).notnilp();
}

CL_DEFMETHOD core::T_sp FFNonbondDb_O::FFNonbond_findType(core::Symbol_sp type)
{
  core::T_sp val = this->_Parameters->gethash(type);
  if ( val.fixnump() ) {
//    printf("%s:%d:%s FFNonbondDb_O::findType type: %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(type).c_str());
    uint index = static_cast<uint>(val.unsafe_fixnum());
    return this->_Terms[index];
  }
  return _Nil<core::T_O>();
}

CL_LISPIFY_NAME("findTypeIndex");
CL_DEFMETHOD uint FFNonbondDb_O::findTypeIndex(core::Symbol_sp type)
{
  core::T_sp index = this->_Parameters->gethash(type);
  if (index.fixnump()) return static_cast<uint>(index.unsafe_fixnum());
  SIMPLE_ERROR(BF("Could not find FFNonbondDb type index for type %s") % _rep_(type));
}

CL_DEFMETHOD size_t FFNonbondDb_O::ffnonbond_find_atom_type_position(core::Symbol_sp type)
{
  size_t ti = this->findTypeIndex(type);
  if ( ti == UndefinedUnsignedInt )
  {
    SIMPLE_ERROR(BF("Unknown type %s") % _rep_(type));
  }
  return ti;
}

uint FFNonbondDb_O::findTypeMajorIndex(core::Symbol_sp type)
{
  uint index = this->findTypeIndex(type);
  if (index == UndefinedUnsignedInt) return index;
  return index*this->_Terms.size();
}

CL_DEFMETHOD FFNonbond_sp FFNonbondDb_O::getFFNonbondUsingTypeIndex(uint typeIdx)
{
  return this->_Terms[typeIdx];
}



CL_LISPIFY_NAME("numberOfTypes");
CL_DEFMETHOD size_t FFNonbondDb_O::numberOfTypes() const
{
  return this->_Terms.size();
}




void FFNonbond_O::initialize()
{
  this->Base::initialize();
  this->_Type = _Nil<core::Symbol_O>();
  this->_Radius_Nanometers = 0.0;
  this->_Epsilon_kJ = 0.0; // Depth of the VDW well
  this->_Apol = 0.0;
  this->_Neff = 0.0;
  this->_Mass = 0.0;
  this->_Polarizability = 0.0;
  this->_DonorAcceptor = daNeither;
  this->_InitialCharge = 0.0;
  this->_Fcadj = 0.0;
  this->_Pbci = 0.0;
  this->_SameParms = "";
}


void FFNonbond_O::setRadius_Nanometers(double n)
{
  this->_Radius_Nanometers = n;
}

void FFNonbond_O::setRadius_Angstroms(double n)
{
  this->_Radius_Nanometers = n/10.0;
}

double FFNonbond_O::getRadius_Nanometers() const
{
  return this->_Radius_Nanometers;
}

double FFNonbond_O::getRadius_Angstroms() const
{
  return this->_Radius_Nanometers*10.0;
}



void FFNonbond_O::setEpsilon_kJ(double kj)
{
  this->_Epsilon_kJ = kj;
}

void FFNonbond_O::setEpsilon_kCal(double kcal)
{
  this->_Epsilon_kJ = kCal_to_kJ(kcal);
}


double FFNonbond_O::getEpsilon_kCal() const
{
  return kJ_to_kCal(this->_Epsilon_kJ);
}

double FFNonbond_O::getEpsilon_kJ() const
{
  return this->_Epsilon_kJ;
}

void FFNonbond_O::setSameParms(string sameparms)
{
  this->_SameParms = sameparms;
}

string FFNonbond_O::getSameParms() const
{
  return this->_SameParms;
}

CL_DEFMETHOD double FFNonbond_O::getMass() const
{
  return this->_Mass;
}


string	FFNonbond_O::levelDescription()
{
  stringstream	desc;
  desc << this->FFParameter_O::levelDescription();
  desc << " nonbond type "<<this->_Type;
  return desc.str();
}

CL_DEFUN FFNonbondDb_sp chem__make_ffnonbond_db() {
  return FFNonbondDb_O::create();
}

};
