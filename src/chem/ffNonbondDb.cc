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
      double epsilonij = sqrt(ffNonbond1->getEpsilon_kcal()*ffNonbond2->getEpsilon_kcal());
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
{
  ASSERT_lt(idx,this->_CrossTerms.size());
  return this->_CrossTerms[idx];
}

uint FFNonbondCrossTermTable_O::nonbondCrossTermIndexFromTypes(uint type1, uint type2 )
{
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


CL_LAMBDA(type &key (radius-nanometers 0.0) (epsilon-kj 0.0) (apol 0.0) (neff 0.0) (mass 0.0) (polarizability 0.0) (initial-charge 0.0) (fcadj 0.0) (pbci 0.0) (donor-acceptor :da-neither));
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
  auto  res = gctools::GC<FFNonbond_O>::allocate();
  res->_Type = type;
  res->_Radius_Nanometers = radius_nanometers;
  res->_Epsilon_kj = epsilon_kj;
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
  node->field(INTERN_(kw,well),this->_Epsilon_kj);
  node->field_if_not_default(INTERN_(kw,apol),this->_Apol,0.0);
  node->field_if_not_default(INTERN_(kw,neff),this->_Neff,0.0);
  node->field(INTERN_(kw,mass),this->_Mass);
  node->field_if_not_default(INTERN_(kw,pol),this->_Polarizability,0.0);
  node->field_if_not_default(INTERN_(kw,initCharge),this->_InitialCharge,0.0);
  node->field_if_not_default(INTERN_(kw,fcadj),this->_Fcadj,0.0);
  node->field_if_not_default(INTERN_(kw,pbci),this->_Pbci,0.0);
  node->field(INTERN_(kw,da),this->_DonorAcceptor );
  node->field(INTERN_(kw,sameParms),this->_SameParms);
  node->field_if_not_nil(INTERN_(kw,type_comment),this->_TypeComment);
  node->field_if_not_nil(INTERN_(kw,nonbond_comment),this->_NonbondComment);
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
    core::T_sp found = this->_Parameters->gethash(other->_Terms[i]->getType());
    if (found.notnilp()) {
      this->_Terms[found.unsafe_fixnum()] = other->_Terms[i];
    } else {
      this->add(other->_Terms[i]);
    }
  }
  this->Base::Base::forceFieldMerge(other);
}


DOCGROUP(cando);
CL_DEFUN void chem__FFNonbondDb_merge(FFNonbondDb_sp dest, FFNonbondDb_sp source)
{
  dest->forceFieldMerge(source);
}


#if 0
CL_LISPIFY_NAME(FFNonbondDb-concatenate);
CL_DOCSTRING(R"dx(Join the other FFNonbondDb_sp to the current one - signal an error if there are any type conflicts)dx");
CL_DEFMETHOD void FFNonbondDb_O::FFNonbondDb_concatenate(FFNonbondDb_sp other, core::Symbol_sp other_name) {
  // What do we do with the nonbond parameters?
  // There needs to be something to control how nonbond interactions work across the whole system
  core::HashTable_sp terms_ht = core::HashTable_O::createEq();
  for ( size_t myi=0; myi<this->_Terms.size(); ++myi ) {
    terms_ht->setf_gethash(this->_Terms[myi]->_Type,_lisp->_true());
  }
  for ( size_t otheri=0; otheri<other->_Terms.size(); ++otheri ) {
    FFNonbond_sp term = this->_Terms[otheri];
    core::T_sp type = term->_Type;
    core::T_sp found = terms_ht->gethash(type,_Nil<core::T_O>());
    if (found.notnilp()) {
      SIMPLE_ERROR("The type {} exists in the force-field {} and another force-field - there can be only one" , _rep_(type) , _rep_(other_name));
    }
    this->_Terms.push_back(term);
  }
}
#endif

CL_DOCSTRING(R"dx(Pass an alist of (nonbond-force-field . force-field-name) and a nonbond-force-field will be returned that combines them all. The global_nonbond defines the global parameters only.)dx");
DOCGROUP(cando);
CL_DEFUN FFNonbondDb_sp chem__combine_nonbond_force_fields(FFNonbondDb_sp global_nonbond, core::List_sp nonbond_force_fields, core::HashTable_sp atomTypes) {
  // What do we do with the nonbond parameters?
  // There needs to be something to control how nonbond interactions work across the whole system
  FFNonbondDb_sp conc = FFNonbondDb_O::create();
  conc->forceFieldMergeGlobalParameters(global_nonbond);
  core::HashTable_sp terms_ht = core::HashTable_O::createEq();
  for ( auto cur : nonbond_force_fields ) {
    core::Cons_sp pair = gc::As<core::Cons_sp>(CONS_CAR(cur));
    FFNonbondDb_sp nb = gc::As<FFNonbondDb_sp>(CONS_CAR(pair));
    core::T_sp nb_name = CONS_CDR(pair);
    for ( size_t otheri=0; otheri<nb->_Terms.size(); ++otheri ) {
      FFNonbond_sp term = nb->_Terms[otheri];
      core::T_sp type = term->_Type;
      core::T_sp found = terms_ht->gethash(type,nil<core::T_O>());
      if (found.notnilp()) {
        SIMPLE_ERROR("The type {} exists in the force-field {} and {} - there can be only one" , _rep_(type) , _rep_(nb_name) , _rep_(found));
      }
      terms_ht->setf_gethash(type,nb_name);
      conc->add(term);
    }
  }
  return conc;
}
  


void	FFNonbondDb_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,terms),this->_Terms);
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
  this->Base::fields(node);
}


CL_LISPIFY_NAME(FFNonbondDb/termVector);
CL_DEFMETHOD core::SimpleVector_sp FFNonbondDb_O::termVector() const
{
  core::SimpleVector_sp terms = core::SimpleVector_O::make(this->_Terms.size());
  for ( size_t ii=0; ii<this->_Terms.size(); ii++ ) {
    terms->rowMajorAset(ii,this->_Terms[ii]);
  }
  return terms;
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

CL_DEFUN core::T_sp chem__FFNonbond_findType(FFNonbondDb_sp ffNonbondDb, core::Symbol_sp type)
{
  core::T_sp val = ffNonbondDb->_Parameters->gethash(type);
  if ( val.fixnump() ) {
//    printf("%s:%d:%s FFNonbondDb_O::findType type: %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(type).c_str());
    uint index = static_cast<uint>(val.unsafe_fixnum());
    return ffNonbondDb->_Terms[index];
  }
  return nil<core::T_O>();
}

CL_LISPIFY_NAME("findTypeIndex");
CL_DEFMETHOD uint FFNonbondDb_O::findTypeIndex(core::Symbol_sp type)
{
  core::T_sp index = this->_Parameters->gethash(type);
  if (index.fixnump()) return static_cast<uint>(index.unsafe_fixnum());
  SIMPLE_ERROR("Could not find FFNonbondDb type index for type {}" , _rep_(type));
}

CL_DEFMETHOD size_t FFNonbondDb_O::ffnonbond_find_atom_type_position(core::Symbol_sp type)
{
  size_t ti = this->findTypeIndex(type);
  if ( ti == UndefinedUnsignedInt )
  {
    SIMPLE_ERROR("Unknown type {}" , _rep_(type));
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
  this->_Type = nil<core::Symbol_O>();
  this->_Radius_Nanometers = 0.0;
  this->_Epsilon_kj = 0.0; // Depth of the VDW well
  this->_Apol = 0.0;
  this->_Neff = 0.0;
  this->_Mass = 0.0;
  this->_Polarizability = 0.0;
  this->_DonorAcceptor = daNeither;
  this->_InitialCharge = 0.0;
  this->_Fcadj = 0.0;
  this->_Pbci = 0.0;
}


CL_LISPIFY_NAME(FFNonbond/setRadius_Nanometers);
CL_DEFMETHOD
void FFNonbond_O::setRadius_Nanometers(double n)
{
  this->_Radius_Nanometers = n;
}

CL_LISPIFY_NAME(FFNonbond/setRadius_Angstroms);
CL_DEFMETHOD
void FFNonbond_O::setRadius_Angstroms(double n)
{
  this->_Radius_Nanometers = n/10.0;
}

CL_LISPIFY_NAME(FFNonbond/getRadius_Nanometers);
CL_DEFMETHOD
double FFNonbond_O::getRadius_Nanometers() const
{
  return this->_Radius_Nanometers;
}

CL_LISPIFY_NAME(FFNonbond/getRadius_Angstroms);
CL_DEFMETHOD
double FFNonbond_O::getRadius_Angstroms() const
{
  return this->_Radius_Nanometers*10.0;
}



CL_LISPIFY_NAME(FFNonbond/setEpsilon_kj);
CL_DEFMETHOD
void FFNonbond_O::setEpsilon_kj(double kj)
{
  this->_Epsilon_kj = kj;
}

CL_LISPIFY_NAME(FFNonbond/setEpsilon_kcal);
CL_DEFMETHOD
void FFNonbond_O::setEpsilon_kcal(double kcal)
{
  this->_Epsilon_kj = kcal_to_kj(kcal);
}


CL_LISPIFY_NAME(FFNonbond/getEpsilon_kcal);
CL_DEFMETHOD
double FFNonbond_O::getEpsilon_kcal() const
{
  return kj_to_kcal(this->_Epsilon_kj);
}

CL_LISPIFY_NAME(FFNonbond/getEpsilon_kj);
CL_DEFMETHOD
double FFNonbond_O::getEpsilon_kj() const
{
  return this->_Epsilon_kj;
}

CL_LISPIFY_NAME(FFNonbond/setSameParms);
CL_DEFMETHOD
void FFNonbond_O::setSameParms(core::T_sp sameparms)
{
  this->_SameParms = sameparms;
}

CL_LISPIFY_NAME(FFNonbond/getSameParms);
CL_DEFMETHOD
core::T_sp FFNonbond_O::getSameParms() const
{
  // Can return unbound<>()
  return this->_SameParms;
}

CL_LISPIFY_NAME(FFNonbond/getMass);
CL_DEFMETHOD double FFNonbond_O::getMass() const
{
  return this->_Mass;
}


CL_LISPIFY_NAME(FFNonbond/getTypeComment);
CL_DEFMETHOD
core::T_sp FFNonbond_O::getTypeComment() const
{
  return this->_TypeComment;
}

CL_LISPIFY_NAME(FFNonbond/setTypeComment);
CL_DEFMETHOD
void FFNonbond_O::setTypeComment(core::T_sp comment)
{
  this->_TypeComment = comment;
}

CL_LISPIFY_NAME(FFNonbond/getNonbondComment);
CL_DEFMETHOD
core::T_sp FFNonbond_O::getNonbondComment() const
{
  return this->_NonbondComment;
}

CL_LISPIFY_NAME(FFNonbond/setNonbondComment);
CL_DEFMETHOD
void FFNonbond_O::setNonbondComment(core::T_sp comment)
{
  this->_NonbondComment = comment;
}


string	FFNonbond_O::levelDescription()
{
  stringstream	desc;
  desc << this->FFParameter_O::levelDescription();
  desc << " nonbond type "<<this->_Type;
  return desc.str();
}

DOCGROUP(cando);
CL_DEFUN FFNonbondDb_sp chem__make_ffnonbond_db() {
  return FFNonbondDb_O::create();
}



///////////////////////////////////////////////////
//
//
void FFLKSolvation_O::initialize() {
  this->Base::initialize();
  this->_Type = nil<core::Symbol_O>();
  this->_TypeComment = nil<core::T_O>();
  this->_LJ_radius = 0.0;
  this->_LJ_wdepth = 0.0;
  this->_LK_dgfree = 0.0;
  this->_LK_lambda = 1.0;
  this->_LK_volume = 0.0;
}

string FFLKSolvation_O::__repr__() const {
  stringstream ss;
  ss << "#<" << cl__class_of(this->asSmartPtr())->_classNameAsString()
     << " :type " << _rep_(this->_Type) << ">";
  return ss.str();
}

void FFLKSolvation_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,type), this->_Type);
  node->field_if_not_nil(INTERN_(kw,type_comment), this->_TypeComment);
  node->field(INTERN_(kw,lj_radius), this->_LJ_radius);
  node->field(INTERN_(kw,lj_wdepth), this->_LJ_wdepth);
  node->field(INTERN_(kw,lk_dgfree), this->_LK_dgfree);
  node->field(INTERN_(kw,lk_lambda), this->_LK_lambda);
  node->field(INTERN_(kw,lk_volume), this->_LK_volume);
  this->Base::fields(node);
}

CL_LAMBDA(type &key (type-comment nil) (lj-radius 0.0) (lj-wdepth 0.0)
          (lk-dgfree 0.0) (lk-lambda 1.0) (lk-volume 0.0));
  CL_DEF_CLASS_METHOD
  FFLKSolvation_sp FFLKSolvation_O::make_FFLKSolvation(core::Symbol_sp type,
                                                        core::T_sp typeComment,
                                                        double LJ_radius,
                                                        double LJ_wdepth,
                                                        double LK_dgfree,
                                                        double LK_lambda,
                                                        double LK_volume)
  {
    auto res = gctools::GC<FFLKSolvation_O>::allocate_with_default_constructor();
    res->_Type = type;
    res->_TypeComment = typeComment;
    res->_LJ_radius = LJ_radius;
    res->_LJ_wdepth = LJ_wdepth;
    res->_LK_dgfree = LK_dgfree;
    res->_LK_lambda = LK_lambda;
    res->_LK_volume = LK_volume;
    return res;
  }







void FFLKSolvationDb_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,lk_terms),this->_LKTerms);
  this->Base::fields(node);
}

CL_LISPIFY_NAME(FFLKSolvationDb_add);
CL_DEFMETHOD void FFLKSolvationDb_O::add(FFLKSolvation_sp param)
{
  uint index = this->_LKTerms.size();
  this->_LKTerms.push_back(param);
  this->_Parameters->setf_gethash(param->getType(), core::clasp_make_fixnum(index));
}

CL_DEFMETHOD bool FFLKSolvationDb_O::hasType(core::Symbol_sp type)
{
  return this->_Parameters->gethash(type).notnilp();
}

CL_DEFMETHOD FFLKSolvation_sp FFLKSolvationDb_O::getLKSolvation(core::Symbol_sp type)
{
  core::T_sp val = this->_Parameters->gethash(type);
  if (val.fixnump()) {
    uint index = static_cast<uint>(val.unsafe_fixnum());
    return this->_LKTerms[index];
  }
  return nil<FFLKSolvation_O>();
}

CL_DEFUN core::T_sp chem__FFLKSolvation_findType(FFLKSolvationDb_sp ffLKSolvationDb, core::Symbol_sp type)
{
  core::T_sp val = ffLKSolvationDb->_Parameters->gethash(type);
  if (val.fixnump()) {
    uint index = static_cast<uint>(val.unsafe_fixnum());
    return ffLKSolvationDb->_LKTerms[index];
  }
  return nil<core::T_O>();
}


};
