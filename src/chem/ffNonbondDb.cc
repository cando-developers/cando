       
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
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/loop.h>
#include <cando/chem/units.h>
#include <clasp/core/wrappers.h>



namespace chem
{

    EXPOSE_CLASS(chem, FFNonbondCrossTermTable_O );
    EXPOSE_CLASS(chem, FFNonbond_O );
    EXPOSE_CLASS(chem, FFNonbondDb_O );


void FFNonbondCrossTermTable_O::exposeCando(core::Lisp_sp e)
{
    core::class_<FFNonbondCrossTermTable_O>()
    ;
}


    void FFNonbondCrossTermTable_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,FFNonbondCrossTermTable,"","",_lisp)
	;
#endif
}
void FFNonbondCrossTermTable_O::fillUsingFFNonbondDb( FFNonbondDb_sp db )
{_G();
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
{_G();
    ASSERT_lt(typeIndex,this->_NumberOfTypes);
    return typeIndex*this->_NumberOfTypes;
}

FFNonbondCrossTerm& FFNonbondCrossTermTable_O::nonbondCrossTermFromTypes(uint type1, uint type2 )
{_G();
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







void FFNonbond_O::exposeCando(core::Lisp_sp e)
{
    core::class_<FFNonbond_O>()
    ;
}

    void FFNonbond_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,FFNonbond,"","",_lisp)
    ;
#endif
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



void FFNonbondDb_O::exposeCando(core::Lisp_sp e)
{
    core::class_<FFNonbondDb_O>()
    .def("findTypeIndex",&FFNonbondDb_O::findTypeIndex)
    .def("numberOfTypes",&FFNonbondDb_O::numberOfTypes)
    ;

    SYMBOL_SC_(ChemPkg,STARDonorAcceptorEnumConverterSTAR);
    DECLARE_SYMBOL_TO_ENUM_CONVERTER(daEnum,"DonorAcceptorEnum",ChemKwPkg,_sym_STARDonorAcceptorEnumConverterSTAR);

    SYMBOL_SC_(ChemPkg,STAREleDielectricEnumConverterSTAR);
    DECLARE_SYMBOL_TO_ENUM_CONVERTER(dielectricEnum,"EleDielectricEnum",ChemKwPkg,_sym_STAREleDielectricEnumConverterSTAR);

    SYMBOL_SC_(ChemPkg,STARVdwMixRadiusEnumConverterSTAR);
    DECLARE_SYMBOL_TO_ENUM_CONVERTER(vmrEnum,"VdwMixRadiusEnum",ChemKwPkg,_sym_STARVdwMixRadiusEnumConverterSTAR);

    SYMBOL_SC_(ChemPkg,STARVdwMixWellEnumConverterSTAR);
    DECLARE_SYMBOL_TO_ENUM_CONVERTER(vmwEnum,"VdwMixWellEnum",ChemKwPkg,_sym_STARVdwMixWellEnumConverterSTAR);





}

    void FFNonbondDb_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,FFNonbondDb,"","",_lisp)
    .def("findTypeIndex",&FFNonbondDb_O::findTypeIndex)
    .def("numberOfTypes",&FFNonbondDb_O::numberOfTypes)
    ;
#endif
    

}





void    FFNonbondDb_O::initialize()
{
    this->Base::initialize();
    this->_Terms.clear();
    this->_Lookup = core::HashTableEq_O::create_default();
    this->_EleDielectricValue = 0.0;
    this->_EleBuffer = 0.0;
    this->_EleScale14 = 0.0;
    this->_EleChargeFcn = _Nil<core::T_O>();
    this->_VdwScale14 = 0.0;
    this->_VdwScaleBufferA = 0.0;
    this->_VdwScaleBufferB = 0.0;
    this->_EleDielectricCode = edConstant;
    this->_VdwMixRadius = vmrAverage;
    this->_VdwMixWell = vmwAverage;
}


void	FFNonbond_O::fields(core::Record_sp node)
{_G();
  this->Base::fields(node);
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
}

void	FFNonbondDb_O::fields(core::Record_sp node)
{_G();
  this->Base::fields(node);
  node->field_if_not_default(INTERN_(kw,eleDielectricValue), this->_EleDielectricValue, 0.0);
  node->field_if_not_default(INTERN_(kw,eleBuffer), this->_EleBuffer, 0.0);
  node->field_if_not_default(INTERN_(kw,eleScale14), this->_EleScale14, 0.0);
  node->field_if_not_nil(INTERN_(kw,eleChargeFcn), this->_EleChargeFcn);
  node->field_if_not_default(INTERN_(kw,vdwScale14), this->_VdwScale14, 0.0);
  node->field_if_not_default(INTERN_(kw,vdwScaleBufferA), this->_VdwScaleBufferA, 0.0);
  node->field_if_not_default(INTERN_(kw,vdwScaleBufferB), this->_VdwScaleBufferB, 0.0);
  node->field(INTERN_(kw,nonbonds),this->_Terms );
  node->field(INTERN_(kw,map),this->_Lookup );
  node->field( INTERN_(kw,eleDielectricCode), this->_EleDielectricCode);
  node->field( INTERN_(kw,vdwMixRadius), this->_VdwMixRadius);
  node->field( INTERN_(kw,vdwMixWell), this->_VdwMixWell);
}


void    FFNonbondDb_O::add(FFNonbond_sp nb)
{_OF();
  if ( this->_Lookup->contains(nb->getType()) ) {  // This was negated
        SIMPLE_ERROR(BF("Adding nonbonded( %s ) to database but it's already there") % nb->getType() );
    }
    uint index = this->_Terms.size();
    this->_Terms.push_back(nb);
    this->_Lookup->setf_gethash(nb->getType(),core::clasp_make_fixnum(index));
}

bool    FFNonbondDb_O::hasType(core::Symbol_sp type)
{
    if ( this->_Lookup->contains(type) ) return true;
    return false;
}

gc::Nilable<FFNonbond_sp>   FFNonbondDb_O::findType(core::Symbol_sp type)
{
  if ( this->_Lookup->contains(type) ) {
//    printf("%s:%d:%s FFNonbondDb_O::findType type: %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(type).c_str());
    core::T_sp val = this->_Lookup->gethash(type);
    uint index = core::clasp_to_fixnum(val);
    return this->_Terms[index];
  }
  return _Nil<core::T_O>();
}

CL_LISPIFY_NAME("findTypeIndex");
CL_DEFMETHOD uint FFNonbondDb_O::findTypeIndex(core::Symbol_sp type)
{
    if ( this->_Lookup->contains(type) ) {
      uint index = core::clasp_to_fixnum(this->_Lookup->gethash(type));
      return index;
    }
    return UndefinedUnsignedInt;
}

uint FFNonbondDb_O::findTypeIndexOrThrow(core::Symbol_sp type)
{_OF();
    uint ti = this->findTypeIndex(type);
    if ( ti == UndefinedUnsignedInt )
    {
	SIMPLE_ERROR(BF("Unknown type %s") % _rep_(type));
    }
    return ti;
}
uint FFNonbondDb_O::findTypeMajorIndex(core::Symbol_sp type)
{
    if ( this->_Lookup->contains(type)) {
      core::T_sp val = this->_Lookup->gethash(type);
      return core::clasp_to_fixnum(val)*this->_Terms.size();
    }
    return UndefinedUnsignedInt;
}

FFNonbond_sp FFNonbondDb_O::getFFNonbondUsingTypeIndex(uint typeIdx)
{_OF();
    return this->_Terms[typeIdx];
}



CL_LISPIFY_NAME("numberOfTypes");
CL_DEFMETHOD uint FFNonbondDb_O::numberOfTypes()
{_G();
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


#if 0 //[
adapt::QDomNode_sp	FFNonbond_O::asXml()
{
    adapt::QDomNode_sp	node;
    node = adapt::QDomNode_O::create(env,"FFNonbond");
    node->addAttributeString("type",this->_Type);
    node->addAttributeDoubleScientific("Radius",this->_Radius);
    node->addAttributeDoubleScientific("Well",this->_Well);
    return node;
}

void	FFNonbond_O::parseFromXml(adapt::QDomNode_sp node)
{
    this->_Type = node->getAttributeString("type");
    this->_Radius = node->getAttributeDouble("Radius");
    this->_Well = node->getAttributeDouble("Well");
}
#endif //]


string	FFNonbond_O::levelDescription()
{
stringstream	desc;
    desc << this->FFParameter_O::levelDescription();
    desc << " nonbond type "<<this->_Type;
    return desc.str();
}


};
