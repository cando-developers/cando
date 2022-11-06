/*
    File: ffNonbondDb.h
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


/*
 *	ffNonbond.h
 *
 *	Maintains a database of nonbonded terms
 */

#ifndef FFNONBONDED_H
#define	FFNONBONDED_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
//#include	"conformation.h"
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/ffBaseDb.h>


#include <cando/chem/chemPackage.h>


namespace chem {


    SMART(FFNonbondDb);

/*! A Nonbond term, this is subclassed in AmberFunction to deal with Amber specific stuff.
 * I should combine energyFunction.cc and amberFunction.cc together.
 */
    class FFNonbondCrossTerm {
    public:
	uint		_Type1Index;
	uint		_Type2Index;
	double		_RStar;
	double		_A;
	double		_C;
    };


    SMART(FFNonbondCrossTermTable);
    class FFNonbondCrossTermTable_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,FFNonbondCrossTermTable_O,"FFNonbondCrossTermTable",core::CxxObject_O);

    public:
    private:
      uint				_NumberOfTypes;
      gctools::Vec0<FFNonbondCrossTerm> _CrossTerms;
    public:
	void fillUsingFFNonbondDb(FFNonbondDb_sp);
    public:
	FFNonbondCrossTerm& nonbondCrossTermFromTypes(uint type1, uint type2 );
	FFNonbondCrossTerm& nonbondCrossTerm(uint idx);

	uint nonbondCrossTermIndexFromTypes(uint type1, uint type2 );	
	uint typeMajorIndex(uint type1);

	DEFAULT_CTOR_DTOR(FFNonbondCrossTermTable_O);
    };




//#define	FFNonbondWildCardCode	-1



    typedef	enum { daDonor, daAcceptor, daNeither } DonorAcceptorEnum ;
    typedef	enum { edConstant, edDistance } EleDielectricEnum;
    typedef	enum { vmrAverage, vmrGeometric, vmrMmff94 } VdwMixRadiusEnum;
    typedef	enum { vmwAverage, vmwMmff94 } VdwMixWellEnum;

    extern core::Symbol_sp& _sym_STARDonorAcceptorEnumConverterSTAR;
    extern core::Symbol_sp& _sym_STAREleDielectricEnumConverterSTAR;
    extern core::Symbol_sp& _sym_STARVdwMixRadiusEnumConverterSTAR;
    extern core::Symbol_sp& _sym_STARVdwMixWellEnumConverterSTAR;
};



ENUM_TRANSLATOR(chem::DonorAcceptorEnum,chem::_sym_STARDonorAcceptorEnumConverterSTAR);
ENUM_TRANSLATOR(chem::EleDielectricEnum,chem::_sym_STAREleDielectricEnumConverterSTAR);
ENUM_TRANSLATOR(chem::VdwMixRadiusEnum,chem::_sym_STARVdwMixRadiusEnumConverterSTAR);
ENUM_TRANSLATOR(chem::VdwMixWellEnum,chem::_sym_STARVdwMixWellEnumConverterSTAR);
namespace chem {

  SMART(FFNonbond);
  class FFNonbond_O : public FFParameter_O
  {
    LISP_CLASS(chem,ChemPkg,FFNonbond_O,"FFNonbond",FFParameter_O);
  public:
    void initialize();
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    core::Symbol_sp 	_Type;
    double		_Radius_Nanometers;
    double		_Epsilon_kj;	// epsilon in AMBER
    double		_Apol;
    double		_Neff;
    double		_Mass;
    double              _Polarizability;
    double		_InitialCharge;
    double		_Fcadj;	// formal charge sharing scale factor
    double		_Pbci;	// partial bond charge increment
    DonorAcceptorEnum	_DonorAcceptor;
    core::T_sp	        _SameParms;
    core::T_sp          _TypeComment;
    core::T_sp          _NonbondComment;

  public:
    static FFNonbond_sp make_FFNonbond(core::Symbol_sp type,
                                       double radius_nanometers,
                                       double epsilon_kj,
                                       double apol,
                                       double neff,
                                       double mass,
                                       double polarizability,
                                       double initial_charge,
                                       double fcadj,
                                       double pbci,
                                       DonorAcceptorEnum donor_acceptor);

  public:
    string __repr__() const;

    CL_LISPIFY_NAME(FFNonbond/setType);
    CL_DEFMETHOD void setType(core::Symbol_sp s) {this->_Type=s;};
    CL_LISPIFY_NAME(FFNonbond/getType);
    CL_DEFMETHOD core::Symbol_sp getType() const { return this->_Type;};
    CL_LISPIFY_NAME(FFNonbond/setPolarizability);
    CL_DEFMETHOD void setPolarizability(double d) { this->_Polarizability = d;};
    CL_LISPIFY_NAME(FFNonbond/getPolarizability);
    CL_DEFMETHOD double getPolarizability() { return this->_Polarizability;};
    CL_LISPIFY_NAME(FFNonbond/setMass);
    CL_DEFMETHOD void setMass(double d) { this->_Mass = d;};

    double getMass() const;
    void setRadius_Angstroms(double val);
    void setRadius_Nanometers(double val);
    double getRadius_Angstroms() const;
    double getRadius_Nanometers() const;
    void setEpsilon_kj(double val);
    void setEpsilon_kcal(double val);
    double getEpsilon_kj() const;
    double getEpsilon_kcal() const;

    void setSameParms(core::T_sp sameparms);
    core::T_sp getSameParms() const;

    core::T_sp getTypeComment() const;
    void setTypeComment(core::T_sp comment);
    core::T_sp getNonbondComment() const;
    void setNonbondComment(core::T_sp comment);

    virtual string	levelDescription();
    virtual ParameterType type() { return nonbond; };
    FFNonbond_O() :
        _SameParms(unbound<core::T_O>()),
        _TypeComment(nil<core::T_O>()),
        _NonbondComment(nil<core::T_O>())
    {};
  };




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
  SMART(FFNonbondDb);
  class FFNonbondDb_O : public FFParameterBaseDb_O
  {
    LISP_CLASS(chem,ChemPkg,FFNonbondDb_O,"FFNonbondDb",FFParameterBaseDb_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    core::T_sp                  _Name;
    core::Symbol_sp		EleChargeFcn;
    double			EleDielectricValue;
    bool                        EleDielectricValueDefined;
    double			EleBuffer;
    bool                        EleBufferDefined;
    double			EleScale14;
    bool                        EleScale14Defined;
    double			VdwScale14;
    bool                        VdwScale14Defined;
    double			VdwScaleBufferA;
    bool                        VdwScaleBufferADefined;
    double			VdwScaleBufferB;
    bool                        VdwScaleBufferBDefined;
    EleDielectricEnum	        EleDielectricCode;
    bool                        EleDielectricCodeDefined;
    VdwMixRadiusEnum	        VdwMixRadius;
    bool                        VdwMixRadiusDefined;
    VdwMixWellEnum		VdwMixWell;
    bool                        VdwMixWellDefined;
    gctools::Vec0<FFNonbond_sp>	_Terms;
  public:

    void set_EleChargeFcn(core::Symbol_sp v) {
      this->EleChargeFcn = v;
    }
    void set_EleDielectricValue(double v) {
      this->EleDielectricValue = v;
      this->EleDielectricValueDefined = true;
    }
  
    void set_EleBuffer(double v) {
      this->EleBuffer = v;
      this->EleBufferDefined = true;
    }
    void set_EleScale14(  double v) {
      this->EleScale14 = v;
      this->EleScale14Defined = true;
    }
    void set_VdwScale14(double v) {
      this->VdwScale14 = v;
      this->VdwScale14Defined = true;
    }
    void set_VdwScaleBufferA(double v) {
      this->VdwScaleBufferA = v;
      this->VdwScaleBufferADefined = true;
    }
    void set_VdwScaleBufferB(double v) {
      this->VdwScaleBufferB = v;
      this->VdwScaleBufferBDefined = true;
    }
    void set_EleDielectricCode(EleDielectricEnum v) {
      this->EleDielectricCode = v;
      this->EleDielectricCodeDefined = true;
    }
    void set_VdwMixRadius(VdwMixRadiusEnum v) {
      this->VdwMixRadius = v;
      this->VdwMixRadiusDefined = true;
    }
    void set_VdwMixWell(VdwMixWellEnum v) {
      this->VdwMixWell = v;
      this->VdwMixWellDefined = true;
    }

  public:
    core::SimpleVector_sp termVector() const;
    void forceFieldMergeGlobalParameters(FFNonbondDb_sp other);
    void forceFieldMerge(FFBaseDb_sp other);
    
    void add( FFNonbond_sp nonbonded );
    bool hasType( core::Symbol_sp type);
    core::T_sp FFNonbond_findType( core::Symbol_sp type );
    FFNonbond_sp getFFNonbondUsingTypeIndex( uint typeIdx );
    uint findTypeIndex( core::Symbol_sp type );
    size_t ffnonbond_find_atom_type_position( core::Symbol_sp type );
    /*! Return the type index multiplied by the number
     * of type indices.  Another type index can be added
     * to this number to get an index into a Nonbond crossterm table
     */
    uint findTypeMajorIndex( core::Symbol_sp type );

    /*! Return the number of Nonbond atom types
     */
    size_t  numberOfTypes() const;

    FFNonbondDb_O() : _Name(nil<core::T_O>()),
                      EleDielectricValueDefined(false),
                      EleBufferDefined(false),
                      EleScale14Defined(false),
                      VdwScale14Defined(false),
                      VdwScaleBufferADefined(false),
                      VdwScaleBufferBDefined(false),
                      EleDielectricCodeDefined(false),
                      VdwMixRadiusDefined(false),
                      VdwMixWellDefined(false) {};
  };



};
#endif
