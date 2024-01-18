/*
    File: energyDihedral.cc
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
#define	DEBUG_LEVEL_NONE


#include <cando/main/extension.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/ql.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffPtorDb.h>
#include <cando/chem/ffItorDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/wrappers.h>

#define VEC8(v) {v,v,v,v,v,v,v,v}
#define VEC4(v) {v,v,v,v}
#define VEC2(v) {v,v}
  
namespace chem {

SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn1);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn2);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn3);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn4);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn5);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn6);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralUnknown);

core::T_sp dihedral_type(int in) {
  switch (in) {
  case 1:
      return chem::_sym_EnergyDihedralIn1;
  case 2:
      return chem::_sym_EnergyDihedralIn2;
  case 3:
      return chem::_sym_EnergyDihedralIn3;
  case 4:
      return chem::_sym_EnergyDihedralIn4;
  case 5:
      return chem::_sym_EnergyDihedralIn5;
  case 6:
      return chem::_sym_EnergyDihedralIn6;
  }
  return chem::_sym_EnergyDihedralUnknown;
}

#undef DIHEDRAL_APPLY_ATOM_MASK
#define DIHEDRAL_APPLY_ATOM_MASK(I1,I2,I3,I4) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         && bitvectorActiveAtomMask->testBit(I4/3)) \
    ) goto SKIP_term;
#define USE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,dihedral_type(IN), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3), core::make_fixnum(I4)); \
    }

#define IGNORE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) {}

core::List_sp EnergyDihedral::encode() const {
  ql::list result;
  result << core::Cons_O::create(INTERN_(kw,v), core::clasp_make_double_float(this->term.V))
         << core::Cons_O::create(INTERN_(kw,in), core::make_fixnum(this->term.IN))
         << core::Cons_O::create(INTERN_(kw,phase), core::clasp_make_double_float(this->_PhaseRad))
         << core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1))
         << core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.I2))
         << core::Cons_O::create(INTERN_(kw,i3), core::make_fixnum(this->term.I3))
         << core::Cons_O::create(INTERN_(kw,i4), core::make_fixnum(this->term.I4))
         << core::Cons_O::create(INTERN_(kw,proper), _lisp->_boolean(this->_Proper))
         << core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1)
         << core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2)
         << core::Cons_O::create(INTERN_(kw,atom3), this->_Atom3)
         << core::Cons_O::create(INTERN_(kw,atom4), this->_Atom4)
    ;
  return result.cons();
}

void EnergyDihedral::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergyDihedral");
}

#ifdef XML_ARCHIVE
void	EnergyDihedral::archive(core::ArchiveP node)
{
//    node->attribute("_Type1",this->_Type1);
//    node->attribute("_Type2",this->_Type2);
//    node->attribute("_Type3",this->_Type3);
//    node->attribute("_Type4",this->_Type4);
  node->attribute("_Proper",this->_Proper);
  node->attribute("_PhaseRad",this->_PhaseRad);
  node->attribute("sinPhase",this->term.sinPhase);
  node->attribute("cosPhase",this->term.cosPhase);
  node->attribute("V",this->term.V);
  node->attribute("DN",this->term.DN);
  node->attribute("IN",this->term.IN);
  node->attribute("I1",this->term.I1);
  node->attribute("I2",this->term.I2);
  node->attribute("I3",this->term.I3);
  node->attribute("I4",this->term.I4);
  node->attribute("a1",this->_Atom1);
  node->attribute("a2",this->_Atom2);
  node->attribute("a3",this->_Atom3);
  node->attribute("a4",this->_Atom4);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
  node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
  node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
  node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include <cando/chem/energy_functions/_Dihedral_debugEvalSerialize.cc>
#endif //]
}
#endif

template <int N,typename Real>
struct SinCos {
  static void sinNPhiCosNPhi(Real& sinNPhi, Real& cosNPhi,
                      Real sinPhi, Real cosPhi ) {
    Real sinNm1Phi;
    Real cosNm1Phi;
    SinCos<N-1,Real>::sinNPhiCosNPhi(sinNm1Phi,cosNm1Phi,sinPhi,cosPhi);
    sinNPhi = cosPhi*sinNm1Phi+sinPhi*cosNm1Phi;
    cosNPhi = cosPhi*cosNm1Phi-sinPhi*sinNm1Phi;
  }
};

template <typename Real>
struct SinCos<1,Real> {
  static void sinNPhiCosNPhi(Real& sinNPhi, Real& cosNPhi,
                      Real sinPhi, Real cosPhi ) {
    sinNPhi = sinPhi;
    cosNPhi = cosPhi;
  }
};


template <typename Real>
void	sinNPhiCosNPhi(int n, Real& sinNPhi, Real& cosNPhi,
                       Real sinPhi, Real cosPhi )
{
  switch (n) {
  case 1:
      SinCos<1,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 2:
      SinCos<2,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 3:
      SinCos<3,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 4:
      SinCos<4,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 5:
      SinCos<5,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 6:
      SinCos<6,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  default:
      SIMPLE_ERROR("N {} should always be 1-6", n);
      break;
  }
}


void EnergyDihedral::defineFrom( int n, FFPtor_sp ffterm , EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale)
{
//    this->_Term = term;
  this->_Proper = true;
//    this->_Type1 = term->_T1;
//    this->_Type2 = term->_T2;
//    this->_Type3 = term->_T3;
//    this->_Type4 = term->_T4;
  this->_Atom1 = ea1->atom();
  this->_Atom2 = ea2->atom();
  this->_Atom3 = ea3->atom();
  this->_Atom4 = ea4->atom();
  this->_PhaseRad = ffterm->getPhaseRadians(n);
  this->term.sinPhase = sin(this->_PhaseRad);
  this->term.cosPhase = cos(this->_PhaseRad);
  this->term.V = ffterm->getV_kcal(n)*scale;
  this->term.DN = n;
  this->term.IN = n;
  this->term.I1 = ea1->coordinateIndexTimes3();
  this->term.I2 = ea2->coordinateIndexTimes3();
  this->term.I3 = ea3->coordinateIndexTimes3();
  this->term.I4 = ea4->coordinateIndexTimes3();
}

double	EnergyDihedral::getDihedral()
{
  Vector3	pos1, pos2, pos3, pos4;
  num_real	angle;
  pos1 = this->_Atom1->getPosition();
  pos2 = this->_Atom2->getPosition();
  pos3 = this->_Atom3->getPosition();
  pos4 = this->_Atom4->getPosition();
  angle = geom::calculateDihedral( pos1, pos2, pos3, pos4);
  return angle;
}

double	EnergyDihedral::getDihedralDeviation()
{
  Vector3	pos1, pos2, pos3, pos4;
  num_real	phi, dev;
  pos1 = this->_Atom1->getPosition();
  pos2 = this->_Atom2->getPosition();
  pos3 = this->_Atom3->getPosition();
  pos4 = this->_Atom4->getPosition();
  phi = geom::calculateDihedral( pos1, pos2, pos3, pos4);
  dev = 1.0+cos(this->term.DN*phi-this->_PhaseRad);
  return dev;
}


void EnergyDihedral::defineFrom( int n, FFItor_sp term , EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale)
{
//    this->_Term = term;
  this->_Proper = false;
//    this->_Type1 = term->_T1;
//    this->_Type2 = term->_T2;
//    this->_Type3 = term->_T3;
//    this->_Type4 = term->_T4;
  this->_Atom1 = ea1->atom();
  this->_Atom2 = ea2->atom();
  this->_Atom3 = ea3->atom();
  this->_Atom4 = ea4->atom();
  this->term.DN = n;
  this->term.IN = n;
  this->term.V = term->getV_kcal(n)*scale;
  this->_PhaseRad = term->getPhaseRadians(n);
  this->term.cosPhase = cos(this->_PhaseRad);
  this->term.sinPhase = sin(this->_PhaseRad);
  this->term.I1 = ea1->coordinateIndexTimes3();
  this->term.I2 = ea2->coordinateIndexTimes3();
  this->term.I3 = ea3->coordinateIndexTimes3();
  this->term.I4 = ea4->coordinateIndexTimes3();
}

#if 0 //[
void	EnergyDihedral::defineMissingProper( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4 )
{
  this->_Proper = true;
//    this->_Type1 = ea1->_Atom->getType(atomTypes);
//    this->_Type2 = ea2->_Atom->getType(atomTypes);
  /    this->_Type3 = ea3->_Atom->getType(atomTypes);
  this->_Type4 = ea4->_Atom->getType(atomTypes);
}
#endif //]


#if 0
adapt::QDomNode_sp	EnergyDihedral::asXml()
{
  adapt::QDomNode_sp	node;
  Vector3	vdiff;

  node = adapt::QDomNode_O::create(env,"EnergyDihedral");
  node->addAttributeBool("proper",this->_Proper );
  node->addAttributeString("atom1Name",this->_Atom1->getName());
  node->addAttributeString("atom2Name",this->_Atom2->getName());
  node->addAttributeString("atom3Name",this->_Atom3->getName());
  node->addAttributeString("atom4Name",this->_Atom4->getName());
//    node->addAttributeString("atom1Type",this->_Type1 );
//    node->addAttributeString("atom2Type",this->_Type2 );
//    node->addAttributeString("atom3Type",this->_Type3 );
//    node->addAttributeString("atom4Type",this->_Type4 );
  node->addAttributeInt("I1",this->term.I1);
  node->addAttributeInt("I2",this->term.I2);
  node->addAttributeInt("I3",this->term.I3);
  node->addAttributeInt("I4",this->term.I4);
  node->addAttributeDoubleScientific("DN",this->term.DN);
  node->addAttributeInt("IN",this->term.IN);
  node->addAttributeDoubleScientific("V",this->term.V);
  node->addAttributeDoubleScientific("PhaseRad",this->_PhaseRad);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
  adapt::QDomNode_sp xml = adapt::QDomNode_O::create(env,"Evaluated");
  xml->addAttributeBool("calcForce",this->_calcForce );
  xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
  xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include <_Dihedral_debugEvalXml.cc>
  node->addChild(xml);
#endif
  return node;
}

void	EnergyDihedral::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
                                                   AtomTable_sp at )
{
  this->term.I1 = xml->getAttributeInt("I1");
  this->term.I2 = xml->getAttributeInt("I2");
  this->term.I3 = xml->getAttributeInt("I3");
  this->term.I4 = xml->getAttributeInt("I4");
  this->_Atom1 = at->findEnergyAtomWithCoordinateIndex(this->term.I1)->atom();
  this->_Atom2 = at->findEnergyAtomWithCoordinateIndex(this->term.I2)->atom();
  this->_Atom3 = at->findEnergyAtomWithCoordinateIndex(this->term.I3)->atom();
  this->_Atom4 = at->findEnergyAtomWithCoordinateIndex(this->term.I4)->atom();
  this->term.DN = xml->getAttributeDouble("DN");
  this->term.IN = xml->getAttributeInt("IN");
  this->term.V = xml->getAttributeDouble("V");
  this->_PhaseRad = xml->getAttributeDouble("PhaseRad");
  this->_Proper = xml->getAttributeBool("proper");
}
#endif


//
// Copy this from implementAmberFunction.cc
//
num_real	_evaluateEnergyOnly_Dihedral(
                                     num_real x1, num_real y1, num_real z1,
                                     num_real x2, num_real y2, num_real z2,
                                     num_real x3, num_real y3, num_real z3,
                                     num_real x4, num_real y4, num_real z4,
                                     num_real V, num_real DN, int IN,
                                     num_real cosPhase, num_real sinPhase )
{
  IMPLEMENT_ME();
#if 0
  num_real	EraseLinearDihedral;
  num_real	SinNPhi, CosNPhi;

#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) {}
#undef	DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(i,o,v) {}
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	DIHEDRAL_CALC_FORCE	// Don't calculate FORCE or HESSIAN


#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
#pragma clang diagnostic pop
  fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
  fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
  fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
  fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
#undef VEC_CONST
#define VEC_CONST(x) (x)
#undef ZERO_SMALL_LEN
#define ZERO_SMALL_LEN(RL,L) {double fabs_ = fabs(L); int cmp = (fabs_ < TENM3); RL = cmp ? 0.0 : RL; }
#undef DO_sinNPhiCosNPhi
#define DO_sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi) sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi)
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>
#undef DIHEDRAL_APPLY_ATOM_MASK

  return Energy;
#endif
}


void EnergyDihedral_O::addTerm(const EnergyDihedral& term)
{
  if (this->_Terms.size() == this->_Terms.capacity()-1) {
    this->_Terms.reserve(this->_Terms.size()*2);
  }
  this->_Terms.push_back(term);
}


string EnergyDihedral_O::beyondThresholdInteractionsAsString()
{
  return component_beyondThresholdInteractionsAsString<EnergyDihedral_O,EnergyDihedral>(*this);
}



void	EnergyDihedral_O::dumpTerms(core::HashTable_sp atomTypes)
{
  gctools::Vec0<EnergyDihedral>::iterator	edi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4, type;
  for (edi=this->_Terms.begin();edi!=this->_Terms.end();edi++)
  {
    if ( edi->term.V == 0.0 ) continue;
    if ( edi->_Proper ) {
      type = "1PROPER  ";
    } else {
      type = "9IMPROPER";
    }
    as1 = atomLabel(edi->_Atom1);
    as2 = atomLabel(edi->_Atom2);
    as3 = atomLabel(edi->_Atom3);
    as4 = atomLabel(edi->_Atom4);
    if ( edi->_Proper ) {
      if ( as1 < as4 ) {
        str1 = as1;
        str2 = as2;
        str3 = as3;
        str4 = as4;
      } else {
        str4 = as1;
        str3 = as2;
        str2 = as3;
        str1 = as4;
      }
    } else {
      stringstream ss;
      ss.str("");
      ss << as3 << "#Cr";
      str1 = ss.str();
      ss.str("");
      ss << as1 << "#A1";
      str2 = ss.str();
      ss.str("");
      ss << as2 << "#A2";
      str3 = ss.str();
      ss.str("");
      ss << as4 << "#A4";
      str4 = ss.str();
    }
    core::clasp_write_string(fmt::format("TERM 3DIH {} {:<9} - {:<9} - {:<9} - {:<9} {:8.2f} {:8.2f} {:2.0f}\n" , type , str1 , str2 , str3 , str4 , edi->term.V , edi->_PhaseRad , edi->term.DN));
  }
}


CL_DEFMETHOD
core::List_sp	EnergyDihedral_O::lookupDihedralTerms(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, core::HashTable_sp atomTypes )
{
  gctools::Vec0<EnergyDihedral>::iterator	edi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4, type;
  ql::list  result;
  core::T_sp tia1 = atomTable->_AtomTableIndices->gethash(a1);
  core::T_sp tia2 = atomTable->_AtomTableIndices->gethash(a2);
  core::T_sp tia3 = atomTable->_AtomTableIndices->gethash(a3);
  core::T_sp tia4 = atomTable->_AtomTableIndices->gethash(a4);
  if (!tia1.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a1));
  if (!tia2.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a2));
  if (!tia3.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a3));
  if (!tia4.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a4));
  for (edi=this->_Terms.begin();edi!=this->_Terms.end();edi++) {
    if ((edi->_Atom1==a1 &&
         edi->_Atom2==a2 &&
         edi->_Atom3==a3 &&
         edi->_Atom4==a4)
        || (edi->_Atom4==a1 &&
            edi->_Atom3==a2 &&
            edi->_Atom2==a3 &&
            edi->_Atom1==a4)) {
      ql::list oneResult;
      oneResult << INTERN_(kw,type1) << edi->_Atom1->getType(atomTypes)
                << INTERN_(kw,type2) << edi->_Atom2->getType(atomTypes)
                << INTERN_(kw,type3) << edi->_Atom3->getType(atomTypes)
                << INTERN_(kw,type4) << edi->_Atom4->getType(atomTypes)
                << INTERN_(kw,proper) << _lisp->_boolean(edi->_Proper)
                << INTERN_(kw,multiplicity) << core::clasp_make_fixnum(edi->term.IN)
                << INTERN_(kw,phase_degrees) << core::clasp_make_double_float(edi->_PhaseRad)
                << INTERN_(kw,v) << core::clasp_make_double_float(edi->term.V);
      result << oneResult.result();
    }
  }
  return result.result();
}



void	EnergyDihedral_O::setupHessianPreconditioner(
                                                     NVector_sp nvPosition,
                                                     AbstractLargeSquareMatrix_sp m,
                                                     core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  core::T_sp debugInteractions = nil<core::T_O>();
  MAYBE_SETUP_DEBUG_INTERACTIONS(false);
  bool		calcForce = true;
  bool		calcDiagonalHessian = true;
  bool		calcOffDiagonalHessian = true;

//
// Copy from implementAmberFunction::setupHessianPreconditioner
//
// -----------------------

#undef ZERO_SMALL_LEN
#define ZERO_SMALL_LEN(RL,L) {double fabs_ = fabs(L); int cmp = (fabs_ < TENM3); RL = cmp ? 0.0 : RL; }
#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x=di->term.x;}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x=nvPosition->element(ii+of);}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) {}
#undef	DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(i,o,v) {}
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
    m->addToElement((i1)+(o1),(i2)+(o2),v);\
  }
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
    m->addToElement((i1)+(o1),(i2)+(o2),v);\
  }
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include	<cando/chem/energy_functions/_Dihedral_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
    fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
    num_real EraseLinearDihedral;
    int	I1, I2, I3, I4, IN;
#define VEC_CONST(v) v
#undef DO_sinNPhiCosNPhi
#define DO_sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi) sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi)
    num_real sinPhase, cosPhase, SinNPhi, CosNPhi;
    for ( gctools::Vec0<EnergyDihedral>::iterator di=this->_Terms.begin();
          di!=this->_Terms.end(); di++ ) {
#define DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) USE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4)
#include	<cando/chem/energy_functions/_Dihedral_termCode.cc>
#undef DIHEDRAL_DEBUG_INTERACTIONS
#undef DIHEDRAL_APPLY_ATOM_MASK
    }
  }
}

num_real	EnergyDihedral_O::evaluateAllComponentSingle(
    gctools::Vec0<EnergyDihedral>::iterator di_start,
    gctools::Vec0<EnergyDihedral>::iterator di_end,
    ScoringFunction_sp score,
    NVector_sp 	pos,
    bool 		calcForce,
    gc::Nilable<NVector_sp> 	force,
    bool		calcDiagonalHessian,
    bool		calcOffDiagonalHessian,
    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
    gc::Nilable<NVector_sp>	hdvec,
    gc::Nilable<NVector_sp> dvec,
    core::T_sp activeAtomMask,
    core::T_sp debugInteractions )
{
#undef DIHEDRAL_APPLY_ATOM_MASK
#define DIHEDRAL_APPLY_ATOM_MASK(I1,I2,I3,I4) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         && bitvectorActiveAtomMask->testBit(I4/3)) \
    ) goto SKIP_term_and_angle_test;
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(("%s {\n") , this->className());
  }
  num_real totalEnergy = 0.0;
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());


//
// Copy from implementAmberFunction::evaluateAll
//
// -----------------------
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x=di->term.x;}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) {ASSERT_NOT_NAN(e); totalEnergy += (e); }
#undef	DIHEDRAL_FORCE_ACCUMULATE
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE 		ForceAcc
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
#pragma clang diagnostic pop
  fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
  fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
  fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
  fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
  num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
  num_real EraseLinearDihedral;
  int	I1, I2, I3, I4, IN;
  num_real sinPhase, cosPhase, SinNPhi, CosNPhi;
  gctools::Vec0<EnergyDihedral>::iterator di;

  int i = 0;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpass-failed"
  
#pragma clang loop vectorize(enable)
#pragma clang loop interleave(enable)

  for ( di=di_start; di!= di_end; di++ ) {
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
    if ( this->_Debug_NumberOfDihedralTermsToCalculate > 0 ) {
      if ( i>= this->_Debug_NumberOfDihedralTermsToCalculate ) {
        break;
      }
    }
#endif
#undef DO_sinNPhiCosNPhi
#define DO_sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi) sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi)
#define DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) USE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4)
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>
#undef DIHEDRAL_DEBUG_INTERACTIONS
#undef DIHEDRAL_APPLY_ATOM_MASK
    if ( EraseLinearDihedral == 0.0 ) {
      ERROR(chem::_sym_LinearDihedralError,core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(di->_Atom1,di->_Atom2,di->_Atom3,di->_Atom4),
                                                                    kw::_sym_coordinates,pos,
                                                                    kw::_sym_indices,core::Cons_O::createList(core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3), core::make_fixnum(I4))));
    }
  SKIP_term_and_angle_test: (void)0;

#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
    di->_calcForce = calcForce;
    di->_calcDiagonalHessian = calcDiagonalHessian;
    di->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ di->eval.var=val;};
#include	<cando/chem/energy_functions/_Dihedral_debugEvalSet.cc>
#endif //]

    if ( this->_DebugEnergy ) 
    {
      LOG_ENERGY(( "MEISTER dihedral %d args cando\n") , (i+1) );
      LOG_ENERGY(( "MEISTER dihedral %d V %lf\n") , (i+1) , V );
      LOG_ENERGY(( "MEISTER dihedral %d DN %lf\n") , (i+1) , DN );
      LOG_ENERGY(( "MEISTER dihedral %d IN %d\n") , (i+1) , IN );
//		    LOG_ENERGY(( "MEISTER dihedral %d phase %lf\n") , (i+1) , phase );
      LOG_ENERGY(( "MEISTER dihedral %d sinPhase %lf\n") , (i+1) , sinPhase );
      LOG_ENERGY(( "MEISTER dihedral %d cosPhase %lf\n") , (i+1) , cosPhase );
      LOG_ENERGY(( "MEISTER dihedral %d x1 %5.3lf %d\n") , (i+1) , x1 , (I1/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d y1 %5.3lf %d\n") , (i+1) , y1 , (I1/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d z1 %5.3lf %d\n") , (i+1) , z1 , (I1/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d x2 %5.3lf %d\n") , (i+1) , x2 , (I2/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d y2 %5.3lf %d\n") , (i+1) , y2 , (I2/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d z2 %5.3lf %d\n") , (i+1) , z2 , (I2/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d x3 %5.3lf %d\n") , (i+1) , x3 , (I3/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d y3 %5.3lf %d\n") , (i+1) , y3 , (I3/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d z3 %5.3lf %d\n") , (i+1) , z3 , (I3/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d x4 %5.3lf %d\n") , (i+1) , x4 , (I4/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d y4 %5.3lf %d\n") , (i+1) , y4 , (I4/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d z4 %5.3lf %d\n") , (i+1) , z4 , (I4/3+1) );
      LOG_ENERGY(( "MEISTER dihedral %d results\n") , (i+1) );
	//	LOG_ENERGY(( "MEISTER dihedral %d Phi %lf (%lf degrees)\n") , (i+1) %
	//		    	Phi, Phi/0.0174533 );
      LOG_ENERGY(( "MEISTER dihedral %d Energy %lf\n") , (i+1) , Energy);
      if ( calcForce ) {
        LOG_ENERGY(( "MEISTER dihedral %d fx1 %8.5lf %d\n") , (i+1) , fx1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fy1 %8.5lf %d\n") , (i+1) , fy1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fz1 %8.5lf %d\n") , (i+1) , fz1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fx2 %8.5lf %d\n") , (i+1) , fx2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fy2 %8.5lf %d\n") , (i+1) , fy2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fz2 %8.5lf %d\n") , (i+1) , fz2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fx3 %8.5lf %d\n") , (i+1) , fx3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fy3 %8.5lf %d\n") , (i+1) , fy3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fz3 %8.5lf %d\n") , (i+1) , fz3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fx4 %8.5lf %d\n") , (i+1) , fx4 , (I4/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fy4 %8.5lf %d\n") , (i+1) , fy4 , (I4/3+1) );
        LOG_ENERGY(( "MEISTER dihedral %d fz4 %8.5lf %d\n") , (i+1) , fz4 , (I4/3+1) );
      }
      LOG_ENERGY(( "MEISTER dihedral %d stop\n") , (i+1) );
    }
			/* Add the forces */

    if ( calcForce ) {
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx3>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy3>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz3>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx4>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy4>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz4>10000.0);
    }
    i++;
  }
#pragma diagnostic pop

  if ( this->_DebugEnergy )
  {
    LOG_ENERGY(("%s }\n") , this->className());
  }
  return totalEnergy;
}






void	EnergyDihedral_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                                 NVector_sp 	pos)
{
  IMPLEMENT_ME();
  #if 0
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm
//
//------------------
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x=di->term.x;}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) {}
#undef	DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(i,o,v) {}
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    fx3 = 0.0; fy3 = 0.0; fz3 = 0.0;
    fx4 = 0.0; fy4 = 0.0; fz4 = 0.0;
#if !USE_EXPLICIT_DECLARES
    num_real dhx1x1 = 0.0;
    num_real ohx1y1 = 0.0;
    num_real ohx1z1 = 0.0;
    num_real ohx1x2 = 0.0;
    num_real ohx1y2 = 0.0;
    num_real ohx1z2 = 0.0;
    num_real ohx1x3 = 0.0;
    num_real ohx1y3 = 0.0;
    num_real ohx1z3 = 0.0;
    num_real ohx1x4 = 0.0;
    num_real ohx1y4 = 0.0;
    num_real ohx1z4 = 0.0;
    num_real dhy1y1 = 0.0;
    num_real ohy1z1 = 0.0;
    num_real ohy1x2 = 0.0;
    num_real ohy1y2 = 0.0;
    num_real ohy1z2 = 0.0;
    num_real ohy1x3 = 0.0;
    num_real ohy1y3 = 0.0;
    num_real ohy1z3 = 0.0;
    num_real ohy1x4 = 0.0;
    num_real ohy1y4 = 0.0;
    num_real ohy1z4 = 0.0;
    num_real dhz1z1 = 0.0;
    num_real ohz1x2 = 0.0;
    num_real ohz1y2 = 0.0;
    num_real ohz1z2 = 0.0;
    num_real ohz1x3 = 0.0;
    num_real ohz1y3 = 0.0;
    num_real ohz1z3 = 0.0;
    num_real ohz1x4 = 0.0;
    num_real ohz1y4 = 0.0;
    num_real ohz1z4 = 0.0;
    num_real dhx2x2 = 0.0;
    num_real ohx2y2 = 0.0;
    num_real ohx2z2 = 0.0;
    num_real ohx2x3 = 0.0;
    num_real ohx2y3 = 0.0;
    num_real ohx2z3 = 0.0;
    num_real ohx2x4 = 0.0;
    num_real ohx2y4 = 0.0;
    num_real ohx2z4 = 0.0;
    num_real dhy2y2 = 0.0;
    num_real ohy2z2 = 0.0;
    num_real ohy2x3 = 0.0;
    num_real ohy2y3 = 0.0;
    num_real ohy2z3 = 0.0;
    num_real ohy2x4 = 0.0;
    num_real ohy2y4 = 0.0;
    num_real ohy2z4 = 0.0;
    num_real dhz2z2 = 0.0;
    num_real ohz2x3 = 0.0;
    num_real ohz2y3 = 0.0;
    num_real ohz2z3 = 0.0;
    num_real ohz2x4 = 0.0;
    num_real ohz2y4 = 0.0;
    num_real ohz2z4 = 0.0;
    num_real dhx3x3 = 0.0;
    num_real ohx3y3 = 0.0;
    num_real ohx3z3 = 0.0;
    num_real ohx3x4 = 0.0;
    num_real ohx3y4 = 0.0;
    num_real ohx3z4 = 0.0;
    num_real dhy3y3 = 0.0;
    num_real ohy3z3 = 0.0;
    num_real ohy3x4 = 0.0;
    num_real ohy3y4 = 0.0;
    num_real ohy3z4 = 0.0;
    num_real dhz3z3 = 0.0;
    num_real ohz3x4 = 0.0;
    num_real ohz3y4 = 0.0;
    num_real ohz3z4 = 0.0;
    num_real dhx4x4 = 0.0;
    num_real ohx4y4 = 0.0;
    num_real ohx4z4 = 0.0;
    num_real dhy4y4 = 0.0;
    num_real ohy4z4 = 0.0;
    num_real dhz4z4 = 0.0;
#endif
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
    num_real EraseLinearDihedral;
    int	I1, I2, I3, I4, IN, i;
#undef ENERGY_FUNCTION
#define ENERGY_FUNCTION _evaluateEnergyOnly_Dihedral
#undef DO_sinNPhiCosNPhi
#define DO_sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi) sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi)
    num_real sinPhase, cosPhase, SinNPhi, CosNPhi;
    gctools::Vec0<EnergyDihedral>::iterator di;
    for ( i=0,di=this->_Terms.begin(); di!=this->_Terms.end(); di++,i++ ) {
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>
#undef DIHEDRAL_APPLY_ATOM_MASK
      int index = i;
#include <cando/chem/energy_functions/_Dihedral_debugFiniteDifference.cc>
    }
  }
#endif
}



void EnergyDihedral_O::initialize()
{
  this->Base::initialize();
  this->setErrorThreshold(3.0);
}

void EnergyDihedral_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

SYMBOL_EXPORT_SC_(KeywordPkg,v);
SYMBOL_EXPORT_SC_(KeywordPkg,in);
SYMBOL_EXPORT_SC_(KeywordPkg,dn);
SYMBOL_EXPORT_SC_(KeywordPkg,phase);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,i3);
SYMBOL_EXPORT_SC_(KeywordPkg,i4);
SYMBOL_EXPORT_SC_(KeywordPkg,proper);
SYMBOL_EXPORT_SC_(KeywordPkg,atom1);
SYMBOL_EXPORT_SC_(KeywordPkg,atom2);
SYMBOL_EXPORT_SC_(KeywordPkg,atom3);
SYMBOL_EXPORT_SC_(KeywordPkg,atom4);
SYMBOL_EXPORT_SC_(KeywordPkg,coordinates);
SYMBOL_EXPORT_SC_(KeywordPkg,indices);
SYMBOL_EXPORT_SC_(KeywordPkg,simd_index);

CL_DEFMETHOD core::List_sp EnergyDihedral_O::extract_vectors_as_alist() const{
  size_t size = this->_Terms.size();
  NVector_sp v_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp in_vec = core::SimpleVector_int32_t_O::make(size);
  NVector_sp phase_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i3_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i4_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_sp proper_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom1_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom2_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom3_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom4_vec = core::SimpleVector_O::make(size);
  for (size_t i=0;i<size;++i){
    const EnergyDihedral& entry = this->_Terms[i];
    (*v_vec)[i] = entry.term.V;
    (*in_vec)[i] = entry.term.IN;
    (*phase_vec)[i] = entry._PhaseRad;
    (*i1_vec)[i] = entry.term.I1;
    (*i2_vec)[i] = entry.term.I2;
    (*i3_vec)[i] = entry.term.I3;
    (*i4_vec)[i] = entry.term.I4;
    (*proper_vec)[i] = _lisp->_boolean(entry._Proper);
    (*atom1_vec)[i] = entry._Atom1;
    (*atom2_vec)[i] = entry._Atom2;
    (*atom3_vec)[i] = entry._Atom3;
    (*atom4_vec)[i] = entry._Atom4;
  }
  ql::list result;
  result
    << core::Cons_O::create(kw::_sym_v, v_vec)
    << core::Cons_O::create(kw::_sym_in, in_vec)
    << core::Cons_O::create(kw::_sym_phase, phase_vec)
    << core::Cons_O::create(kw::_sym_i1, i1_vec)
    << core::Cons_O::create(kw::_sym_i2, i2_vec)
    << core::Cons_O::create(kw::_sym_i3, i3_vec)
    << core::Cons_O::create(kw::_sym_i4, i4_vec)
    << core::Cons_O::create(kw::_sym_proper, proper_vec)
    << core::Cons_O::create(kw::_sym_atom1, atom1_vec)
    << core::Cons_O::create(kw::_sym_atom2, atom2_vec)
    << core::Cons_O::create(kw::_sym_atom3, atom3_vec)
    << core::Cons_O::create(kw::_sym_atom4, atom4_vec);
  return result.cons();
}

CL_DEFMETHOD void EnergyDihedral_O::fill_from_vectors_in_alist(core::List_sp vectors)
{
  NVector_sp v_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_v));
  core::SimpleVector_int32_t_sp in_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_in));
  NVector_sp dn_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_dn));
  NVector_sp phase_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_phase));
  core::SimpleVector_int32_t_sp i1_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i1));
  core::SimpleVector_int32_t_sp i2_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i2));
  core::SimpleVector_int32_t_sp i3_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i3));
  core::SimpleVector_int32_t_sp i4_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i4));
  core::SimpleVector_sp proper_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_proper));
  core::SimpleVector_sp atom1_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom1));
  core::SimpleVector_sp atom2_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom2));
  core::SimpleVector_sp atom3_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom3));
  core::SimpleVector_sp atom4_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom4));
  this->_Terms.resize(v_vec->length());
  for (size_t i=0, iEnd(v_vec->length());i<iEnd;++i){
    EnergyDihedral& entry = this->_Terms[i];
    entry.term.V = (*v_vec)[i];
    entry.term.IN = (*in_vec)[i];
    entry.term.DN = (*dn_vec)[i];
    entry._PhaseRad = (*phase_vec)[i];
    entry.term.cosPhase = cos((*phase_vec)[i]);
    entry.term.sinPhase = sin((*phase_vec)[i]);
    entry.term.I1 = (*i1_vec)[i];
    entry.term.I2 = (*i2_vec)[i];
    entry.term.I3 = (*i3_vec)[i];
    entry.term.I4 = (*i4_vec)[i];
    entry._Proper =  ((*proper_vec)[i]).isTrue();
    entry._Atom1 = gc::As_unsafe<Atom_sp>((*atom1_vec)[i]);
    entry._Atom2 = gc::As_unsafe<Atom_sp>((*atom2_vec)[i]);
    entry._Atom3 = gc::As_unsafe<Atom_sp>((*atom3_vec)[i]);
    entry._Atom4 = gc::As_unsafe<Atom_sp>((*atom4_vec)[i]);
  }
}

CL_DEFMETHOD void EnergyDihedral_O::addDihedralTerm(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4,
                                                    double phase, bool proper, double v, int multiplicity)
{
  EnergyAtom* ea1 = atomTable->getEnergyAtomPointer(a1);
  EnergyAtom* ea2 = atomTable->getEnergyAtomPointer(a2);
  EnergyAtom* ea3 = atomTable->getEnergyAtomPointer(a3);
  EnergyAtom* ea4 = atomTable->getEnergyAtomPointer(a4);
  EnergyDihedral energyDihedral(a1,a2,a3,a4,phase, proper, 
                                ea1->coordinateIndexTimes3(),
                                ea2->coordinateIndexTimes3(),
                                ea3->coordinateIndexTimes3(),
                                ea4->coordinateIndexTimes3(),
                                std::sin(phase),
                                std::cos(phase),
                                v,
                                1.0*multiplicity,
                                multiplicity);
  this->addTerm(energyDihedral);
}

EnergyDihedral_sp EnergyDihedral_O::copyFilter(core::T_sp keepInteraction) {
  EnergyDihedral_sp copy = EnergyDihedral_O::create();
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    Atom_sp a3 = edi->_Atom3;
    Atom_sp a4 = edi->_Atom4;
    if ( skipInteraction( keepInteraction, EnergyDihedral_O::staticClass(), a1, a2, a3, a4 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}


#ifndef _TARGET_OS_DARWIN
#define VREAL8_WIDTH 8
typedef double real;
typedef real real8 __attribute__((vector_size(8*VREAL8_WIDTH))) __attribute__((aligned(4)));
typedef int64_t int8 __attribute__((vector_size(8*VREAL8_WIDTH))) __attribute__((aligned(4)));

num_real EnergyDihedral_O::evaluateAllComponentSimd8(
    gctools::Vec0<EnergyDihedral>::iterator di_start8,
    gctools::Vec0<EnergyDihedral>::iterator di_end8,
    ScoringFunction_sp          score,
    NVector_sp 	                pos,
    bool 		        calcForce,
    gc::Nilable<NVector_sp> 	force,
    bool		        calcDiagonalHessian,
    bool		        calcOffDiagonalHessian,
    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
    gc::Nilable<NVector_sp>	hdvec,
    gc::Nilable<NVector_sp>     dvec,
    core::T_sp activeAtomMask )
{
  
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  if ( this->_DebugEnergy ) {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(("%s {\n") , this->className());
  }
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());

  num_real totalEnergy = 0.0;
//
// Copy from implementAmberFunction::evaluateAll
//
// -----------------------
#undef mysqrt
#define mysqrt(xx) { sqrt(xx[0]), sqrt(xx[1]), sqrt(xx[2]), sqrt(xx[3]), sqrt(xx[4]), sqrt(xx[5]), sqrt(xx[6]), sqrt(xx[7]) }
#undef ZERO_SMALL_LEN
#define ZERO_SMALL_LEN(RL,L) {real8 fabs = __builtin_elementwise_abs(L); int8 cmp = (fabs < TENM3); RL = cmp ? 0.0 : RL; }
#define SIMD 1
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN

#undef VEC_CONST
#define VEC_CONST(x) VEC8(x)
  int8 bmask;
  real8 mask1 = VEC_CONST(1.0);
  real8 mask;
  int8 I1_;
  int8 I2_;
  int8 I3_;
  int8 I4_;

#undef DIHEDRAL_APPLY_ATOM_MASK
#define DIHEDRAL_APPLY_ATOM_MASK(I1,I2,I3,I4) \
  if (hasActiveAtomMask) { \
    I1_ = I1/3; \
    I2_ = I2/3; \
    I3_ = I3/3; \
    I4_ = I4/3; \
    bmask = { \
        (bitvectorActiveAtomMask->testBit(I1_[0]) \
         && bitvectorActiveAtomMask->testBit(I2_[0]) \
         && bitvectorActiveAtomMask->testBit(I3_[0]) \
         && bitvectorActiveAtomMask->testBit(I4_[0])), \
        (bitvectorActiveAtomMask->testBit(I1_[1]) \
         && bitvectorActiveAtomMask->testBit(I2_[1]) \
         && bitvectorActiveAtomMask->testBit(I3_[1]) \
         && bitvectorActiveAtomMask->testBit(I4_[1])), \
        (bitvectorActiveAtomMask->testBit(I1_[2]) \
         && bitvectorActiveAtomMask->testBit(I2_[2]) \
         && bitvectorActiveAtomMask->testBit(I3_[2]) \
         && bitvectorActiveAtomMask->testBit(I4_[2])),\
        (bitvectorActiveAtomMask->testBit(I1_[3]) \
         && bitvectorActiveAtomMask->testBit(I2_[3]) \
         && bitvectorActiveAtomMask->testBit(I3_[3]) \
         && bitvectorActiveAtomMask->testBit(I4_[3])), \
        (bitvectorActiveAtomMask->testBit(I1_[4]) \
         && bitvectorActiveAtomMask->testBit(I2_[4]) \
         && bitvectorActiveAtomMask->testBit(I3_[4]) \
         && bitvectorActiveAtomMask->testBit(I4_[4])), \
        (bitvectorActiveAtomMask->testBit(I1_[5]) \
         && bitvectorActiveAtomMask->testBit(I2_[5]) \
         && bitvectorActiveAtomMask->testBit(I3_[5]) \
         && bitvectorActiveAtomMask->testBit(I4_[5])), \
        (bitvectorActiveAtomMask->testBit(I1_[6]) \
         && bitvectorActiveAtomMask->testBit(I2_[6]) \
         && bitvectorActiveAtomMask->testBit(I3_[6]) \
         && bitvectorActiveAtomMask->testBit(I4_[6])), \
        (bitvectorActiveAtomMask->testBit(I1_[7]) \
         && bitvectorActiveAtomMask->testBit(I2_[7]) \
         && bitvectorActiveAtomMask->testBit(I3_[7]) \
         && bitvectorActiveAtomMask->testBit(I4_[7])) \
    };\
    if (!(bmask[0]|| bmask[1]|| bmask[2]|| bmask[3]|| bmask[4]|| bmask[5]|| bmask[6]|| bmask[7])) goto SKIP_term_past_linear_dihedral_test; \
    mask = (bmask ? 1.0 : 0.0);  \
    /*printf("%s:%d:%s mask = %f,%f,%f,%f,%f,%f,%f,%f\n", __FILE__, __LINE__, __FUNCTION__, mask[0], mask[1], mask[2], mask[3], mask[4], mask[5], mask[6], mask[7] ); */\
  } \
    else { \
      mask = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }; \
    }



#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x={ di[0].term.x, di[1].term.x, di[2].term.x, di[3].term.x, di[4].term.x, di[5].term.x, di[6].term.x, di[7].term.x };}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x={ pos->element(ii[0]+of), pos->element(ii[1]+of), pos->element(ii[2]+of), pos->element(ii[3]+of), pos->element(ii[4]+of), pos->element(ii[5]+of), pos->element(ii[6]+of), pos->element(ii[7]+of) }; }
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) totalEnergy += (e[0]*mask[0]+e[1]*mask[1]+e[2]*mask[2]+e[3]*mask[3]+e[4]*mask[4]+e[5]*mask[5]+e[6]*mask[6]+e[7]*mask[7]);
#undef	DIHEDRAL_FORCE_ACCUMULATE
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#undef ONE_ForceAcc
#define ONE_ForceAcc(IDX,i,o,v) force->setElement((i)+(o),(v)*mask[IDX]+force->getElement((i)+(o)));
#undef DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(I,O,V) { \
    if (hasForce) { \
      ONE_ForceAcc(0,I[0],O,V[0]); \
      ONE_ForceAcc(1,I[1],O,V[1]); \
      ONE_ForceAcc(2,I[2],O,V[2]); \
      ONE_ForceAcc(3,I[3],O,V[3]); \
      ONE_ForceAcc(4,I[4],O,V[4]); \
      ONE_ForceAcc(5,I[5],O,V[5]); \
      ONE_ForceAcc(6,I[6],O,V[6]); \
      ONE_ForceAcc(7,I[7],O,V[7]); \
    }}
#undef ONE_OffDiagHessAcc_hasHessian
#define ONE_OffDiagHessAcc_hasHessian(vi,i1,o1,i2,o2,v) hessian->addToElement((i1[vi])+(o1),(i2[vi])+(o2),v[vi]*mask[vi]);
#undef ONE_OffDiagHessAcc_hasHdAndD
#define ONE_OffDiagHessAcc_hasHdAndD(vi,i1,o1,i2,o2,v) {\
    hdvec->addToElement((i1[vi])+(o1),v[vi]*mask[vi]*dvec->element((i2[vi])+(o2)));\
    hdvec->addToElement((i2[vi])+(o2),v[vi]*mask[vi]*dvec->element((i1[vi])+(o1)));\
  }
#undef ONE_DiagHessAcc_hasHessian
#define ONE_DiagHessAcc_hasHessian(vi,i1,o1,i2,o2,v)  hessian->addToElement((i1[vi])+(o1),(i2[vi])+(o2),v[vi]*mask[vi]);
#undef ONE_DiagHessAcc_hasHdAndD
#define ONE_DiagHessAcc_hasHdAndD(vi,i1,o1,i2,o2,v)   hdvec->addToElement((i1[vi])+(o1),v[vi]*mask[vi]*dvec->element((i1[vi])+(o1)));

#undef DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) { \
    if (hasHessian) { \
      ONE_DiagHessAcc_hasHessian(0,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(1,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(2,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(3,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(4,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(5,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(6,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(7,i1,o1,i2,o2,v); \
    } \
    if (hasHdAndD) { \
      ONE_DiagHessAcc_hasHdAndD(0,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(1,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(2,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(3,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(4,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(5,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(6,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(7,i1,o1,i2,o2,v); \
    } \
  }
#undef DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) { \
    if (hasHessian) { \
      ONE_OffDiagHessAcc_hasHessian(0,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(1,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(2,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(3,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(4,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(5,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(6,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(7,i1,o1,i2,o2,v); \
    } \
    if (hasHdAndD) { \
      ONE_OffDiagHessAcc_hasHdAndD(0,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(1,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(2,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(3,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(4,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(5,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(6,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(7,i1,o1,i2,o2,v); \
    } \
  }
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#undef DECLARE_FLOAT
#define DECLARE_FLOAT(xx) real8 xx
#include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
#pragma clang diagnostic pop
#undef VEC_CONST
#define VEC_CONST(x) VEC8(x)
  
  fx1 = VEC_CONST(0.0);
  fy1 = VEC_CONST(0.0);
  fz1 = VEC_CONST(0.0);
  fx2 = VEC_CONST(0.0);
  fy2 = VEC_CONST(0.0);
  fz2 = VEC_CONST(0.0);
  fx3 = VEC_CONST(0.0);
  fy3 = VEC_CONST(0.0);
  fz3 = VEC_CONST(0.0);
  fx4 = VEC_CONST(0.0);
  fy4 = VEC_CONST(0.0);
  fz4 = VEC_CONST(0.0);
  real8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
  real8 EraseLinearDihedral;
  int8	I1, I2, I3, I4, IN;
  real8 sinPhase, cosPhase, SinNPhi, CosNPhi;

#undef ONE_sinNPhiCosNPhi
#define ONE_sinNPhiCosNPhi(INI,IN,SinNPhi,CosNPhi,SinPhi,CosPhi) {\
    real tSinNPhi, tCosNPhi; \
    sinNPhiCosNPhi(IN[INI],tSinNPhi,tCosNPhi,SinPhi[INI],CosPhi[INI]); \
    SinNPhi[INI] = tSinNPhi; CosNPhi[INI] = tCosNPhi; \
  }

#undef DO_sinNPhiCosNPhi
#define DO_sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi) { \
    ONE_sinNPhiCosNPhi(0,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(1,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(2,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(3,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(4,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(5,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(6,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(7,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
  }

    
  gctools::Vec0<EnergyDihedral>::iterator di;
  int i=0;
  for ( di=di_start8; di<di_end8; di += VREAL8_WIDTH ) {
#define DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) IGNORE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4)
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>
#undef DIHEDRAL_DEBUG_INTERACTIONS
#undef DIHEDRAL_APPLY_ATOM_MASK
#undef SIMD
#undef MAYBE_ERROR_LINEAR_DIHEDRAL
#define MAYBE_ERROR_LINEAR_DIHEDRAL(IDX,di,I1,I2,I3,I4) { \
      if ( EraseLinearDihedral[IDX] == 0.0 && bmask[IDX] ) { \
          ERROR(chem::_sym_LinearDihedralError, \
                core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(di[IDX]._Atom1,\
                                                                                 di[IDX]._Atom2,\
                                                                                 di[IDX]._Atom3,\
                                                                                 di[IDX]._Atom4 ), \
                                         kw::_sym_simd_index, core::make_fixnum(IDX), \
                                         kw::_sym_coordinates,pos, \
                                         kw::_sym_indices,core::Cons_O::createList(core::make_fixnum(I1[IDX]), core::make_fixnum(I2[IDX]), core::make_fixnum(I3[IDX]), core::make_fixnum(I4[IDX])))); \
      }}

    MAYBE_ERROR_LINEAR_DIHEDRAL(0,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(1,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(2,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(3,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(4,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(5,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(6,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(7,di,I1,I2,I3,I4);

    // If no dihedrals were calculated then we skip all linear dihedral tests
  SKIP_term_past_linear_dihedral_test: (void)0;
    
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
    di->_calcForce = calcForce;
    di->_calcDiagonalHessian = calcDiagonalHessian;
    di->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ di->eval.var=val;};
#include	<cando/chem/energy_functions/_Dihedral_debugEvalSet.cc>
#endif //]
    i += VREAL8_WIDTH;
  }
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY(("%s }\n") , this->className());
  }
  return totalEnergy;
}

//
//
//  -------------------------------------------------------------
//
//  Simd4 function
//
//

#define VREAL4_WIDTH 4
typedef double real;
typedef real real4 __attribute__((vector_size(8*VREAL4_WIDTH)));
typedef int64_t int4 __attribute__((vector_size(8*VREAL4_WIDTH)));

num_real EnergyDihedral_O::evaluateAllComponentSimd4(
    gctools::Vec0<EnergyDihedral>::iterator di_start4,
    gctools::Vec0<EnergyDihedral>::iterator di_end4,
    ScoringFunction_sp          score,
    NVector_sp 	                pos,
    bool 		        calcForce,
    gc::Nilable<NVector_sp> 	force,
    bool		        calcDiagonalHessian,
    bool		        calcOffDiagonalHessian,
    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
    gc::Nilable<NVector_sp>	hdvec,
    gc::Nilable<NVector_sp>     dvec,
    core::T_sp activeAtomMask )
{
  IMPLEMENT_ME();
#if 0
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  if ( this->_DebugEnergy ) {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(("%s {\n") , this->className());
  }
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());

  double totalEnergy = 0.0;
//
// Copy from implementAmberFunction::evaluateAll
//
// -----------------------
#undef mysqrt
#define mysqrt(xx) { sqrt(xx[0]), sqrt(xx[1]), sqrt(xx[2]), sqrt(xx[3]) }
#undef ZERO_SMALL_LEN
#define ZERO_SMALL_LEN(RL,L) {real4 fabs = __builtin_elementwise_abs(L); int4 cmp = (fabs < TENM3); RL = cmp ? 0.0 : RL; }
#define SIMD 1
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x={ di[0].term.x, di[1].term.x, di[2].term.x, di[3].term.x};}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x={ pos->element(ii[0]+of), pos->element(ii[1]+of), pos->element(ii[2]+of), pos->element(ii[3]+of)};}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) totalEnergy += (e[0]+e[1]+e[2]+e[3]);
#undef	DIHEDRAL_FORCE_ACCUMULATE
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#undef ONE_ForceAcc
#define ONE_ForceAcc(i,o,v) force->setElement((i)+(o),(v)+force->getElement((i)+(o)));
#undef DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(I,O,V) { \
    if (hasForce) { \
      ONE_ForceAcc(I[0],O,V[0]); \
      ONE_ForceAcc(I[1],O,V[1]); \
      ONE_ForceAcc(I[2],O,V[2]); \
      ONE_ForceAcc(I[3],O,V[3]); \
    }}
#undef ONE_OffDiagHessAcc_hasHessian
#define ONE_OffDiagHessAcc_hasHessian(vi,i1,o1,i2,o2,v) hessian->addToElement((i1[vi])+(o1),(i2[vi])+(o2),v[vi]);
#undef ONE_OffDiagHessAcc_hasHdAndD
#define ONE_OffDiagHessAcc_hasHdAndD(vi,i1,o1,i2,o2,v) {\
    hdvec->addToElement((i1[vi])+(o1),v[vi]*dvec->element((i2[vi])+(o2)));\
    hdvec->addToElement((i2[vi])+(o2),v[vi]*dvec->element((i1[vi])+(o1)));\
  }
#undef ONE_DiagHessAcc_hasHessian
#define ONE_DiagHessAcc_hasHessian(vi,i1,o1,i2,o2,v)  hessian->addToElement((i1[vi])+(o1),(i2[vi])+(o2),v[vi]);
#undef ONE_DiagHessAcc_hasHdAndD
#define ONE_DiagHessAcc_hasHdAndD(vi,i1,o1,i2,o2,v)   hdvec->addToElement((i1[vi])+(o1),v[vi]*dvec->element((i1[vi])+(o1)));

#undef DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) { \
    if (hasHessian) { \
      ONE_DiagHessAcc_hasHessian(0,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(1,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(2,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(3,i1,o1,i2,o2,v); \
    } \
    if (hasHdAndD) { \
      ONE_DiagHessAcc_hasHdAndD(0,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(1,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(2,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(3,i1,o1,i2,o2,v); \
    } \
  }
#undef DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) { \
    if (hasHessian) { \
      ONE_OffDiagHessAcc_hasHessian(0,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(1,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(2,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(3,i1,o1,i2,o2,v); \
    } \
    if (hasHdAndD) { \
      ONE_OffDiagHessAcc_hasHdAndD(0,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(1,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(2,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(3,i1,o1,i2,o2,v); \
    } \
  }
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#undef DECLARE_FLOAT
#define DECLARE_FLOAT(xx) real4 xx
#include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
#pragma clang diagnostic pop
#undef VEC_CONST
#define VEC_CONST(x) VEC4(x)
  
  fx1 = VEC_CONST(0.0);
  fy1 = VEC_CONST(0.0);
  fz1 = VEC_CONST(0.0);
  fx2 = VEC_CONST(0.0);
  fy2 = VEC_CONST(0.0);
  fz2 = VEC_CONST(0.0);
  fx3 = VEC_CONST(0.0);
  fy3 = VEC_CONST(0.0);
  fz3 = VEC_CONST(0.0);
  fx4 = VEC_CONST(0.0);
  fy4 = VEC_CONST(0.0);
  fz4 = VEC_CONST(0.0);
  real4 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
  real4 EraseLinearDihedral;
  int4	I1, I2, I3, I4, IN;
  real4 sinPhase, cosPhase, SinNPhi, CosNPhi;

#undef ONE_sinNPhiCosNPhi
#define ONE_sinNPhiCosNPhi(INI,IN,SinNPhi,CosNPhi,SinPhi,CosPhi) {\
    real tSinNPhi, tCosNPhi; \
    sinNPhiCosNPhi(IN[INI],tSinNPhi,tCosNPhi,SinPhi[INI],CosPhi[INI]); \
    SinNPhi[INI] = tSinNPhi; CosNPhi[INI] = tCosNPhi; \
  }

#undef DO_sinNPhiCosNPhi
#define DO_sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi) { \
    ONE_sinNPhiCosNPhi(0,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(1,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(2,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(3,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
  }

  gctools::Vec0<EnergyDihedral>::iterator di;
  int i = 0;
  for ( di=di_start4; di<di_end4; di += VREAL4_WIDTH ) {
#define DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) IGNORE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4)
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>
#undef DIHEDRAL_DEBUG_INTERACTIONS
#undef DIHEDRAL_APPLY_ATOM_MASK
#undef SIMD
#undef MAYBE_ERROR_LINEAR_DIHEDRAL
#define MAYBE_ERROR_LINEAR_DIHEDRAL(IDX,di,I1,I2,I3,I4) { \
      if ( EraseLinearDihedral[IDX] == 0.0 ) { \
          ERROR(chem::_sym_LinearDihedralError, \
                core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(di[IDX]._Atom1,\
                                                                                 di[IDX]._Atom2,\
                                                                                 di[IDX]._Atom3,\
                                                                                 di[IDX]._Atom4 ), \
                                         kw::_sym_coordinates,pos, \
                                         kw::_sym_indices,core::Cons_O::createList(core::make_fixnum(I1[IDX]), core::make_fixnum(I2[IDX]), core::make_fixnum(I3[IDX]), core::make_fixnum(I4[IDX])))); \
      }}

    MAYBE_ERROR_LINEAR_DIHEDRAL(0,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(1,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(2,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(3,di,I1,I2,I3,I4);

#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
    di->_calcForce = calcForce;
    di->_calcDiagonalHessian = calcDiagonalHessian;
    di->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ di->eval.var=val;};
#include	<cando/chem/energy_functions/_Dihedral_debugEvalSet.cc>
#endif //]
    i += VREAL4_WIDTH;
  }
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY(("%s }\n") , this->className());
  }
  return totalEnergy;
#endif
}

//
//
//  -------------------------------------------------------------
//
//  Simd2 function
//
//

#define VREAL2_WIDTH 2
typedef double real;
typedef real real2 __attribute__((vector_size(8*VREAL2_WIDTH)));
typedef int64_t int2 __attribute__((vector_size(8*VREAL2_WIDTH)));

num_real EnergyDihedral_O::evaluateAllComponentSimd2(
    gctools::Vec0<EnergyDihedral>::iterator di_start2,
    gctools::Vec0<EnergyDihedral>::iterator di_end2,
    ScoringFunction_sp          score,
    NVector_sp 	                pos,
    bool 		        calcForce,
    gc::Nilable<NVector_sp> 	force,
    bool		        calcDiagonalHessian,
    bool		        calcOffDiagonalHessian,
    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
    gc::Nilable<NVector_sp>	hdvec,
    gc::Nilable<NVector_sp>     dvec,
    core::T_sp activeAtomMask )
{
  IMPLEMENT_ME();
#if 0
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  if ( this->_DebugEnergy ) {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(("%s {\n") , this->className());
  }
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());

  double totalEnergy = 0.0;
//
// Copy from implementAmberFunction::evaluateAll
//
// -----------------------
#undef mysqrt
#define mysqrt(xx) { sqrt(xx[0]), sqrt(xx[1]) }
#undef ZERO_SMALL_LEN
#define ZERO_SMALL_LEN(RL,L) {real2 fabs = __builtin_elementwise_abs(L); int2 cmp = (fabs < TENM3); RL = cmp ? 0.0 : RL; }
#define SIMD 1
#define DIHEDRAL_CALC_FORCE
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
#undef	DIHEDRAL_SET_PARAMETER
#define	DIHEDRAL_SET_PARAMETER(x)	{x={ di[0].term.x, di[1].term.x };}
#undef	DIHEDRAL_SET_POSITION
#define	DIHEDRAL_SET_POSITION(x,ii,of)	{x={ pos->element(ii[0]+of), pos->element(ii[1]+of)};}
#undef	DIHEDRAL_ENERGY_ACCUMULATE
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) totalEnergy += (e[0]+e[1])
#undef	DIHEDRAL_FORCE_ACCUMULATE
#undef	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#undef	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#undef ONE_ForceAcc
#define ONE_ForceAcc(i,o,v) force->setElement((i)+(o),(v)+force->getElement((i)+(o)));
#undef DIHEDRAL_FORCE_ACCUMULATE
#define	DIHEDRAL_FORCE_ACCUMULATE(I,O,V) { \
    if (hasForce) { \
      ONE_ForceAcc(I[0],O,V[0]); \
      ONE_ForceAcc(I[1],O,V[1]); \
    }}
#undef ONE_OffDiagHessAcc_hasHessian
#define ONE_OffDiagHessAcc_hasHessian(vi,i1,o1,i2,o2,v) hessian->addToElement((i1[vi])+(o1),(i2[vi])+(o2),v[vi]);
#undef ONE_OffDiagHessAcc_hasHdAndD
#define ONE_OffDiagHessAcc_hasHdAndD(vi,i1,o1,i2,o2,v) {\
    hdvec->addToElement((i1[vi])+(o1),v[vi]*dvec->element((i2[vi])+(o2)));\
    hdvec->addToElement((i2[vi])+(o2),v[vi]*dvec->element((i1[vi])+(o1)));\
  }
#undef ONE_DiagHessAcc_hasHessian
#define ONE_DiagHessAcc_hasHessian(vi,i1,o1,i2,o2,v)  hessian->addToElement((i1[vi])+(o1),(i2[vi])+(o2),v[vi]);
#undef ONE_DiagHessAcc_hasHdAndD
#define ONE_DiagHessAcc_hasHdAndD(vi,i1,o1,i2,o2,v)   hdvec->addToElement((i1[vi])+(o1),v[vi]*dvec->element((i1[vi])+(o1)));

#undef DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) { \
    if (hasHessian) { \
      ONE_DiagHessAcc_hasHessian(0,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHessian(1,i1,o1,i2,o2,v); \
    } \
    if (hasHdAndD) { \
      ONE_DiagHessAcc_hasHdAndD(0,i1,o1,i2,o2,v); \
      ONE_DiagHessAcc_hasHdAndD(1,i1,o1,i2,o2,v); \
    } \
  }
#undef DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) { \
    if (hasHessian) { \
      ONE_OffDiagHessAcc_hasHessian(0,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHessian(1,i1,o1,i2,o2,v); \
    } \
    if (hasHdAndD) { \
      ONE_OffDiagHessAcc_hasHdAndD(0,i1,o1,i2,o2,v); \
      ONE_OffDiagHessAcc_hasHdAndD(1,i1,o1,i2,o2,v); \
    } \
  }
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#undef DECLARE_FLOAT
#define DECLARE_FLOAT(xx) real2 xx
#include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
#pragma clang diagnostic pop
#undef VEC_CONST
#define VEC_CONST(x) VEC2(x)
  
  fx1 = VEC_CONST(0.0);
  fy1 = VEC_CONST(0.0);
  fz1 = VEC_CONST(0.0);
  fx2 = VEC_CONST(0.0);
  fy2 = VEC_CONST(0.0);
  fz2 = VEC_CONST(0.0);
  fx3 = VEC_CONST(0.0);
  fy3 = VEC_CONST(0.0);
  fz3 = VEC_CONST(0.0);
  fx4 = VEC_CONST(0.0);
  fy4 = VEC_CONST(0.0);
  fz4 = VEC_CONST(0.0);
  real2 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
  real2 EraseLinearDihedral;
  int2	I1, I2, I3, I4, IN;
  real2 sinPhase, cosPhase, SinNPhi, CosNPhi;

#undef ONE_sinNPhiCosNPhi
#define ONE_sinNPhiCosNPhi(INI,IN,SinNPhi,CosNPhi,SinPhi,CosPhi) {\
    real tSinNPhi, tCosNPhi; \
    sinNPhiCosNPhi(IN[INI],tSinNPhi,tCosNPhi,SinPhi[INI],CosPhi[INI]); \
    SinNPhi[INI] = tSinNPhi; CosNPhi[INI] = tCosNPhi; \
  }

#undef DO_sinNPhiCosNPhi
#define DO_sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi) { \
    ONE_sinNPhiCosNPhi(0,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
    ONE_sinNPhiCosNPhi(1,IN,SinNPhi,CosNPhi,SinPhi,CosPhi); \
  }

  gctools::Vec0<EnergyDihedral>::iterator di;
  int i = 0;
  for ( di=di_start2; di<di_end2; di += VREAL2_WIDTH ) {
#define DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) IGNORE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4)
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>
#undef DIHEDRAL_DEBUG_INTERACTIONS
#undef DIHEDRAL_APPLY_ATOM_MASK
#undef SIMD
#undef MAYBE_ERROR_LINEAR_DIHEDRAL
#define MAYBE_ERROR_LINEAR_DIHEDRAL(IDX,di,I1,I2,I3,I4) { \
      if ( EraseLinearDihedral[IDX] == 0.0 ) { \
          ERROR(chem::_sym_LinearDihedralError, \
                core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(di[IDX]._Atom1,\
                                                                                 di[IDX]._Atom2,\
                                                                                 di[IDX]._Atom3,\
                                                                                 di[IDX]._Atom4 ), \
                                         kw::_sym_coordinates,pos, \
                                         kw::_sym_indices,core::Cons_O::createList(core::make_fixnum(I1[IDX]), core::make_fixnum(I2[IDX]), core::make_fixnum(I3[IDX]), core::make_fixnum(I4[IDX])))); \
      }}

    MAYBE_ERROR_LINEAR_DIHEDRAL(0,di,I1,I2,I3,I4);
    MAYBE_ERROR_LINEAR_DIHEDRAL(1,di,I1,I2,I3,I4);

#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
    di->_calcForce = calcForce;
    di->_calcDiagonalHessian = calcDiagonalHessian;
    di->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ di->eval.var=val;};
#include	<cando/chem/energy_functions/_Dihedral_debugEvalSet.cc>
#endif //]
    i += VREAL2_WIDTH;
  }
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY(("%s }\n") , this->className());
  }
  return totalEnergy;
#endif
}
#endif // !_TARGET_OS_DARWIN

num_real EnergyDihedral_O::evaluateAllComponent(ScoringFunction_sp          score,
                                                NVector_sp 	                pos,
                                                core::T_sp               componentEnergy,
                                                bool 		        calcForce,
                                                gc::Nilable<NVector_sp> 	force,
                                                bool		        calcDiagonalHessian,
                                                bool		        calcOffDiagonalHessian,
                                                gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                gc::Nilable<NVector_sp>	hdvec,
                                                gc::Nilable<NVector_sp>     dvec,
                                                core::T_sp activeAtomMask,
                                                core::T_sp debugInteractions )
{
  num_real  energy = 0.0;
  this->_Evaluations++;
#ifdef _TARGET_OS_DARWIN
    energy += this->evaluateAllComponentSingle(this->_Terms.begin(),
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               debugInteractions );
#else
  if (cando::global_simd_width ==1 || debugInteractions.notnilp()) {
    energy += this->evaluateAllComponentSingle(this->_Terms.begin(),
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               debugInteractions );
  } else if (cando::global_simd_width == 8 ) {
    gctools::Vec0<EnergyDihedral>::iterator di_vector_end8 = this->_Terms.begin()+((int)(std::distance(this->_Terms.begin(),this->_Terms.end())/VREAL8_WIDTH))*VREAL8_WIDTH;
    energy += this->evaluateAllComponentSimd8(this->_Terms.begin(),
                                              di_vector_end8,
                                              score,
                                              pos,
                                              calcForce,
                                              force,
                                              calcDiagonalHessian,
                                              calcOffDiagonalHessian,
                                              hessian,
                                              hdvec,
                                              dvec,
                                              activeAtomMask );
    energy += this->evaluateAllComponentSingle(di_vector_end8,
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               nil<core::T_O>());
  } else if (cando::global_simd_width == 4 ) {
    gctools::Vec0<EnergyDihedral>::iterator di_vector_end4 = this->_Terms.begin()+((int)(std::distance(this->_Terms.begin(),this->_Terms.end())/VREAL4_WIDTH))*VREAL4_WIDTH;
    energy += this->evaluateAllComponentSimd4(this->_Terms.begin(),
                                              di_vector_end4,
                                              score,
                                              pos,
                                              calcForce,
                                              force,
                                              calcDiagonalHessian,
                                              calcOffDiagonalHessian,
                                              hessian,
                                              hdvec,
                                              dvec,
                                              activeAtomMask );
    energy += this->evaluateAllComponentSingle(di_vector_end4,
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               nil<core::T_O>());
  } else if (cando::global_simd_width == 2 ) {
    gctools::Vec0<EnergyDihedral>::iterator di_vector_end2 = this->_Terms.begin()+((int)(std::distance(this->_Terms.begin(),this->_Terms.end())/VREAL2_WIDTH))*VREAL2_WIDTH;
    energy += this->evaluateAllComponentSimd2(this->_Terms.begin(),
                                              di_vector_end2,
                                              score,
                                              pos,
                                              calcForce,
                                              force,
                                              calcDiagonalHessian,
                                              calcOffDiagonalHessian,
                                              hessian,
                                              hdvec,
                                              dvec,
                                              activeAtomMask );
    energy += this->evaluateAllComponentSingle(di_vector_end2,
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               nil<core::T_O>() );
  }
#endif
  maybeSetEnergy( componentEnergy, EnergyDihedral_O::static_classSymbol(), energy );
  return energy;
};

};
