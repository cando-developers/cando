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
#include <clasp/core/wrappers.h>

namespace chem {


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
  SIMPLE_ERROR(("Implement decode of EnergyDihedral"));
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


void	sinNPhiCosNPhi(int n, double* sinNPhi, double* cosNPhi,
                       double sinPhi, double cosPhi )
{
  double	sinNm1Phi, cosNm1Phi;
  if ( n==1 ) {
    *sinNPhi = sinPhi;
    *cosNPhi = cosPhi;
    return;
  };
  sinNPhiCosNPhi(n-1,&sinNm1Phi,&cosNm1Phi,sinPhi,cosPhi);
  *sinNPhi = cosPhi*sinNm1Phi+sinPhi*cosNm1Phi;
  *cosNPhi = cosPhi*cosNm1Phi-sinPhi*sinNm1Phi;
  return;
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
  this->_PhaseRad = ffterm->getPhaseRad(n);
  this->term.sinPhase = sin(this->_PhaseRad);
  this->term.cosPhase = cos(this->_PhaseRad);
  this->term.V = ffterm->getV_kCal(n)*scale;
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
  double	angle;
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
  double	phi, dev;
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
  this->term.V = term->getV_kCal(n)*scale;
  this->_PhaseRad = term->getPhaseRad(n);
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
//    this->_Type1 = ea1->_Atom->getType();
//    this->_Type2 = ea2->_Atom->getType();
  /    this->_Type3 = ea3->_Atom->getType();
  this->_Type4 = ea4->_Atom->getType();
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
double	_evaluateEnergyOnly_Dihedral(
                                     double x1, double y1, double z1,
                                     double x2, double y2, double z2,
                                     double x3, double y3, double z3,
                                     double x4, double y4, double z4,
                                     double V, double DN, int IN,
                                     double cosPhase, double sinPhase )
{
  double	EraseLinearDihedral;
  double	SinNPhi, CosNPhi;

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
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>

  return Energy;
}


void EnergyDihedral_O::addTerm(const EnergyDihedral& term)
{
  this->_Terms.push_back(term);
}


string EnergyDihedral_O::beyondThresholdInteractionsAsString()
{
  return component_beyondThresholdInteractionsAsString<EnergyDihedral_O,EnergyDihedral>(*this);
}



void	EnergyDihedral_O::dumpTerms()
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
    core::writeln_bf_stream(fmt::sprintf("TERM 3DIH %s %-9s - %-9s - %-9s - %-9s %8.2lf %8.2lf %2.0lf" , type , str1 , str2 , str3 , str4 , edi->term.V , edi->_PhaseRad , edi->term.DN));
  }
}


CL_DEFMETHOD
core::List_sp	EnergyDihedral_O::lookupDihedralTerms(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4)
{
  gctools::Vec0<EnergyDihedral>::iterator	edi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4, type;
  ql::list  result;
  core::T_sp tia1 = atomTable->_AtomTableIndices->gethash(a1);
  core::T_sp tia2 = atomTable->_AtomTableIndices->gethash(a2);
  core::T_sp tia3 = atomTable->_AtomTableIndices->gethash(a3);
  core::T_sp tia4 = atomTable->_AtomTableIndices->gethash(a4);
  if (!tia1.fixnump()) SIMPLE_ERROR(("Could not find %s in energy function") , _rep_(a1));
  if (!tia2.fixnump()) SIMPLE_ERROR(("Could not find %s in energy function") , _rep_(a2));
  if (!tia3.fixnump()) SIMPLE_ERROR(("Could not find %s in energy function") , _rep_(a3));
  if (!tia4.fixnump()) SIMPLE_ERROR(("Could not find %s in energy function") , _rep_(a4));
  int ia1 = tia1.unsafe_fixnum();
  int ia2 = tia2.unsafe_fixnum();
  int ia3 = tia3.unsafe_fixnum();
  int ia4 = tia4.unsafe_fixnum();
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
      oneResult << INTERN_(kw,type1) << edi->_Atom1->getType()
                << INTERN_(kw,type2) << edi->_Atom2->getType()
                << INTERN_(kw,type3) << edi->_Atom3->getType()
                << INTERN_(kw,type4) << edi->_Atom4->getType()
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
                                                     AbstractLargeSquareMatrix_sp m )
{
  bool		calcForce = true;
  bool		calcDiagonalHessian = true;
  bool		calcOffDiagonalHessian = true;

//
// Copy from implementAmberFunction::setupHessianPreconditioner
//
// -----------------------

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
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
    double EraseLinearDihedral;
    int	I1, I2, I3, I4, IN;
    double sinPhase, cosPhase, SinNPhi, CosNPhi;
    for ( gctools::Vec0<EnergyDihedral>::iterator di=this->_Terms.begin();
          di!=this->_Terms.end(); di++ ) {
#include	<cando/chem/energy_functions/_Dihedral_termCode.cc>
    }
  }


}





double	EnergyDihedral_O::evaluateAllComponent( ScoringFunction_sp score,
                                       NVector_sp 	pos,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec,
                                       gc::Nilable<NVector_sp> dvec)
{ 
  this->_Evaluations++;
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(BF("%s {\n")% this->className());
  }

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
#define	DIHEDRAL_ENERGY_ACCUMULATE(e) this->_TotalEnergy += (e);
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
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
    double EraseLinearDihedral;
    int	I1, I2, I3, I4, IN, i;
    double sinPhase, cosPhase, SinNPhi, CosNPhi;
    gctools::Vec0<EnergyDihedral>::iterator di;
    for ( i=0,di=this->_Terms.begin();
          di!=this->_Terms.end(); di++,i++ ) {
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfDihedralTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfDihedralTermsToCalculate ) {
          break;
        }
      }
#endif
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>
      if ( EraseLinearDihedral == 0.0 ) {
        ERROR(chem::_sym_LinearDihedralError,core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(di->_Atom1,di->_Atom2,di->_Atom3,di->_Atom4),
                                                                      kw::_sym_coordinates,pos,
                                                                      kw::_sym_indices,core::Cons_O::createList(core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3), core::make_fixnum(I4))));
      }
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
        LOG_ENERGY(BF( "MEISTER dihedral %d args cando\n")% (i+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d V %lf\n")% (i+1) % V );
        LOG_ENERGY(BF( "MEISTER dihedral %d DN %lf\n")% (i+1) % DN );
        LOG_ENERGY(BF( "MEISTER dihedral %d IN %d\n")% (i+1) % IN );
//		    LOG_ENERGY(BF( "MEISTER dihedral %d phase %lf\n")% (i+1) % phase );
        LOG_ENERGY(BF( "MEISTER dihedral %d sinPhase %lf\n")% (i+1) % sinPhase );
        LOG_ENERGY(BF( "MEISTER dihedral %d cosPhase %lf\n")% (i+1) % cosPhase );
        LOG_ENERGY(BF( "MEISTER dihedral %d x1 %5.3lf %d\n")%(i+1) % x1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d y1 %5.3lf %d\n")%(i+1) % y1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d z1 %5.3lf %d\n")%(i+1) % z1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d x2 %5.3lf %d\n")%(i+1) % x2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d y2 %5.3lf %d\n")%(i+1) % y2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d z2 %5.3lf %d\n")%(i+1) % z2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d x3 %5.3lf %d\n")%(i+1) % x3 % (I3/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d y3 %5.3lf %d\n")%(i+1) % y3 % (I3/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d z3 %5.3lf %d\n")%(i+1) % z3 % (I3/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d x4 %5.3lf %d\n")%(i+1) % x4 % (I4/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d y4 %5.3lf %d\n")%(i+1) % y4 % (I4/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d z4 %5.3lf %d\n")%(i+1) % z4 % (I4/3+1) );
        LOG_ENERGY(BF( "MEISTER dihedral %d results\n")% (i+1) );
	//	LOG_ENERGY(BF( "MEISTER dihedral %d Phi %lf (%lf degrees)\n")% (i+1) %
	//		    	Phi, Phi/0.0174533 );
        LOG_ENERGY(BF( "MEISTER dihedral %d Energy %lf\n")% (i+1) % Energy);
        if ( calcForce ) {
          LOG_ENERGY(BF( "MEISTER dihedral %d fx1 %8.5lf %d\n")%(i+1) % fx1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fy1 %8.5lf %d\n")%(i+1) % fy1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fz1 %8.5lf %d\n")%(i+1) % fz1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fx2 %8.5lf %d\n")%(i+1) % fx2 % (I2/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fy2 %8.5lf %d\n")%(i+1) % fy2 % (I2/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fz2 %8.5lf %d\n")%(i+1) % fz2 % (I2/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fx3 %8.5lf %d\n")%(i+1) % fx3 % (I3/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fy3 %8.5lf %d\n")%(i+1) % fy3 % (I3/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fz3 %8.5lf %d\n")%(i+1) % fz3 % (I3/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fx4 %8.5lf %d\n")%(i+1) % fx4 % (I4/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fy4 %8.5lf %d\n")%(i+1) % fy4 % (I4/3+1) );
          LOG_ENERGY(BF( "MEISTER dihedral %d fz4 %8.5lf %d\n")%(i+1) % fz4 % (I4/3+1) );
        }
        LOG_ENERGY(BF( "MEISTER dihedral %d stop\n")% (i+1) );
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
    }
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY(BF("%s }\n")% this->className());
  }
  return this->_TotalEnergy;
}







void	EnergyDihedral_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                                 NVector_sp 	pos)
{_OF();
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
    double dhx1x1 = 0.0;
    double ohx1y1 = 0.0;
    double ohx1z1 = 0.0;
    double ohx1x2 = 0.0;
    double ohx1y2 = 0.0;
    double ohx1z2 = 0.0;
    double ohx1x3 = 0.0;
    double ohx1y3 = 0.0;
    double ohx1z3 = 0.0;
    double ohx1x4 = 0.0;
    double ohx1y4 = 0.0;
    double ohx1z4 = 0.0;
    double dhy1y1 = 0.0;
    double ohy1z1 = 0.0;
    double ohy1x2 = 0.0;
    double ohy1y2 = 0.0;
    double ohy1z2 = 0.0;
    double ohy1x3 = 0.0;
    double ohy1y3 = 0.0;
    double ohy1z3 = 0.0;
    double ohy1x4 = 0.0;
    double ohy1y4 = 0.0;
    double ohy1z4 = 0.0;
    double dhz1z1 = 0.0;
    double ohz1x2 = 0.0;
    double ohz1y2 = 0.0;
    double ohz1z2 = 0.0;
    double ohz1x3 = 0.0;
    double ohz1y3 = 0.0;
    double ohz1z3 = 0.0;
    double ohz1x4 = 0.0;
    double ohz1y4 = 0.0;
    double ohz1z4 = 0.0;
    double dhx2x2 = 0.0;
    double ohx2y2 = 0.0;
    double ohx2z2 = 0.0;
    double ohx2x3 = 0.0;
    double ohx2y3 = 0.0;
    double ohx2z3 = 0.0;
    double ohx2x4 = 0.0;
    double ohx2y4 = 0.0;
    double ohx2z4 = 0.0;
    double dhy2y2 = 0.0;
    double ohy2z2 = 0.0;
    double ohy2x3 = 0.0;
    double ohy2y3 = 0.0;
    double ohy2z3 = 0.0;
    double ohy2x4 = 0.0;
    double ohy2y4 = 0.0;
    double ohy2z4 = 0.0;
    double dhz2z2 = 0.0;
    double ohz2x3 = 0.0;
    double ohz2y3 = 0.0;
    double ohz2z3 = 0.0;
    double ohz2x4 = 0.0;
    double ohz2y4 = 0.0;
    double ohz2z4 = 0.0;
    double dhx3x3 = 0.0;
    double ohx3y3 = 0.0;
    double ohx3z3 = 0.0;
    double ohx3x4 = 0.0;
    double ohx3y4 = 0.0;
    double ohx3z4 = 0.0;
    double dhy3y3 = 0.0;
    double ohy3z3 = 0.0;
    double ohy3x4 = 0.0;
    double ohy3y4 = 0.0;
    double ohy3z4 = 0.0;
    double dhz3z3 = 0.0;
    double ohz3x4 = 0.0;
    double ohz3y4 = 0.0;
    double ohz3z4 = 0.0;
    double dhx4x4 = 0.0;
    double ohx4y4 = 0.0;
    double ohx4z4 = 0.0;
    double dhy4y4 = 0.0;
    double ohy4z4 = 0.0;
    double dhz4z4 = 0.0;
#endif
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,V,DN;
    double EraseLinearDihedral;
    int	I1, I2, I3, I4, IN, i;
    double sinPhase, cosPhase, SinNPhi, CosNPhi;
    gctools::Vec0<EnergyDihedral>::iterator di;
    for ( i=0,di=this->_Terms.begin();
          di!=this->_Terms.end(); di++,i++ ) {
#include <cando/chem/energy_functions/_Dihedral_termCode.cc>
      int index = i;
#include <cando/chem/energy_functions/_Dihedral_debugFiniteDifference.cc>
    }
  }

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

CL_DEFMETHOD core::List_sp EnergyDihedral_O::extract_vectors_as_alist() const{
  size_t size = this->_Terms.size();
  core::SimpleVector_double_sp v_vec = core::SimpleVector_double_O::make(size);
  core::SimpleVector_int32_t_sp in_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_double_sp phase_vec = core::SimpleVector_double_O::make(size);
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
  core::SimpleVector_double_sp v_vec = (safe_alist_lookup<core::SimpleVector_double_sp>(vectors,kw::_sym_v));
  core::SimpleVector_int32_t_sp in_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_in));
  core::SimpleVector_double_sp dn_vec = (safe_alist_lookup<core::SimpleVector_double_sp>(vectors,kw::_sym_dn));
  core::SimpleVector_double_sp phase_vec = (safe_alist_lookup<core::SimpleVector_double_sp>(vectors,kw::_sym_phase));
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

};
