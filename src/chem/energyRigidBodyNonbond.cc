/*
    File: energyNonbond.cc
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
#define	DEBUG_LEVEL_FULL

#include <clasp/core/foundation.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energyRigidBodyNonbond.h>
#include <cando/chem/rigidBodyEnergyFunction.h>
#include <clasp/core/array.h>
#include <clasp/core/evaluator.h>
#include <cando/geom/color.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/nativeVector.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/matter.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/nVector.h>
#include <cando/chem/loop.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
//#define LOG_ENERGY BFORMAT_T

namespace chem
{



EnergyRigidBodyNonbond_sp EnergyRigidBodyNonbond_O::make(core::Array_sp end_atoms) {
  if (end_atoms->length()<1) {
    SIMPLE_ERROR(BF("You must provide a vector of end atom indices with at least one end atom"));
  }
  core::SimpleVector_byte32_t_sp sv = core::SimpleVector_byte32_t_O::make(end_atoms->length());
  // Copy the atom indices
  size_t istart = 0;
  for ( size_t i(0), iEnd(end_atoms->length()); i<iEnd; ++i ) {
    Fixnum f = core::clasp_to_fixnum(end_atoms->rowMajorAref(i));
    if (istart>=f) {
      SIMPLE_ERROR(BF("The list of atom indices must be strictly increasing"));
    }
    (*sv)[i] = f;
  }
  GC_ALLOCATE_VARIADIC(EnergyRigidBodyNonbond_O,nb,sv);
  return nb;
}

void	EnergyRigidBodyNonbond_O::zeroEnergy()
{
  this->Base::zeroEnergy();
  this->_EnergyElectrostatic = 0.0;
  this->_EnergyVdw = 0.0;
}

double	EnergyRigidBodyNonbond_O::getEnergy()
{
  double	e;
  e = this->getVdwEnergy();
  e += this->getElectrostaticEnergy();
  return e;
}

void EnergyRigidBodyNonbond_O::energyRigidBodyNonbondSetTerm(gc::Fixnum index, Atom_sp atom, double radius, double epsilon, double charge, const Vector3& position) {
  if (index < 0 || index >= this->_AtomInfoTable.size()) {
    SIMPLE_ERROR(BF("Index out of range %d - max is %d") % index % this->_AtomInfoTable.size());
  }
  RigidBodyAtomInfo info(atom,radius,epsilon,charge,position);
  this->_AtomInfoTable[index] = info;
}



void	EnergyRigidBodyNonbond_O::dumpTerms()
{
  size_t istart = 0;
  for ( size_t rb = 0; rb<this->_RigidBodyEndAtom->length(); ++rb) {
    for ( size_t i = istart; i<(*this->_RigidBodyEndAtom)[rb]; ++i ) {
      RigidBodyAtomInfo& ai = this->_AtomInfoTable[i];
      BFORMAT_T(BF("I1 = %3d  %s  %lf %lf %lf %s\n") % (rb*7) % _rep_(ai._Atom) % ai._Radius % ai._Epsilon % ai._Charge % ai._Position.asString());
    }
    istart = (*this->_RigidBodyEndAtom)[rb];
  }
}

CL_DEFMETHOD core::List_sp EnergyRigidBodyNonbond_O::parts_as_list(NVector_sp pos)
{
  ql::list result;
  size_t istart = 0;
#undef	NONBOND_POSITION_RB_SET_PARAMETER
#define	NONBOND_POSITION_RB_SET_PARAMETER(x)	{}
#undef	NONBOND_POSITION_RB_SET_POSITION
#define	NONBOND_POSITION_RB_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_POSITION_RB_SET_POINT
#define	NONBOND_POSITION_RB_SET_POINT(x,ii,of)	{x=ii._Position.of;}
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_NONBOND_POSITIONS_termDeclares.cc>
#pragma clang diagnostic pop
  double am, bm, cm, dm, xm, ym, zm;
  double pxm, pym, pzm;
  int	I1;
  size_t I1start = 0;
  printf("%s:%d   _RigidBodyEndAtom->length() -> %lu\n", __FILE__, __LINE__, this->_RigidBodyEndAtom->length());
  for ( size_t iI1 = 0; iI1<this->_RigidBodyEndAtom->length(); ++iI1 ) {
    printf("%s:%d    iI1 = %lu\n", __FILE__, __LINE__, iI1);
    ql::list helix;
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
    printf("%s:%d    iI1 = %lu  I1start-> %lu I1end -> %lu\n", __FILE__, __LINE__, iI1, I1start, I1end);
    for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
      RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
      I1 = iI1*7;
#include <cando/energy-functions/_NONBOND_POSITIONS_termCode.cc>
      helix << core::Cons_O::createList(core::DoubleFloat_O::create(plabmx),
                                        core::DoubleFloat_O::create(plabmy),
                                        core::DoubleFloat_O::create(plabmz),
                                        core::DoubleFloat_O::create(ea1._Radius));
    }
    I1start = I1end;
    result << helix.cons();
  }
  return result.cons();
}


void	EnergyRigidBodyNonbond_O::setupHessianPreconditioner(
                                                    NVector_sp nvPosition,
                                                    AbstractLargeSquareMatrix_sp m )
{
  SIMPLE_ERROR(BF("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!!"));
}

double	EnergyRigidBodyNonbond_O::evaluateAll(NVector_sp 	pos,
                                     bool 		calcForce,
                                     gc::Nilable<NVector_sp> 	force,
                                     bool		calcDiagonalHessian,
                                     bool		calcOffDiagonalHessian,
                                     gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                     gc::Nilable<NVector_sp>	hdvec, 
                                     gc::Nilable<NVector_sp> 	dvec )
{
  double vdwScale = this->getVdwScale();
  double electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
  this->_EnergyElectrostatic = 0.0;
  this->_EnergyVdw = 0.0;
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
#define NONBONDRB_CALC_FORCE
#define NONBONDRB_CALC_DIAGONAL_HESSIAN
#define NONBONDRB_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBONDRB_SET_PARAMETER
#define	NONBONDRB_SET_PARAMETER(x)	{}
#undef	NONBONDRB_SET_POSITION
#define	NONBONDRB_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBONDRB_SET_POINT
#define	NONBONDRB_SET_POINT(x,ii,of)	{x=ii._Position.of;}
#undef	NONBONDRB_EEEL_ENERGY_ACCUMULATE
#define	NONBONDRB_EEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);}
#undef	NONBONDRB_EVDW_ENERGY_ACCUMULATE
#define	NONBONDRB_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);}
#undef	NONBONDRB_ENERGY_ACCUMULATE
#define	NONBONDRB_ENERGY_ACCUMULATE(e) {};
#undef	NONBONDRB_FORCE_ACCUMULATE
#undef	NONBONDRB_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBONDRB_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBONDRB_FORCE_ACCUMULATE 		ForceAcc
#define	NONBONDRB_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBONDRB_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
  if ( !this->isEnabled() ) return 0.0;
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
  LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_NONBONDRB_termDeclares.cc>
#pragma clang diagnostic pop
  double dA,dC,dQ1Q2;
  double am, bm, cm, dm, xm, ym, zm;
  double pxm, pym, pzm;
  double an, bn, cn, dn, xn, yn, zn;
  double pxn, pyn, pzn;
  int	I1, I2;
  size_t I1start = 0;
  for ( size_t iI1 = 0; iI1<(this->_RigidBodyEndAtom->length()-1); ++iI1 ) {
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
    size_t I2start = I1end;
    for (size_t iI2 = iI1+1; iI2<this->_RigidBodyEndAtom->length(); ++iI2 ) {
      size_t I2end = (*this->_RigidBodyEndAtom)[iI2];
      for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
        RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
        double charge1 = ea1._Charge;
        double electrostatic_scaled_charge1 = charge1*electrostaticScale;
        for ( size_t I2cur = I2start; I2cur<I2end; ++I2cur ) {
          RigidBodyAtomInfo& ea2 = this->_AtomInfoTable[I2cur];
          double rStar = ea1._Radius+ea2._Radius;
          double epsilonij = sqrt(ea1._Epsilon*ea2._Epsilon);
          double rStar2 = rStar*rStar;
          double rStar6 = rStar2*rStar2*rStar2;
          double rStar12 = rStar6*rStar6;
          dA = epsilonij*rStar12*vdwScale;
          dC = 2.0*epsilonij*rStar6*vdwScale;
          double charge2 = ea2._Charge;
          dQ1Q2 = electrostatic_scaled_charge1*charge2;
          I1 = iI1*7;
          I2 = iI2*7;
#include <cando/energy-functions/_NONBONDRB_termCode.cc>
#if 0
          printf("    %s:%d coordinate offsets %d - %d\n", __FILE__, __LINE__, I1, I2 );
          printf("      dA -> %lf   dC -> %lf  dQ1Q2 -> %lf  Energy -> %lf\n", dA, dC, dQ1Q2, Energy );
#endif
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
          nbi->_calcForce = calcForce;
          nbi->_calcDiagonalHessian = calcDiagonalHessian;
          nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/energy-functions/_NONBONDRB_debugEvalSet.cc>
#endif //]
#ifdef DEBUG_NONBOND_TERM
          if ( this->_DebugEnergy ) {
            std::string key;
            if (gc::IsA<Atom_sp>(ea1._Atom) && gc::IsA<Atom_sp>(ea2._Atom)) {
              Atom_sp ea1_atom = gc::As_unsafe<Atom_sp>(ea1._Atom);
              Atom_sp ea2_atom = gc::As_unsafe<Atom_sp>(ea2._Atom);
              if ( ea1_atom->getName()->symbolNameAsString()
                   < ea2_atom->getName()->symbolNameAsString() ) {
                key = ea1_atom->getName()->symbolNameAsString()+"-"+ea2_atom->getName()->symbolNameAsString();
              } else {
                key = ea2_atom->getName()->symbolNameAsString()+"-"+ea1_atom->getName()->symbolNameAsString();
              }
            }
            LOG_ENERGY(BF( "MEISTER nonbond %s args cando\n")% key );
            LOG_ENERGY(BF( "MEISTER nonbond %s dA %5.3lf\n")% key % dA );
            LOG_ENERGY(BF( "MEISTER nonbond %s dC %5.3lf\n")% key % dC );
            LOG_ENERGY(BF( "MEISTER nonbond %s dQ1Q2 %5.3lf\n")% key % dQ1Q2 );
            LOG_ENERGY(BF( "MEISTER nonbond %s x1 %5.3lf %d\n")% key % x1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s y1 %5.3lf %d\n")% key % y1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s z1 %5.3lf %d\n")% key % z1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s x2 %5.3lf %d\n")% key % x2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s y2 %5.3lf %d\n")% key % y2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s z2 %5.3lf %d\n")% key % z2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s results\n")% key );
            LOG_ENERGY(BF( "MEISTER nonbond %s evdw %lf\n")% key % Evdw);
            LOG_ENERGY(BF( "MEISTER nonbond %s eeel %lf\n")% key % Eeel);
            LOG_ENERGY(BF( "MEISTER nonbond %s Enonbond(evdw+eeel) %lf\n")% key % (Evdw+Eeel) );
            if ( calcForce ) {
              LOG_ENERGY(BF( "MEISTER nonbond %s fx1 %lf %d\n")% key % fx1 % (I1/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fy1 %lf %d\n")% key % fy1 % (I1/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fz1 %lf %d\n")% key % fz1 % (I1/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fx2 %lf %d\n")% key % fx2 % (I2/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fy2 %lf %d\n")% key % fy2 % (I2/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fz2 %lf %d\n")% key % fz2 % (I2/3+1) );
            }
            LOG_ENERGY(BF( "MEISTER nonbond %s stop\n")% key );
          }
#endif
        }
      }
    }
  }
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy }\n"));
  this->_TotalEnergy = this->_EnergyVdw+this->_EnergyElectrostatic;
#if 0
  printf("%s:%d NonbondRigidBody    EnergyVdw -> %lf\n", __FILE__, __LINE__, this->_EnergyVdw );
  printf("%s:%d NonbondRigidBody    EnergyElectrostatic -> %lf\n", __FILE__, __LINE__, this->_EnergyElectrostatic );
  printf("%s:%d NonbondRigidBody    TotalEnergy -> %lf\n", __FILE__, __LINE__, this->_TotalEnergy );
#endif
  return this->_TotalEnergy;
}
    
    



#if 0
ALWAYS_INLINE void _calculate_nonbond_term(EnergyRigidBodyNonbond_O& me, double vdwScale, double electrostatic_scaled_charge1, FFNonbond_sp ffNonbond1, FFNonbond_sp ffNonbond2, AtomTable_O::iterator iea1, AtomTable_O::iterator iea2,
                                           NVector_sp 	pos,
                                           bool 		calcForce,
                                           gc::Nilable<NVector_sp> 	force,
                                           bool		calcDiagonalHessian,
                                           bool		calcOffDiagonalHessian,
                                           gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                           gc::Nilable<NVector_sp>	hdvec, 
                                           gc::Nilable<NVector_sp> 	dvec )
{
#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {me._EnergyElectrostatic +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {me._EnergyVdw+=(e);}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_NONBONDRB_termDeclares.cc>
#pragma clang diagnostic pop
  double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
  double rStar = ffNonbond1->getRadius_Angstroms()+ffNonbond2->getRadius_Angstroms();
  // OPTIMIZE: Can I precalculate the square-roots!!!!!!!!
  double epsilonij = sqrt(ffNonbond1->getEpsilon_kCal()*ffNonbond2->getEpsilon_kCal()); 
  double rStar6 = pow(rStar,6.0);
  double rStar12 = rStar*rStar;
  dA = epsilonij*rStar12*vdwScale;
  dC = 2.0*epsilonij*rStar6*vdwScale;
  double charge2 = iea2->atom()->getCharge();
  dQ1Q2 = electrostatic_scaled_charge1*charge2;
  LOG(BF( "Calc dQ1Q2 electrostaticScale= %lf")% (double)(electrostaticScale));
  LOG(BF( "Calc dQ1Q2 Dielectric constant = %lf")% (double)(dielectricConstant));
  LOG(BF( "Calc dQ1Q2 Charge1 = %lf")% (double)(me._Charge1));
  LOG(BF( "Calc dQ1Q2 Charge2 = %lf")% (double)(me._Charge2));
  LOG(BF( "dQ1Q2 = %lf")% (double)(me.term.dQ1Q2));
  int I1 = iea1->coordinateIndexTimes3();
  int I2 = iea2->coordinateIndexTimes3();
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
  if ( me._Debug_NumberOfNonbondTermsToCalculate > 0 ) {
    if ( i>= me._Debug_NumberOfNonbondTermsToCalculate ) {
      break;
    }
  }
#endif
#include <cando/energy-functions/_NONBONDRB_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
  nbi->_calcForce = calcForce;
  nbi->_calcDiagonalHessian = calcDiagonalHessian;
  nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/energy-functions/_NONBONDRB_debugEvalSet.cc>
#endif //]
#ifdef DEBUG_NONBOND_TERM
  if ( me._DebugEnergy ) {
    LOG_ENERGY(BF( "MEISTER nonbond %d args cando\n")% (i+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d dA %5.3lf\n")% (i+1) % (nbi->term.dA) );
    LOG_ENERGY(BF( "MEISTER nonbond %d dC %5.3lf\n")% (i+1) % (nbi->term.dC) );
    LOG_ENERGY(BF( "MEISTER nonbond %d dQ1Q2 %5.3lf\n")% (i+1) % (nbi->term.dQ1Q2) );
    LOG_ENERGY(BF( "MEISTER nonbond %d x1 %5.3lf %d\n")% (i+1) % x1 % (I1/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d y1 %5.3lf %d\n")% (i+1) % y1 % (I1/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d z1 %5.3lf %d\n")% (i+1) % z1 % (I1/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d x2 %5.3lf %d\n")% (i+1) % x2 % (I2/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d y2 %5.3lf %d\n")% (i+1) % y2 % (I2/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d z2 %5.3lf %d\n")% (i+1) % z2 % (I2/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d results\n")% (i+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d Enonbond %lf\n")% (i+1) % (Evdw+Eeel) );
    LOG_ENERGY(BF( "MEISTER nonbond %d evdw %lf\n")% (i+1) % Evdw);
    LOG_ENERGY(BF( "MEISTER nonbond %d eeel %lf\n")% (i+1) % Eeel);
    if ( calcForce ) {
      LOG_ENERGY(BF( "MEISTER nonbond %d fx1 %lf %d\n")% (i+1) % fx1 % (I1/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fy1 %lf %d\n")% (i+1) % fy1 % (I1/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fz1 %lf %d\n")% (i+1) % fz1 % (I1/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fx2 %lf %d\n")% (i+1) % fx2 % (I2/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fy2 %lf %d\n")% (i+1) % fy2 % (I2/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fz2 %lf %d\n")% (i+1) % fz2 % (I2/3+1) );
    }
    LOG(BF( "MEISTER nonbond %d stop\n")% (i+1) );
  }
#endif
}
#endif






void	EnergyRigidBodyNonbond_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                                NVector_sp 	pos)
{_OF();
  IMPLEMENT_ME();
  // Use the code below
#if 0
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//------------------


#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{x=nbi->term.x;}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE(i,o,v) {}
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

  if ( this->isEnabled() ) {
    _BLOCK_TRACE("NonbondEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_NONBONDRB_termDeclares.cc>
#pragma clang diagnostic pop
    double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
    int	I1, I2,i;
    gctools::Vec0<EnergyRigidBodyNonbond>::iterator nbi;
    for ( i=0,nbi=this->_Terms.begin();
          nbi!=this->_Terms.end(); nbi++,i++ ) {
#include <cando/energy-functions/_NONBONDRB_termCode.cc>
      int index = i;
#include <cando/energy-functions/_NONBONDRB_debugFiniteDifference.cc>

    }
  }
#endif
}




void EnergyRigidBodyNonbond_O::initialize()
{
  this->Base::initialize();
  this->setErrorThreshold(1.0);
  this->setDielectricConstant(80.0);
  this->setVdwScale(1.0);
  this->setElectrostaticScale(1.0);
}






};
