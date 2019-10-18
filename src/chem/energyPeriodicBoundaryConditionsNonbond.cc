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
#include <cando/chem/energyPeriodicBoundaryConditionsNonbond.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/array.h>
#include <clasp/core/evaluator.h>
#include <cando/geom/color.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/nativeVector.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/matter.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/loop.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
//#define LOG_ENERGY core::write_bf_stream

namespace chem
{

inline double periodic_boundary_adjust(const double& delta, const double& rsize, const double& size)
{
  double fdelta = fabs(delta);
  fdelta -= static_cast<int>(fdelta*rsize + 0.5)*size;
  return fdelta;
}

/* Periodic boundary conditons/bounding-box are used in this code.
   So make these macros take care of that.
   See energyNonbond.cc for counter-example.
*/

#define PBX(_delta_) periodic_boundary_adjust(_delta_,x_rsize,x_size)
#define PBY(_delta_) periodic_boundary_adjust(_delta_,y_rsize,y_size)
#define PBZ(_delta_) periodic_boundary_adjust(_delta_,z_rsize,z_size)




double	_evaluateEnergyOnly_PeriodicBoundaryConditionsNonbond(double x1, double y1, double z1,
                                                              double x2, double y2, double z2,
                                                              double dA, double dC, double dQ1Q2)
{
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {}
#undef	NONBOND_FORCE_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE(i,o,v) {}
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(x) {}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(x) {}

#undef	NONBOND_CALC_FORCE	// Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>

  return Energy;
}



void	EnergyPeriodicBoundaryConditionsNonbond_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                                              AbstractLargeSquareMatrix_sp m )
{
  SIMPLE_ERROR(BF("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!!"));
}


void	EnergyPeriodicBoundaryConditionsNonbond_O::evaluateTerms(NVector_sp 	pos,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec, 
                                       gc::Nilable<NVector_sp> 	dvec )
{
//  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{x=nbi->term.x;}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);}

#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
  if ( this->isEnabled() ) 
  {
    gctools::Vec0<EnergyNonbond>::iterator firstElement = this->_Terms.begin();
    int nonBondTerms = this->_Terms.size();
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
    {
      LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>

#pragma clang diagnostic pop

      double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
      int	I1, I2,i;
      gctools::Vec0<EnergyNonbond>::iterator nbi;
      for ( i=0;i<nonBondTerms; i++ )
      {
        nbi = firstElement+i;
//        printf("%s:%d:%s i %d \n", __FILE__, __LINE__, __FUNCTION__,i );

#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
        if ( this->_Debug_NumberOfNonbondTermsToCalculate > 0 ) {
  printf("%s:%d:%s \n", __FILE__, __LINE__, __FUNCTION__ );
          if ( i>= this->_Debug_NumberOfNonbondTermsToCalculate ) {
            break;
          }
        }
#endif
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>

#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
        nbi->_calcForce = calcForce;
        nbi->_calcDiagonalHessian = calcDiagonalHessian;
        nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET

#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/chem/energy_functions/_Nonbond_debugEvalSet.cc>
#endif //]

#ifdef DEBUG_NONBOND_TERM
        if ( this->_DebugEnergy ) {
          std::string key;
          std::string key1 = nbi->_Atom1->getName()->symbolNameAsString();
          std::string key2 = nbi->_Atom2->getName()->symbolNameAsString();
          if ( key1<key2 ) {
            key = key1+"-"+key2;
          } else {
            key = key2+"-"+key1;
          }
          if (nbi->_Is14 ) {
            LOG_ENERGY(BF( "MEISTER nonbond %s 1-4 args cando\n")% key );
          } else {
            LOG_ENERGY(BF( "MEISTER nonbond %s args cando\n")% key );
          }
          LOG_ENERGY(BF( "MEISTER nonbond %s dA %5.3lf\n")% key % (nbi->term.dA) );
          LOG_ENERGY(BF( "MEISTER nonbond %s dC %5.3lf\n")% key % (nbi->term.dC) );
          LOG_ENERGY(BF( "MEISTER nonbond %s dQ1Q2 %5.3lf\n")% key % (nbi->term.dQ1Q2) );
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
  } else {
    LOG(BF("Nonbond component is not enabled"));
  }
//  printf( "Nonbond energy vdw(%lf) electrostatic(%lf)\n", (double)this->_EnergyVdw, this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF("Nonbond energy }\n"));
}



void	EnergyPeriodicBoundaryConditionsNonbond_O::evaluateUsingExcludedAtoms(NVector_sp 	pos,
                                                    bool 		calcForce,
                                                    gc::Nilable<NVector_sp> 	force,
                                                    bool		calcDiagonalHessian,
                                                    bool		calcOffDiagonalHessian,
                                                    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                    gc::Nilable<NVector_sp>	hdvec, 
                                                    gc::Nilable<NVector_sp> 	dvec )
{
//  printf("%s:%d In evaluateUsingExcludedAtoms starting this->_DebugEnergy -> %d\n", __FILE__, __LINE__, this->_DebugEnergy );
  if (!this->_iac_vec) {
    SIMPLE_ERROR(BF("The nonbonded excluded atoms parameters have not been set up"));
  }
  core::SimpleVector_int32_t_sp numberOfExcludedAtoms = this->_NumberOfExcludedAtomIndices;
  core::SimpleVector_int32_t_sp excludedAtomIndices = this->_ExcludedAtomIndices;
  double vdwScale = this->getVdwScale();
  double electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
//  printf("%s:%d electrostaticcharge %lf\n", __FILE__, __LINE__, electrostaticScale );
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
//  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  if ( !this->isEnabled() ) return;
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
  LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
  // printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2,dA_old,dC_old,dQ1Q2_old;
  int	I1, I2;
  int i = 0;
  int maxIndex = pos->length()/3;
//  int maxIndex = this->_AtomTable->getNumberOfAtoms();
  int excludedAtomIndex = 0;
  int nlocaltype = 0;
  // Find the max local type
  for (i=0; i<this->_iac_vec->length(); ++i){
    if (nlocaltype < (*this->_iac_vec)[i]){
      nlocaltype = (*this->_iac_vec)[i];
    }
  }

  int index1_end = maxIndex-1;
  for ( int index1 = 0; index1 <index1_end; ++index1 ) {
    LOG(BF("%s ====== top of outer loop - index1 = %d\n") % __FUNCTION__ % index1 );
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    LOG(BF("Read numberOfExcludedAtomsRemaining from numberOfExcludedAtoms[%d]= %d\n") % index1 % numberOfExcludedAtomsRemaining);
          // Skip 0 in excluded atom list that amber requires
    int maybe_excluded_atom = (*excludedAtomIndices)[excludedAtomIndex];
    if (maybe_excluded_atom<0) {
      ++excludedAtomIndex;
    } else {
      double charge11 = (*this->_charge_vector)[index1];
      double electrostatic_scaled_charge11 = charge11*electrostaticScale;
      int number_excluded = 0;
      for ( int index2 = index1+1, index2_end(maxIndex); index2 < index2_end; ++index2 ) {
        LOG(BF("    --- top of inner loop   numberOfExcludedAtomsRemaining -> %d    index2 -> %d\n") % numberOfExcludedAtomsRemaining % index2 );
        maybe_excluded_atom = (*excludedAtomIndices)[excludedAtomIndex];
        if ( numberOfExcludedAtomsRemaining>0 && maybe_excluded_atom == index2) {
          LOG(BF("    Excluding atom %d\n") % index2);
          ++excludedAtomIndex;
          --numberOfExcludedAtomsRemaining;
          number_excluded++;
          continue;
        }
      //   cn1i    (*this->_cn1_vec)[cn1i]
        int localindex1 = (*this->_iac_vec)[index1];
        int localindex2 = (*this->_iac_vec)[index2];
        dA = (*this->_cn1_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
        dC = (*this->_cn2_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
//      printf("%s:%d localindex1 %d and localindex2 %d\n", __FILE__, __LINE__, localindex1, localindex2);
//      printf("%s:%d dA     %lf and dC     %lf\n", __FILE__, __LINE__, dA, dC);
        double charge22 = (*this->_charge_vector)[index2];
//      dQ1Q2 = electrostatic_scaled_charge11*charge22;
        dQ1Q2 = charge11*charge22;
//      printf("%s:%d charge1     %lf and charge2     %lf\n", __FILE__, __LINE__, charge11, charge22);
//      printf("%s:%d electrostaticScale     %lf and dQ1Q2     %lf\n", __FILE__, __LINE__, electrostaticScale, dQ1Q2);
#ifdef DEBUG_NONBOND_TERM
        if ( this->_DebugEnergy ) {
          LOG_ENERGY(BF( "nonbond localindex1 %d\n")% localindex1 );
          LOG_ENERGY(BF( "nonbond localindex2 %d\n")% localindex2 );
          LOG_ENERGY(BF( "nonbond dA %5.3lf\n")% dA );
          LOG_ENERGY(BF( "nonbond dC %5.3lf\n")% dC );
          LOG_ENERGY(BF( "nonbond dQ1Q2 %5.3lf\n")% dQ1Q2 );
        }
#endif
      ////////////////////////////////////////////////////////////
      //
      // To here
      //
      ////////////////////////////////////////////////////////////
        I1 = index1*3; 
        I2 = index2*3; 
#ifdef DEBUG_NONBOND_TERM
        if ( this->_DebugEnergy ) {
          LOG_ENERGY(BF( "nonbond I1 %d\n")% I1 );
          LOG_ENERGY(BF( "nonbond I2 %d\n")% I2 );
        }
#endif
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
        if ( this->_Debug_NumberOfNonbondTermsToCalculate > 0 ) {
          if ( i>= this->_Debug_NumberOfNonbondTermsToCalculate ) {
            break;
          }
        }
#endif
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
        nbi->_calcForce = calcForce;
        nbi->_calcDiagonalHessian = calcDiagonalHessian;
        nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/chem/energy_functions/_Nonbond_debugEvalSet.cc>
#endif //]
#ifdef DEBUG_NONBOND_TERM
        if ( this->_DebugEnergy ) {
          std::string key;
          std::string atom1Name = gc::As<core::Symbol_sp>((*this->_atom_name_vector)[index1])->symbolNameAsString();
          std::string atom2Name = gc::As<core::Symbol_sp>((*this->_atom_name_vector)[index2])->symbolNameAsString();
          if ( atom1Name < atom2Name ) {
            key = atom1Name+"-"+atom2Name;
          } else {
            key = atom2Name+"-"+atom1Name;
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
//  printf( "Nonbond energy vdw(%lf) electrostatic(%lf)\n", (double)this->_EnergyVdw,  this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy }\n"));
}


void	EnergyPeriodicBoundaryConditionsNonbond_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp 	pos)
{_OF();
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
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
    double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
    int	I1, I2,i;
    gctools::Vec0<EnergyNonbond>::iterator nbi;
    for ( i=0,nbi=this->_Terms.begin();
          nbi!=this->_Terms.end(); nbi++,i++ ) {
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
      int index = i;
#include <cando/chem/energy_functions/_Nonbond_debugFiniteDifference.cc>

    }
  }
}

SYMBOL_EXPORT_SC_(KeywordPkg,nonbond);
SYMBOL_EXPORT_SC_(KeywordPkg,vdwRadii);
SYMBOL_EXPORT_SC_(KeywordPkg,distance);
SYMBOL_EXPORT_SC_(KeywordPkg,force);

core::List_sp	EnergyPeriodicBoundaryConditionsNonbond_O::checkForBeyondThresholdInteractionsWithPosition(NVector_sp pos, double threshold )
{
  ql::list result;
  bool calcForce = true;
  printf("%s:%d:%s   number of entries pos -> %lu\n", __FILE__, __LINE__, __FUNCTION__, pos->length()/3);
  
//  printf("%s:%d In evaluateUsingExcludedAtoms starting this->_DebugEnergy -> %d\n", __FILE__, __LINE__, this->_DebugEnergy );
  if (!this->_iac_vec) {
    SIMPLE_ERROR(BF("The nonbonded excluded atoms parameters have not been set up"));
  }
  core::SimpleVector_int32_t_sp numberOfExcludedAtoms = this->_NumberOfExcludedAtomIndices;
  core::SimpleVector_int32_t_sp excludedAtomIndices = this->_ExcludedAtomIndices;
  double vdwScale = this->getVdwScale();
  double electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
//  printf("%s:%d electrostaticcharge %lf\n", __FILE__, __LINE__, electrostaticScale );
#define NONBOND_CALC_FORCE
#undef  NONBOND_CALC_DIAGONAL_HESSIAN
#undef  NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE(ix,iy,iz) 		{}
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	{}
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE {}
//  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
  LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
  // printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2,dA_old,dC_old,dQ1Q2_old;
  int	I1, I2;
  int i = 0;
  int maxIndex = pos->length()/3;
//  int maxIndex = this->_AtomTable->getNumberOfAtoms();
  int excludedAtomIndex = 0;
  int nlocaltype = 0;
  // Find the max local type
  for (i=0; i<this->_iac_vec->length(); ++i){
    if (nlocaltype < (*this->_iac_vec)[i]){
      nlocaltype = (*this->_iac_vec)[i];
    }
  }
  int index1_end = maxIndex-1;
  for ( int index1 = 0; index1<index1_end; ++index1 ) {
    LOG(BF("%s ====== top of outer loop - index1 = %d\n") % __FUNCTION__ % index1 );
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
//    printf("%s:%d:%s index1 = %d numberOfExcludedAtomsRemaining = %d  excludedAtomIndex = %d\n", __FILE__, __LINE__, __FUNCTION__, index1, numberOfExcludedAtomsRemaining, excludedAtomIndex);
    LOG(BF("Read numberOfExcludedAtomsRemaining from numberOfExcludedAtoms[%d]= %d\n") % index1 % numberOfExcludedAtomsRemaining);
          // Skip 0 in excluded atom list that amber requires
    int maybe_excluded_atom = (*excludedAtomIndices)[excludedAtomIndex];
    if (maybe_excluded_atom<0) {
//      printf("%s:%d:%s SKIPPING numberOfExcludedAtomsRemaining<0\n", __FILE__, __LINE__, __FUNCTION__);
      ++excludedAtomIndex;
    } else {
      double charge11 = (*this->_charge_vector)[index1];
      double electrostatic_scaled_charge11 = charge11*electrostaticScale;
      int number_excluded = 0;
      for ( int index2 = index1+1, index2_end(maxIndex); index2 < index2_end; ++index2 ) {
        LOG(BF("    --- top of inner loop   numberOfExcludedAtomsRemaining -> %d    index2 -> %d\n") % numberOfExcludedAtomsRemaining % index2 );
        maybe_excluded_atom = (*excludedAtomIndices)[excludedAtomIndex];
        if ( numberOfExcludedAtomsRemaining>0 && maybe_excluded_atom == index2) {
          LOG(BF("    Excluding atom %d\n") % index2);
          ++excludedAtomIndex;
          --numberOfExcludedAtomsRemaining;
          number_excluded++;
//          printf("%s:%d:%s Excluding atom %d\n", __FILE__, __LINE__, __FUNCTION__, index2);
          continue;
        }
      //   cn1i    (*this->_cn1_vec)[cn1i]
      // printf("%s:%d:%s Interaction between atom %d = %d\n", __FILE__, __LINE__, __FUNCTION__, index1, index2);
        int localindex1 = (*this->_iac_vec)[index1];
        int localindex2 = (*this->_iac_vec)[index2];
        dA = (*this->_cn1_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
        dC = (*this->_cn2_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
//      printf("%s:%d localindex1 %d and localindex2 %d\n", __FILE__, __LINE__, localindex1, localindex2);
//      printf("%s:%d dA     %lf and dC     %lf\n", __FILE__, __LINE__, dA, dC);
        double charge22 = (*this->_charge_vector)[index2];
//      dQ1Q2 = electrostatic_scaled_charge11*charge22;
        dQ1Q2 = charge11*charge22;
//      printf("%s:%d charge1     %lf and charge2     %lf\n", __FILE__, __LINE__, charge11, charge22);
//      printf("%s:%d electrostaticScale     %lf and dQ1Q2     %lf\n", __FILE__, __LINE__, electrostaticScale, dQ1Q2);
#ifdef DEBUG_NONBOND_TERM
        if ( this->_DebugEnergy ) {
          LOG_ENERGY(BF( "nonbond localindex1 %d\n")% localindex1 );
          LOG_ENERGY(BF( "nonbond localindex2 %d\n")% localindex2 );
          LOG_ENERGY(BF( "nonbond dA %5.3lf\n")% dA );
          LOG_ENERGY(BF( "nonbond dC %5.3lf\n")% dC );
          LOG_ENERGY(BF( "nonbond dQ1Q2 %5.3lf\n")% dQ1Q2 );
        }
#endif
      ////////////////////////////////////////////////////////////
      //
      // To here
      //
      ////////////////////////////////////////////////////////////
        I1 = index1*3; 
        I2 = index2*3; 
#ifdef DEBUG_NONBOND_TERM
        if ( this->_DebugEnergy ) {
          LOG_ENERGY(BF( "nonbond I1 %d\n")% I1 );
          LOG_ENERGY(BF( "nonbond I2 %d\n")% I2 );
        }
#endif
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
        if ( this->_Debug_NumberOfNonbondTermsToCalculate > 0 ) {
          if ( i>= this->_Debug_NumberOfNonbondTermsToCalculate ) {
            break;
          }
        }
#endif
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
        EnergyAtom& ea1 = this->_AtomTable->energyAtomEntry(index1);
        EnergyAtom& ea2 = this->_AtomTable->energyAtomEntry(index2);
        if ( NonbondDistance < this->_ErrorThreshold ||
             (ea1._Flag && ea2._Flag)) {
          Atom_sp a1, a2;
          a1 = ea1.atom();
          a2 = ea2.atom();
          ql::list one;
          one << kw::_sym_nonbond
              << kw::_sym_atoms << core::Cons_O::createList(a1,a2)
              << kw::_sym_distance << core::clasp_make_double_float(NonbondDistance)
              << kw::_sym_force << (ql::list()
                                    << geom::OVector3_O::make(fx1,fy1,fz1)
                                    << geom::OVector3_O::make(fx2,fy2,fz2)).result();
          result << one.result();
        }
      }
      if (numberOfExcludedAtoms->operator[](index1)!=number_excluded) {
        printf("%s:%d:%s ERROR the number of excluded atoms should be %d but we excluded %d\n", __FILE__, __LINE__, __FUNCTION__, numberOfExcludedAtoms->operator[](index1), number_excluded);
      }
    }
  }
  return result.result();
}

};

