/*
    File: energyPeriodicBoundaryConditionsNonbond.cc
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
#include <cando/chem/energyNonbond.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

#define BAIL_OUT_IF_CUTOFF(deltaSquared) if (deltaSquared>CUTOFF_SQUARED) goto SKIP_term;

//#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
//#define LOG_ENERGY core::clasp_write_string

namespace chem
{

core::T_sp pbnonbond_type(bool is14) {
  if (is14) return chem::_sym_EnergyNonbond14;
  return chem::_sym_EnergyNonbond;
}

#define NONBOND_APPLY_ATOM_MASK(I1,I2) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         ) \
    ) goto SKIP_term;

#define NONBOND_DEBUG_INTERACTIONS(I1,I2) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,pbnonbond_type(InteractionIs14), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2)); \
    }


inline double periodic_boundary_adjust(const num_real& delta, const num_real& rsize, const num_real& size)
{
  double result = delta;
  result -= size*std::nearbyint(delta*rsize);
  return result;
}

/* Periodic boundary conditons/bounding-box are used in this code.
   So make these macros take care of that.
   See energyNonbond.cc for counter-example.
*/

#define PBX(_delta_) periodic_boundary_adjust(_delta_,x_rwidth,x_width)
#define PBY(_delta_) periodic_boundary_adjust(_delta_,y_rwidth,y_width)
#define PBZ(_delta_) periodic_boundary_adjust(_delta_,z_rwidth,z_width)




double	_evaluateEnergyOnly_PeriodicBoundaryConditionsNonbond(ScoringFunction_sp score,
                                                              num_real x1, num_real y1, num_real z1,
                                                              num_real x2, num_real y2, num_real z2,
                                                              num_real dA, num_real dC, num_real dQ1Q2)
{
  SIMPLE_ERROR("Support {}" , __FUNCTION__);
#if 0
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

  EnergyFunction_sp energyFunction = gc::As<EnergyFunction_sp>(score);
  num_real x_rwidth = energyFunction->boundingBox()->get_x_rwidth();
  num_real y_rwidth = energyFunction->boundingBox()->get_y_rwidth();
  num_real z_rwidth = energyFunction->boundingBox()->get_z_rwidth();
  num_real x_width = energyFunction->boundingBox()->get_x_width();
  num_real y_width = energyFunction->boundingBox()->get_y_width();
  num_real z_width = energyFunction->boundingBox()->get_z_width();
  num_real half_x_size = x_width*0.5;
  num_real half_y_size = y_width*0.5;
  num_real half_z_size = z_width*0.5;
#undef	NONBOND_CALC_FORCE	// Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
  double vdwScale = 1.0;
  double eelScale = 1.0;
  double DIELECTRIC = 1.0;
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#undef CUTOFF_SQUARED
#undef DIELECTRIC
 DONE:
  return Energy;
#endif

}


void	EnergyPeriodicBoundaryConditionsNonbond_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                                              AbstractLargeSquareMatrix_sp m,
                                                                              core::T_sp activeAtomMask )
{
  SIMPLE_ERROR("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!!");
}


void	EnergyPeriodicBoundaryConditionsNonbond_O::evaluateTerms(ScoringFunction_sp score,
                                                                 NVector_sp 	pos,
                                                                 bool 		calcForce,
                                                                 gc::Nilable<NVector_sp> 	force,
                                                                 bool		calcDiagonalHessian,
                                                                 bool		calcOffDiagonalHessian,
                                                                 gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                                 gc::Nilable<NVector_sp>	hdvec, 
                                                                 gc::Nilable<NVector_sp> 	dvec,
                                                                 core::T_sp activeAtomMask )
{
  double dielectricConstant;
  double dQ1Q2Scale;
  double Evdw = 0.0;
  double Eeel = 0.0;
  double cutoff;
  energyFunctionNonbondParameters(score,dielectricConstant,dQ1Q2Scale,cutoff);
#define CUTOFF_SQUARED (cutoff*cutoff)
  double nonbondCutoffSquared = cutoff*cutoff;
  double nonbondCutoffReciprocal6 = 1.0/(cutoff*cutoff*cutoff*cutoff*cutoff*cutoff);
  double nonbondCutoffReciprocal12 = nonbondCutoffReciprocal6*nonbondCutoffReciprocal6;

  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  core::T_sp debugInteractions = nil<core::T_O>();
  MAYBE_SETUP_DEBUG_INTERACTIONS(false);
  
//  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  this->_Evaluations++;
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
  EnergyFunction_sp energyFunction = gc::As<EnergyFunction_sp>(score);
  num_real x_rwidth = energyFunction->boundingBox()->get_x_rwidth();
  num_real y_rwidth = energyFunction->boundingBox()->get_y_rwidth();
  num_real z_rwidth = energyFunction->boundingBox()->get_z_rwidth();
  num_real x_width = energyFunction->boundingBox()->get_x_width();
  num_real y_width = energyFunction->boundingBox()->get_y_width();
  num_real z_width = energyFunction->boundingBox()->get_z_width();
  [[maybe_unused]]num_real half_x_size = x_width*0.5;
  [[maybe_unused]]num_real half_y_size = y_width*0.5;
  [[maybe_unused]]num_real half_z_size = z_width*0.5;

#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{x=nbi->term.x;}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {Eeel +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {Evdw +=(e)-vljrc;}

#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
  {
    gctools::Vec0<EnergyNonbond>::iterator firstElement = this->_Terms.begin();
    int nonBondTerms = this->_Terms.size();
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
    {
      LOG("Nonbond component is enabled" );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>

#pragma clang diagnostic pop

      num_real x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
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
        num_real vljrc = nbi->term.dA*nonbondCutoffReciprocal12-nbi->term.dC*nonbondCutoffReciprocal6;
        bool InteractionIs14 = nbi->_Is14;
  double vdwScale = 1.0;
  double eelScale = 1.0;
  double DIELECTRIC = 1.0;
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#undef CUTOFF_SQUARED
#undef DIELECTRIC
      DONE:

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
            LOG_ENERGY(( "MEISTER nonbond %s 1-4 args cando\n") , key );
          } else {
            LOG_ENERGY(( "MEISTER nonbond %s args cando\n") , key );
          }
          LOG_ENERGY(( "MEISTER nonbond %s dA %5.3lf\n") , key , (nbi->term.dA) );
          LOG_ENERGY(( "MEISTER nonbond %s dC %5.3lf\n") , key , (nbi->term.dC) );
          LOG_ENERGY(( "MEISTER nonbond %s dQ1Q2 %5.3lf\n") , key , (nbi->term.dQ1Q2) );
          LOG_ENERGY(( "MEISTER nonbond %s x1 %5.3lf %d\n") , key , x1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s y1 %5.3lf %d\n") , key , y1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s z1 %5.3lf %d\n") , key , z1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s x2 %5.3lf %d\n") , key , x2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s y2 %5.3lf %d\n") , key , y2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s z2 %5.3lf %d\n") , key , z2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s results\n") , key );
          LOG_ENERGY(( "MEISTER nonbond %s evdw %lf\n") , key , Evdw);
          LOG_ENERGY(( "MEISTER nonbond %s eeel %lf\n") , key , Eeel);
          LOG_ENERGY(( "MEISTER nonbond %s Enonbond(evdw+eeel) %lf\n") , key , (Evdw+Eeel) );
          if ( calcForce ) {
            LOG_ENERGY(( "MEISTER nonbond %s fx1 %lf %d\n") , key , fx1 , (I1/3+1) );
            LOG_ENERGY(( "MEISTER nonbond %s fy1 %lf %d\n") , key , fy1 , (I1/3+1) );
            LOG_ENERGY(( "MEISTER nonbond %s fz1 %lf %d\n") , key , fz1 , (I1/3+1) );
            LOG_ENERGY(( "MEISTER nonbond %s fx2 %lf %d\n") , key , fx2 , (I2/3+1) );
            LOG_ENERGY(( "MEISTER nonbond %s fy2 %lf %d\n") , key , fy2 , (I2/3+1) );
            LOG_ENERGY(( "MEISTER nonbond %s fz2 %lf %d\n") , key , fz2 , (I2/3+1) );
          }
          LOG_ENERGY(( "MEISTER nonbond %s stop\n") , key );
        }
#endif
        (void)0;
      }
    }
  } 
//  printf( "Nonbond energy vdw(%lf) electrostatic(%lf)\n", (double)this->_EnergyVdw, this->_EnergyElectrostatic );
//  LOG( "Nonbond energy vdw({}) electrostatic({})\n" , (double)this->_EnergyVdw , this->_EnergyElectrostatic );
}



void	EnergyPeriodicBoundaryConditionsNonbond_O::evaluateUsingExcludedAtoms(ScoringFunction_sp score,
                                                                              NVector_sp 	pos,
                                                                              bool 		calcForce,
                                                                              gc::Nilable<NVector_sp> 	force,
                                                                              bool		calcDiagonalHessian,
                                                                              bool		calcOffDiagonalHessian,
                                                                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                                              gc::Nilable<NVector_sp>	hdvec, 
                                                                              gc::Nilable<NVector_sp> 	dvec,
                                                                              core::T_sp activeAtomMask )
{
  double dielectricConstant;
  double dQ1Q2Scale;
  double cutoff;
  EnergyFunction_sp energyFunction = energyFunctionNonbondParameters(score,dielectricConstant,dQ1Q2Scale,cutoff);
  double nonbondCutoffSquared = cutoff*cutoff;
#define CUTOFF_SQUARED nonbondCutoffSquared
  double nonbondCutoffReciprocal6 = 1.0/(cutoff*cutoff*cutoff*cutoff*cutoff*cutoff);
  double nonbondCutoffReciprocal12 = nonbondCutoffReciprocal6*nonbondCutoffReciprocal6;
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  core::T_sp debugInteractions = nil<core::T_O>();
  MAYBE_SETUP_DEBUG_INTERACTIONS(false);
  this->_Evaluations++;
//  printf("%s:%d In evaluateUsingExcludedAtoms starting this->_DebugEnergy -> %d\n", __FILE__, __LINE__, this->_DebugEnergy );
  if (!this->_iac_vec) {
    SIMPLE_ERROR("The nonbonded excluded atoms parameters have not been set up");
  }
  core::SimpleVector_int32_t_sp numberOfExcludedAtoms = this->_NumberOfExcludedAtomIndices;
  core::SimpleVector_int32_t_sp excludedAtomIndices = this->_ExcludedAtomIndices;
  [[maybe_unused]]num_real vdwScale = 1.0;
  num_real electrostaticScale = dQ1Q2Scale;
#define DIELECTRIC dielectricConstant
//  printf("%s:%d electrostaticcharge %lf\n", __FILE__, __LINE__, electrostaticScale );
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
  num_real x_rwidth = energyFunction->boundingBox()->get_x_rwidth();
  num_real y_rwidth = energyFunction->boundingBox()->get_y_rwidth();
  num_real z_rwidth = energyFunction->boundingBox()->get_z_rwidth();
  num_real x_width = energyFunction->boundingBox()->get_x_width();
  num_real y_width = energyFunction->boundingBox()->get_y_width();
  num_real z_width = energyFunction->boundingBox()->get_z_width();
  [[maybe_unused]]num_real half_x_size = x_width*0.5;
  [[maybe_unused]]num_real half_y_size = y_width*0.5;
  [[maybe_unused]]num_real half_z_size = z_width*0.5;
  
#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {Eeel +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {Evdw +=(e)-vljrc;}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
  // If you are going to use openmp here, you need to control access to the force and hessian
  // arrays so that only one thread updates each element at a time.
  LOG("Nonbond component is enabled" );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
  // printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  num_real x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
  int	I1, I2;
  int i = 0;
  int endIndex = pos->length()/3;
//  int endIndex = this->_AtomTable->getNumberOfAtoms();
  int excludedAtomIndex = 0;
  int nlocaltype = 0;
  // Find the max local type
  for (i=0; i<this->_iac_vec->length(); ++i){
    if (nlocaltype < (*this->_iac_vec)[i]){
      nlocaltype = (*this->_iac_vec)[i];
    }
  }
//  printf("%s:%d About to do nonbond calc endIndex -> %d\n", __FILE__, __LINE__, endIndex );
  int index1_end = endIndex-1;
  for ( int index1 = 0; index1 <index1_end; ++index1 ) {
    LOG("{} ====== top of outer loop - index1 = {}\n" , __FUNCTION__ , index1 );
          // Skip 0 in excluded atom list that amber requires
    bool has_excluded_atoms = ((*excludedAtomIndices)[excludedAtomIndex] >= 0);
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    num_real charge11 = (*this->_charge_vector)[index1];
    num_real electrostatic_scaled_charge11 = charge11*electrostaticScale;
    for ( int index2 = index1+1, index2_end(endIndex); index2 < index2_end; ++index2 ) {
      LOG("    --- top of inner loop   numberOfExcludedAtomsRemaining -> {}    index2 -> {}\n" , numberOfExcludedAtomsRemaining , index2 );
      int maybe_excluded_atom = (*excludedAtomIndices)[excludedAtomIndex];
      if (numberOfExcludedAtomsRemaining>0 && maybe_excluded_atom == index2) {
        LOG("    Excluding atom {}\n" , index2);
        ++excludedAtomIndex;
        --numberOfExcludedAtomsRemaining;
        continue;
      }
      //   cn1i    (*this->_cn1_vec)[cn1i]
      int localindex1 = (*this->_iac_vec)[index1];
      int localindex2 = (*this->_iac_vec)[index2];
      dA = (*this->_cn1_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
      dC = (*this->_cn2_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
      num_real vljrc = dA*nonbondCutoffReciprocal12-dC*nonbondCutoffReciprocal6;
//        printf("%s:%d localindex1 %d and localindex2 %d\n", __FILE__, __LINE__, localindex1, localindex2);
//        printf("%s:%d dA     %lf and dC     %lf\n", __FILE__, __LINE__, dA, dC);
      num_real charge22 = (*this->_charge_vector)[index2];
      dQ1Q2 = electrostatic_scaled_charge11*charge22;
      dQ1Q2 = charge11*charge22;
#if 0
      core::clasp_write_string(fmt::format("{}:{}  dA -> {}   dC -> {}  dQ1Q2 -> {}\n" , __FILE__ , __LINE__ , dA , dC , dQ1Q2));
#endif
//        printf("%s:%d charge1     %lf and charge2     %lf\n", __FILE__, __LINE__, charge11, charge22);
//        printf("%s:%d electrostaticScale     %lf and dQ1Q2     %lf\n", __FILE__, __LINE__, electrostaticScale, dQ1Q2);
#ifdef DEBUG_NONBOND_TERM
      if ( this->_DebugEnergy ) {
        LOG_ENERGY(( "nonbond localindex1 %d\n") , localindex1 );
        LOG_ENERGY(( "nonbond localindex2 %d\n") , localindex2 );
        LOG_ENERGY(( "nonbond dA %5.3lf\n") , dA );
        LOG_ENERGY(( "nonbond dC %5.3lf\n") , dC );
        LOG_ENERGY(( "nonbond dQ1Q2 %5.3lf\n") , dQ1Q2 );
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
        LOG_ENERGY(( "nonbond I1 %d\n") , I1 );
        LOG_ENERGY(( "nonbond I2 %d\n") , I2 );
      }
#endif
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfNonbondTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfNonbondTermsToCalculate ) {
          break;
        }
      }
#endif
      bool InteractionIs14 = false; // Always false for excluded atoms
  double vdwScale = 1.0;
  double eelScale = 1.0;
  double DIELECTRIC = 1.0;
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#undef CUTOFF_SQUARED
#undef DIELECTRIC
    DONE:
#if 0
      core::clasp_write_string(fmt::format("{}:{}  DistanceSquared -> {}   Evdw -> {}  Evdw_rc -> {}   Eeel -> {}\n" , __FILE__ , __LINE__ , DistanceSquared , Evdw , vljrc , Eeel));
      core::clasp_write_string(fmt::format("{}:{}  x1, y1, z1 -> {}, {}, {}\n" , __FILE__ , __LINE__ , x1 , y1 , z1 ));
      core::clasp_write_string(fmt::format("{}:{}  x2, y2, z2 -> {}, {}, {}\n" , __FILE__ , __LINE__ , x2 , y2 , z2 ));
      core::clasp_write_string(fmt::format("{}:{}  fx1, fy1, fz1 -> {}, {}, {}\n" , __FILE__ , __LINE__ , fx1 , fy1 , fz1 ));
      core::clasp_write_string(fmt::format("{}:{}  fx2, fy2, fz2 -> {}, {}, {}\n" , __FILE__ , __LINE__ , fx2 , fy2 , fz2 ));
#endif
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
        LOG_ENERGY(( "MEISTER nonbond %s args cando\n") , key );
        LOG_ENERGY(( "MEISTER nonbond %s dA %5.3lf\n") , key , dA );
        LOG_ENERGY(( "MEISTER nonbond %s dC %5.3lf\n") , key , dC );
        LOG_ENERGY(( "MEISTER nonbond %s dQ1Q2 %5.3lf\n") , key , dQ1Q2 );
        LOG_ENERGY(( "MEISTER nonbond %s x1 %5.3lf %d\n") , key , x1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER nonbond %s y1 %5.3lf %d\n") , key , y1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER nonbond %s z1 %5.3lf %d\n") , key , z1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER nonbond %s x2 %5.3lf %d\n") , key , x2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER nonbond %s y2 %5.3lf %d\n") , key , y2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER nonbond %s z2 %5.3lf %d\n") , key , z2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER nonbond %s results\n") , key );
        LOG_ENERGY(( "MEISTER nonbond %s evdw %lf\n") , key , Evdw);
        LOG_ENERGY(( "MEISTER nonbond %s eeel %lf\n") , key , Eeel);
        LOG_ENERGY(( "MEISTER nonbond %s Enonbond(evdw+eeel) %lf\n") , key , (Evdw+Eeel) );
        if ( calcForce ) {
          LOG_ENERGY(( "MEISTER nonbond %s fx1 %lf %d\n") , key , fx1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s fy1 %lf %d\n") , key , fy1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s fz1 %lf %d\n") , key , fz1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s fx2 %lf %d\n") , key , fx2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s fy2 %lf %d\n") , key , fy2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER nonbond %s fz2 %lf %d\n") , key , fz2 , (I2/3+1) );
        }
        LOG_ENERGY(( "MEISTER nonbond %s stop\n") , key );
      }
#endif
      (void)0;
    }
    if (!has_excluded_atoms) {
      // No excluded atoms for the current index2 so increment excludedAtomIndex because inner loop didn't do it.
      // see Swails document http://ambermd.org/prmtop.pdf
      ++excludedAtomIndex;
    }
  }
//  printf( "Nonbond energy vdw(%lf) electrostatic(%lf)\n", (double)this->_EnergyVdw,  this->_EnergyElectrostatic );
  LOG( "Nonbond energy vdw({}) electrostatic({})\n" , (double)this->_EnergyVdw , this->_EnergyElectrostatic );
  LOG( "Nonbond energy }\n");
}


void	EnergyPeriodicBoundaryConditionsNonbond_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(ScoringFunction_sp score, NVector_sp 	pos)
{
  IMPLEMENT_ME();
#if 0
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;

  EnergyFunction_sp energyFunction = gc::As<EnergyFunction_sp>(score);
  num_real x_rwidth = energyFunction->boundingBox()->get_x_rwidth();
  num_real y_rwidth = energyFunction->boundingBox()->get_y_rwidth();
  num_real z_rwidth = energyFunction->boundingBox()->get_z_rwidth();
  num_real x_width = energyFunction->boundingBox()->get_x_width();
  num_real y_width = energyFunction->boundingBox()->get_y_width();
  num_real z_width = energyFunction->boundingBox()->get_z_width();
  [[maybe_unused]]num_real half_x_size = x_width*0.5;
  [[maybe_unused]]num_real half_y_size = y_width*0.5;
  [[maybe_unused]]num_real half_z_size = z_width*0.5;
  

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

  {
    
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
    int	I1, I2,i;
    gctools::Vec0<EnergyNonbond>::iterator nbi;
    for ( i=0,nbi=this->_Terms.begin();
          nbi!=this->_Terms.end(); nbi++,i++ ) {
  double vdwScale = 1.0;
  double eelScale = 1.0;
  double DIELECTRIC = 1.0;
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#undef CUTOFF_SQUARED
#undef DIELECTRIC
    DONE:
      int index = i;
#define ENERGY_FUNCTION score
#include <cando/chem/energy_functions/_Nonbond_debugFiniteDifference.cc>
#undef ENERGY_FUNCTION

    }
  }
#endif
}

SYMBOL_EXPORT_SC_(KeywordPkg,nonbond);
SYMBOL_EXPORT_SC_(KeywordPkg,vdwRadii);
SYMBOL_EXPORT_SC_(KeywordPkg,distance);
SYMBOL_EXPORT_SC_(KeywordPkg,force);

core::List_sp	EnergyPeriodicBoundaryConditionsNonbond_O::checkForBeyondThresholdInteractionsWithPosition(ScoringFunction_sp score, NVector_sp pos, double threshold )
{
  IMPLEMENT_ME();
#if 0
  EnergyFunction_sp energyFunction = gc::As<EnergyFunction_sp>(score);
  num_real x_rwidth = energyFunction->boundingBox()->get_x_rwidth();
  num_real y_rwidth = energyFunction->boundingBox()->get_y_rwidth();
  num_real z_rwidth = energyFunction->boundingBox()->get_z_rwidth();
  num_real x_width = energyFunction->boundingBox()->get_x_width();
  num_real y_width = energyFunction->boundingBox()->get_y_width();
  num_real z_width = energyFunction->boundingBox()->get_z_width();
  [[maybe_unused]]num_real half_x_size = x_width*0.5;
  [[maybe_unused]]num_real half_y_size = y_width*0.5;
  [[maybe_unused]]num_real half_z_size = z_width*0.5;
  ql::list result;
  bool calcForce = true;
//  printf("%s:%d:%s   number of entries pos -> %lu\n", __FILE__, __LINE__, __FUNCTION__, pos->length()/3);
  
//  printf("%s:%d In evaluateUsingExcludedAtoms starting this->_DebugEnergy -> %d\n", __FILE__, __LINE__, this->_DebugEnergy );
  if (!this->_iac_vec) {
    SIMPLE_ERROR("The nonbonded excluded atoms parameters have not been set up");
  }
  core::SimpleVector_int32_t_sp numberOfExcludedAtoms = this->_NumberOfExcludedAtomIndices;
  core::SimpleVector_int32_t_sp excludedAtomIndices = this->_ExcludedAtomIndices;
  [[maybe_unused]]num_real vdwScale = this->getVdwScale();
  num_real electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
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
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e)-vljrc;}
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
  LOG("Nonbond component is enabled" );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
  // printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  num_real x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2; //,dA_old,dC_old,dQ1Q2_old;
  int	I1, I2;
  int i = 0;
  int endIndex = pos->length()/3;
//  int endIndex = this->_AtomTable->getNumberOfAtoms();
  int excludedAtomIndex = 0;
  int nlocaltype = 0;
  // Find the max local type
  for (i=0; i<this->_iac_vec->length(); ++i){
    if (nlocaltype < (*this->_iac_vec)[i]){
      nlocaltype = (*this->_iac_vec)[i];
    }
  }
  int index1_end = endIndex-1;
  for ( int index1 = 0; index1<index1_end; ++index1 ) {
    LOG("{} ====== top of outer loop - index1 = {}\n" , __FUNCTION__ , index1 );
          // Skip 0 in excluded atom list that amber requires
    bool has_excluded_atoms = ((*excludedAtomIndices)[excludedAtomIndex] >= 0);
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    num_real charge11 = (*this->_charge_vector)[index1];
    [[maybe_unused]]num_real electrostatic_scaled_charge11 = charge11*electrostaticScale;
    for ( int index2 = index1+1, index2_end(endIndex); index2 < index2_end; ++index2 ) {
      LOG("    --- top of inner loop   numberOfExcludedAtomsRemaining -> {}    index2 -> {}\n" , numberOfExcludedAtomsRemaining , index2 );
      int maybe_excluded_atom = (*excludedAtomIndices)[excludedAtomIndex];
      if (numberOfExcludedAtomsRemaining>0 && maybe_excluded_atom == index2) {
        LOG("    Excluding atom {}\n" , index2);
        ++excludedAtomIndex;
        --numberOfExcludedAtomsRemaining;
        continue;
      }
      //   cn1i    (*this->_cn1_vec)[cn1i]
      // printf("%s:%d:%s Interaction between atom %d = %d\n", __FILE__, __LINE__, __FUNCTION__, index1, index2);
      int localindex1 = (*this->_iac_vec)[index1];
      int localindex2 = (*this->_iac_vec)[index2];
      dA = (*this->_cn1_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
      dC = (*this->_cn2_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
      num_real vljrc = dA*this->_NonbondCutoffReciprocal12-dC*this->_NonbondCutoffReciprocal6;
//      printf("%s:%d localindex1 %d and localindex2 %d\n", __FILE__, __LINE__, localindex1, localindex2);
//      printf("%s:%d dA     %lf and dC     %lf\n", __FILE__, __LINE__, dA, dC);
      num_real charge22 = (*this->_charge_vector)[index2];
//      dQ1Q2 = electrostatic_scaled_charge11*charge22;
      dQ1Q2 = charge11*charge22;
//      printf("%s:%d charge1     %lf and charge2     %lf\n", __FILE__, __LINE__, charge11, charge22);
//      printf("%s:%d electrostaticScale     %lf and dQ1Q2     %lf\n", __FILE__, __LINE__, electrostaticScale, dQ1Q2);
#ifdef DEBUG_NONBOND_TERM
      if ( this->_DebugEnergy ) {
        LOG_ENERGY(( "nonbond localindex1 %d\n") , localindex1 );
        LOG_ENERGY(( "nonbond localindex2 %d\n") , localindex2 );
        LOG_ENERGY(( "nonbond dA %5.3lf\n") , dA );
        LOG_ENERGY(( "nonbond dC %5.3lf\n") , dC );
        LOG_ENERGY(( "nonbond dQ1Q2 %5.3lf\n") , dQ1Q2 );
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
        LOG_ENERGY(( "nonbond I1 %d\n") , I1 );
        LOG_ENERGY(( "nonbond I2 %d\n") , I2 );
      }
#endif
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfNonbondTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfNonbondTermsToCalculate ) {
          break;
        }
      }
#endif
  double vdwScale = 1.0;
  double eelScale = 1.0;
  double DIELECTRIC = 1.0;
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#undef CUTOFF_SQUARED
#undef DIELECTRIC
    DONE:
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
    if (!has_excluded_atoms) {
      // No excluded atoms for the current index2 so increment excludedAtomIndex because inner loop didn't do it.
      // see Swails document http://ambermd.org/prmtop.pdf
      ++excludedAtomIndex;
    }

  }
  return result.result();
#endif
}

};

