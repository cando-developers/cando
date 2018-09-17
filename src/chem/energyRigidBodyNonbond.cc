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
//#define	DEBUG_LEVEL_FULL

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

#if 0
#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
#define LOG_ENERGY BFORMAT_T
#endif

namespace chem
{

struct TypeParameters {
  double _Radius;
  double _Epsilon;
  TypeParameters(double radius, double epsilon) : _Radius(radius), _Epsilon(epsilon) {};
};

bool operator<(const TypeParameters& l, const TypeParameters& r) {
  return (l._Radius<r._Radius || (l._Radius==r._Radius && l._Epsilon<r._Epsilon));
}

CL_DEFMETHOD void EnergyRigidBodyNonbond_O::initializeCrossTerms(bool verbose)
{
  typedef std::map<TypeParameters,size_t> double_pair_map;
  double_pair_map types;
  std::vector<TypeParameters> typeParameters;
  size_t next_type_index = 0;
  for ( size_t i=0; i<this->_AtomInfoTable.size(); ++i ) {
    RigidBodyAtomInfo& info = this->_AtomInfoTable[i];
    TypeParameters key(info._Radius,info._Epsilon);
    double_pair_map::iterator it = types.find(key);
    if (it == types.end()) {
      types[key] = next_type_index;
      this->_AtomInfoTable[i]._TypeIndex = next_type_index;
      typeParameters.push_back(key);
      next_type_index++;
    } else {
      this->_AtomInfoTable[i]._TypeIndex = it->second;
    }
  }
  this->_NumberOfTypes = next_type_index;
  if (verbose) {
    BFORMAT_T(BF("There are %d types\n") % this->_NumberOfTypes);
  }
  this->_CrossTerms.resize(next_type_index*next_type_index);
  double vdwScale = this->getVdwScale();
  for ( size_t xi=0; xi<next_type_index; xi++ ) {
    for (size_t yi=xi; yi<next_type_index; yi++ ) {
      TypeParameters& ea1 = typeParameters[xi];
      TypeParameters& ea2 = typeParameters[yi];
      double rStar = ea1._Radius+ea2._Radius;
      double epsilonij = sqrt(ea1._Epsilon*ea2._Epsilon);
      double rStar2 = rStar*rStar;
      double rStar6 = rStar2*rStar2*rStar2;
      double rStar12 = rStar6*rStar6;
      double dA = epsilonij*rStar12*vdwScale;
      double dC = 2.0*epsilonij*rStar6*vdwScale;
      RigidBodyNonbondCrossTerm crossTerm(dA,dC);
      if (verbose) {
        BFORMAT_T(BF("Cross term for types %d - %d\n") % xi % yi );
        BFORMAT_T(BF("    type: %d   radius-> %f   epsilon-> %f\n") % xi % ea1._Radius % ea1._Epsilon);
        BFORMAT_T(BF("    type: %d   radius-> %f   epsilon-> %f\n") % yi % ea2._Radius % ea2._Epsilon);
        BFORMAT_T(BF("    dA -> %f      dC -> %f\n") % dA % dC );
      }
      this->crossTerm(xi,yi) = crossTerm;
      this->crossTerm(yi,xi) = crossTerm;
    }
  }
}
          
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

 void EnergyRigidBodyNonbond_O::energyRigidBodyNonbondSetTerm(gc::Fixnum index, core::T_sp atom, double radius, double epsilon, double charge, const Vector3& position) {
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
      BFORMAT_T(BF("I1 = %3d  %s  %lf %lf %lf %s\n") % (rb*7) % _rep_(ai._Object) % ai._Radius % ai._Epsilon % ai._Charge % ai._Position.asString());
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
  for ( size_t iI1 = 0; iI1<this->_RigidBodyEndAtom->length(); ++iI1 ) {
    ql::list helix;
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
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
  if (this->_CrossTerms.size() == 0 ) this->initializeCrossTerms(false);
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
  size_t interactions = 0;
  size_t I1start = 0;
  for ( size_t iHelix1 = 0; iHelix1<(this->_RigidBodyEndAtom->length()-1); ++iHelix1 ) {
    size_t I1end = (*this->_RigidBodyEndAtom)[iHelix1];
    for (size_t iHelix2 = iHelix1+1; iHelix2<this->_RigidBodyEndAtom->length(); ++iHelix2 ) {
      size_t I2start = (*this->_RigidBodyEndAtom)[iHelix2-1];
      size_t I2end = (*this->_RigidBodyEndAtom)[iHelix2];
      for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
        RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
        double charge1 = ea1._Charge;
        double electrostatic_scaled_charge1 = charge1*electrostaticScale;
        for ( size_t I2cur = I2start; I2cur<I2end; ++I2cur ) {
          ++interactions;
          RigidBodyAtomInfo& ea2 = this->_AtomInfoTable[I2cur];
          RigidBodyNonbondCrossTerm& crossTerm = this->crossTerm(ea1._TypeIndex,ea2._TypeIndex);
          dA = crossTerm.dA;
          dC = crossTerm.dC;
#if 0          
          double rStar = ea1._Radius+ea2._Radius;
          double epsilonij = sqrt(ea1._Epsilon*ea2._Epsilon);
          double rStar2 = rStar*rStar;
          double rStar6 = rStar2*rStar2*rStar2;
          double rStar12 = rStar6*rStar6;
          dA = epsilonij*rStar12*vdwScale;
          dC = 2.0*epsilonij*rStar6*vdwScale;
#endif
          double charge2 = ea2._Charge;
          dQ1Q2 = electrostatic_scaled_charge1*charge2;
          I1 = iHelix1*7;
          I2 = iHelix2*7;
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
            key = _rep_(ea1._Object)+"-"+_rep_(ea2._Object);
            LOG_ENERGY(BF( "MEISTER nonbond %s args cando\n")% key );
            LOG_ENERGY(BF( "MEISTER nonbond %s iHelix1->%3d iHelix2->%3d  I1Cur->%3d I2Cur->%3d\n")% key % iHelix1 % iHelix2 % I1cur % I2cur );
            LOG_ENERGY(BF( "MEISTER nonbond %s dA %5.3lf\n")% key % dA );
            LOG_ENERGY(BF( "MEISTER nonbond %s dC %5.3lf\n")% key % dC );
            LOG_ENERGY(BF( "MEISTER nonbond %s dQ1Q2 %5.3lf\n")% key % dQ1Q2 );
#if 0
            LOG_ENERGY(BF( "MEISTER nonbond %s x1 %5.3lf %d\n")% key % x1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s y1 %5.3lf %d\n")% key % y1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s z1 %5.3lf %d\n")% key % z1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s x2 %5.3lf %d\n")% key % x2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s y2 %5.3lf %d\n")% key % y2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s z2 %5.3lf %d\n")% key % z2 % (I2/3+1) );
#endif
            LOG_ENERGY(BF( "MEISTER nonbond %s results\n")% key );
            LOG_ENERGY(BF( "MEISTER nonbond %s EnergyVdw %lf\n")% key % EnergyVdw);
            LOG_ENERGY(BF( "MEISTER nonbond %s EnergyElectrostatic %lf\n")% key % EnergyElectrostatic);
            LOG_ENERGY(BF( "MEISTER nonbond %s Enonbond(EnergyVdw+EnergyElectrostatic) %lf\n")% key % (EnergyVdw+EnergyElectrostatic) );
            if ( calcForce ) {
              LOG_ENERGY(BF( "MEISTER nonbond %s fam %lf %d\n")% key % fam % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fbm %lf %d\n")% key % fbm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fcm %lf %d\n")% key % fcm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fdm %lf %d\n")% key % fdm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fxm %lf %d\n")% key % fxm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fym %lf %d\n")% key % fym % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fzm %lf %d\n")% key % fzm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fan %lf %d\n")% key % fan % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fbn %lf %d\n")% key % fbn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fcn %lf %d\n")% key % fcn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fdn %lf %d\n")% key % fdn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fxn %lf %d\n")% key % fxn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fyn %lf %d\n")% key % fyn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fzn %lf %d\n")% key % fzn % (I2cur/3+1) );
            }
            LOG_ENERGY(BF( "MEISTER nonbond %s stop\n")% key );
          }
#endif
        }
      }
    }
    I1start = I1end;
  }
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy }\n"));
  this->_TotalEnergy = this->_EnergyVdw+this->_EnergyElectrostatic;
#ifdef DEBUG_NONBOND_TERM
  if ( this->_DebugEnergy ) {
    LOG_ENERGY(BF( "MEISTER nonbond interactions -> %d\n") % interactions);
  }
#endif
  
#if 0
  printf("%s:%d NonbondRigidBody    EnergyVdw -> %lf\n", __FILE__, __LINE__, this->_EnergyVdw );
  printf("%s:%d NonbondRigidBody    EnergyElectrostatic -> %lf\n", __FILE__, __LINE__, this->_EnergyElectrostatic );
  printf("%s:%d NonbondRigidBody    TotalEnergy -> %lf\n", __FILE__, __LINE__, this->_TotalEnergy );
#endif
  return this->_TotalEnergy;
}
    
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
