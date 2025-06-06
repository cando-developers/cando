/*
    File: energyStretch.cc
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

#include <cando/chem/energyRigidBodyStaple.h>
#include <clasp/core/numerics.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/ql.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/bformat.h>
#include <cando/chem/nVector.h>
#include <cando/geom/ovector3.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


namespace chem {

core::List_sp EnergyRigidBodyStaple::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,ks),mk_double_float(this->ks)), 
                                  core::Cons_O::create(INTERN_(kw,r0),mk_double_float(this->r0)),
                                  core::Cons_O::create(INTERN_(kw,rigidBodyK),core::make_fixnum(this->rigidBodyK)),
                                  core::Cons_O::create(INTERN_(kw,pointK), translate::to_object<Vector3>::convert(this->pointK) ),
                                  core::Cons_O::create(INTERN_(kw,rigidBodyL),core::make_fixnum(this->rigidBodyL)),
                                  core::Cons_O::create(INTERN_(kw,pointL), translate::to_object<Vector3>::convert(this->pointL) )
                                  );
}

void EnergyRigidBodyStaple::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergyRigidBodyStaple");
}


void EnergyRigidBodyStaple_O::energy_rigid_body_staple_add_term( double ks, double r0, size_t rba, const Vector3& pointa, size_t rbb, const Vector3& pointb)
{
  LOG("Defining EnergyRigidBodyStaple with ks={} r0={} rba={} pointa={}   rbb={}  pointb={}\n"
      , ks
      , r0
      , rba
      , pointa.asString()
      , rbb
      , pointb.asString());
  this->_Terms.emplace_back(ks,r0,rba,pointa,rbb,pointb);
}


#if 0
num_real	_evaluateEnergyOnly_STAPLE (num_real ak,
                                    num_real bk,
                                    num_real ck,
                                    num_real dk,
                                    num_real xk,
                                    num_real yk,
                                    num_real zk,
                                    num_real al,
                                    num_real bl,
                                    num_real cl,
                                    num_real dl,
                                    num_real xl,
                                    num_real yl,
                                    num_real zl,
                                    num_real ks )
{
#undef	STAPLE_SET_PARAMETER
#define	STAPLE_SET_PARAMETER(x,y)	{}
#undef STAPLE_SET_POINT
#define	STAPLE_SET_POINT(x,p,i)         {x = p.i;}
#undef	STAPLE_SET_POSITION
#define	STAPLE_SET_POSITION(x,ii,of)	{}
#undef	STAPLE_ENERGY_ACCUMULATE
#define	STAPLE_ENERGY_ACCUMULATE(e) {}
#undef	STAPLE_FORCE_ACCUMULATE
#define	STAPLE_FORCE_ACCUMULATE(i,o,v) {}
#undef	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	STAPLE_CALC_FORCE	// Don't calculate FORCE or HESSIAN


#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_STAPLE_termDeclares.cc>
#pragma clang diagnostic pop
#pragma clang diagnostic pop
 
#if !USE_EXPLICIT_DECLARES
  num_real fx1 = 0.0; 
  num_real fy1 = 0.0; 
  num_real fz1 = 0.0;
  num_real fx2 = 0.0;
  num_real fy2 = 0.0;
  num_real fz2 = 0.0;
#endif
#include <cando/energy-functions/_STAPLE_termCode.cc>

  return Energy;
}
#endif




void	EnergyRigidBodyStaple_O::setupHessianPreconditioner(
                                                            NVector_sp nvPosition,
                                                            AbstractLargeSquareMatrix_sp m,
                                                            core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  core::T_sp debugInteractions = nil<core::T_O>();
  MAYBE_SETUP_DEBUG_INTERACTIONS(false);
  bool		calcForce = true;
  num_real        pxk, pyk, pzk;
  num_real        pxl, pyl, pzl;
  num_real        ak, bk, ck, dk, xk, yk, zk;
  num_real        al, bl, cl, dl, xl, yl, zl;
#undef	STAPLE_SET_PARAMETER
#define	STAPLE_SET_PARAMETER(x,y)	{x = si->y;}
#undef STAPLE_SET_POINT
#define	STAPLE_SET_POINT(x,p,i)         {x = si->p.i;}
#undef	STAPLE_SET_POSITION
#define	STAPLE_SET_POSITION(x,ii,of)	{x = nvPosition->getElement(ii+of);}
#undef	STAPLE_ENERGY_ACCUMULATE
#define	STAPLE_ENERGY_ACCUMULATE(e) {}
#undef	STAPLE_FORCE_ACCUMULATE
#define	STAPLE_FORCE_ACCUMULATE(i,o,v) {}
#undef	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {	\
    m->addToElement((i1)+(o1),(i2)+(o2),v);		\
  }
#undef	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {	\
    m->addToElement((i1)+(o1),(i2)+(o2),v);			\
  }
#define STAPLE_CALC_FORCE
#define STAPLE_CALC_DIAGONAL_HESSIAN
#define STAPLE_CALC_OFF_DIAGONAL_HESSIAN

  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_STAPLE_termDeclares.cc>
#pragma clang diagnostic pop
    //    double x1,y1,z1,x2,y2,z2,kxb,r0;
    double ks, r0;
    int I1, I2;
//	stretchScale = this->getScale();
    for ( gctools::Vec0<EnergyRigidBodyStaple>::iterator si=this->_Terms.begin();
          si!=this->_Terms.end(); si++ ) {
#include <cando/energy-functions/_STAPLE_termCode.cc>
    }
  }
}


  
double EnergyRigidBodyStaple_O::evaluateAllComponent( ScoringFunction_sp score,
                                                      NVector_sp 	pos,
                                                      core::T_sp energyScale,
                                                      core::T_sp componentEnergy,
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
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
//  SIMPLE_WARN("FIXactiveAtomMask How do I deal with activeAtomMask");
  this->_Evaluations++;
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  num_real        ks, r0;
  num_real        pxk, pyk, pzk;
  num_real        pxl, pyl, pzl;
  num_real        ak, bk, ck, dk, xk, yk, zk;
  num_real        al, bl, cl, dl, xl, yl, zl;
  double totalEnergy = 0.0;
#define STAPLE_CALC_FORCE
#define STAPLE_CALC_DIAGONAL_HESSIAN
#define STAPLE_CALC_OFF_DIAGONAL_HESSIAN
#undef	STAPLE_SET_PARAMETER
#define	STAPLE_SET_PARAMETER(x,y)	{x = si->y;}
#undef STAPLE_SET_POINT
#define	STAPLE_SET_POINT(x,p,i)         {x = si->p.i;}
#undef	STAPLE_SET_POSITION
#define	STAPLE_SET_POSITION(x,ii,of)	{x = pos->getElement(ii+of);}
#undef	STAPLE_ENERGY_ACCUMULATE
#define	STAPLE_ENERGY_ACCUMULATE(e) {totalEnergy += (e);};
#undef	STAPLE_FORCE_ACCUMULATE
#undef	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE
#undef	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_FORCE_ACCUMULATE 		ForceAcc
#define	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc

  { 
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_STAPLE_termDeclares.cc>
#pragma clang diagnostic pop
    int I1, I2,i;
    gctools::Vec0<EnergyRigidBodyStaple>::iterator si;
    for ( i=0,si=this->_Terms.begin();
          si!=this->_Terms.end(); si++,i++ ) {
#ifdef DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfTermsToCalculate ) {
          break;
        }
      }
#endif

		// ************* Evaluate the stretch energy/force/hessian
		// using code generated by Mathematica ***************
#include <cando/energy-functions/_STAPLE_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
      si->_calcForce = calcForce;
      si->_calcDiagonalHessian = calcDiagonalHessian;
      si->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define EVAL_SET(var,val) {si->eval.var=val;}
#include <cando/energy-functions/_STAPLE_debugEvalSet.cc>
#endif //]
    }
  }
  maybeSetEnergy( componentEnergy, EnergyRigidBodyStaple_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}


void	EnergyRigidBodyStaple_O::dumpTerms(core::HashTable_sp atomTypes)
{
  gctools::Vec0<EnergyRigidBodyStaple>::iterator	esi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4;
  LOG("Dumping EnergyRigidBodyStaple terms");
  uint idx;
  for ( idx = 0, esi=this->_Terms.begin(); esi!=this->_Terms.end(); esi++, idx++ )
  {
    core::clasp_write_string(fmt::format("TERM 1RBST {:<8.2f} {:<8.2f} {:<9} {:<8.2f} {:<8.2f} {:<8.2f} - {:<9} {:<8.2f} {:<8.2f} {:<8.2f}\n" , esi->ks , esi->r0 , esi->rigidBodyK , esi->pointK.getX() , esi->pointK.getY() , esi->pointK.getZ() , esi->rigidBodyL , esi->pointL.getX() , esi->pointL.getY() , esi->pointL.getZ() ));
  }
}


void EnergyRigidBodyStaple_O::initialize()
{
  this->Base::initialize();
}


void EnergyRigidBodyStaple_O::addTerm(const EnergyRigidBodyStaple& term)
{
  this->_Terms.push_back(term);
}


CL_DEFMETHOD core::List_sp EnergyRigidBodyStaple_O::parts_as_list(NVector_sp 	pos) {
  ql::list result;
  num_real        ks, r0;
  num_real        pxk, pyk, pzk;
  num_real        pxl, pyl, pzl;
  num_real        ak, bk, ck, dk, xk, yk, zk;
  num_real        al, bl, cl, dl, xl, yl, zl;
  int I1, I2;
#undef	STAPLE_POSITION_SET_PARAMETER
#define	STAPLE_POSITION_SET_PARAMETER(x,y)	{x = si->y;}
#undef  STAPLE_POSITION_SET_POINT
#define	STAPLE_POSITION_SET_POINT(x,p,i)         {x = si->p.i;}
#undef	STAPLE_POSITION_SET_POSITION
#define	STAPLE_POSITION_SET_POSITION(x,ii,of)	{x = pos->getElement(ii+of);}

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_STAPLE_POSITIONS_termDeclares.cc>
#pragma clang diagnostic pop
  int i;
  gctools::Vec0<EnergyRigidBodyStaple>::iterator si;
  for ( i=0,si=this->_Terms.begin();
        si!=this->_Terms.end(); si++,i++ ) {
		// ************* Evaluate the stretch energy/force/hessian
		// using code generated by Mathematica ***************
#include <cando/energy-functions/_STAPLE_POSITIONS_termCode.cc>
#if 0
    printf("%s:%d  I1 -> %d   I2 -> %d\n", __FILE__, __LINE__, I1, I2);
    printf("      ak -> %lf bk -> %lf ck -> %lf dk -> %lf xk -> %lf yk -> %lf zk -> %lf\n", ak, bk, ck, dk, xk, yk, zk);
    printf("      al -> %lf bl -> %lf cl -> %lf dl -> %lf xl -> %lf yl -> %lf zl -> %lf\n", al, bl, cl, dl, xl, yl, zl);
#endif
    result << core::Cons_O::createList(core::DoubleFloat_O::create(plabkx),
                                       core::DoubleFloat_O::create(plabky),
                                       core::DoubleFloat_O::create(plabkz),
                                       core::DoubleFloat_O::create(plablx),
                                       core::DoubleFloat_O::create(plably),
                                       core::DoubleFloat_O::create(plablz),
                                       core::DoubleFloat_O::create(si->r0));
  }
  return result.cons();
}

size_t EnergyRigidBodyStaple_O::partsCoordinates(NVector_sp 	pos, size_t idx, core::SimpleVector_float_sp coords) {
  double        ks, r0;
  num_real        pxk, pyk, pzk;
  num_real        pxl, pyl, pzl;
  num_real        ak, bk, ck, dk, xk, yk, zk;
  num_real        al, bl, cl, dl, xl, yl, zl;
  int I1, I2;
#undef	STAPLE_POSITION_SET_PARAMETER
#define	STAPLE_POSITION_SET_PARAMETER(x,y)	{x = si->y;}
#undef  STAPLE_POSITION_SET_POINT
#define	STAPLE_POSITION_SET_POINT(x,p,i)         {x = si->p.i;}
#undef	STAPLE_POSITION_SET_POSITION
#define	STAPLE_POSITION_SET_POSITION(x,ii,of)	{x = pos->getElement(ii+of);}

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_STAPLE_POSITIONS_termDeclares.cc>
#pragma clang diagnostic pop
  int i;
  gctools::Vec0<EnergyRigidBodyStaple>::iterator si;
  for ( i=0,si=this->_Terms.begin();
        si!=this->_Terms.end(); si++,i++ ) {
		// ************* Evaluate the stretch energy/force/hessian
		// using code generated by Mathematica ***************
#include <cando/energy-functions/_STAPLE_POSITIONS_termCode.cc>
#if 0
    printf("%s:%d  I1 -> %d   I2 -> %d\n", __FILE__, __LINE__, I1, I2);
    printf("      ak -> %lf bk -> %lf ck -> %lf dk -> %lf xk -> %lf yk -> %lf zk -> %lf\n", ak, bk, ck, dk, xk, yk, zk);
    printf("      al -> %lf bl -> %lf cl -> %lf dl -> %lf xl -> %lf yl -> %lf zl -> %lf\n", al, bl, cl, dl, xl, yl, zl);
#endif
    (*coords)[idx++] = plabkx;
    (*coords)[idx++] = plabky;
    (*coords)[idx++] = plabkz;
    (*coords)[idx++] = plablx;
    (*coords)[idx++] = plably;
    (*coords)[idx++] = plablz;
  }
  return idx;
}


void EnergyRigidBodyStaple_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms) ,this->_Terms);
};

};
