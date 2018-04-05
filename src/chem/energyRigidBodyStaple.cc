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
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/bformat.h>
#include <cando/chem/nVector.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


namespace chem {

void EnergyRigidBodyStaple_O::energy_rigid_body_staple_add_term( double ks, double r0, size_t rba, const Vector3& pointa, size_t rbb, const Vector3& pointb)
{
  LOG(BF("Defining EnergyRigidBodyStaple with ks=%lf r0=%lf rba=%lu pointa=%s   rbb=%lu  pointb=%s\n")
      % ks
      % r0
      % rba
      % pointa.asString()
      % rbb
      % pointb.asString());
  this->_Terms.emplace_back(ks,r0,rba,pointa,rbb,pointb);
}


#if 0
double	_evaluateEnergyOnly_STAPLE (double ak,
                                    double bk,
                                    double ck,
                                    double dk,
                                    double xk,
                                    double yk,
                                    double zk,
                                    double al,
                                    double bl,
                                    double cl,
                                    double dl,
                                    double xl,
                                    double yl,
                                    double zl,
                                    double ks )
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
  double fx1 = 0.0; 
  double fy1 = 0.0; 
  double fz1 = 0.0;
  double fx2 = 0.0;
  double fy2 = 0.0;
  double fz2 = 0.0;
#endif
#include <cando/energy-functions/_STAPLE_termCode.cc>

  return Energy;
}
#endif




void	EnergyRigidBodyStaple_O::setupHessianPreconditioner(
                                                            NVector_sp nvPosition,
                                                            AbstractLargeSquareMatrix_sp m )
{
  bool		calcForce = true;
  bool		calcDiagonalHessian = true;
  bool		calcOffDiagonalHessian = true;
  double        ks;
  size_t        rba, rbb;
  double        pxk, pyk, pzk;
  double        pxl, pyl, pzl;
  double        ak, bk, ck, dk, xk, yk, zk;
  double        al, bl, cl, dl, xl, yl, zl;
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

  if ( this->isEnabled() ) {
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


  
double EnergyRigidBodyStaple_O::evaluateAll(
                                          NVector_sp 	pos,
                                          bool 		calcForce,
                                          gc::Nilable<NVector_sp> 	force,
                                          bool		calcDiagonalHessian,
                                          bool		calcOffDiagonalHessian,
                                          gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                          gc::Nilable<NVector_sp>	hdvec,
                                          gc::Nilable<NVector_sp> dvec)
{ 
  if ( this->_DebugEnergy ) 
  {
    //BFORMAT_T_CLEAR();
    BFORMAT_T(BF("%s {\n")% this->className());
  }

  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
  double        ks, r0;
  size_t        rba, rbb;
  double        pxk, pyk, pzk;
  double        pxl, pyl, pzl;
  double        ak, bk, ck, dk, xk, yk, zk;
  double        al, bl, cl, dl, xl, yl, zl;
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

  if ( this->isEnabled() ) 
  { _BLOCK_TRACE("Stretch");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_STAPLE_termDeclares.cc>
#pragma clang diagnostic pop
    fxk = 0.0; fyk = 0.0; fzk = 0.0;
    fxl = 0.0; fyl = 0.0; fzl = 0.0;
//    double x1,y1,z1,x2,y2,z2,kb,r0;
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
      if ( this->_DebugEnergy ) {
        BFORMAT_T(BF( "MEISTER staple %d args cando\n")% (i+1) );
        BFORMAT_T(BF( "MEISTER staple %d address %x \n")% (i+1) % (&(*si)) );
        BFORMAT_T(BF( "MEISTER staple %d r0  %5.3lf\n")% (i+1) % r0 );
        BFORMAT_T(BF( "MEISTER staple %d ks  %5.1lf\n")% (i+1) % ks );
        BFORMAT_T(BF( "MEISTER staple %d I1  %4d\n")% (i+1) % I1 );
        BFORMAT_T(BF( "MEISTER staple %d I2  %4d\n")% (i+1) % I2 );
        BFORMAT_T(BF( "MEISTER staple %d ak  %5.3lf %d\n")% (i+1) % ak % (I1/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d bk  %5.3lf %d\n")% (i+1) % bk % (I1/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d ck  %5.3lf %d\n")% (i+1) % ck % (I1/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d dk  %5.3lf %d\n")% (i+1) % dk % (I1/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d xk  %5.3lf %d\n")% (i+1) % xk % (I1/7+1) ); 
        BFORMAT_T(BF( "MEISTER staple %d yk  %5.3lf %d\n")% (i+1) % yk % (I1/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d zk  %5.3lf %d\n")% (i+1) % zk % (I1/7+1) ); 
        BFORMAT_T(BF( "MEISTER staple %d pxk %5.3lf %d\n")% (i+1) % pxk % (I1/7+1) ); 
        BFORMAT_T(BF( "MEISTER staple %d pyk %5.3lf %d\n")% (i+1) % pyk % (I1/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d pzk %5.3lf %d\n")% (i+1) % pzk % (I1/7+1) ); 
        BFORMAT_T(BF( "MEISTER staple %d al  %5.3lf %d\n")% (i+1) % al % (I2/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d bl  %5.3lf %d\n")% (i+1) % bl % (I2/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d cl  %5.3lf %d\n")% (i+1) % cl % (I2/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d dl  %5.3lf %d\n")% (i+1) % dl % (I2/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d xl  %5.3lf %d\n")% (i+1) % xl % (I2/7+1) ); 
        BFORMAT_T(BF( "MEISTER staple %d yl  %5.3lf %d\n")% (i+1) % yl % (I2/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d zl  %5.3lf %d\n")% (i+1) % zl % (I2/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d pxl %5.3lf %d\n")% (i+1) % pxl % (I2/7+1) ); 
        BFORMAT_T(BF( "MEISTER staple %d pyl %5.3lf %d\n")% (i+1) % pyl % (I2/7+1) );
        BFORMAT_T(BF( "MEISTER staple %d pzl %5.3lf %d\n")% (i+1) % pzl % (I2/7+1) ); 
        BFORMAT_T(BF( "MEISTER staple %d results\n")% (i+1) );
        BFORMAT_T(BF( "MEISTER staple %d Energy %lf\n")% (i+1) % Energy );
        if ( calcForce ) {
          BFORMAT_T(BF( "MEISTER staple %d fak  %5.3lf %d\n")% (i+1) % fak % (I1/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fbk  %5.3lf %d\n")% (i+1) % fbk % (I1/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fck  %5.3lf %d\n")% (i+1) % fck % (I1/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fdk  %5.3lf %d\n")% (i+1) % fdk % (I1/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fxk  %5.3lf %d\n")% (i+1) % fxk % (I1/7+1) ); 
          BFORMAT_T(BF( "MEISTER staple %d fyk  %5.3lf %d\n")% (i+1) % fyk % (I1/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fzk  %5.3lf %d\n")% (i+1) % fzk % (I1/7+1) ); 
          BFORMAT_T(BF( "MEISTER staple %d fal  %5.3lf %d\n")% (i+1) % fal % (I2/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fbl  %5.3lf %d\n")% (i+1) % fbl % (I2/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fcl  %5.3lf %d\n")% (i+1) % fcl % (I2/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fdl  %5.3lf %d\n")% (i+1) % fdl % (I2/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fxl  %5.3lf %d\n")% (i+1) % fxl % (I2/7+1) ); 
          BFORMAT_T(BF( "MEISTER staple %d fyl  %5.3lf %d\n")% (i+1) % fyl % (I2/7+1) );
          BFORMAT_T(BF( "MEISTER staple %d fzl  %5.3lf %d\n")% (i+1) % fzl % (I2/7+1) );
        }
        BFORMAT_T(BF( "MEISTER staple %d stop\n")% (i+1) );
      }
		/* Add the forces */

//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz2>10000.0);
    }
  }
  if ( this->_DebugEnergy ) 
  {
    BFORMAT_T(BF("%s }")% this->className());
  }
  return totalEnergy;
}






#if 0
// THis is going to take some more work to implement
void	EnergyRigidBodyStaple_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm( NVector_sp 	pos)
{_OF();
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;
  double        ks;
  size_t        rba, rbb;
  double        xh1, yh1, zh1;
  double        xh2, yh2, zh2;
  double        ak, bk, ck, dk, xk, yk, zk;
  double        al, bl, cl, dl, xl, yl, zl;
  
  {

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
#define	STAPLE_ENERGY_ACCUMULATE(e) {}
#undef	STAPLE_FORCE_ACCUMULATE
#define	STAPLE_FORCE_ACCUMULATE(i,o,v) {}
#undef	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

    if ( this->isEnabled() ) {
      _BLOCK_TRACE("StapleEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_STAPLE_termDeclares.cc>
#pragma clang diagnostic pop
      fxk = 0.0; fyk = 0.0; fzk = 0.0;
      fxl = 0.0; fyl = 0.0; fzl = 0.0;
//      double x1,y1,z1,x2,y2,z2,kb,r0;
      int I1, I2,i;
//	    stretchScale = this->getScale();
      gctools::Vec0<EnergyRigidBodyStaple>::iterator si;
      for ( i=0,si=this->_Terms.begin();
            si!=this->_Terms.end(); si++,i++ ) {
		    // ************* Evaluate the stretch energy/force/hessian
		    // using code generated by Mathematica ***************
#include <cando/energy-functions/_STAPLE_termCode.cc>
        LOG(BF("Energy = %le") % Energy  );
        int index = i;
        LOG(BF("x1 = %le") % x1  );
        LOG(BF("y1 = %le") % y1  );
        LOG(BF("z1 = %le") % z1  );
        LOG(BF("x2 = %le") % x2  );
        LOG(BF("y2 = %le") % y2  );
        LOG(BF("z2 = %le") % z2  );
#include <cando/energy-functions/_STAPLE_debugFiniteDifference.cc>
      }
    }
  }
}
#endif

#if 0
int	EnergyRigidBodyStaple_O::checkForBeyondThresholdInteractions(
                                                                     stringstream& info, NVector_sp pos )
{_OF();
  int	fails = 0;
//	bool	calcForce = false;
//	bool	calcDiagonalHessian = false;
//	bool	calcOffDiagonalHessian = false;

  this->_BeyondThresholdTerms.clear();

#undef STAPLE_CALC_FORCE
#undef STAPLE_CALC_DIAGONAL_HESSIAN
#undef STAPLE_CALC_OFF_DIAGONAL_HESSIAN
#undef	STAPLE_SET_PARAMETER
#define	STAPLE_SET_PARAMETER(x,y)	{x = si->y;}
#undef STAPLE_SET_POINT
#define	STAPLE_SET_POINT(x,p,i)         {x = si->p.i;}
#undef	STAPLE_SET_POSITION
#define	STAPLE_SET_POSITION(x,ii,of)	{x = pos->getElement(ii+of);}
#undef	STAPLE_ENERGY_ACCUMULATE
#define	STAPLE_ENERGY_ACCUMULATE(e) {}
#undef	STAPLE_FORCE_ACCUMULATE
#define	STAPLE_FORCE_ACCUMULATE(i,o,v) {}
#undef	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	STAPLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

  if ( this->isEnabled() ) {
    _BLOCK_TRACE("StretchEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_STAPLE_termDeclares.cc>
#pragma clang diagnostic pop
    fx1 = 0.0; fy1 = 0.0; fz1 = 0.0;
    fx2 = 0.0; fy2 = 0.0; fz2 = 0.0;
    double x1,y1,z1,x2,y2,z2,r0;
    int I1, I2,i;
//	    stretchScale = this->getScale();
    gctools::Vec0<EnergyRigidBodyStaple>::iterator si;
    for ( i=0,si=this->_Terms.begin();
          si!=this->_Terms.end(); si++,i++ ) {
		// ************* Evaluate the stretch energy/force/hessian
		// using code generated by Mathematica ***************
#include <cando/energy-functions/_STAPLE_termCode.cc>
      ASSERT(r0 != 0.0);
      if ( fabs(StretchDeviation)/r0 > this->_ErrorThreshold ) {
        Atom_sp a1, a2;
        a1 = (*si)._Atom1;
        a2 = (*si)._Atom2;
        info<< "StretchDeviation ";
        info<< "value " << fabs(StretchDeviation)/r0 << " ";
        info<<"threshold " << this->_ErrorThreshold << " Atoms(";
        info << a1->getName() << " ";
        info << a2->getName() << ")";
        info << std::endl;
        this->_BeyondThresholdTerms.push_back(*si);
        fails++;
      }
    }
  }
  return fails;
}

#endif




void	EnergyRigidBodyStaple_O::dumpTerms()
{_OF();
  gctools::Vec0<EnergyRigidBodyStaple>::iterator	esi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4;
  LOG(BF("Dumping EnergyRigidBodyStaple terms"));
  uint idx;
  for ( idx = 0, esi=this->_Terms.begin(); esi!=this->_Terms.end(); esi++, idx++ )
  {
    _lisp->print(BF("TERM 1RBST %-8.2lf %-8.2lf %-9lu %-8.2lf %-8.2lf %-8.2lf - %-9lu %-8.2lf %-8.2lf %-8.2lf") 
                 % esi->ks
                 % esi->r0
                 % esi->rigidBodyK
                 % esi->pointK.getX()
                 % esi->pointK.getY()
                 % esi->pointK.getZ()
                 % esi->rigidBodyL
                 % esi->pointL.getX()
                 % esi->pointL.getY()
                 % esi->pointL.getZ() );
  }
}


void EnergyRigidBodyStaple_O::initialize()
{
  this->Base::initialize();
  this->setErrorThreshold(0.05);
}


void EnergyRigidBodyStaple_O::addTerm(const EnergyRigidBodyStaple& term)
{
  this->_Terms.push_back(term);
}


};
