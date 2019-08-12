/*
    File: energyPointToLineRestraint.cc
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

#include <cando/chem/energyPointToLineRestraint.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <clasp/core/wrappers.h>

namespace chem {

    EnergyPointToLineRestraint::EnergyPointToLineRestraint()
{
  this->_Atom1 = _Nil<Atom_O>();
}

EnergyPointToLineRestraint::~EnergyPointToLineRestraint()
{
}


//
// Copy this from implementAmberFunction.cc
//
double	_evaluateEnergyOnly_PointToLineRestraint(
		double x1, double y1, double z1,
		double xa, double ya, double za,
		double xb, double yb, double zb,
		double ka, double ra)
{
#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of)	{}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANCHOR_RESTRAINT_CALC_FORCE	// Don't calculate FORCE or HESSIAN


#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_PointToLineRestraint_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_PointToLineRestraint_termCode.cc>

    return Energy;
}


void	EnergyPointToLineRestraint_O::addTerm(const EnergyPointToLineRestraint& r)
{
    this->_Terms.push_back(r);
}

void	EnergyPointToLineRestraint_O::dumpTerms()
{
}



string EnergyPointToLineRestraint_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyPointToLineRestraint_O,EnergyPointToLineRestraint>(*this);
}







void	EnergyPointToLineRestraint_O::setupHessianPreconditioner(
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

#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{x=cri->term.x;}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of) {x=nvPosition->element(ii+of);}
#undef	ANCHOR_RESTRAINT_PHI_SET
#define	ANCHOR_RESTRAINT_PHI_SET(x) {}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define ANCHOR_RESTRAINT_CALC_FORCE
#define ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN

    if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_PointToLineRestraint_termDeclares.cc>
#pragma clang diagnostic pop
	;
	double x1,y1,z1,xa,ya,za,xb,yb,zb,ka,ra;
	int	I1;
	for ( gctools::Vec0<EnergyPointToLineRestraint>::iterator cri=this->_Terms.begin();
		    cri!=this->_Terms.end(); cri++ ) {
#include <cando/chem/energy_functions/_PointToLineRestraint_termCode.cc>
	}
    }





}





double EnergyPointToLineRestraint_O::evaluateAll( NVector_sp 	pos,
                                             bool 		calcForce,
                                             gc::Nilable<NVector_sp> 	force,
                                             bool		calcDiagonalHessian,
                                             bool		calcOffDiagonalHessian,
                                             gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                             gc::Nilable<NVector_sp>	hdvec,
                                             gc::Nilable<NVector_sp> dvec)
{
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
//
// Copy from implementAmberFunction::evaluateAll
//
// -----------------------

#define ANCHOR_RESTRAINT_CALC_FORCE
#define ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	ANCHOR_RESTRAINT_PHI_SET
#define	ANCHOR_RESTRAINT_PHI_SET(x) {}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) this->_TotalEnergy += (e);
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE 		ForceAcc
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc


  if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_PointToLineRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    double x1,y1,z1,xa,ya,za,xb,yb,zb,ka,ra;
    int	I1, i;
    gctools::Vec0<EnergyPointToLineRestraint>::iterator cri;
    for ( i=0,cri=this->_Terms.begin();
          cri!=this->_Terms.end(); cri++,i++ ) {
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfPointToLineRestraintTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfPointToLineRestraintTermsToCalculate ) {
          break;
        }
      }
#endif

			    /* Obtain all the parameters necessary to calculate */
			    /* the amber and forces */
#include <cando/chem/energy_functions/_PointToLineRestraint_termCode.cc>

#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
      cri->_calcForce = calcForce;
      cri->_calcDiagonalHessian = calcDiagonalHessian;
      cri ->_calcOffDiagonalHessian = calcOffDiagonalHessian;
			// Now write all of the calculated values into the eval structure
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ cri->eval.var=val;};
#include <cando/chem/energy_functions/_PointToLineRestraint_debugEvalSet.cc>
#endif //]

			    /* Add the forces */

      if ( calcForce ) {
#if !USE_EXPLICIT_DECLARES
        double fx1 = 0.0;
        double fy1 = 0.0;
        double fz1 = 0.0;
#endif
//		    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
//		    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
//		    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
      }
    }
  }
  return this->_TotalEnergy;
}







void	EnergyPointToLineRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
		NVector_sp 	pos)
{
int	fails = 0;
bool	calcForce = true;
bool	calcDiagonalHessian = true;
bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//------------------

#define ANCHOR_RESTRAINT_CALC_FORCE
#define ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	ANCHOR_RESTRAINT_PHI_SET
#define	ANCHOR_RESTRAINT_PHI_SET(x) {}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


	if ( this->isEnabled() ) {
		_BLOCK_TRACE("PointToLineRestraintEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_PointToLineRestraint_termDeclares.cc>
#pragma clang diagnostic pop
                double x1,y1,z1,xa,ya,za,xb,yb,zb,ka,ra;

	    int	I1, i;
	    gctools::Vec0<EnergyPointToLineRestraint>::iterator cri;
	    for ( i=0,cri=this->_Terms.begin();
			cri!=this->_Terms.end(); cri++,i++ ) {
			  /* Obtain all the parameters necessary to calculate */
			  /* the amber and forces */
#include <cando/chem/energy_functions/_PointToLineRestraint_termCode.cc>
		LOG(BF("fx1 = %le") % fx1 );
		LOG(BF("fy1 = %le") % fy1 );
		LOG(BF("fz1 = %le") % fz1 );
		int index = i;
#if !USE_EXPLICIT_DECLARES
		double fx1 = 0.0;
		double fy1 = 0.0;
		double fz1 = 0.0;
		double dhx1x1 = 0.0;
		double ohx1y1 = 0.0;
		double dhy1y1 = 0.0;
		double ohy1z1 = 0.0;
		double dhz1z1 = 0.0;
		double ohx1z1 = 0.0;
#endif
#include <cando/chem/energy_functions/_PointToLineRestraint_debugFiniteDifference.cc>

	    }
	}
    IMPLEMENT_ME(); // must return some sort of integer value
}


#if 0
int	EnergyPointToLineRestraint_O::checkForBeyondThresholdInteractions(
			stringstream& info, NVector_sp pos )
{
int	fails = 0;

    this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------

#undef ANCHOR_RESTRAINT_CALC_FORCE
#undef ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN
#undef ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	ANCHOR_RESTRAINT_PHI_SET
#define	ANCHOR_RESTRAINT_PHI_SET(x) {}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


	if ( this->isEnabled() ) {
		_BLOCK_TRACE("PointToLineRestraintEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_PointToLineRestraint_termDeclares.cc>
#pragma clang diagnostic pop
double x1,y1,z1,xa,ya,za,xb,yb,zb,ka,ra;
	    int	I1, i;
	    gctools::Vec0<EnergyPointToLineRestraint>::iterator cri;
	    for ( i=0,cri=this->_Terms.begin();
			cri!=this->_Terms.end(); cri++,i++ ) {
			  /* Obtain all the parameters necessary to calculate */
			  /* the amber and forces */
#include <cando/chem/energy_functions/_PointToLineRestraint_termCode.cc>
		if ( PointToLineDeviation>this->_ErrorThreshold ) {
		    Atom_sp a1, a2, a3, a4;
		    a1 = (*cri)._Atom1;
		    info<< "PointToLineRestraintDeviation ";
//		    info<< a1->getAbsoluteIdPath() << " ";
		    info<< "value " << PointToLineDeviation << " ";
		    info << a1->getName() << " ";
		    info << std::endl;
		    this->_BeyondThresholdTerms.push_back(*cri);
		    fails++;
		}
	    }
	}

    return fails;
}
#endif










void EnergyPointToLineRestraint_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(0.2);
}


#ifdef XML_ARCHIVE
void EnergyPointToLineRestraint_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    IMPLEMENT_ME();
}
#endif


};
