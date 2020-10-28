/*
    File: energyOutOfZPlane.cc
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

#include <cando/chem/energyOutOfZPlane.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <clasp/core/wrappers.h>

namespace chem {


EnergyOutOfZPlane::~EnergyOutOfZPlane()
{
}

#ifdef XML_ARCHIVE
    void	EnergyOutOfZPlane::archive(core::ArchiveP node)
{
    node->attribute("kb",this->term.kb);
    node->attribute("za",this->term.za);
    node->attribute("I1",this->term.I1);
    node->attribute("a1",this->_Atom1);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
    node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
    node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
    node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include <cando/chem/energy_functions/_Oozp_debugEvalSerialize.cc>
#endif //]
}
#endif

//
// Copy this from implementAmberFunction.cc
//
double	_evaluateEnergyOnly_Oozp(
		double x1, double y1, double z1,
		double za,
		double kb)
{
#undef	OOZP_SET_PARAMETER
#define	OOZP_SET_PARAMETER(x)	{}
#undef	OOZP_SET_POSITION
#define	OOZP_SET_POSITION(x,ii,of)	{}
#undef	OOZP_ENERGY_ACCUMULATE
#define	OOZP_ENERGY_ACCUMULATE(e) {}
#undef	OOZP_FORCE_ACCUMULATE
#define	OOZP_FORCE_ACCUMULATE(i,o,v) {}
#undef	OOZP_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	OOZP_CALC_FORCE	// Don't calculate FORCE or HESSIAN


#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Oozp_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_Oozp_termCode.cc>

    return Energy;
}


void	EnergyOutOfZPlane_O::addTerm(const EnergyOutOfZPlane& r)
{
    this->_Terms.push_back(r);
}



core::List_sp EnergyOutOfZPlane::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kb),core::clasp_make_double_float(this->term.kb)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1)),
                                  core::Cons_O::create(INTERN_(kw,za), core::make_fixnum(this->term.za)));
}

CL_DEFMETHOD
void EnergyOutOfZPlane_O::addOutOfZPlaneTerm(size_t index, double kb, double r0) {
  EnergyOutOfZPlane term(kb,index,r0);
  this->addTerm(term);
}

void	EnergyOutOfZPlane_O::dumpTerms()
{
}

void EnergyOutOfZPlane_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}


string EnergyOutOfZPlane_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyOutOfZPlane_O,EnergyOutOfZPlane>(*this);
}







void	EnergyOutOfZPlane_O::setupHessianPreconditioner(
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

#undef	OOZP_SET_PARAMETER
#define	OOZP_SET_PARAMETER(x)	{x=cri->term.x;}
#undef	OOZP_SET_POSITION
#define	OOZP_SET_POSITION(x,ii,of) {x=nvPosition->element(ii+of);}
#undef	OOZP_PHI_SET
#define	OOZP_PHI_SET(x) {}
#undef	OOZP_ENERGY_ACCUMULATE
#define	OOZP_ENERGY_ACCUMULATE(e) {}
#undef	OOZP_FORCE_ACCUMULATE
#define	OOZP_FORCE_ACCUMULATE(i,o,v) {}
#undef	OOZP_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define OOZP_CALC_FORCE
#define OOZP_CALC_DIAGONAL_HESSIAN
#define OOZP_CALC_OFF_DIAGONAL_HESSIAN

    if ( this->isEnabled() ) {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Oozp_termDeclares.cc>
#pragma clang diagnostic pop
	;
	double x1,y1,z1,za,kb;
	int	I1;
	for ( gctools::Vec0<EnergyOutOfZPlane>::iterator cri=this->_Terms.begin();
		    cri!=this->_Terms.end(); cri++ ) {
#include <cando/chem/energy_functions/_Oozp_termCode.cc>
	}
    }





}





double EnergyOutOfZPlane_O::evaluateAllComponent( ScoringFunction_sp score,
                                         NVector_sp 	pos,
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

#define OOZP_CALC_FORCE
#define OOZP_CALC_DIAGONAL_HESSIAN
#define OOZP_CALC_OFF_DIAGONAL_HESSIAN
#undef	OOZP_SET_PARAMETER
#define	OOZP_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	OOZP_SET_POSITION
#define	OOZP_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	OOZP_PHI_SET
#define	OOZP_PHI_SET(x) {}
#undef	OOZP_ENERGY_ACCUMULATE
#define	OOZP_ENERGY_ACCUMULATE(e) this->_TotalEnergy += (e);
#undef	OOZP_FORCE_ACCUMULATE
#undef	OOZP_DIAGONAL_HESSIAN_ACCUMULATE
#undef	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_FORCE_ACCUMULATE 		ForceAcc
#define	OOZP_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Oozp_termDeclares.cc>
#pragma clang diagnostic pop


  if ( this->isEnabled() ) {
    if (chem__verbose(1)) {
      core::write_bf_stream(BF("Evaluating Oozp component %d terms\n") % this->_Terms.size());
    }
    double x1,y1,z1,za,kb;
    int	I1, i;
    gctools::Vec0<EnergyOutOfZPlane>::iterator cri;
    for ( i=0,cri=this->_Terms.begin();
          cri!=this->_Terms.end(); cri++,i++ ) {
			    /* Obtain all the parameters necessary to calculate */
			    /* the amber and forces */
      OOZP_SET_PARAMETER(kb);
      OOZP_SET_PARAMETER(I1);
      OOZP_SET_POSITION(x1,I1,0);
      OOZP_SET_POSITION(y1,I1,1);
      OOZP_SET_POSITION(z1,I1,2);
      if (z1 > 0.0) {
        OOZP_FORCE_ACCUMULATE(I1,2,-kb);
      } else {
        OOZP_FORCE_ACCUMULATE(I1,2,kb);
      }
//#include <cando/chem/energy_functions/_Oozp_termCode.cc>
#if 1
      if (chem__verbose(2)) {
        core::write_bf_stream(BF("I1,kb,za -> %d, %f, %f\n") % I1 % kb % za );
        core::write_bf_stream(BF("x1,y1,z1 = %f, %f, %f\n") % x1 % y1 % z1 );
        core::write_bf_stream(BF("Oopz Energy = %f\n") % Energy );
        core::write_bf_stream(BF("fx1,fy1,fz1 = %f, %f, %f\n") % force->getElement(I1+0) % force->getElement(I1+1) % force->getElement(I1+2));
      }
#endif
    }
  }
  if (chem__verbose(1)) {
    core::write_bf_stream(BF("Oopz Total Energy = %f\n") % this->_TotalEnergy );
  }
  return this->_TotalEnergy;
}






void	EnergyOutOfZPlane_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
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

#define OOZP_CALC_FORCE
#define OOZP_CALC_DIAGONAL_HESSIAN
#define OOZP_CALC_OFF_DIAGONAL_HESSIAN
#undef	OOZP_SET_PARAMETER
#define	OOZP_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	OOZP_SET_POSITION
#define	OOZP_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	OOZP_PHI_SET
#define	OOZP_PHI_SET(x) {}
#undef	OOZP_ENERGY_ACCUMULATE
#define	OOZP_ENERGY_ACCUMULATE(e) {}
#undef	OOZP_FORCE_ACCUMULATE
#define	OOZP_FORCE_ACCUMULATE(i,o,v) {}
#undef	OOZP_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


	if ( this->isEnabled() ) {
		_BLOCK_TRACE("OutOfZPlaneEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Oozp_termDeclares.cc>
#pragma clang diagnostic pop
	    double x1,y1,z1,za,kb;
	    int	I1, i;
	    gctools::Vec0<EnergyOutOfZPlane>::iterator cri;
	    for ( i=0,cri=this->_Terms.begin();
			cri!=this->_Terms.end(); cri++,i++ ) {
			  /* Obtain all the parameters necessary to calculate */
			  /* the amber and forces */
#include <cando/chem/energy_functions/_Oozp_termCode.cc>
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
#include <cando/chem/energy_functions/_Oozp_debugFiniteDifference.cc>

	    }
	}
    IMPLEMENT_ME(); // must return some sort of integer value
}

#if 0
int	EnergyOutOfZPlane_O::checkForBeyondThresholdInteractions(
			stringstream& info, NVector_sp pos )
{
int	fails = 0;

    this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------

#undef OOZP_CALC_FORCE
#undef OOZP_CALC_DIAGONAL_HESSIAN
#undef OOZP_CALC_OFF_DIAGONAL_HESSIAN
#undef	OOZP_SET_PARAMETER
#define	OOZP_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	OOZP_SET_POSITION
#define	OOZP_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	OOZP_PHI_SET
#define	OOZP_PHI_SET(x) {}
#undef	OOZP_ENERGY_ACCUMULATE
#define	OOZP_ENERGY_ACCUMULATE(e) {}
#undef	OOZP_FORCE_ACCUMULATE
#define	OOZP_FORCE_ACCUMULATE(i,o,v) {}
#undef	OOZP_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


	if ( this->isEnabled() ) {
		_BLOCK_TRACE("OutOfZPlaneEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Oozp_termDeclares.cc>
#pragma clang diagnostic pop
	    double x1,y1,z1,za,kb;
	    int	I1, i;
	    gctools::Vec0<EnergyOutOfZPlane>::iterator cri;
	    for ( i=0,cri=this->_Terms.begin();
			cri!=this->_Terms.end(); cri++,i++ ) {
			  /* Obtain all the parameters necessary to calculate */
			  /* the amber and forces */
#include <cando/chem/energy_functions/_Oozp_termCode.cc>
		if ( AnchorDeviation>this->_ErrorThreshold ) {
                  Atom_sp a1;
		    a1 = (*cri)._Atom1;
		    info<< "OutOfZPlaneDeviation ";
//		    info<< a1->getAbsoluteIdPath() << " ";
		    info<< "value " << AnchorDeviation << " ";
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

void EnergyOutOfZPlane_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(0.2);
}

void EnergyOutOfZPlane_O::reset()
{
}


};
