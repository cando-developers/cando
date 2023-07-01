/*
    File: energyDihedralRestraint.cc
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

#include <clasp/core/foundation.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/geom/vector3.fwd.h>
#include <cando/chem/energyDihedralRestraint.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <clasp/core/wrappers.h>


namespace chem {

core::List_sp EnergyDihedralRestraint::encode() const {
  return core::Cons_O::create(core::Cons_O::create(INTERN_(kw,k),core::clasp_make_double_float(this->term.K)),
                              core::Cons_O::create(core::Cons_O::create(INTERN_(kw,U),core::clasp_make_double_float(this->term.U)),
                                                   core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,L),core::clasp_make_double_float(this->term.L)),
                                                                            core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1)),
                                                                            core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.I2)),
                                                                            core::Cons_O::create(INTERN_(kw,i3), core::make_fixnum(this->term.I3)),
                                                                            core::Cons_O::create(INTERN_(kw,i4), core::make_fixnum(this->term.I4)),
                                                                            core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                                                            core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2),
                                                                            core::Cons_O::create(INTERN_(kw,atom3), this->_Atom3),
                                                                            core::Cons_O::create(INTERN_(kw,atom4), this->_Atom4)
                                                                            )));
}


double EnergyDihedralRestraint::getAngle()
{
    Vector3	pos1, pos2, pos3, pos4;
    pos1 = this->_Atom1->getPosition();
    pos2 = this->_Atom2->getPosition();
    pos3 = this->_Atom3->getPosition();
    pos4 = this->_Atom4->getPosition();
    return geom::calculateDihedral(pos1,pos2,pos3,pos4);
}

//
// Copy this from implementAmberFunction.cc
//
num_real	_evaluateEnergyOnly_ImproperRestraint(
		num_real x1, num_real y1, num_real z1,
		num_real x2, num_real y2, num_real z2,
		num_real x3, num_real y3, num_real z3,
		num_real x4, num_real y4, num_real z4,
		num_real K, num_real L, num_real U )
{
num_real	EraseLinearDihedral;
num_real	UShift, PhiShift;
bool	RestraintActive;

#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of)	{}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x)	{}
#undef	IMPROPER_RESTRAINT_CALC_FORCE	// Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ImproperRestraint_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_ImproperRestraint_termCode.cc>

    return Energy;
}


void EnergyDihedralRestraint_O::addTerm(const EnergyDihedralRestraint& e)
{
    this->_Terms.push_back(e);
}


CL_DEFMETHOD
void EnergyDihedralRestraint_O::addDihedralRestraint(EnergyFunction_sp energyFunction, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double minRadians, double maxRadians, double weight)
{
  EnergyDihedralRestraint energyTerm;
  energyTerm._Atom1 = a1;
  energyTerm._Atom2 = a2;
  energyTerm._Atom3 = a3;
  energyTerm._Atom4 = a4;
  EnergyAtom* ea1 = energyFunction->getEnergyAtomPointer(energyTerm._Atom1);
  EnergyAtom* ea2 = energyFunction->getEnergyAtomPointer(energyTerm._Atom2);
  EnergyAtom* ea3 = energyFunction->getEnergyAtomPointer(energyTerm._Atom3);
  EnergyAtom* ea4 = energyFunction->getEnergyAtomPointer(energyTerm._Atom4);
  energyTerm.term.I1 = ea1->coordinateIndexTimes3();
  energyTerm.term.I2 = ea2->coordinateIndexTimes3();
  energyTerm.term.I3 = ea3->coordinateIndexTimes3();
  energyTerm.term.I4 = ea4->coordinateIndexTimes3();
  energyTerm.term.U = minRadians;
  energyTerm.term.L = maxRadians;
  energyTerm.term.K = weight;
  this->addTerm(energyTerm);
}





void	EnergyDihedralRestraint_O::dumpTerms(core::HashTable_sp atomTypes)
{
}

void EnergyDihedralRestraint_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}


string EnergyDihedralRestraint_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyDihedralRestraint_O,EnergyDihedralRestraint>(*this);
}




void	EnergyDihedralRestraint_O::setupHessianPreconditioner(
					chem::NVector_sp nvPosition,
					chem::AbstractLargeSquareMatrix_sp m )
{
bool		calcForce = true;
bool		calcDiagonalHessian = true;
bool		calcOffDiagonalHessian = true;

//
// Copy from implementAmberFunction::setupHessianPreconditioner
//
// -----------------------

#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{x=iri->term.x;}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of) {x=nvPosition->element(ii+of);}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x) {}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define IMPROPER_RESTRAINT_CALC_FORCE
#define IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN

 {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include	<cando/chem/energy_functions/_ImproperRestraint_termDeclares.cc>
#pragma clang diagnostic pop
	num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
	num_real K,U,L;
	num_real EraseLinearDihedral, UShift, PhiShift;
	bool RestraintActive;
	int	I1, I2, I3, I4;
	for ( gctools::Vec0<EnergyDihedralRestraint>::iterator iri=this->_Terms.begin();
		    iri!=this->_Terms.end(); iri++ ) {
#define DEBUG_IMPROPER_RESTRAINT 1
#include	<cando/chem/energy_functions/_ImproperRestraint_termCode.cc>
#undef DEBUG_IMPROPER_RESTRAINT
	}
    }
}





num_real EnergyDihedralRestraint_O::evaluateAllComponent( ScoringFunction_sp score,
                                                        chem::NVector_sp 	pos,
                                                        core::T_sp componentEnergy,
                                                        bool 		calcForce,
                                                        gc::Nilable<chem::NVector_sp> 	force,
                                                        bool		calcDiagonalHessian,
                                                        bool		calcOffDiagonalHessian,
                                                        gc::Nilable<chem::AbstractLargeSquareMatrix_sp>	hessian,
                                                        gc::Nilable<chem::NVector_sp>	hdvec,
                                                        gc::Nilable<chem::NVector_sp> dvec)
{
  num_real totalEnergy = 0.0;
  this->_Evaluations++;
  if ( this->_DebugEnergy ) {
    core::clasp_write_string(fmt::format("{}\n" , __FUNCTION__ ));
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(("%s {") , this->className());
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
#define IMPROPER_RESTRAINT_CALC_FORCE
#define IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{x = iri->term.x;}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x) {}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) totalEnergy += (e);
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE 		ForceAcc
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc

  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ImproperRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    num_real K,U,L;
    num_real EraseLinearDihedral, UShift, PhiShift;
    bool RestraintActive;
    int	I1, I2, I3, I4,i;
    gctools::Vec0<EnergyDihedralRestraint>::iterator iri;
    for ( i=0,iri =this->_Terms.begin();
          iri!=this->_Terms.end(); iri++,i++ ) {
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfDihedralRestraintTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfDihedralRestraintTermsToCalculate ) {
          break;
        }
      }
#endif
  if ( this->_DebugEnergy ) {
    core::clasp_write_string(fmt::format("Evaluating term: {}\n" , i ));
  }
    
#define DEBUG_IMPROPER_RESTRAINT 1
#include	<cando/chem/energy_functions/_ImproperRestraint_termCode.cc>
#undef DEBUG_IMPROPER_RESTRAINT
      if ( EraseLinearDihedral == 0.0 ) {
        ERROR(chem::_sym_LinearDihedralError,core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(iri->_Atom1,iri->_Atom2,iri->_Atom3,iri->_Atom4),
                                                                      kw::_sym_coordinates,pos,
                                                                      kw::_sym_indices,core::Cons_O::createList(core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3), core::make_fixnum(I4))));
      }
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
      iri->_calcForce = calcForce;
      iri->_calcDiagonalHessian = calcDiagonalHessian;
      iri->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define EVAL_SET(var,val)	{ iri->eval.var=val;};
#include <cando/chem/energy_functions/_ImproperRestraint_debugEvalSet.cc>
#endif //]

      if ( this->_DebugEnergy ) 
      {
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} args cando\n" , (i+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} K {:f}\n" , (i+1) , K ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} U {:f}\n" , (i+1) , U ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} L {:f}\n" , (i+1) , L ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} x1 {:5.3f} {}\n" , (i+1) , x1 , (I1/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} y1 {:5.3f} {}\n" , (i+1) , y1 , (I1/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} z1 {:5.3f} {}\n" , (i+1) , z1 , (I1/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} x2 {:5.3f} {}\n" , (i+1) , x2 , (I2/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} y2 {:5.3f} {}\n" , (i+1) , y2 , (I2/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} z2 {:5.3f} {}\n" , (i+1) , z2 , (I2/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} x3 {:5.3f} {}\n" , (i+1) , x3 , (I3/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} y3 {:5.3f} {}\n" , (i+1) , y3 , (I3/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} z3 {:5.3f} {}\n" , (i+1) , z3 , (I3/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} x4 {:5.3f} {}\n" , (i+1) , x4 , (I4/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} y4 {:5.3f} {}\n" , (i+1) , y4 , (I4/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} z4 {:5.3f} {}\n" , (i+1) , z4 , (I4/3+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} results\n" , (i+1) ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} CosPhi {:f}\n" , (i+1) , CosPhi ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} SinPhi {:f}\n" , (i+1) , SinPhi ));
        if ( CosPhi>0.1 ) {
          Phi = safe_asin(SinPhi);
        } else {
          Phi = safe_acos(CosPhi)*SIGN(SinPhi);
        }
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} Phi {:f}\n" , (i+1) , Phi ));
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} Energy {:f}\n" , (i+1) , Energy));
        if ( calcForce ) 
        {
//			LOG_ENERGY(( "MEISTER improperRestraint %d DePhi %lf\n") , (i+1) , DePhi);
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fx1 {:8.5f} {}\n", (i+1) , fx1 , (I1/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fy1 {:8.5f} {}\n", (i+1) , fy1 , (I1/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fz1 {:8.5f} {}\n", (i+1) , fz1 , (I1/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fx2 {:8.5f} {}\n", (i+1) , fx2 , (I2/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fy2 {:8.5f} {}\n", (i+1) , fy2 , (I2/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fz2 {:8.5f} {}\n", (i+1) , fz2 , (I2/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fx3 {:8.5f} {}\n", (i+1) , fx3 , (I3/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fy3 {:8.5f} {}\n", (i+1) , fy3 , (I3/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fz3 {:8.5f} {}\n", (i+1) , fz3 , (I3/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fx4 {:8.5f} {}\n", (i+1) , fx4 , (I4/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fy4 {:8.5f} {}\n", (i+1) , fy4 , (I4/3+1) ));
          core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} fz4 {:8.5f} {}\n", (i+1) , fz4 , (I4/3+1) ));
        }
        core::clasp_write_string(fmt::format( "MEISTER improperRestraint {} stop\n" , (i+1) ));
        LOG_ENERGY(( "MEISTER improperRestraint %d args cando\n") , (i+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d K %lf\n") , (i+1) , K );
        LOG_ENERGY(( "MEISTER improperRestraint %d U %lf\n") , (i+1) , U );
        LOG_ENERGY(( "MEISTER improperRestraint %d L %lf\n") , (i+1) , L );
        LOG_ENERGY(( "MEISTER improperRestraint %d x1 %5.3lf %d\n") , (i+1) , x1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d y1 %5.3lf %d\n") , (i+1) , y1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d z1 %5.3lf %d\n") , (i+1) , z1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d x2 %5.3lf %d\n") , (i+1) , x2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d y2 %5.3lf %d\n") , (i+1) , y2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d z2 %5.3lf %d\n") , (i+1) , z2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d x3 %5.3lf %d\n") , (i+1) , x3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d y3 %5.3lf %d\n") , (i+1) , y3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d z3 %5.3lf %d\n") , (i+1) , z3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d x4 %5.3lf %d\n") , (i+1) , x4 , (I4/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d y4 %5.3lf %d\n") , (i+1) , y4 , (I4/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d z4 %5.3lf %d\n") , (i+1) , z4 , (I4/3+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d results\n") , (i+1) );
        LOG_ENERGY(( "MEISTER improperRestraint %d CosPhi %lf\n") , (i+1) , CosPhi );
        LOG_ENERGY(( "MEISTER improperRestraint %d SinPhi %lf\n") , (i+1) , SinPhi );
        if ( CosPhi>0.1 ) {
          Phi = safe_asin(SinPhi);
        } else {
          Phi = safe_acos(CosPhi)*SIGN(SinPhi);
        }
        LOG_ENERGY(( "MEISTER improperRestraint %d Phi %lf\n") , (i+1) , Phi );
        LOG_ENERGY(( "MEISTER improperRestraint %d Energy %lf\n") , (i+1) , Energy);
        if ( calcForce ) {
//			LOG_ENERGY(( "MEISTER improperRestraint %d DePhi %lf\n") , (i+1) , DePhi);
          LOG_ENERGY(( "MEISTER improperRestraint %d fx1 %8.5lf %d\n") , (i+1) , fx1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fy1 %8.5lf %d\n") , (i+1) , fy1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fz1 %8.5lf %d\n") , (i+1) , fz1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fx2 %8.5lf %d\n") , (i+1) , fx2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fy2 %8.5lf %d\n") , (i+1) , fy2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fz2 %8.5lf %d\n") , (i+1) , fz2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fx3 %8.5lf %d\n") , (i+1) , fx3 , (I3/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fy3 %8.5lf %d\n") , (i+1) , fy3 , (I3/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fz3 %8.5lf %d\n") , (i+1) , fz3 , (I3/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fx4 %8.5lf %d\n") , (i+1) , fx4 , (I4/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fy4 %8.5lf %d\n") , (i+1) , fy4 , (I4/3+1) );
          LOG_ENERGY(( "MEISTER improperRestraint %d fz4 %8.5lf %d\n") , (i+1) , fz4 , (I4/3+1) );
        }
        LOG_ENERGY(( "MEISTER improperRestraint %d stop\n") , (i+1) );
      }
			/* Add the forces */

    }
  }
  LOG_ENERGY(( "DihedralRestraint energy = %lf\n") , (double)(this->_TotalEnergy) );
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY(("%s }\n") , this->className());
  }
  maybeSetEnergy( componentEnergy, EnergyDihedralRestraint_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}

void	EnergyDihedralRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(chem::NVector_sp 	pos)
{
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//------------------
#define IMPROPER_RESTRAINT_CALC_FORCE
#define IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{x = iri->term.x;}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of) {x = pos->element(ii+of);}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x) {}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


  {
    
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ImproperRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    num_real K,U,L;
    num_real EraseLinearDihedral, UShift, PhiShift;
    bool RestraintActive;
    int	I1, I2, I3, I4,i;
    gctools::Vec0<EnergyDihedralRestraint>::iterator iri;
    for ( i=0,iri =this->_Terms.begin();
          iri!=this->_Terms.end(); iri++,i++ ) {
#define DEBUG_IMPROPER_RESTRAINT 1
#include	<cando/chem/energy_functions/_ImproperRestraint_termCode.cc>
#undef DEBUG_IMPROPER_RESTRAINT
      int index = i;
#include <cando/chem/energy_functions/_ImproperRestraint_debugFiniteDifference.cc>
    }
  }


}



int	EnergyDihedralRestraint_O::checkForBeyondThresholdInteractions(
			stringstream& info, chem::NVector_sp pos )
{
  int	fails = 0;

  this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------
#undef IMPROPER_RESTRAINT_CALC_FORCE
#undef IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN
#undef IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	IMPROPER_RESTRAINT_SET_PARAMETER
#define	IMPROPER_RESTRAINT_SET_PARAMETER(x)	{x = iri->term.x;}
#undef	IMPROPER_RESTRAINT_SET_POSITION
#define	IMPROPER_RESTRAINT_SET_POSITION(x,ii,of) {x = pos->element(ii+of);}
#undef	IMPROPER_RESTRAINT_PHI_SET
#define	IMPROPER_RESTRAINT_PHI_SET(x) {}
#undef	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE
#define	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	IMPROPER_RESTRAINT_FORCE_ACCUMULATE
#define	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


  {
    
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ImproperRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    num_real K,U,L;
    num_real EraseLinearDihedral, UShift, PhiShift;
    bool RestraintActive;
    int	I1, I2, I3, I4,i;
    gctools::Vec0<EnergyDihedralRestraint>::iterator iri;
    for ( i=0,iri =this->_Terms.begin();
          iri!=this->_Terms.end(); iri++,i++ ) {
#define DEBUG_IMPROPER_RESTRAINT 1
#include	<cando/chem/energy_functions/_ImproperRestraint_termCode.cc>
#undef DEBUG_IMPROPER_RESTRAINT
      if ( RestraintActive ) {
        chem::Atom_sp a1, a2, a3, a4;
        a1 = (*iri)._Atom1;
        a2 = (*iri)._Atom2;
        a3 = (*iri)._Atom3;
        a4 = (*iri)._Atom4;
        info<< "DihedralRestraintDeviation ";
        info<< "Phi[degrees](" << Phi/0.0174533 << ") ";
        info<< "U[degrees](" << U/0.0174533 << ") ";
        info<< "L[degrees](" << L/0.0174533 << ") ";
        info << a1->description() << " ";
        info << a2->description() << " ";
        info << a3->description() << " ";
        info << a4->description() << " ";
        info << std::endl;
        this->_BeyondThresholdTerms.push_back(*iri);
        fails++;
      }
    }
  }
  return fails;
}

EnergyDihedralRestraint_sp EnergyDihedralRestraint_O::copyFilter(core::T_sp keepInteraction) {
  EnergyDihedralRestraint_sp copy = EnergyDihedralRestraint_O::create();
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    Atom_sp a3 = edi->_Atom3;
    Atom_sp a4 = edi->_Atom4;
    if ( skipInteraction( keepInteraction, EnergyDihedralRestraint_O::staticClass(), a1, a2, a3, a4 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}


void EnergyDihedralRestraint_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(0.001);
}

#ifdef XML_ARCHIVE
void EnergyDihedralRestraint_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    IMPLEMENT_ME();
}
#endif


};
