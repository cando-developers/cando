/*
    File: energyFixedNonbond.cc
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
#include <clasp/core/evaluator.h>
#include <cando/chem/energyFixedNonbond.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

namespace chem {


core::List_sp FixedNonbondRestraint::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,fixed_atom),this->_FixedAtom),
                                  core::Cons_O::create(INTERN_(kw,fixed_charge),core::clasp_make_double_float(this->_FixedCharge)),
                                  core::Cons_O::create(INTERN_(kw,fixed_type),core::make_fixnum(this->_FixedType)),
                                  core::Cons_O::create(INTERN_(kw,fixed_pos_x),core::clasp_make_double_float(this->_FixedPosition.getX())),
                                  core::Cons_O::create(INTERN_(kw,fixed_pos_y),core::clasp_make_double_float(this->_FixedPosition.getY())),
                                  core::Cons_O::create(INTERN_(kw,fixed_pos_z),core::clasp_make_double_float(this->_FixedPosition.getZ())) );
}



void EnergyFixedNonbondRestraint_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

//
// Copy this from implementAmberFunction.cc
//
double	_evaluateEnergyOnly_FixedNonbond(
		double x1, double y1, double z1,
		double xf, double yf, double zf,
		double dA, double dC, double dQ1Q2)
{
    IMPLEMENT_ME();
#if 0
#undef	FNONBOND_SET_PARAMETER
#define	FNONBOND_SET_PARAMETER(x)	{}
#undef	FNONBOND_SET_POSITION
#define	FNONBOND_SET_POSITION(x,ii,of)	{}
#undef	FNONBOND_ENERGY_ACCUMULATE
#define	FNONBOND_ENERGY_ACCUMULATE(e) {}
#undef	FNONBOND_FORCE_ACCUMULATE
#define	FNONBOND_FORCE_ACCUMULATE(i,o,v) {}
#undef	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	FNONBOND_EFVDW_ENERGY_ACCUMULATE
#define	FNONBOND_EFVDW_ENERGY_ACCUMULATE(x) {}
#undef	FNONBOND_EFEEL_ENERGY_ACCUMULATE
#define	FNONBOND_EFEEL_ENERGY_ACCUMULATE(x) {}

#undef	FNONBOND_CALC_FORCE	// Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_FixedNonbond_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_FixedNonbond_termCode.cc>

    return Energy;
#endif
}








void	EnergyFixedNonbondRestraint_O::setupForEvaluation(AtomTable_sp atomTable, FFNonbondCrossTermTable_sp crossTerms)
{
    this->_MobileAtomTable = atomTable;
    this->_NonbondCrossTermTable = crossTerms;
}

string EnergyFixedNonbondRestraint_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyFixedNonbondRestraint_O,FixedNonbondRestraint>(*this);
}


void	EnergyFixedNonbondRestraint_O::zeroEnergy()
{
    this->Base::zeroEnergy();
    this->_EnergyElectrostatic = 0.0;
    this->_EnergyVdw = 0.0;
}


double	EnergyFixedNonbondRestraint_O::getEnergy()
{
double	e;
    e = this->getVdwEnergy();
    e += this->getElectrostaticEnergy();
    return e;
}


SYMBOL_EXPORT_SC_(ChemPkg,find_atom_type_position)

void	EnergyFixedNonbondRestraint_O::addFixedAtom(core::T_sp nonbondDb, Atom_sp fa)
{
    FixedNonbondRestraint entry;
    entry._FixedAtom = fa;
    entry._FixedCharge = fa->getCharge();
    try
    {
      if (fa->getType().nilp()) {
        SIMPLE_ERROR(BF("The atom type of %s is NIL!") % _rep_(fa));
      }
      core::T_sp pos = core::eval::funcall(_sym_find_atom_type_position,nonbondDb,fa->getType());
      if (pos.fixnump()) {
        entry._FixedType = pos.unsafe_fixnum();
      } else {
        TYPE_ERROR(pos,cl::_sym_fixnum);
      }
    } catch (UnknownType& err)
    {
	stringstream serr;
	serr << "Unknown type("<<fa->getType()<<") for fixed atom: " << fa->description() << std::endl;
	serr << "Were types assigned for the fixed molecule?" << std::endl;
	SIMPLE_ERROR(BF(serr.str()));
    }
    entry._FixedPosition = fa->getPosition();
    this->_Terms.push_back(entry);
}



void	EnergyFixedNonbondRestraint_O::dumpTerms()
{
  gc::Vec0<FixedNonbondRestraint>::iterator	eni;
string					as1;
string					str1, str2, str3, str4;
    for ( eni=this->_Terms.begin();
	    eni!=this->_Terms.end(); eni++ ) 
    {
	as1 = atomLabel(eni->_FixedAtom);
	_lisp->print(BF( "TERM 6FIXED_NONBOND %-9s") % as1.c_str() );
    }
}

#ifdef RENDER
int EnergyFixedNonbondRestraint_O::countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn)
{
IMPLEMENT_ME();
#if 0
vector<FixedNonbondRestraint>::iterator	eni;
string				as1,as2,as3,as4;
string				str1, str2, str3, str4;
int				overlapCount = 0;
int				ia1, ia2;
double				x1, y1, z1;
double				x2, y2, z2;
double				dx,dy,dz;
double				cutoff, distSquared;
bool				render;
RPGrColor			color;
geom::GrLine_sp			line;
Vector3				v1,v2;
    render = false;
    ANN(displayIn);
    if ( displayIn.notnilp() )
    {
        render = true;
	displayIn->clear();
	color = O_GrColor::systemColor(_lisp->symbol(_kw_yellow));
	displayIn->add(color);
    }
    for ( eni=this->_Terms.begin();
	    eni!=this->_Terms.end(); eni++ )
    {
	ia1 = eni->I1;
	ia2 = eni->I2;
	x1 = pos->element(ia1+0);
	y1 = pos->element(ia1+1);
	z1 = pos->element(ia1+2);
	x2 = pos->element(ia2+0);
	y2 = pos->element(ia2+1);
	z2 = pos->element(ia2+2);
	cutoff = eni->_RStar*scaleSumOfVdwRadii;
	cutoff = cutoff*cutoff;
	dx = x1-x2;
	dy = y1-y2;
	dz = z1-z2;
	distSquared = dx*dx+dy*dy+dz*dz;
	if ( distSquared<cutoff )
	{
	    LOG(BF("Found a close contact with distance = %lf and cutoff %lf") % sqrt(distSquared) % sqrt(cutoff)  );
	    LOG(BF("Atom1 = %s") % eni->_Atom1->description().c_str()  );
	    LOG(BF("Atom2 = %s") % eni->_Atom2->description().c_str()  );
	    if ( render )
	    {
	        v1.set(x1,y1,z1);
	        v2.set(x2,y2,z2);
	        line = geom::GrLine_O::create(v1,v2);
		displayIn->add(line);
	    }
	    overlapCount++;
	}
    }
    return overlapCount;
#endif
}
#endif





double EnergyFixedNonbondRestraint_O::evaluateAllComponent( ScoringFunction_sp score,
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
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(BF("%s {\n")% this->className());
  }
  if ( this->_Terms.size() == 0 ) return 0.0;
  ASSERTNOTNULL(this->_MobileAtomTable);
  ASSERT(this->_MobileAtomTable.notnilp());
  ASSERTNOTNULL(this->_NonbondCrossTermTable);
  ASSERT(this->_NonbondCrossTermTable.notnilp());

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

#define FNONBOND_CALC_FORCE
#define FNONBOND_CALC_DIAGONAL_HESSIAN
#define FNONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	FNONBOND_SET_PARAMETER
#define	FNONBOND_SET_PARAMETER(x)	{x=nbi->x;}
#undef	FNONBOND_SET_POSITION
#define	FNONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	FNONBOND_EFEEL_ENERGY_ACCUMULATE
#define	FNONBOND_EFEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);}
#undef	FNONBOND_EFVDW_ENERGY_ACCUMULATE
#define	FNONBOND_EFVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);}
#undef	FNONBOND_ENERGY_ACCUMULATE
#define	FNONBOND_ENERGY_ACCUMULATE(e) {};
#undef	FNONBOND_FORCE_ACCUMULATE
#undef	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	FNONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc


  if ( this->isEnabled() ) 
  {
	// If you are going to use openmp here, you need to control access to the force and hessian
	// arrays so that only one thread updates each element at a time.
    {
      LOG(BF("FixedNonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_FixedNonbond_termDeclares.cc>
#pragma clang diagnostic pop
      double x1,y1,z1,xf,yf,zf, dQ1Q2;
      uint I1;
      {   
		//
		// From here to termCode is the fixedNonbond interaction loop
		//
        uint ifixed, imobile;
        uint fixedNonbondAtoms = this->_Terms.size();
        uint mobileNonbondAtoms = this->_MobileAtomTable->getNumberOfAtoms();
        uint i = 0;
        for ( ifixed = 0; ifixed < fixedNonbondAtoms; ifixed++ )
        {
          FixedNonbondRestraint& fixedAtomEntry = this->_Terms[ifixed];
          xf = fixedAtomEntry._FixedPosition.getX();
          yf = fixedAtomEntry._FixedPosition.getY();
          zf = fixedAtomEntry._FixedPosition.getZ();
          double fixedChargeMultiplier = fixedAtomEntry._FixedCharge * this->getElectrostaticScale() / this->_DielectricConstant * ELECTROSTATIC_MODIFIER;
          uint fixedTypeMajorIndex = this->_NonbondCrossTermTable->typeMajorIndex(fixedAtomEntry._FixedType);
          for ( imobile=0; imobile< mobileNonbondAtoms; imobile++ )
          {
            EnergyAtom& mobileAtomEntry = this->_MobileAtomTable->energyAtomEntry(imobile);
            I1 = mobileAtomEntry.coordinateIndexTimes3();
            dQ1Q2 = mobileAtomEntry._Charge * fixedChargeMultiplier;
            uint crossTermIndex = fixedTypeMajorIndex + mobileAtomEntry._TypeIndex;
            FFNonbondCrossTerm crossTerm = this->_NonbondCrossTermTable->nonbondCrossTerm(crossTermIndex);
            double dA = crossTerm._A*this->getVdwScale();
            double dC = crossTerm._C*this->getVdwScale();
#include	<cando/chem/energy_functions/_FixedNonbond_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
            fixedAtomEntry._calcForce = calcForce;
            fixedAtomEntry._calcDiagonalHessian = calcDiagonalHessian;
            fixedAtomEntry._calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ fixedAtomEntry.eval.var=val;};
#include	<cando/chem/energy_functions/_FixedNonbond_debugEvalSet.cc>
#endif //]

            if ( calcForce ) {
//			    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
//			    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
//			    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
            }

            if ( this->_DebugEnergy ) 
            {
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d args cando\n")% (i+1) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d dA %5.3lf\n")% (i+1) % dA );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d dC %5.3lf\n")% (i+1) % dC );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d dQ1Q2 %5.3lf\n")% (i+1) % dQ1Q2 );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d x1 %5.3lf %d\n")% (i+1) % x1 % (I1/3+1) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d y1 %5.3lf %d\n")% (i+1) % y1 % (I1/3+1) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d z1 %5.3lf %d\n")% (i+1) % z1 % (I1/3+1) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d xf %5.3lf %d\n")% (i+1) % xf % (ifixed) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d yf %5.3lf %d\n")% (i+1) % yf % (ifixed) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d zf %5.3lf %d\n")% (i+1) % zf % (ifixed) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d results\n")% (i+1) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d Efixed nonbond %lf\n")% (i+1) % (Efvdw+Efeel) );
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d efvdw %lf\n")% (i+1) % Efvdw);
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d efeel %lf\n")% (i+1) % Efeel);
              if ( calcForce ) {
                LOG_ENERGY(BF( "MEISTER fixed nonbond %d fx1 %lf %d\n")% (i+1) % fx1 % (I1/3+1) );
                LOG_ENERGY(BF( "MEISTER fixed nonbond %d fy1 %lf %d\n")% (i+1) % fy1 % (I1/3+1) );
                LOG_ENERGY(BF( "MEISTER fixed nonbond %d fz1 %lf %d\n")% (i+1) % fz1 % (I1/3+1) );
              }
              LOG_ENERGY(BF( "MEISTER fixed nonbond %d stop\n")% (i+1) );
            }
            i++;
          }
        }
      }
    }
  } else {
    LOG_ENERGY( "FixedNonbond component is not enabled" );
  }
  LOG_ENERGY(BF( "          Vdw energy = %lf\n")% (double)this->_EnergyVdw);
  LOG_ENERGY(BF( "Electrostatic energy = %lf\n")% (double)this->_EnergyElectrostatic);
  LOG_ENERGY(BF("%s }\n")% this->className() );
  return this->_EnergyVdw+this->_EnergyElectrostatic;
}







void	EnergyFixedNonbondRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
		NVector_sp 	pos)
{
    IMPLEMENT_ME();
#if 0
int	fails = 0;
bool	calcForce = true;
bool	calcDiagonalHessian = true;
bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//------------------


#define FNONBOND_CALC_FORCE
#define FNONBOND_CALC_DIAGONAL_HESSIAN
#define FNONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	FNONBOND_SET_PARAMETER
#define	FNONBOND_SET_PARAMETER(x)	{x=nbi->x;}
#undef	FNONBOND_SET_POSITION
#define	FNONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	FNONBOND_EFEEL_ENERGY_ACCUMULATE
#define	FNONBOND_EFEEL_ENERGY_ACCUMULATE(e) {}
#undef	FNONBOND_EFVDW_ENERGY_ACCUMULATE
#define	FNONBOND_EFVDW_ENERGY_ACCUMULATE(e) {}
#undef	FNONBOND_ENERGY_ACCUMULATE
#define	FNONBOND_ENERGY_ACCUMULATE(e) {};
#undef	FNONBOND_FORCE_ACCUMULATE
#define	FNONBOND_FORCE_ACCUMULATE(i,o,v) {}
#undef	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

	if ( this->isEnabled() ) {
		_BLOCK_TRACE("FixedNonbondRestraint finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_FixedNonbond_termDeclares.cc>
#pragma clang diagnostic pop
	    double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
	    int	I1, I2,i;
	    vector<FixedNonbondRestraint>::iterator nbi;
	    for ( i=0,nbi=this->_Terms.begin();
			nbi!=this->_Terms.end(); nbi++,i++ ) {
#include	<cando/chem/energy_functions/_FixedNonbond_termCode.cc>
		int index = i;
#include <cando/chem/energy_functions/_FixedNonbond_debugFiniteDifference.cc>

	    }
	}
#endif
}



#if 0
int	EnergyFixedNonbondRestraint_O::checkForBeyondThresholdInteractions(
			stringstream& info, NVector_sp pos )
{
int	fails = 0;

    this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------
#undef FNONBOND_CALC_FORCE
#undef FNONBOND_CALC_DIAGONAL_HESSIAN
#undef FNONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	FNONBOND_SET_PARAMETER
#define	FNONBOND_SET_PARAMETER(x)	{x=nbi->x;}
#undef	FNONBOND_SET_POSITION
#define	FNONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	FNONBOND_EFEEL_ENERGY_ACCUMULATE
#define	FNONBOND_EFEEL_ENERGY_ACCUMULATE(e) {}
#undef	FNONBOND_EFVDW_ENERGY_ACCUMULATE
#define	FNONBOND_EFVDW_ENERGY_ACCUMULATE(e) {}
#undef	FNONBOND_ENERGY_ACCUMULATE
#define	FNONBOND_ENERGY_ACCUMULATE(e) {};
#undef	FNONBOND_FORCE_ACCUMULATE
#define	FNONBOND_FORCE_ACCUMULATE(i,o,v) {}
#undef	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

    if ( this->isEnabled() ) 
    { _BLOCK_TRACE("FixedNonbondRestraint finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_FixedNonbond_termDeclares.cc>
#pragma clang diagnostic pop
	double x1,y1,z1,xf,yf,zf,dQ1Q2;
	int	I1;
	{
	    //
	    // From here to termCode is the fixedNonbond interaction loop
	    //
	    uint ifixed, imobile;
	    uint fixedNonbondAtoms = this->_Terms.size();
	    uint mobileNonbondAtoms = this->_MobileAtomTable->getNumberOfAtoms();
	    for ( ifixed = 0; ifixed < fixedNonbondAtoms; ifixed++ )
	    {
		FixedNonbondRestraint& fixedAtomEntry = this->_Terms[ifixed];
		xf = fixedAtomEntry._FixedPosition.getX();
		yf = fixedAtomEntry._FixedPosition.getY();
		zf = fixedAtomEntry._FixedPosition.getZ();
		double fixedChargeMultiplier = fixedAtomEntry._FixedCharge * this->getElectrostaticScale() / this->_DielectricConstant * ELECTROSTATIC_MODIFIER;
		uint fixedTypeMajorIndex = this->_NonbondCrossTermTable->typeMajorIndex(fixedAtomEntry._FixedType);
		for ( imobile=0; imobile< mobileNonbondAtoms; imobile++ )
		{
		    EnergyAtom& mobileAtomEntry = this->_MobileAtomTable->energyAtomEntry(imobile);
		    I1 = mobileAtomEntry.coordinateIndexTimes3();
		    dQ1Q2 = mobileAtomEntry._Charge * fixedChargeMultiplier;
		    uint crossTermIndex = fixedTypeMajorIndex + mobileAtomEntry._TypeIndex;
		    FFNonbondCrossTerm crossTerm = this->_NonbondCrossTermTable->nonbondCrossTerm(crossTermIndex);
		    double dA = crossTerm._A*this->getVdwScale();
		    double dC = crossTerm._C*this->getVdwScale();
#include	<cando/chem/energy_functions/_FixedNonbond_termCode.cc>
		    if ( NonbondDistance < this->_ErrorThreshold ) 
		    {
			Atom_sp a1, af;
			a1 = mobileAtomEntry.atom();
			af = fixedAtomEntry._FixedAtom;
			info<< "FixedNonbondDeviation";
	//		info<< a1->getAbsoluteIdPath() << " ";
	//		info<< af->getAbsoluteIdPath() << " ";
			info<< "value " << NonbondDistance << " ";
			info<<"threshold " << this->_ErrorThreshold;
			info << a1->getName() << " ";
			info << af->getName() << " ";
			info << std::endl;
			EnergyFixedNonbondRestraint_O::BeyondThresholdTermType bt;
			bt._MobileAtom = a1;
			bt._FixedAtom = af;
			bt._NonbondDistance = NonbondDistance;
			this->_BeyondThresholdTerms.push_back(bt);
			fails++;
		    }
		}
	    }
	}
    }
    return fails;
}
#endif














void EnergyFixedNonbondRestraint_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(1.0);
    this->setDielectricConstant(80.0);
    this->setVdwScale(1.0);
    this->setElectrostaticScale(1.0);
    this->_MobileAtomTable = _Nil<AtomTable_O>();
    this->_NonbondCrossTermTable = _Nil<FFNonbondCrossTermTable_O>();
}


#ifdef XML_ARCHIVE
void EnergyFixedNonbondRestraint_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    IMPLEMENT_ME();
}
#endif


};
