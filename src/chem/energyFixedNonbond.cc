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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */

#define DEBUG_LEVEL_NONE

#include <clasp/core/foundation.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/energyFixedNonbond.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

namespace chem {

core::List_sp FixedNonbondRestraint::encode() const {
  return core::Cons_O::createList(
      core::Cons_O::create(INTERN_(kw, fixed_atom), this->_FixedAtom),
      core::Cons_O::create(INTERN_(kw, fixed_charge), core::clasp_make_double_float(this->_FixedCharge)),
      core::Cons_O::create(INTERN_(kw, fixed_type), core::make_fixnum(this->_FixedType)),
      core::Cons_O::create(INTERN_(kw, fixed_pos_x), core::clasp_make_double_float(this->_FixedPosition.getX())),
      core::Cons_O::create(INTERN_(kw, fixed_pos_y), core::clasp_make_double_float(this->_FixedPosition.getY())),
      core::Cons_O::create(INTERN_(kw, fixed_pos_z), core::clasp_make_double_float(this->_FixedPosition.getZ())));
}

void EnergyFixedNonbondRestraint_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms), this->_Terms);
  this->Base::fields(node);
}

//
// Copy this from implementAmberFunction.cc
//
double _evaluateEnergyOnly_FixedNonbond(num_real x1, num_real y1, num_real z1, num_real xf, num_real yf, num_real zf, num_real dA,
                                        num_real dC, num_real dQ1Q2) {
  IMPLEMENT_ME();
#if 0
#undef FNONBOND_SET_PARAMETER
#define FNONBOND_SET_PARAMETER(x)                                                                                                  \
  {}
#undef FNONBOND_SET_POSITION
#define FNONBOND_SET_POSITION(x, ii, of)                                                                                           \
  {}
#undef FNONBOND_ENERGY_ACCUMULATE
#define FNONBOND_ENERGY_ACCUMULATE(e)                                                                                              \
  {}
#undef FNONBOND_FORCE_ACCUMULATE
#define FNONBOND_FORCE_ACCUMULATE(i, o, v)                                                                                         \
  {}
#undef FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1, o1, i2, o2, v)                                                                    \
  {}
#undef FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1, o1, i2, o2, v)                                                                \
  {}
#undef FNONBOND_EFVDW_ENERGY_ACCUMULATE
#define FNONBOND_EFVDW_ENERGY_ACCUMULATE(x)                                                                                        \
  {}
#undef FNONBOND_EFEEL_ENERGY_ACCUMULATE
#define FNONBOND_EFEEL_ENERGY_ACCUMULATE(x)                                                                                        \
  {}

#undef FNONBOND_CALC_FORCE // Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_FixedNonbond_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_FixedNonbond_termCode.cc>

    return Energy;
#endif
}

void EnergyFixedNonbondRestraint_O::setupForEvaluation(AtomTable_sp atomTable, FFNonbondCrossTermTable_sp crossTerms) {
  this->_MobileAtomTable = atomTable;
  this->_NonbondCrossTermTable = crossTerms;
}

string EnergyFixedNonbondRestraint_O::beyondThresholdInteractionsAsString() {
  return component_beyondThresholdInteractionsAsString<EnergyFixedNonbondRestraint_O, FixedNonbondRestraint>(*this);
}

SYMBOL_EXPORT_SC_(ChemPkg, find_atom_type_position)

void EnergyFixedNonbondRestraint_O::addFixedAtom(core::T_sp nonbondDb, Atom_sp fa, core::HashTable_sp atomTypes) {
  FixedNonbondRestraint entry;
  entry._FixedAtom = fa;
  entry._FixedCharge = fa->getCharge();
  try {
    if (fa->getType(atomTypes).nilp()) {
      SIMPLE_ERROR("In EnergyFixedNonbondRestraint_O::addFixedAtom - the atom type of {} is NIL!", _rep_(fa));
    }
    core::T_sp pos = core::eval::funcall(_sym_find_atom_type_position, nonbondDb, fa->getType(atomTypes));
    if (pos.fixnump()) {
      entry._FixedType = pos.unsafe_fixnum();
    } else {
      TYPE_ERROR(pos, cl::_sym_fixnum);
    }
  } catch (UnknownType &err) {
    SIMPLE_ERROR("Unknown type({}) for fixed atom: {}\nWere types assigned for the fixed molecule?", fa->getType(atomTypes),
                 fa->description());
  }
  entry._FixedPosition = fa->getPosition();
  this->_Terms.push_back(entry);
}

void EnergyFixedNonbondRestraint_O::dumpTerms(core::HashTable_sp atomTypes) {
  gc::Vec0<FixedNonbondRestraint>::iterator eni;
  string as1;
  string str1, str2, str3, str4;
  for (eni = this->_Terms.begin(); eni != this->_Terms.end(); eni++) {
    as1 = atomLabel(eni->_FixedAtom);
    core::clasp_write_string(fmt::format("TERM 6FIXED_NONBOND {:<9}\n", as1.c_str()));
  }
}

#ifdef RENDER
int EnergyFixedNonbondRestraint_O::countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn) {
  IMPLEMENT_ME();
#if 0
vector<FixedNonbondRestraint>::iterator	eni;
string				as1,as2,as3,as4;
string				str1, str2, str3, str4;
int				overlapCount = 0;
int				ia1, ia2;
num_real				x1, y1, z1;
num_real				x2, y2, z2;
num_real				dx,dy,dz;
num_real				cutoff, distSquared;
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
	    LOG("Found a close contact with distance = {} and cutoff {}" , sqrt(distSquared) , sqrt(cutoff)  );
	    LOG("Atom1 = {}" , eni->_Atom1->description().c_str()  );
	    LOG("Atom2 = {}" , eni->_Atom2->description().c_str()  );
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

SYMBOL_EXPORT_SC_(ChemPkg, energyVdwExcludedAtoms);
SYMBOL_EXPORT_SC_(ChemPkg, energyElectrostaticExcludedAtoms);
SYMBOL_EXPORT_SC_(ChemPkg, energyVdw);
SYMBOL_EXPORT_SC_(ChemPkg, energyElectrostatic);
SYMBOL_EXPORT_SC_(ChemPkg, energyVdw14);
SYMBOL_EXPORT_SC_(ChemPkg, energyElectrostatic14);

double EnergyFixedNonbondRestraint_O::evaluateAllComponent(ScoringFunction_sp score,
                                                           NVector_sp pos,
                                                           core::T_sp energyScale,
                                                           core::T_sp componentEnergy,
                                                           bool calcForce, gc::Nilable<NVector_sp> force, bool calcDiagonalHessian,
                                                           bool calcOffDiagonalHessian,
                                                           gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                                                           gc::Nilable<NVector_sp> hdvec, gc::Nilable<NVector_sp> dvec,
                                                           core::T_sp activeAtomMask, core::T_sp debugInteractions) {

  double dielectricConstant = energyScaleDielectricConstant(energyScale);
  double amber_charge_conversion_18dot2223 =
      core::Number_O::as_double_float(gc::As<core::Number_sp>(_sym_STARamber_charge_conversion_18_DOT_2223STAR->symbolValue()));
  double dQ1Q2Scale = amber_charge_conversion_18dot2223 * amber_charge_conversion_18dot2223;
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  this->_Evaluations++;
  if (this->_Terms.size() == 0)
    return 0.0;
  ASSERTNOTNULL(this->_MobileAtomTable);
  ASSERT(this->_MobileAtomTable.notnilp());
  ASSERTNOTNULL(this->_NonbondCrossTermTable);
  ASSERT(this->_NonbondCrossTermTable.notnilp());

  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool hasForce = force.notnilp();
  bool hasHessian = hessian.notnilp();
  bool hasHdAndD = (hdvec.notnilp()) && (dvec.notnilp());

  //
  // Copy from implementAmberFunction::evaluateAll
  //
  // -----------------------

  num_real energyElectrostatic = 0.0;
  num_real energyVdw = 0.0;
#define FNONBOND_CALC_FORCE
#define FNONBOND_CALC_DIAGONAL_HESSIAN
#define FNONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef FNONBOND_SET_PARAMETER
#define FNONBOND_SET_PARAMETER(x)                                                                                                  \
  { x = nbi->x; }
#undef FNONBOND_SET_POSITION
#define FNONBOND_SET_POSITION(x, ii, of)                                                                                           \
  { x = pos->element(ii + of); }
#undef FNONBOND_EFEEL_ENERGY_ACCUMULATE
#define FNONBOND_EFEEL_ENERGY_ACCUMULATE(e)                                                                                        \
  { energyElectrostatic += (e); }
#undef FNONBOND_EFVDW_ENERGY_ACCUMULATE
#define FNONBOND_EFVDW_ENERGY_ACCUMULATE(e)                                                                                        \
  { energyVdw += (e); }
#undef FNONBOND_ENERGY_ACCUMULATE
#define FNONBOND_ENERGY_ACCUMULATE(e) {};
#undef FNONBOND_FORCE_ACCUMULATE
#undef FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define FNONBOND_FORCE_ACCUMULATE ForceAcc
#define FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE DiagHessAcc
#define FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc

  {
    // If you are going to use openmp here, you need to control access to the force and hessian
    // arrays so that only one thread updates each element at a time.
    {
      LOG("FixedNonbond component is enabled");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_FixedNonbond_termDeclares.cc>
#pragma clang diagnostic pop
      num_real x1, y1, z1, xf, yf, zf, dQ1Q2;
      uint I1;
      {
        //
        // From here to termCode is the fixedNonbond interaction loop
        //
        uint ifixed, imobile;
        uint fixedNonbondAtoms = this->_Terms.size();
        uint mobileNonbondAtoms = this->_MobileAtomTable->getNumberOfAtoms();
        uint i = 0;
        for (ifixed = 0; ifixed < fixedNonbondAtoms; ifixed++) {
          FixedNonbondRestraint &fixedAtomEntry = this->_Terms[ifixed];
          xf = fixedAtomEntry._FixedPosition.getX();
          yf = fixedAtomEntry._FixedPosition.getY();
          zf = fixedAtomEntry._FixedPosition.getZ();
          num_real fixedChargeMultiplier =
              fixedAtomEntry._FixedCharge * this->getElectrostaticScale() / dielectricConstant * dQ1Q2Scale;
          uint fixedTypeMajorIndex = this->_NonbondCrossTermTable->typeMajorIndex(fixedAtomEntry._FixedType);
          for (imobile = 0; imobile < mobileNonbondAtoms; imobile++) {
            EnergyAtom &mobileAtomEntry = this->_MobileAtomTable->energyAtomEntry(imobile);
            I1 = mobileAtomEntry.coordinateIndexTimes3();
            dQ1Q2 = mobileAtomEntry._Charge * fixedChargeMultiplier;
            uint crossTermIndex = fixedTypeMajorIndex + mobileAtomEntry._TypeIndex;
            FFNonbondCrossTerm crossTerm = this->_NonbondCrossTermTable->nonbondCrossTerm(crossTermIndex);
            num_real dA = crossTerm._A * this->getVdwScale();
            num_real dC = crossTerm._C * this->getVdwScale();
#include <cando/chem/energy_functions/_FixedNonbond_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
            fixedAtomEntry._calcForce = calcForce;
            fixedAtomEntry._calcDiagonalHessian = calcDiagonalHessian;
            fixedAtomEntry._calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define EVAL_SET(var, val)                                                                                                         \
  { fixedAtomEntry.eval.var = val; };
#include <cando/chem/energy_functions/_FixedNonbond_debugEvalSet.cc>
#endif //]

            if (calcForce) {
              //			    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
              //			    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
              //			    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
            }
            i++;
          }
        }
      }
    }
  }
  LOG_ENERGY(("          Vdw energy = %lf\n"), (double)this->_EnergyVdw);
  LOG_ENERGY(("Electrostatic energy = %lf\n"), (double)this->_EnergyElectrostatic);
  LOG_ENERGY(("%s }\n"), this->className());
  maybeSetEnergy(componentEnergy, _sym_energyVdw, energyVdw);
  maybeSetEnergy(componentEnergy, _sym_energyElectrostatic, energyElectrostatic);
  return energyVdw + energyElectrostatic;
}

void EnergyFixedNonbondRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos) {
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
#undef FNONBOND_SET_PARAMETER
#define FNONBOND_SET_PARAMETER(x)                                                                                                  \
  { x = nbi->x; }
#undef FNONBOND_SET_POSITION
#define FNONBOND_SET_POSITION(x, ii, of)                                                                                           \
  { x = pos->element(ii + of); }
#undef FNONBOND_EFEEL_ENERGY_ACCUMULATE
#define FNONBOND_EFEEL_ENERGY_ACCUMULATE(e)                                                                                        \
  {}
#undef FNONBOND_EFVDW_ENERGY_ACCUMULATE
#define FNONBOND_EFVDW_ENERGY_ACCUMULATE(e)                                                                                        \
  {}
#undef FNONBOND_ENERGY_ACCUMULATE
#define FNONBOND_ENERGY_ACCUMULATE(e) {};
#undef FNONBOND_FORCE_ACCUMULATE
#define FNONBOND_FORCE_ACCUMULATE(i, o, v)                                                                                         \
  {}
#undef FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1, o1, i2, o2, v)                                                                    \
  {}
#undef FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1, o1, i2, o2, v)                                                                \
  {}

 {

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_FixedNonbond_termDeclares.cc>
#pragma clang diagnostic pop
	    num_real x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
	    int	I1, I2,i;
	    vector<FixedNonbondRestraint>::iterator nbi;
	    for ( i=0,nbi=this->_Terms.begin();
			nbi!=this->_Terms.end(); nbi++,i++ ) {
#include <cando/chem/energy_functions/_FixedNonbond_termCode.cc>
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
#undef FNONBOND_SET_PARAMETER
#define FNONBOND_SET_PARAMETER(x)                                                                                                  \
  { x = nbi->x; }
#undef FNONBOND_SET_POSITION
#define FNONBOND_SET_POSITION(x, ii, of)                                                                                           \
  { x = pos->element(ii + of); }
#undef FNONBOND_EFEEL_ENERGY_ACCUMULATE
#define FNONBOND_EFEEL_ENERGY_ACCUMULATE(e)                                                                                        \
  {}
#undef FNONBOND_EFVDW_ENERGY_ACCUMULATE
#define FNONBOND_EFVDW_ENERGY_ACCUMULATE(e)                                                                                        \
  {}
#undef FNONBOND_ENERGY_ACCUMULATE
#define FNONBOND_ENERGY_ACCUMULATE(e) {};
#undef FNONBOND_FORCE_ACCUMULATE
#define FNONBOND_FORCE_ACCUMULATE(i, o, v)                                                                                         \
  {}
#undef FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1, o1, i2, o2, v)                                                                    \
  {}
#undef FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1, o1, i2, o2, v)                                                                \
  {}

    {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_FixedNonbond_termDeclares.cc>
#pragma clang diagnostic pop
	num_real x1,y1,z1,xf,yf,zf,dQ1Q2;
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
		num_real fixedChargeMultiplier = fixedAtomEntry._FixedCharge * this->getElectrostaticScale() / this->_DielectricConstant * ELECTROSTATIC_MODIFIER;
		uint fixedTypeMajorIndex = this->_NonbondCrossTermTable->typeMajorIndex(fixedAtomEntry._FixedType);
		for ( imobile=0; imobile< mobileNonbondAtoms; imobile++ )
		{
		    EnergyAtom& mobileAtomEntry = this->_MobileAtomTable->energyAtomEntry(imobile);
		    I1 = mobileAtomEntry.coordinateIndexTimes3();
		    dQ1Q2 = mobileAtomEntry._Charge * fixedChargeMultiplier;
		    uint crossTermIndex = fixedTypeMajorIndex + mobileAtomEntry._TypeIndex;
		    FFNonbondCrossTerm crossTerm = this->_NonbondCrossTermTable->nonbondCrossTerm(crossTermIndex);
		    num_real dA = crossTerm._A*this->getVdwScale();
		    num_real  dC = crossTerm._C*this->getVdwScale();
#include <cando/chem/energy_functions/_FixedNonbond_termCode.cc>
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

void EnergyFixedNonbondRestraint_O::initialize() {
  this->Base::initialize();
  this->setDielectricConstant(80.0);
  this->setVdwScale(1.0);
  this->setElectrostaticScale(1.0);
  this->_MobileAtomTable = nil<AtomTable_O>();
  this->_NonbondCrossTermTable = nil<FFNonbondCrossTermTable_O>();
}

EnergyFixedNonbondRestraint_sp EnergyFixedNonbondRestraint_O::copyFilter(core::T_sp keepInteractionFactory) {
  EnergyFixedNonbondRestraint_sp copy = EnergyFixedNonbondRestraint_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyFixedNonbondRestraint_O::staticClass() );
  for (auto edi = this->_Terms.begin(); edi != this->_Terms.end(); edi++) {
    Atom_sp a1 = edi->_FixedAtom;
    if (skipInteraction(keepInteraction, a1)) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}

#ifdef XML_ARCHIVE
void EnergyFixedNonbondRestraint_O::archiveBase(core::ArchiveP node) {
  this->Base::archiveBase(node);
  IMPLEMENT_ME();
}
#endif

}; // namespace chem
