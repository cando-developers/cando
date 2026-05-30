/*
    File: energyChiralRestraint.cc
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


#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/wrappers.h>



//#define	CHIRAL_RESTRAINT_WEIGHT 1000000.0	// 10000.0 is used in MOE


namespace chem {

SYMBOL_EXPORT_SC_(KeywordPkg,ChiralRestraint)

#define CHIRAL_RESTRAINT_APPLY_ATOM_MASK(I1,I2,I3,I4) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         && bitvectorActiveAtomMask->testBit(I4/3) \
         ) \
    ) goto SKIP_term;
#define CHIRAL_RESTRAINT_DEBUG_INTERACTIONS(I1,I2,I3,I4) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,kw::_sym_ChiralRestraint, \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3), core::make_fixnum(I4)); \
    }

#include "cando/chem/energyKernels/chiral_restraint.c"

string	EnergyChiralRestraint::description()
{
    stringstream ss;
    ss << "EnergyChiralRestraint[";
    ss << "atoms: ";
    ss << this->_Atom1->description() << "-";
    ss << this->_Atom2->description() << "-";
    ss << this->_Atom3->description() << "-";
    ss << this->_Atom4->description() ;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
    ss << " eval.Energy=" << this->eval.Energy;
#endif
    return ss.str();
}





double	_evaluateEnergyOnly_ChiralRestraint(
		num_real x1, num_real y1, num_real z1,
		num_real x2, num_real y2, num_real z2,
		num_real x3, num_real y3, num_real z3,
		num_real x4, num_real y4, num_real z4,
		num_real K, num_real CO )
{
  double localPos[12] = { (double)x1, (double)y1, (double)z1,
                          (double)x2, (double)y2, (double)z2,
                          (double)x3, (double)y3, (double)z3,
                          (double)x4, (double)y4, (double)z4 };
  double energy_accum = 0.0;
  Chiral_Restraint<NoHessian> chiral;
  chiral_term lt((double)K, (double)CO, 0, 3, 6, 9);
  return chiral.energy(lt, localPos, &energy_accum);
}

SYMBOL_EXPORT_SC_(KeywordPkg,central_atom3);

core::List_sp EnergyChiralRestraint::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,k),core::clasp_make_double_float(this->term.k)),
                                  core::Cons_O::create(INTERN_(kw,CO),core::clasp_make_double_float(this->term.co)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.i3x1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.i3x2)),
                                  core::Cons_O::create(INTERN_(kw,i3), core::make_fixnum(this->term.i3x3)),
                                  core::Cons_O::create(INTERN_(kw,i4), core::make_fixnum(this->term.i3x4)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                  core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2),
                                  core::Cons_O::create(INTERN_(kw,central_atom3), this->_Atom3),
                                  core::Cons_O::create(INTERN_(kw,atom4), this->_Atom4));
}

std::string EnergyChiralRestraint_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}




EnergyComponent_sp EnergyChiralRestraint_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyChiralRestraint_O::staticClass() );
  EnergyChiralRestraint_sp copy = EnergyChiralRestraint_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    Atom_sp a3 = edi->_Atom3;
    Atom_sp a4 = edi->_Atom4;
    if ( skipInteraction( keepInteraction, a1, a2, a3, a4 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}


void EnergyChiralRestraint_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}


CL_LISPIFY_NAME(make-energy-chiral-restraint);
CL_DEF_CLASS_METHOD
EnergyChiralRestraint_sp EnergyChiralRestraint_O::make(EnergyFunction_sp energyFunction) {
  return ensureComponent<EnergyChiralRestraint_O>(energyFunction);
}

CL_DOCSTRING(R"dx(Add one chiral restraint term over atoms a1-a2-a3(center)-a4 with
force constant k and offset co.  The per-term energy is the one-sided cubic
E = k*(co+Q)^3 when that value is positive, else 0, where Q is the normalized
signed volume of the four atoms.  Flipping the sign of k selects the opposite
chirality domain.)dx")
CL_DEFMETHOD
size_t EnergyChiralRestraint_O::addChiralRestraintTerm(EnergyFunction_sp energyFunction, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double k, double co) {
  AtomTable_sp atomTable = energyFunction->atomTable();
  EnergyChiralRestraint t;
  t._Atom1 = a1;
  t._Atom2 = a2;
  t._Atom3 = a3;
  t._Atom4 = a4;
  t.term.k = k;
  t.term.co = co;
  t.term.i3x1 = atomTable->getEnergyAtomPointer(a1)->coordinateIndexTimes3();
  t.term.i3x2 = atomTable->getEnergyAtomPointer(a2)->coordinateIndexTimes3();
  t.term.i3x3 = atomTable->getEnergyAtomPointer(a3)->coordinateIndexTimes3();
  t.term.i3x4 = atomTable->getEnergyAtomPointer(a4)->coordinateIndexTimes3();
  this->_Terms.push_back(t);
  return this->_Terms.size()-1;
}

void EnergyChiralRestraint_O::addTerm(const EnergyChiralRestraint& e)
{
    this->_Terms.push_back(e);
}

void EnergyChiralRestraint_O::dumpTerms(core::HashTable_sp atomTypes)
{
    gctools::Vec0<EnergyChiralRestraint>::iterator cri;
    int i;
    string as1, as2, as3, as4;
    for ( i=0,cri=this->_Terms.begin();
	  cri!=this->_Terms.end(); cri++,i++ ) {
	as1 = atomLabel(cri->_Atom1);
	as2 = atomLabel(cri->_Atom2);
	as3 = atomLabel(cri->_Atom3);
	as4 = atomLabel(cri->_Atom4);
        core::clasp_write_string(fmt::format("TERM 7CHIRAL {:<9} {:<9} {:<9} {:<9} {:8.2f} {:8.2f} ; I1={} I2={} I3={} I4={}\n" , as1 , as2 , as3 , as4 , cri->term.k , cri->term.co , cri->term.i3x1 , cri->term.i3x2 , cri->term.i3x3 , cri->term.i3x4 ));
    }
}



void	EnergyChiralRestraint_O::setupHessianPreconditioner(chem::NVector_sp nvPosition,
                                                            chem::AbstractLargeSquareMatrix_sp m,
                                                            core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  Chiral_Restraint<double*> chiral;
  constexpr size_t PS = Chiral_Restraint<double*>::PositionSize;  // 12
  double localPos[PS];
  double localForce[PS];
  double localHess[PS*PS];
  double localDvec[PS];
  double localHdvec[PS];
  for ( gctools::Vec0<EnergyChiralRestraint>::iterator cri=this->_Terms.begin();
        cri!=this->_Terms.end(); cri++ ) {
    int I1 = cri->term.i3x1;
    int I2 = cri->term.i3x2;
    int I3 = cri->term.i3x3;
    int I4 = cri->term.i3x4;
    if (hasActiveAtomMask &&
        !(bitvectorActiveAtomMask->testBit(I1/3) &&
          bitvectorActiveAtomMask->testBit(I2/3) &&
          bitvectorActiveAtomMask->testBit(I3/3) &&
          bitvectorActiveAtomMask->testBit(I4/3))) continue;
    localPos[0]  = (*nvPosition)[I1+0];
    localPos[1]  = (*nvPosition)[I1+1];
    localPos[2]  = (*nvPosition)[I1+2];
    localPos[3]  = (*nvPosition)[I2+0];
    localPos[4]  = (*nvPosition)[I2+1];
    localPos[5]  = (*nvPosition)[I2+2];
    localPos[6]  = (*nvPosition)[I3+0];
    localPos[7]  = (*nvPosition)[I3+1];
    localPos[8]  = (*nvPosition)[I3+2];
    localPos[9]  = (*nvPosition)[I4+0];
    localPos[10] = (*nvPosition)[I4+1];
    localPos[11] = (*nvPosition)[I4+2];
    // Pack the four atoms' 3 coords into a contiguous 12-vector at offsets 0, 3, 6 and 9.
    // Zero local accumulator buffers (kernel uses += into them).
    for (size_t k = 0; k < PS; ++k) {
      localForce[k] = 0.0; localDvec[k] = 0.0; localHdvec[k] = 0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) localHess[k] = 0.0;
    double energyAccum = 0.0;
    chiral_term lt(cri->term.k, cri->term.co, 0, 3, 6, 9);
    chiral.hessian(lt, localPos, &energyAccum, localForce, localHess,
                   localDvec, localHdvec);
    // Scatter the 12x12 local Hessian into the global sparse matrix.
    // Local offsets 0..2 map to global I1+0..I1+2, 3..5 to I2, 6..8 to I3 and 9..11 to I4.
    // The kernel writes BOTH (r,c) and (c,r) for off-diagonal entries
    // (energyComponent.h:269-270), so accumulate only the lower triangle here.
    int globalIdx[PS] = { I1+0, I1+1, I1+2, I2+0, I2+1, I2+2,
                          I3+0, I3+1, I3+2, I4+0, I4+1, I4+2 };
    for (size_t r = 0; r < PS; ++r) {
      for (size_t c = 0; c <= r; ++c) {
        double v = localHess[r*PS + c];
        if (v != 0.0) {
          m->addToElement(globalIdx[r], globalIdx[c], v);
        }
      }
    }
  }
}





double EnergyChiralRestraint_O::evaluateAllComponent( ScoringFunction_sp score,
                                                      chem::NVector_sp 	pos,
                                                      core::T_sp energyScale,
                                                      core::T_sp energyComponents,
                                                      bool 		calcForce,
                                                      gc::Nilable<chem::NVector_sp> 	force,
                                                      bool		calcDiagonalHessian,
                                                      bool		calcOffDiagonalHessian,
                                                      gc::Nilable<chem::AbstractLargeSquareMatrix_sp>	hessian,
                                                      gc::Nilable<chem::NVector_sp>	hdvec,
                                                      gc::Nilable<chem::NVector_sp> dvec,
                                                      core::T_sp activeAtomMask,
                                                      core::T_sp debugInteractions )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  this->_Evaluations++;
  double totalEnergy = 0.0;
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);

  auto evalType = determineEnergyComponentEvalType(force,hdvec,dvec);

  double Energy = 0.0;  // referenced by CHIRAL_RESTRAINT_DEBUG_INTERACTIONS macro
  double* position = &(*pos)[0];
  double* rforce = NULL;
  double* rdvec = NULL;
  double* rhdvec = NULL;
  Chiral_Restraint<NoHessian> chiral;
  gctools::Vec0<EnergyChiralRestraint>::iterator cri;
  if (evalType==energyEval) {
    for ( cri=this->_Terms.begin(); cri!=this->_Terms.end(); cri++ ) {
      if (hasActiveAtomMask &&
          !(bitvectorActiveAtomMask->testBit(cri->term.i3x1/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x2/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x3/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x4/3))) continue;
      Energy = chiral.energy(cri->term, position, &totalEnergy);
      CHIRAL_RESTRAINT_DEBUG_INTERACTIONS(cri->term.i3x1, cri->term.i3x2, cri->term.i3x3, cri->term.i3x4);
    }
  } else if (evalType==gradientEval) {
    rforce = &(*force)[0];
    for ( cri=this->_Terms.begin(); cri!=this->_Terms.end(); cri++ ) {
      if (hasActiveAtomMask &&
          !(bitvectorActiveAtomMask->testBit(cri->term.i3x1/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x2/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x3/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x4/3))) continue;
      Energy = chiral.gradient(cri->term, position, &totalEnergy, rforce);
      CHIRAL_RESTRAINT_DEBUG_INTERACTIONS(cri->term.i3x1, cri->term.i3x2, cri->term.i3x3, cri->term.i3x4);
    }
  } else { // hessianEval
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for ( cri=this->_Terms.begin(); cri!=this->_Terms.end(); cri++ ) {
      if (hasActiveAtomMask &&
          !(bitvectorActiveAtomMask->testBit(cri->term.i3x1/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x2/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x3/3) &&
            bitvectorActiveAtomMask->testBit(cri->term.i3x4/3))) continue;
      Energy = chiral.hessian(cri->term, position, &totalEnergy, rforce, NoHessian(), rdvec, rhdvec);
      CHIRAL_RESTRAINT_DEBUG_INTERACTIONS(cri->term.i3x1, cri->term.i3x2, cri->term.i3x3, cri->term.i3x4);
    }
  }
  maybeSetEnergy( energyComponents, EnergyChiralRestraint_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}







void	EnergyChiralRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
		chem::NVector_sp 	pos)
{
  constexpr size_t PS = Chiral_Restraint<double*>::PositionSize;  // 12
  const double energyTol = 1.0e-7;
  const double forceTol  = 1.0e-5;
  const double hessTol   = 1.0e-3;
  double localPos[PS];
  double force_a[PS], force_fd[PS];
  double hess_a[PS*PS], hess_fd[PS*PS];
  double dvec_a[PS], dvec_fd[PS];
  double hdvec_a[PS], hdvec_fd[PS];
  Chiral_Restraint<double*> chiral;
  int fails = 0;
  int idx = 0;
  for ( auto cri=this->_Terms.begin(); cri!=this->_Terms.end(); cri++, idx++ ) {
    int I1 = cri->term.i3x1;
    int I2 = cri->term.i3x2;
    int I3 = cri->term.i3x3;
    int I4 = cri->term.i3x4;
    localPos[0]  = pos->getElement(I1+0);
    localPos[1]  = pos->getElement(I1+1);
    localPos[2]  = pos->getElement(I1+2);
    localPos[3]  = pos->getElement(I2+0);
    localPos[4]  = pos->getElement(I2+1);
    localPos[5]  = pos->getElement(I2+2);
    localPos[6]  = pos->getElement(I3+0);
    localPos[7]  = pos->getElement(I3+1);
    localPos[8]  = pos->getElement(I3+2);
    localPos[9]  = pos->getElement(I4+0);
    localPos[10] = pos->getElement(I4+1);
    localPos[11] = pos->getElement(I4+2);
    for (size_t k = 0; k < PS; ++k) {
      force_a[k]=0.0; force_fd[k]=0.0;
      dvec_a[k]=0.0;  dvec_fd[k]=0.0;
      hdvec_a[k]=0.0; hdvec_fd[k]=0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) {
      hess_a[k]=0.0; hess_fd[k]=0.0;
    }
    double e_a = 0.0, e_fd = 0.0;
    chiral_term lt(cri->term.k, cri->term.co, 0, 3, 6, 9);
    chiral.hessian(   lt, localPos, &e_a,  force_a,  hess_a,  dvec_a,  hdvec_a);
    chiral.hessian_fd(lt, localPos, &e_fd, force_fd, hess_fd, dvec_fd, hdvec_fd);
    bool termFails = false;
    if (std::fabs(e_a - e_fd) > energyTol) {
      termFails = true;
      core::lisp_write(fmt::format(
        "  energy mismatch term {}: ana={} fd={} diff={}\n",
        idx, e_a, e_fd, e_a - e_fd));
    }
    for (size_t k = 0; k < PS; ++k) {
      if (std::fabs(force_a[k] - force_fd[k]) > forceTol) {
        termFails = true;
        core::lisp_write(fmt::format(
          "  force mismatch term {} dof {}: ana={} fd={} diff={}\n",
          idx, k, force_a[k], force_fd[k], force_a[k] - force_fd[k]));
      }
    }
    for (size_t r = 0; r < PS; ++r) {
      for (size_t c = 0; c < PS; ++c) {
        double a = hess_a[r*PS+c];
        double f = hess_fd[r*PS+c];
        if (std::fabs(a - f) > hessTol) {
          termFails = true;
          core::lisp_write(fmt::format(
            "  hessian mismatch term {} ({},{}): ana={} fd={} diff={}\n",
            idx, r, c, a, f, a - f));
        }
      }
    }
    if (termFails) {
      fails++;
      core::lisp_write(fmt::format(
        "  -> term {} I1={} I2={} I3={} I4={} K={} CO={}\n",
        idx, I1, I2, I3, I4, cri->term.k, cri->term.co));
    }
  }
  core::lisp_write(fmt::format(
    "compareAnalyticalAndNumericalForceAndHessianTermByTerm: {} fails out of {} terms\n",
    fails, idx));
}





void EnergyChiralRestraint_O::initialize()
{
    this->Base::initialize();
}

size_t EnergyChiralRestraint_O::runTestCalls(core::T_sp stream, chem::NVector_sp coords) const {
#define POS_SIZE 12
  double energy_new;
  double energy_ground;
  double position[POS_SIZE];
  double force_new[POS_SIZE];
  double force_ground[POS_SIZE];
  double hessian_new[POS_SIZE*POS_SIZE];
  double hessian_ground[POS_SIZE*POS_SIZE];
  double dvec_new[POS_SIZE];
  double dvec_ground[POS_SIZE];
  double hdvec_new[POS_SIZE];
  double hdvec_ground[POS_SIZE];
  size_t idx=0;
  size_t errs = 0;
  Chiral_Restraint<double*> chiral;
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    position[0]  = coords[si->term.i3x1];
    position[1]  = coords[si->term.i3x1+1];
    position[2]  = coords[si->term.i3x1+2];
    position[3]  = coords[si->term.i3x2];
    position[4]  = coords[si->term.i3x2+1];
    position[5]  = coords[si->term.i3x2+2];
    position[6]  = coords[si->term.i3x3];
    position[7]  = coords[si->term.i3x3+1];
    position[8]  = coords[si->term.i3x3+2];
    position[9]  = coords[si->term.i3x4];
    position[10] = coords[si->term.i3x4+1];
    position[11] = coords[si->term.i3x4+2];
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    chiral_term lt(si->term.k, si->term.co, 0, 3, 6, 9);
    chiral.gradient(    lt, position, &energy_new,    force_new );
    chiral.gradient_fd( lt, position, &energy_ground, force_ground );
    if (!test_match( stream, "chiral_gradient", POS_SIZE,
                     force_new, force_ground,
                     0, 0,
                     0, 0 )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH chiral_gradient #{} K = {}  CO = {}\n",
                              idx, si->term.k, si->term.co ), stream );
    }
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    chiral_term lt2(si->term.k, si->term.co, 0, 3, 6, 9);
    chiral.hessian(    lt2, position, &energy_new,    force_new,    hessian_new,    dvec_new,    hdvec_new );
    chiral.hessian_fd( lt2, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    if (!test_match( stream, "chiral_hessian", POS_SIZE,
                     force_new, force_ground,
                     hessian_new, hessian_ground,
                     hdvec_new, hdvec_ground )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH chiral_hessian #{} K = {}  CO = {}\n",
                              idx, si->term.k, si->term.co ), stream );
    }
    idx++;
  }
  core::print(fmt::format("chiral errors = {}\n", errs), stream);
  return errs;
#undef POS_SIZE
}

};
