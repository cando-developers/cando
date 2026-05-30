/*
    File: energyAnchorRestraint.cc
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

#include <cando/chem/energyAnchorRestraint.h>
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
#include <clasp/core/evaluator.h>
#include <clasp/core/wrappers.h>

namespace chem {

#include "cando/chem/energyKernels/anchor_restraint.c"

#define ANCHOR_RESTRAINT_DEBUG_INTERACTIONS(I1) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,EnergyAnchorRestraint_O::static_classSymbol(), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1)); \
    }

EnergyAnchorRestraint::EnergyAnchorRestraint()
{
  this->_Atom1 = nil<Atom_O>();
}

EnergyAnchorRestraint::~EnergyAnchorRestraint()
{
}

core::List_sp EnergyAnchorRestraint::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,ka),core::clasp_make_double_float(this->term.ka)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.i3x1)),
                                  core::Cons_O::create(INTERN_(kw,xa), core::clasp_make_double_float(this->term.xa)),
                                  core::Cons_O::create(INTERN_(kw,ya), core::clasp_make_double_float(this->term.ya)),
                                  core::Cons_O::create(INTERN_(kw,za), core::clasp_make_double_float(this->term.za)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1));
}


void EnergyAnchorRestraint_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

#ifdef XML_ARCHIVE
void	EnergyAnchorRestraint::archive(core::ArchiveP node)
{
  node->attribute("ka",this->term.ka);
  node->attribute("xa",this->term.xa);
  node->attribute("ya",this->term.ya);
  node->attribute("za",this->term.za);
  node->attribute("I1",this->term.i3x1);
  node->attribute("a1",this->_Atom1);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
  node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
  node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
  node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include <cando/chem/energy_functions/_AnchorRestraint_debugEvalSerialize.cc>
#endif //]
}
#endif

#if 0
adapt::QDomNode_sp	EnergyAnchorRestraint::asXml()
{
  adapt::QDomNode_sp	node,child;
  Vector3	vdiff;

  node = adapt::QDomNode_O::create(env,"EnergyAnchorRestraint");
  node->addAttributeString("atomName",this->_Atom1->getName());
  node->addAttributeInt("I1",this->term.i3x1);
  node->addAttributeDoubleScientific("ka",this->term.ka);
  node->addAttributeDoubleScientific("xa",this->term.xa);
  node->addAttributeDoubleScientific("ya",this->term.ya);
  node->addAttributeDoubleScientific("za",this->term.za);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
  adapt::QDomNode_sp xml = adapt::QDomNode_O::create(env,"Evaluated");
  xml->addAttributeBool("calcForce",this->_calcForce );
  xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
  xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include <_AnchorRestraint_debugEvalXml.cc>
  node->addChild(xml);
#endif
  return node;
}

void	EnergyAnchorRestraint::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
                                                          AtomTable_sp at )
{
  this->term.ka = xml->getAttributeDouble("ka");
  this->term.i3x1 = xml->getAttributeInt("I1");
  this->term.xa = xml->getAttributeDouble("xa");
  this->term.ya = xml->getAttributeDouble("ya");
  this->term.za = xml->getAttributeDouble("za");
  this->_Atom1 = at->findEnergyAtomWithCoordinateIndex(this->term.i3x1)->atom();
}
#endif


double	_evaluateEnergyOnly_AnchorRestraint(
                                                    num_real x1, num_real y1, num_real z1,
                                                    num_real xa, num_real ya, num_real za,
                                                    num_real ka, core::T_sp activeAtomMask )
{
  // activeAtomMask is part of the public API for symmetry with the
  // gradient/hessian forms; it has no effect on a single-term energy
  // calculation, so suppress unused-parameter warnings.
  (void)activeAtomMask;
  double localPos[3] = { (double)x1, (double)y1, (double)z1 };
  double energy_accum = 0.0;
  Anchor_Restraint<NoHessian> anchor;
  anchor_term lt((double)ka, (double)xa, (double)ya, (double)za, 0);
  return anchor.energy(lt, localPos, &energy_accum);
}

std::string EnergyAnchorRestraint_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}


CL_LISPIFY_NAME(make-energy-anchor-restraint);
CL_DEF_CLASS_METHOD
EnergyAnchorRestraint_sp EnergyAnchorRestraint_O::make(EnergyFunction_sp energyFunction) {
  return ensureComponent<EnergyAnchorRestraint_O>(energyFunction);
}

void	EnergyAnchorRestraint_O::addTerm(const EnergyAnchorRestraint& r)
{
  this->_Terms.push_back(r);
}

CL_DEFMETHOD
size_t EnergyAnchorRestraint_O::addAnchorRestraintTerm(EnergyFunction_sp energyFunction, Atom_sp a1, Vector3 position, double forceConstant ) {
  AtomTable_sp atomTable = energyFunction->atomTable();
  EnergyAtom* ea = atomTable->getEnergyAtomPointer(a1);
  //EnergyAnchorRestraint restraint(a1,position,ea->coordinateIndexTimes3(),forceConstant);
  this->_Terms.emplace_back(a1,position,ea->coordinateIndexTimes3(),forceConstant);
  return this->_Terms.size()-1;
}



void	EnergyAnchorRestraint_O::dumpTerms(core::HashTable_sp atomTypes)
{
}



void	EnergyAnchorRestraint_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                            AbstractLargeSquareMatrix_sp m,
                                                            core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  Anchor_Restraint<double*> anchor;
  constexpr size_t PS = Anchor_Restraint<double*>::PositionSize;  // 3
  double localPos[PS];
  double localForce[PS];
  double localHess[PS*PS];
  double localDvec[PS];
  double localHdvec[PS];
  for ( gctools::Vec0<EnergyAnchorRestraint>::iterator cri=this->_Terms.begin();
        cri!=this->_Terms.end(); cri++ ) {
    int I1 = cri->term.i3x1;
    if (hasActiveAtomMask && !bitvectorActiveAtomMask->testBit(I1/3)) continue;
    localPos[0] = (*nvPosition)[I1+0];
    localPos[1] = (*nvPosition)[I1+1];
    localPos[2] = (*nvPosition)[I1+2];
    // The single anchored atom's 3 coords occupy a contiguous 3-vector at offset 0.
    // Zero local accumulator buffers (kernel uses += into them).
    for (size_t k = 0; k < PS; ++k) {
      localForce[k] = 0.0; localDvec[k] = 0.0; localHdvec[k] = 0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) localHess[k] = 0.0;
    double energyAccum = 0.0;
    anchor_term lt(cri->term.ka, cri->term.xa, cri->term.ya, cri->term.za, 0);
    anchor.hessian(lt, localPos, &energyAccum, localForce, localHess,
                   localDvec, localHdvec);
    // Scatter the 3x3 local Hessian into the global sparse matrix.
    // Local offsets 0..2 map to global I1+0..I1+2.
    // The kernel writes BOTH (r,c) and (c,r) for off-diagonal entries
    // (energyComponent.h:269-270), so accumulate only the lower triangle here.
    int globalIdx[PS] = { I1+0, I1+1, I1+2 };
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





double EnergyAnchorRestraint_O::evaluateAllComponent( ScoringFunction_sp score,
                                                      NVector_sp 	pos,
                                                      core::T_sp energyScale,
                                                      core::T_sp energyComponents,
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
  double totalEnergy = 0.0;
  this->_Evaluations++;
  auto evalType = determineEnergyComponentEvalType(force,hdvec,dvec);

  double Energy = 0.0;  // referenced by ANCHOR_RESTRAINT_DEBUG_INTERACTIONS macro
  double* position = &(*pos)[0];
  double* rforce = NULL;
  double* rdvec = NULL;
  double* rhdvec = NULL;
  Anchor_Restraint<NoHessian> anchor;
  gctools::Vec0<EnergyAnchorRestraint>::iterator cri;
  if (evalType==energyEval) {
    for ( cri=this->_Terms.begin(); cri!=this->_Terms.end(); cri++ ) {
      if (hasActiveAtomMask && !bitvectorActiveAtomMask->testBit(cri->term.i3x1/3)) continue;
      Energy = anchor.energy(cri->term, position, &totalEnergy);
      ANCHOR_RESTRAINT_DEBUG_INTERACTIONS(cri->term.i3x1);
    }
  } else if (evalType==gradientEval) {
    rforce = &(*force)[0];
    for ( cri=this->_Terms.begin(); cri!=this->_Terms.end(); cri++ ) {
      if (hasActiveAtomMask && !bitvectorActiveAtomMask->testBit(cri->term.i3x1/3)) continue;
      Energy = anchor.gradient(cri->term, position, &totalEnergy, rforce);
      ANCHOR_RESTRAINT_DEBUG_INTERACTIONS(cri->term.i3x1);
    }
  } else { // hessianEval
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for ( cri=this->_Terms.begin(); cri!=this->_Terms.end(); cri++ ) {
      if (hasActiveAtomMask && !bitvectorActiveAtomMask->testBit(cri->term.i3x1/3)) continue;
      Energy = anchor.hessian(cri->term, position, &totalEnergy, rforce, NoHessian(), rdvec, rhdvec);
      ANCHOR_RESTRAINT_DEBUG_INTERACTIONS(cri->term.i3x1);
    }
  }
  maybeSetEnergy( energyComponents, EnergyAnchorRestraint_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}







void	EnergyAnchorRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp 	pos,
                                                                                        core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  constexpr size_t PS = Anchor_Restraint<double*>::PositionSize;  // 3
  const double energyTol = 1.0e-7;
  const double forceTol  = 1.0e-5;
  const double hessTol   = 1.0e-3;
  double localPos[PS];
  double force_a[PS], force_fd[PS];
  double hess_a[PS*PS], hess_fd[PS*PS];
  double dvec_a[PS], dvec_fd[PS];
  double hdvec_a[PS], hdvec_fd[PS];
  Anchor_Restraint<double*> anchor;
  int fails = 0;
  int idx = 0;
  for ( auto cri=this->_Terms.begin(); cri!=this->_Terms.end(); cri++, idx++ ) {
    int I1 = cri->term.i3x1;
    if (hasActiveAtomMask && !bitvectorActiveAtomMask->testBit(I1/3)) continue;
    localPos[0] = pos->getElement(I1+0);
    localPos[1] = pos->getElement(I1+1);
    localPos[2] = pos->getElement(I1+2);
    for (size_t k = 0; k < PS; ++k) {
      force_a[k]=0.0; force_fd[k]=0.0;
      dvec_a[k]=0.0;  dvec_fd[k]=0.0;
      hdvec_a[k]=0.0; hdvec_fd[k]=0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) {
      hess_a[k]=0.0; hess_fd[k]=0.0;
    }
    double e_a = 0.0, e_fd = 0.0;
    anchor_term lt(cri->term.ka, cri->term.xa, cri->term.ya, cri->term.za, 0);
    anchor.hessian(   lt, localPos, &e_a,  force_a,  hess_a,  dvec_a,  hdvec_a);
    anchor.hessian_fd(lt, localPos, &e_fd, force_fd, hess_fd, dvec_fd, hdvec_fd);
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
        "  -> term {} I1={} ka={} xa={} ya={} za={}\n",
        idx, I1, cri->term.ka, cri->term.xa, cri->term.ya, cri->term.za));
    }
  }
  core::lisp_write(fmt::format(
    "compareAnalyticalAndNumericalForceAndHessianTermByTerm: {} fails out of {} terms\n",
    fails, idx));
}



void EnergyAnchorRestraint_O::initialize()
{
  this->Base::initialize();
}

EnergyComponent_sp EnergyAnchorRestraint_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  EnergyAnchorRestraint_sp copy = EnergyAnchorRestraint_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory,EnergyAnchorRestraint_O::staticClass());
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    if ( skipInteraction( keepInteraction, a1 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}

core::T_sp test_allocate_EnergyAnchorRestraint_O() {
#if 1
  // succeeds
  auto thing = EnergyAnchorRestraint_O::create();
#elif 1
  // succeeds ---- I got it to fail once I think but it may have been that I didn't save 
  auto thing = gctools::GC<EnergyAnchorRestraint_O>::allocate();
#elif 0
  // FAILS
  auto kind = gctools::Header_s::StampWtagMtag::make_StampWtagMtag(EnergyAnchorRestraint_O::static_ValueStampWtagMtag);
  auto thing = gctools::GCObjectAllocator<EnergyAnchorRestraint_O>::template allocate_kind<gctools::RuntimeStage>(kind,gctools::sizeof_with_header<EnergyAnchorRestraint_O>());
#else
  // FAILS
  gctools::Header_s::BadgeStampWtagMtag the_header;
  auto thing = gctools::GCObjectAppropriatePoolAllocator<EnergyAnchorRestraint_O,gctools::normal>::allocate_in_appropriate_pool_kind<gctools::RuntimeStage>(the_header,9999);
#endif
  return thing;
}

size_t EnergyAnchorRestraint_O::runTestCalls(core::T_sp stream, chem::NVector_sp coords) const {
#define POS_SIZE 3
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
  Anchor_Restraint<double*> anchor;
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    position[0] = coords[si->term.i3x1];
    position[1] = coords[si->term.i3x1+1];
    position[2] = coords[si->term.i3x1+2];
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    anchor_term lt(si->term.ka, si->term.xa, si->term.ya, si->term.za, 0);
    anchor.gradient(    lt, position, &energy_new,    force_new );
    anchor.gradient_fd( lt, position, &energy_ground, force_ground );
    if (!test_match( stream, "anchor_gradient", POS_SIZE,
                     force_new, force_ground,
                     0, 0,
                     0, 0 )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH anchor_gradient #{} ka = {}\n", idx, si->term.ka ), stream );
    }
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    anchor_term lt2(si->term.ka, si->term.xa, si->term.ya, si->term.za, 0);
    anchor.hessian(    lt2, position, &energy_new,    force_new,    hessian_new,    dvec_new,    hdvec_new );
    anchor.hessian_fd( lt2, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    if (!test_match( stream, "anchor_hessian", POS_SIZE,
                     force_new, force_ground,
                     hessian_new, hessian_ground,
                     hdvec_new, hdvec_ground )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH anchor_hessian #{} ka = {}\n", idx, si->term.ka ), stream );
    }
    idx++;
  }
  core::print(fmt::format("anchor errors = {}\n", errs), stream);
  return errs;
#undef POS_SIZE
}

};
