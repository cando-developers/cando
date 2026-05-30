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

#define ENERGY_FUNCTION I1, I2, activeAtomMask

#include <cando/chem/energyStretch.h>
#include <clasp/core/numerics.h>
#include <clasp/core/ql.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/array.h>
#include <cando/geom/ovector3.h>
#include <cando/chem/matter.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


namespace chem {
CL_LISPIFY_NAME(make-energy-stretch);
CL_DEF_CLASS_METHOD
EnergyStretch_sp EnergyStretch_O::make(EnergyFunction_sp energyFunction) {
  auto component = ensureComponent<EnergyStretch_O>(energyFunction);
  return component;
}


#if 1
# include "cando/chem/energyKernels/stretch.c"
#else
# include "cando/chem/energyKernels/mm_stretch.c"
#endif

}

namespace chem {
void EnergyStretch_O::emitTestCalls(core::T_sp stream, chem::NVector_sp coords) const
{
  size_t idx=0;
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    Vector3 p1(coords,si->term.i3x1,Safe());
    Vector3 p2(coords,si->term.i3x2,Safe());
    core::print(fmt::format("POSITION stretch 6 {} {} {} {} {} {}\n",
                            p1.getX(), p1.getY(), p1.getZ(),
                            p2.getX(), p2.getY(), p2.getZ() ), stream );
    core::print(fmt::format("TEST stretch 6 {} {} {} 0 3\n",
                            idx,
                            si->term.kb, si->term.r0 ), stream );
    idx++;
  }
}

size_t EnergyStretch_O::runTestCalls(core::T_sp stream, chem::NVector_sp coords) const
{
  #define POS_SIZE 6
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
  Stretch<double*> stretch;
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    position[0] = coords[si->term.i3x1];
    position[1] = coords[si->term.i3x1+1];
    position[2] = coords[si->term.i3x1+2];
    position[3] = coords[si->term.i3x2];
    position[4] = coords[si->term.i3x2+1];
    position[5] = coords[si->term.i3x2+2];
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    TermStretch lt(si->term.kb, si->term.r0, 0, 3);
    stretch.gradient(    lt, position, &energy_new,    force_new );
    stretch.gradient_fd( lt, position, &energy_ground, force_ground );
    if (!test_match( stream, "stretch_force", POS_SIZE,
                     force_new, force_ground,
                     0, 0,
                     0, 0 )) {
      errs++;
      if (stream.notnilp()) {
        test_position( stream, POS_SIZE, position );
        core::print(fmt::format("MISMATCH stretch_gradient #{} kb = {}  r0 = {}\n",
                                idx, si->term.kb, si->term.r0 ), stream );
      }
    }
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    TermStretch lt2(si->term.kb, si->term.r0, 0, 3);
    stretch.hessian(    lt2, position, &energy_new,    force_new,    hessian_new,    dvec_new,    hdvec_new );
    stretch.hessian_fd( lt2, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    if (!test_match( stream, "stretch_hessian", POS_SIZE,
                     force_new, force_ground,
                     hessian_new, hessian_ground,
                     hdvec_new, hdvec_ground )) {
      errs++;
      if (stream.notnilp()) {
        test_position( stream, POS_SIZE, position );
        core::print(fmt::format("MISMATCH stretch_hessian #{} kb = {}  r0 = {}\n",
                                idx, si->term.kb, si->term.r0 ), stream );
      }
    }
    idx++;
  }
  return errs;
}

};

namespace chem {

#define STRETCH_APPLY_ATOM_MASK(I1,I2) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3))) goto SKIP_term;

double _evaluateEnergyOnly_Stretch ( int I1,
                                     int I2,
                                     core::T_sp activeAtomMask,
                                     num_real x1,
                                     num_real y1,
                                     num_real z1,
                                     num_real x2,
                                     num_real y2,
                                     num_real z2,
                                     num_real r0,
                                     num_real kb )
{
  // I1, I2, and activeAtomMask are part of the public API for symmetry with
  // the gradient/hessian forms; they have no effect on a single-term energy
  // calculation, so suppress unused-parameter warnings.
  (void)I1; (void)I2; (void)activeAtomMask;
  double localPos[6] = { (double)x1, (double)y1, (double)z1,
                         (double)x2, (double)y2, (double)z2 };
  double energy_accum = 0.0;
  Stretch<NoHessian> stretch;
  TermStretch lt((double)kb, (double)r0, 0, 3);
  return stretch.energy(lt, localPos, &energy_accum);
}

CL_LAMBDA((energy-stretch chem:energy-stretch) pos &optional active-atom-mask);
CL_DEFMETHOD void	EnergyStretch_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm( NVector_sp pos, core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  constexpr size_t PS = Stretch<double*>::PositionSize;  // 6
  const double energyTol = 1.0e-7;
  const double forceTol  = 1.0e-5;
  const double hessTol   = 1.0e-3;
  double localPos[PS];
  double force_a[PS], force_fd[PS];
  double hess_a[PS*PS], hess_fd[PS*PS];
  double dvec_a[PS], dvec_fd[PS];
  double hdvec_a[PS], hdvec_fd[PS];
  Stretch<double*> stretch;
  int fails = 0;
  int idx = 0;
  for ( auto si=this->_Terms.begin(); si!=this->_Terms.end(); si++, idx++ ) {
    int I1 = si->term.i3x1;
    int I2 = si->term.i3x2;
    if (hasActiveAtomMask &&
        !(bitvectorActiveAtomMask->testBit(I1/3) &&
          bitvectorActiveAtomMask->testBit(I2/3))) continue;
    localPos[0] = pos->getElement(I1+0);
    localPos[1] = pos->getElement(I1+1);
    localPos[2] = pos->getElement(I1+2);
    localPos[3] = pos->getElement(I2+0);
    localPos[4] = pos->getElement(I2+1);
    localPos[5] = pos->getElement(I2+2);
    for (size_t k = 0; k < PS; ++k) {
      force_a[k]=0.0; force_fd[k]=0.0;
      dvec_a[k]=0.0;  dvec_fd[k]=0.0;
      hdvec_a[k]=0.0; hdvec_fd[k]=0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) {
      hess_a[k]=0.0; hess_fd[k]=0.0;
    }
    double e_a = 0.0, e_fd = 0.0;
    TermStretch lt(si->term.kb, si->term.r0, 0, 3);
    stretch.hessian(   lt, localPos, &e_a,  force_a,  hess_a,  dvec_a,  hdvec_a);
    stretch.hessian_fd(lt, localPos, &e_fd, force_fd, hess_fd, dvec_fd, hdvec_fd);
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
        "  -> term {} I1={} I2={} kb={} r0={}\n",
        idx, I1, I2, si->term.kb, si->term.r0));
    }
  }
  core::lisp_write(fmt::format(
    "compareAnalyticalAndNumericalForceAndHessianTermByTerm: {} fails out of {} terms\n",
    fails, idx));
}

#define STRETCH_DEBUG_INTERACTIONS(I1,I2) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,EnergyStretch_O::static_classSymbol(), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2)); \
    }

core::List_sp EnergyStretch::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kb),core::clasp_make_double_float(this->term.kb)),
                                  core::Cons_O::create(INTERN_(kw,r0),core::clasp_make_double_float(this->term.r0)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.i3x1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.i3x2)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                  core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2));
}

void EnergyStretch::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergyStretch");
}



void EnergyStretch::defineFrom( FFStretch_sp stretch, EnergyAtom *ea1, EnergyAtom *ea2, double scale)
{
  LOG("Defining EnergyStretch with kb={} r0={}"
      , this->term.kb
      , this->term.r0  );
//    this->_K3 = stretch->_K3;
//    this->_K4 = stretch->_K4;
//    this->_Bci = stretch->_Bci;
  this->_Atom1 = ea1->atom();
  this->_Atom2 = ea2->atom();
  LOG(" Adding stretch between {} - {}"
      , this->_Atom1->description() , this->_Atom2->description() );
  this->term.i3x1 = ea1->coordinateIndexTimes3();
  this->term.i3x2 = ea2->coordinateIndexTimes3();
  this->term.kb = stretch->getKb_kcalPerAngstromSquared();
  this->term.r0 = stretch->getR0_Angstrom();
}

std::string EnergyStretch_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}


 
void	EnergyStretch::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml, AtomTable_sp atomTable )
{
  int	i1, i2;
  i1 = xml->getAttributeInt("idx1");
  i2 = xml->getAttributeInt("idx2");
  LOG("parse i1={} i2={}" , i1 , i2  );
  this->_Atom1 = atomTable->findEnergyAtomWithCoordinateIndex(i1)->atom();
  this->_Atom2 = atomTable->findEnergyAtomWithCoordinateIndex(i2)->atom();
  LOG("_Atom1 = {}  _Atom2 = {}" , i1 , i2  );
  this->term.r0 = xml->getAttributeDouble("r0");
  this->term.kb = xml->getAttributeDouble("kb");
  this->term.i3x1 = xml->getAttributeInt("I1");
  this->term.i3x2 = xml->getAttributeInt("I2");
}



double	EnergyStretch::getR()
{
  Vector3	pos1, pos2, diff;
  pos1 = this->getAtom1()->getPosition();
  pos2 = this->getAtom2()->getPosition();
  diff = pos1.sub(pos2);
  return diff.length();
}


void	EnergyStretch_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                    AbstractLargeSquareMatrix_sp m,
                                                    core::T_sp activeAtomMask)
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  Stretch<double*> stretch;
  constexpr size_t PS = Stretch<double*>::PositionSize;  // 6
  double localPos[PS];
  double localForce[PS];
  double localHess[PS*PS];
  double localDvec[PS];
  double localHdvec[PS];
  for ( gctools::Vec0<EnergyStretch>::iterator si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    {
      int I1 = si->term.i3x1;
      int I2 = si->term.i3x2;
      STRETCH_APPLY_ATOM_MASK(I1, I2);
      localPos[0] = (*nvPosition)[I1+0];
      localPos[1] = (*nvPosition)[I1+1];
      localPos[2] = (*nvPosition)[I1+2];
      localPos[3] = (*nvPosition)[I2+0];
      localPos[4] = (*nvPosition)[I2+1];
      localPos[5] = (*nvPosition)[I2+2];
      // Pack the two atoms' 3 coords into a contiguous 6-vector at offsets 0 and 3.
      // Zero local accumulator buffers (kernel uses += into them).
      for (size_t k = 0; k < PS; ++k) {
        localForce[k] = 0.0; localDvec[k] = 0.0; localHdvec[k] = 0.0;
      }
      for (size_t k = 0; k < PS*PS; ++k) localHess[k] = 0.0;
      double energyAccum = 0.0;
      TermStretch lt(si->term.kb, si->term.r0, 0, 3);
      stretch.hessian(lt, localPos, &energyAccum, localForce, localHess,
                      localDvec, localHdvec);
      // Scatter the 6x6 local Hessian into the global sparse matrix.
      // Local offsets 0..2 map to global I1+0..I1+2 and 3..5 to I2+0..I2+2.
      // The kernel writes BOTH (r,c) and (c,r) for off-diagonal entries
      // (energyComponent.h:269-270), so accumulate only the lower triangle here.
      int globalIdx[PS] = { I1+0, I1+1, I1+2, I2+0, I2+1, I2+2 };
      for (size_t r = 0; r < PS; ++r) {
        for (size_t c = 0; c <= r; ++c) {
          double v = localHess[r*PS + c];
          if (v != 0.0) {
            m->addToElement(globalIdx[r], globalIdx[c], v);
          }
        }
      }
    }
    SKIP_term:;
  }
}


CL_DEFMETHOD core::T_sp EnergyStretch_O::stretchTermBetweenAtoms(Atom_sp x, Atom_sp y)
{
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    if ( (si->_Atom1 == x && si->_Atom2 == y) ||
         (si->_Atom1 == y && si->_Atom2 == x ) ) {
      return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kb), core::clasp_make_double_float(si->term.kb)),
                                      core::Cons_O::create(INTERN_(kw,r0), core::clasp_make_double_float(si->term.r0)));
    }
  }
  return nil<core::T_O>();
}


double EnergyStretch_O::evaluateAllComponent( ScoringFunction_sp score,
                                              NVector_sp 	pos,
                                              core::T_sp energyScale,
                                              core::T_sp energyComponents,
                                              bool 		calcForce,
                                              gc::Nilable<NVector_sp> 	force,
                                              bool		calcDiagonalHessian,
                                              bool		calcOffDiagonalHessian,
                                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                              gc::Nilable<NVector_sp>	hdvec,
                                              gc::Nilable<NVector_sp>   dvec,
                                              core::T_sp activeAtomMask ,
                                              core::T_sp debugInteractions )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  double totalEnergy = 0.0;
  this->_Evaluations++;
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);

  auto evalType = determineEnergyComponentEvalType(force,hdvec,dvec);

  int i;
  double Energy = 0.0;  // referenced by STRETCH_DEBUG_INTERACTIONS macro
  gctools::Vec0<EnergyStretch>::iterator si;
  DOUBLE* position = &(*pos)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  Stretch<NoHessian> stretch;
  if (evalType==energyEval) {
    // energy only
    for ( i=0,si=this->_Terms.begin(); si!=this->_Terms.end(); si++,i++ ) {
      if (hasActiveAtomMask &&
          !(bitvectorActiveAtomMask->testBit(si->term.i3x1/3) &&
            bitvectorActiveAtomMask->testBit(si->term.i3x2/3))) continue;
      Energy = stretch.energy(si->term, position, &totalEnergy);
      STRETCH_DEBUG_INTERACTIONS(si->term.i3x1,si->term.i3x2);
    }
  } else if (evalType==gradientEval) {
    rforce = &(*force)[0];
    for ( i=0,si=this->_Terms.begin(); si!=this->_Terms.end(); si++,i++ ) {
      if (hasActiveAtomMask &&
          !(bitvectorActiveAtomMask->testBit(si->term.i3x1/3) &&
            bitvectorActiveAtomMask->testBit(si->term.i3x2/3))) continue;
      Energy = stretch.gradient(si->term, position, &totalEnergy, rforce);
      STRETCH_DEBUG_INTERACTIONS(si->term.i3x1,si->term.i3x2);
    }
  } else { // hessianEval
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for ( i=0,si=this->_Terms.begin(); si!=this->_Terms.end(); si++,i++ ) {
      if (hasActiveAtomMask &&
          !(bitvectorActiveAtomMask->testBit(si->term.i3x1/3) &&
            bitvectorActiveAtomMask->testBit(si->term.i3x2/3))) continue;
      Energy = stretch.hessian(si->term, position, &totalEnergy, rforce, NoHessian(), rdvec, rhdvec);
      STRETCH_DEBUG_INTERACTIONS(si->term.i3x1,si->term.i3x2);
    }
  }
  maybeSetEnergy( energyComponents, EnergyStretch_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}


SYMBOL_EXPORT_SC_(KeywordPkg,stretch);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);
SYMBOL_EXPORT_SC_(KeywordPkg,r0);
SYMBOL_EXPORT_SC_(KeywordPkg,indexes);
SYMBOL_EXPORT_SC_(KeywordPkg,stretch_deviation);
SYMBOL_EXPORT_SC_(KeywordPkg,force);





void	EnergyStretch_O::dumpTerms(core::HashTable_sp atomTypes)
{
  gctools::Vec0<EnergyStretch>::iterator	esi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4;
  LOG("Dumping EnergyStretch terms");
  uint idx;
  for ( idx = 0, esi=this->_Terms.begin(); esi!=this->_Terms.end(); esi++, idx++ )
  {
    as1 = atomLabel(esi->_Atom1);
    as2 = atomLabel(esi->_Atom2);
    if ( as1 < as2 ) {
      str1 = as1;
      str2 = as2;
    } else {
      str2 = as1;
      str1 = as2;
    }
    core::clasp_write_string(fmt::format("TERM 1BND {:<9} - {:<9} {:8.2f} {:8.2f}\n" , str1 , str2 , esi->term.kb , esi->term.r0 ));
    LOG("TERM[#{}] 1BND %-9s - %-9s %8.2lf %8.2lf ; I1={} I2={}"
        , idx
        , str1 
        , str2
        , esi->term.kb
        , esi->term.r0
        , esi->term.i3x1
        , esi->term.i3x2 );
  }
}










void EnergyStretch_O::initialize()
{
  this->Base::initialize();
}

void EnergyStretch_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

void EnergyStretch_O::addTerm(const EnergyStretch& term)
{
  this->_Terms.push_back(term);
}

SYMBOL_EXPORT_SC_(KeywordPkg,kb);
SYMBOL_EXPORT_SC_(KeywordPkg,r0);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,atom1);
SYMBOL_EXPORT_SC_(KeywordPkg,atom2);

CL_DEFMETHOD core::List_sp EnergyStretch_O::extract_vectors_as_alist() const {
  size_t size = this->_Terms.size();
  printf("%s:%d The number of EnergyStretch terms -> %lu\n", __FILE__, __LINE__, size);
  NVector_sp kb_vec = NVector_O::make(size);
  NVector_sp r0_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_sp atom1_vec    = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom2_vec    = core::SimpleVector_O::make(size);
  for (size_t i=0;i<size;++i) {
    const EnergyStretch& entry = this->_Terms[i];
    (*kb_vec)[i] = entry.term.kb;
    (*r0_vec)[i] = entry.term.r0;
    (*i1_vec)[i] = entry.term.i3x1;
    (*i2_vec)[i] = entry.term.i3x2;
    (*atom1_vec)[i] = entry._Atom1 ;
    (*atom2_vec)[i] = entry._Atom2;
  }
  return core::Cons_O::createList(core::Cons_O::create(kw::_sym_kb,kb_vec),
                                  core::Cons_O::create(kw::_sym_r0,r0_vec),
                                  core::Cons_O::create(kw::_sym_i1,i1_vec),
                                  core::Cons_O::create(kw::_sym_i2,i2_vec),
                                  core::Cons_O::create(kw::_sym_atom1,atom1_vec),
                                  core::Cons_O::create(kw::_sym_atom2,atom2_vec));
};

CL_DEFMETHOD void EnergyStretch_O::fill_from_vectors_in_alist(core::List_sp vectors)
{
  NVector_sp kb_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_kb));
  NVector_sp r0_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_r0));
  core::SimpleVector_int32_t_sp i1_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i1));
  core::SimpleVector_int32_t_sp i2_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i2));
  core::SimpleVector_sp atom1_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom1));
  core::SimpleVector_sp atom2_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom2));
  this->_Terms.resize(kb_vec->length());
  for (size_t i=0,iEnd(kb_vec->length());i<iEnd;++i) {
    EnergyStretch& entry = this->_Terms[i];
    entry.term.kb = (*kb_vec)[i];
    entry.term.r0 = (*r0_vec)[i];
    entry.term.i3x1 = (*i1_vec)[i];
    entry.term.i3x2 = (*i2_vec)[i];
    entry._Atom1 = gc::As_unsafe<Atom_sp>((*atom1_vec)[i]);
    entry._Atom2 = gc::As_unsafe<Atom_sp>((*atom2_vec)[i]);
  }
}

CL_DEFMETHOD void EnergyStretch_O::addStretchTerm(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, double kb, double r0) {
  EnergyAtom* ea1 = atomTable->getEnergyAtomPointer(a1);
  EnergyAtom* ea2 = atomTable->getEnergyAtomPointer(a2);
  EnergyStretch energyStretch(a1,a2,ea1->coordinateIndexTimes3(),ea2->coordinateIndexTimes3(),kb,r0);
  this->addTerm(energyStretch);
}

CL_DOCSTRING(R"dx(Invoke a callback for every term in the energy-stretch.
The callback takes the lambda-list (index atom1 atom2 atom1-index*3 atom2-index*3 kb r0).
Both kb and r0 are single-floats to avoid consing.
Use the modify-stretch-term-kb and modify-stretch-term-r0 to modify the parameters.)dx")
CL_DEFMETHOD void EnergyStretch_O::walkStretchTerms(core::T_sp callback)
{
  for (size_t i=0;i<this->_Terms.size();++i) {
    const EnergyStretch& entry = this->_Terms[i];
    core::eval::funcall(callback,core::make_fixnum(i),
                        entry._Atom1,
                        entry._Atom2,
                        core::make_fixnum(entry.term.i3x1),
                        core::make_fixnum(entry.term.i3x2),
                        core::make_single_float(entry.term.kb),
                        core::make_single_float(entry.term.r0));
  };
};

CL_DEFMETHOD void EnergyStretch_O::modifyStretchTermKb(size_t index, float kb) {
  if (index<this->_Terms.size()) {
    EnergyStretch& entry = this->_Terms[index];
    entry.term.kb = kb;
    return;
  }
  SIMPLE_ERROR("index {} is out of bounds as a energy-stretch term index (#entries {})" , index , this->_Terms.size() );
}

CL_DEFMETHOD void EnergyStretch_O::modifyStretchTermR0(size_t index, float r0) {
  if (index<this->_Terms.size()) {
    EnergyStretch& entry = this->_Terms[index];
    entry.term.r0 = r0;
    return;
  }
  SIMPLE_ERROR("index {} is out of bounds as a energy-stretch term index (#entries {})" , index , this->_Terms.size() );
}

void EnergyStretch_O::reset()
{
  this->_Terms.clear();
}

CL_DEFMETHOD
core::List_sp	EnergyStretch_O::lookupStretchTerms(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, core::HashTable_sp atomTypes)
{
  ql::list  result;
  core::T_sp tia1 = atomTable->_AtomTableIndexes->gethash(a1);
  core::T_sp tia2 = atomTable->_AtomTableIndexes->gethash(a2);
  if (!tia1.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a1));
  if (!tia2.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a2));
  for (auto edi=this->_Terms.begin();edi!=this->_Terms.end();edi++) {
    if ((edi->_Atom1==a1 &&
         edi->_Atom2==a2)
        || (edi->_Atom2==a1 &&
            edi->_Atom1==a2)) {
      ql::list oneResult;
      oneResult << INTERN_(kw,type1) << edi->_Atom1->getType(atomTypes)
                << INTERN_(kw,type2) << edi->_Atom2->getType(atomTypes)
                << INTERN_(kw,kb) << core::clasp_make_double_float(edi->term.kb)
                << INTERN_(kw,r0) << core::clasp_make_double_float(edi->term.r0);
      result << oneResult.result();
    }
  }
  return result.result();
}

EnergyComponent_sp EnergyStretch_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyStretch_O::staticClass() );
  EnergyStretch_sp copy = EnergyStretch_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    if ( skipInteraction( keepInteraction, a1, a2 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}


};




