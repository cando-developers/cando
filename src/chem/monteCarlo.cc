/*
    File: monteCarlo.cc
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

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <random>

#include <clasp/core/common.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/array_int32.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/evaluator.h>
#include <clasp/llvmo/intrinsics.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/monteCarlo.h>

namespace chem {
#include <iostream>

struct KahanAccumulator {
public:
  KahanAccumulator() : sum(0.0), c(0.0) {}

  void Add(double value) {
    double y = value - c;    // Reduced precision in value
    double t = sum + y;      // Sum is accumulated in t
    c = (t - sum) - y;       // (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
    sum = t;                 // Algebraically, c should always be zero. Beware overly aggressive optimizing compilers!
  }

  double GetSum() const {
    return sum;
  }

public:
  double sum; // Accumulated sum
  double c;   // A running compensation for lost low-order bits.
};


SYMBOL_EXPORT_SC_( ChemPkg, intramolecular_backbone_energy );
SYMBOL_EXPORT_SC_( ChemPkg, intermolecular_backbone_energy );
SYMBOL_EXPORT_SC_( ChemPkg, monomer_vector );
SYMBOL_EXPORT_SC_( ChemPkg, monomer_corrections );
SYMBOL_EXPORT_SC_( ChemPkg, mrkey_index_to_lmkey_index );
SYMBOL_EXPORT_SC_( ChemPkg, monomer_locus_max_mrkindex );
SYMBOL_EXPORT_SC_( ChemPkg, intramolecular_single_scan_energy );
SYMBOL_EXPORT_SC_( ChemPkg, intermolecular_single_scan_energy );
SYMBOL_EXPORT_SC_( ChemPkg, pair_scan_energy_lower_triangular_matrix );
SYMBOL_EXPORT_SC_( ChemPkg, intermolecular_p_pair_scan_energy_lower_triangular_matrix );
SYMBOL_EXPORT_SC_( ChemPkg, lmkey_index_to_lmkey_info );

struct Energies;

struct State {
  typedef std::vector<uint32_t> StateType;
  StateType _State;

  State() {};

  // Use a defined state
  State(core::SimpleVector_byte32_t_sp mcstate) {
    size_t len = mcstate->length();
    this->resize(len);
    for ( size_t ii=0; ii<len; ii++ ) {
      this->_State[ii] = (*mcstate)[ii];
    }
  }

  void resize(size_t sz) {
    this->_State.resize(sz);
  }

  void randomState(const Energies& energies );

  void verifyState(core::SimpleVector_byte32_t_sp sv_monomerLocusMaxMrkindex ) {
    int prevMax = -1;
    uint32_t* monomerLocusMaxMrkindex = &(*sv_monomerLocusMaxMrkindex)[0]; 
    for ( size_t ii=0; ii<sv_monomerLocusMaxMrkindex->length(); ii++ ) {
      int max = (*sv_monomerLocusMaxMrkindex)[ii];
      int val = this->_State[ii];
      if (prevMax < val && val <= max ) continue;
      SIMPLE_ERROR("Invalidated state: {}", this->stateAsString());
    }
  }

  std::string stateAsString() {
    stringstream ss;
    ss << "#<mcstate ";
    for ( size_t ii = 0; ii<this->_State.size()-1; ii++ ) {
      ss << this->_State[ii];
      ss << " ";
    }
    ss << ">";
    return ss.str();
  }
};



struct Energies {
  std::mt19937                   _Rng;
  std::uniform_real_distribution<double> _U01;
  std::uniform_int_distribution<size_t>  _PickSlot;
  std::uniform_int_distribution<size_t>  _Rand;
  vector<size_t>                 _SlotIndexes;
  vector<size_t>                 _MonomerIndexes;
  vector<size_t>                 _LmkeyIndexToLocus;
  vector<std::string>            _LmkeyIndexToMonomerName;
  vector<size_t>                 _LocusMonomerCount;
  size_t                         _NumberOfSlots;
  double                         _TemperatureScale;
  double                         _IntramolecularBackboneEnergy;
  double                         _IntermolecularBackboneEnergy;
  core::SimpleVector_sp          _MonomerVector;
  core::SimpleVector_byte32_t_sp _MonomerLocusMaxMrkindex;
  core::SimpleVector_byte32_t_sp _MrkeyIndexToLmkeyIndex;
  core::SimpleVector_double_sp   _MonomerCorrectionsOriginal;
  core::SimpleVector_double_sp   _MonomerCorrections;
  size_t                         _MonomerCorrectionsLength;
  core::SimpleVector_double_sp   _IntramolecularSingleTerms;
  core::SimpleVector_double_sp   _IntermolecularSingleTerms;
  core::SimpleVector_double_sp   _PairTerms;
  core::SimpleVector_byte8_t_sp    _IntermolecularPPairTerms;
  Energies(core::T_sp energies,double temperatureScale = 0.95) : _Rng(12345u), _U01(0.0,1.0), _Rand(RAND_MAX) {
    this->_TemperatureScale = temperatureScale;
    core::T_sp val;
    val = core::eval::funcall( _sym_intramolecular_backbone_energy, energies );
    this->_IntramolecularBackboneEnergy = gc::As<core::DoubleFloat_sp>(val)->get();
    val = core::eval::funcall( _sym_intermolecular_backbone_energy, energies );
    this->_IntermolecularBackboneEnergy = gc::As<core::DoubleFloat_sp>(val)->get();
    val = core::eval::funcall( _sym_monomer_vector, energies );
    this->_MonomerVector = gc::As<core::SimpleVector_sp>(val);
    val = core::eval::funcall( _sym_monomer_locus_max_mrkindex, energies );
    this->_MonomerLocusMaxMrkindex = gc::As<core::SimpleVector_byte32_t_sp>(val);
    val = core::eval::funcall( _sym_mrkey_index_to_lmkey_index, energies );
    this->_MrkeyIndexToLmkeyIndex = gc::As<core::SimpleVector_byte32_t_sp>(val);
    this->_NumberOfSlots = core::cl__length(this->_MonomerLocusMaxMrkindex);
    val = core::eval::funcall( _sym_monomer_corrections, energies );
    this->_MonomerCorrectionsOriginal = gc::As<core::SimpleVector_double_sp>(val);
    this->_MonomerCorrections = gc::As<core::SimpleVector_double_sp>(core::cl__copy_seq(this->_MonomerCorrectionsOriginal));
    this->_MonomerCorrectionsLength = core::cl__length(this->_MonomerCorrections);
    val = core::eval::funcall( _sym_intramolecular_single_scan_energy, energies );
    this->_IntramolecularSingleTerms = gc::As<core::SimpleVector_double_sp>(val);
    val = core::eval::funcall( _sym_intermolecular_single_scan_energy, energies );
    this->_IntermolecularSingleTerms = gc::As<core::SimpleVector_double_sp>(val);
    val = core::eval::funcall( _sym_pair_scan_energy_lower_triangular_matrix, energies );
    this->_PairTerms = gc::As<core::SimpleVector_double_sp>(val);
    val = core::eval::funcall( _sym_intermolecular_p_pair_scan_energy_lower_triangular_matrix, energies );
    this->_IntermolecularPPairTerms = gc::As<core::SimpleVector_byte8_t_sp>(val);
    std::uniform_int_distribution<size_t> pickSlot(0,this->_NumberOfSlots-1);
    this->_PickSlot = pickSlot;
    this->_LocusMonomerCount.resize(this->_NumberOfSlots,0);
    for ( size_t ii=0; ii<this->_MonomerCorrectionsLength; ii++ ) {
      core::T_mv tinfo = core::eval::funcall(chem::_sym_lmkey_index_to_lmkey_info,energies,core::make_fixnum(ii));
      core::MultipleValues& mvn = core::lisp_multipleValues();
      core::T_sp tlocus = mvn.valueGet(1,tinfo.number_of_values());
      core::Symbol_sp sym = gc::As<core::Symbol_sp>(tinfo);
      std::string monomerName = sym->symbolNameAsString();
      size_t locus = core::clasp_to_fixnum(tlocus);
      this->_LocusMonomerCount[locus]++;
      this->_LmkeyIndexToLocus.push_back(locus);
      this->_LmkeyIndexToMonomerName.push_back(monomerName);
    }
  }

  int lowerTriangularIndex(int xx, int yy ) {
    if (yy < xx) {
      int tt = xx;
      xx = yy;
      yy = tt;
    }
    int index = (yy*(yy-1)/2)+xx;
    return index;
  }
#if 0
  __attribute__((optnone))
  double oldEnergyFunction(State& state) {
    double* singleTerms = &(*this->_SingleTerms)[0];
    double* pairTerms = &(*this->_PairTerms)[0];
    KahanAccumulator singleSum;
    for ( int ii=0; ii<state._State.size(); ii++ ) {
      int index = state._State[ii];
      double energyTerm = singleTerms[index];
      singleSum.Add(energyTerm);
      uint32_t lmkey_index = (*this->_MrkeyIndexToLmkeyIndex)[index];
      double correction = this->_MonomerCorrections[lmkey_index];
      singleSum.Add(correction);
    }
    KahanAccumulator pairSum;
    for ( int xii=0; xii<state._State.size()-1; xii++ ) {
      int xxindex = state._State[xii];
      for ( int yii=xii+1; yii<state._State.size(); yii++ ) {
        int yyindex = state._State[yii];
        int ltmIndex = this->lowerTriangularIndex(xxindex,yyindex);
        double energyTerm = pairTerms[ltmIndex];
        pairSum.Add(energyTerm);
      }
    }
    KahanAccumulator total;
    total.Add(this->_IntramolecularBackboneEnergy);
    total.Add(this->_IntermolecularBackboneEnergy);
    total.Add(singleSum.sum);
    total.Add(pairSum.sum);
    return total.sum;
  }
#endif

  double physicalEnergy(const State& state, double lambda=1.0);
  double reducedEnergy(State& state,double lambda=1.0);

  State randomStep( size_t& slotIndex, const State& state) {
    slotIndex = this->_PickSlot(this->_Rng);
    int max = (*this->_MonomerLocusMaxMrkindex)[slotIndex];
    int prevMax = (slotIndex==0) ? -1 : (*this->_MonomerLocusMaxMrkindex)[slotIndex-1];
    int range = max - prevMax;
    size_t mrkeyIndex = (this->_Rand(this->_Rng) % range) + prevMax + 1;
    State tempState(state);
    tempState._State[slotIndex] = mrkeyIndex;
    return tempState;
  }

// Temperature decay function
  double tempDecay(const double& temperature) {
      // Replace with your temperature decay logic
    return temperature * this->_TemperatureScale;
  }

  double U01() {
    return this->_U01(this->_Rng);
  }
};

void State::randomState(const Energies& energies )
{
  this->_State.resize(energies._MonomerLocusMaxMrkindex->length());
  int prevMax = -1;
  uint32_t* monomerLocusMaxMrkindex = &(*energies._MonomerLocusMaxMrkindex)[0]; 
  for ( size_t ii=0; ii<this->_State.size(); ii++ ) {
    int max = monomerLocusMaxMrkindex[ii];
    int range = max - prevMax;
    this->_State[ii] = (rand() % range) + prevMax + 1;
    prevMax = max;
  }
}
CL_DEFUN core::DoubleFloat_sp chem__mcstate_energy(core::T_sp tmcstate, core::T_sp tenergies, double lambda)
{
  Energies energies(tenergies);
  State state(gc::As<core::SimpleVector_byte32_t_sp>(tmcstate));
  double testEnergy = energies.reducedEnergy(state,lambda);
  return core::DoubleFloat_O::create(testEnergy);
}

SYMBOL_EXPORT_SC_( ChemPkg, max_iterations );
SYMBOL_EXPORT_SC_( ChemPkg, wl_increment_stop );
SYMBOL_EXPORT_SC_( ChemPkg, min_temperature );
SYMBOL_EXPORT_SC_( ChemPkg, max_failed_improvements );
SYMBOL_EXPORT_SC_( ChemPkg, max_failed_accepts );

CL_LAMBDA(energies &key seed (start-temperature 1000.0) (max-iterations 1000000) (epoch-steps 100) (temperature-scale 0.95) initial-state accept-callback temperature-drop-callback debug (max-failed-improvements 100000) (max-failed-accepts 100000) (min-temperature 0.001) );
CL_DEFUN core::T_mv chem__simulatedAnnealing(core::T_sp tenergies, core::T_sp seed, double startTemperature, size_t max_iterations, size_t epoch_steps, double temperatureScale, core::T_sp tinitial_state, core::T_sp acceptCallback, core::T_sp temperatureDropCallback, core::T_sp debug, size_t maxFailedImprovements, size_t maxFailedAccepts, double minTemperature ) {
  //  srand(time(NULL)); // Initialize random seed
  Energies energies(tenergies,temperatureScale);
  if (seed.notnilp()) {
    uint32_t s = core::clasp_to_uint32_t(seed);
    energies._Rng.seed(s);
    srand(s);
  }
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  size_t failedAccepts = 0;
  size_t failedImprovements = 0;
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  bool useTemperatureDropCallback = temperatureDropCallback.notnilp();
  if (useAcceptCallback||useTemperatureDropCallback) {
    saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  }
  State initial_state ;
  initial_state.resize(energies._NumberOfSlots);
  if (tinitial_state.nilp()) {
    initial_state.randomState(energies);
  } else {
    auto passed = gc::As<core::SimpleVector_byte32_t_sp>(tinitial_state);
    for ( size_t ii=0; ii<energies._NumberOfSlots; ii++ ) {
      initial_state._State[ii] = (*passed)[ii];
    }
    initial_state.verifyState(energies._MonomerLocusMaxMrkindex);
  }
  State currentState = initial_state;
  double temperature = startTemperature;
  double beta = 1.0/(0.0019872*temperature);
  double curEnergy = energies.reducedEnergy(currentState);
  State testState;
  size_t steps = 0;
  size_t accepts = 0;
  core::T_sp stopReason = nil<core::T_O>();
  for (int i = 0; i < max_iterations; ++i) {
    beta = 1.0/(0.0019872*temperature);
    for (size_t epochi=0; epochi<epoch_steps; epochi++ ) {
      size_t slotIndex;
      testState = energies.randomStep( slotIndex, currentState );
      double testEnergy = energies.reducedEnergy(testState);
      double delta = beta * (testEnergy - curEnergy);
      double x = - delta;
      double p = (x >= 0.0) ? 1.0 : std::exp(x);   // never exp(large +)
      if (debug.notnilp()) core::clasp_write_string(fmt::format("acc_prob test: testE={:.6} prevE={:.6} beta={:.3} delta={:.6} p={:.3}\n", testEnergy, curEnergy, beta, delta, p ),debug);
      double rnd = energies.U01();
      steps++;
      if ( rnd <= p) {
        accepts++;
        failedAccepts = 0;
        if (testEnergy<curEnergy) {
          failedImprovements = 0;
        } else {
          failedImprovements++;
        }
        currentState = testState;
        curEnergy = testEnergy;
        if (debug.notnilp()) core::clasp_write_string(fmt::format(" +++ accepted step\n"),debug);
        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentState._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, saveState, mk_double_float(curEnergy) );
        }
      } else {
        failedAccepts++;
        failedImprovements++;
      }
      if (failedAccepts > maxFailedAccepts) {
        stopReason = chem::_sym_max_failed_accepts;
        goto DONE;
      }
      if (failedImprovements > maxFailedImprovements ) {
        stopReason = chem::_sym_max_failed_improvements;
        goto DONE;
      }
    }
    temperature = energies.tempDecay(temperature);
    if (useTemperatureDropCallback) {
      memcpy(&(*saveState)[0],&currentState._State[0],sizeof(int32_t)*energies._NumberOfSlots);
      core::eval::funcall( temperatureDropCallback, saveState, core::make_fixnum(steps), mk_double_float(curEnergy), mk_double_float(temperature) );
    }
    if ( temperature < minTemperature ) {
      stopReason = chem::_sym_min_temperature;
      goto DONE;
    }
  }
  stopReason = chem::_sym_max_iterations;
 DONE:
  return Values(core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots,
                                                    0, false,
                                                    energies._NumberOfSlots, // initialContentsSize
                                                    &currentState._State[0] ), // initialContents
                core::make_fixnum(steps),
                mk_double_float(curEnergy),
                mk_double_float(temperature),
                stopReason);
}


CL_LAMBDA(energies &key (temperature 100.0) (max-iterations 1000) initial-state step-callback debug);
CL_DEFUN core::T_mv chem__constantTemperatureMonteCarlo(core::T_sp tenergies, double temperature, size_t max_iterations, core::T_sp tinitial_state, core::T_sp stepCallback, core::T_sp debug ) {
  srand(time(NULL)); // Initialize random seed
  double beta = 1.0/(0.0019872*temperature);
  size_t accepts = 0;
  Energies energies(tenergies);
  core::SimpleVector_byte32_t_sp saveState;
  bool useStepCallback = stepCallback.notnilp();
  if (useStepCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  State initial_state ;
  initial_state.resize(energies._NumberOfSlots);
  if (tinitial_state.nilp()) {
    initial_state.randomState(energies);
  } else {
    auto passed = gc::As<core::SimpleVector_byte32_t_sp>(tinitial_state);
    for ( size_t ii=0; ii<energies._NumberOfSlots; ii++ ) {
      initial_state._State[ii] = (*passed)[ii];
    }
    initial_state.verifyState(energies._MonomerLocusMaxMrkindex);
  }
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  State currentState = initial_state;
  double curEnergy = energies.reducedEnergy(currentState);
  State testState;
  for (int i = 0; i < max_iterations; ++i) {
    size_t slotIndex;
    testState = energies.randomStep( slotIndex, currentState );
    double testEnergy = energies.reducedEnergy(testState);
    double delta = beta * (testEnergy - curEnergy);
    double x = - delta;
    double p = (x >= 0.0) ? 1.0 : std::exp(x);   // never exp(large +)
    if (debug.notnilp()) core::clasp_write_string(fmt::format("acc_prob test: testE={:.6} prevE={:.6} beta={:.3} delta={:.6} p={:.3}\n", testEnergy, curEnergy, beta, delta, p ),debug);
    double rnd = energies.U01();
    if ( rnd <= p) {
      accepts++;
      currentState = testState;
      curEnergy = testEnergy;
      if (debug.notnilp()) core::clasp_write_string(fmt::format(" +++ accepted step\n"),debug);
    }
    if (useStepCallback) {
      memcpy(&(*saveState)[0],&currentState._State[0],sizeof(int32_t)*energies._NumberOfSlots);
      core::eval::funcall( stepCallback, saveState, core::make_fixnum(slotIndex), mk_double_float(curEnergy) );
    }
  }
  return Values(core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots,
                                                    0, false,
                                                    energies._NumberOfSlots, // initialContentsSize
                                                    &currentState._State[0] ), // initialContents
                mk_double_float(curEnergy) );
}

CL_LAMBDA(energies beta mcstate)
CL_DEFUN core::T_sp chem__mcstate_energy(core::T_sp tenergies, double beta, core::SimpleVector_byte32_t_sp mcstate ) {
  Energies energies(tenergies);
  State state(mcstate);
  double energy = energies.reducedEnergy(state);
  return core::DoubleFloat_O::create(energy);
}

inline bool flat_enough(core::T_sp debug,
                        const std::vector<size_t>& histogram,
                        const std::vector<size_t>& rawHistogram,
                        const std::vector<size_t>& reachable,
                        size_t bi, double flatness ) {
  size_t hmin = SIZE_MAX;
  size_t hmax = 0.0;
  size_t visits = 0;
  size_t rawVisits = 0;
  size_t nonzero = 0;
  size_t rc = 0;
  for (size_t ii = 0; ii < histogram.size(); ++ii ) {
    if (!reachable[ii]) continue;
    rc++;
    size_t h = histogram[ii];
    visits += h;
    size_t rh = rawHistogram[ii];
    rawVisits += rh;
    if (h > 0) nonzero++;
    hmin = std::min(hmin,h);
    hmax = std::max(hmax,h);
  }
  if (rc == 0 || nonzero < rc) return false; // ensure all reachable bins seen
  size_t min_visits = 20 * rc;               // or scale with rc
  if (rawVisits >= min_visits) {
    double avg = double(visits) / rc;
    // SIMPLE_ERROR("Also check the hmax that hmax/avg < 1/flatness");
    bool flat = (double(hmin)/avg > flatness &&
                 double(hmax)/avg<(1.0/flatness));
    if (debug.notnilp()) core::clasp_write_string(fmt::format("Enough rawVisits for bi={} hmin = {} hmax = {} avg = {} hmin/avg={:.4} hmax/avg={:.4}  flatness = {:.4} 1.0/flatness = {:.4}  {}\n",
                                                              bi, hmin, hmax, avg, hmin/avg, hmax/avg, flatness, 1.0/flatness, flat ? "FLAT" : "not-flat" ),debug);
    if (flat) return true;
  } else {
    if (debug.notnilp()) core::clasp_write_string(fmt::format("Too early for flatness test for betai={} rawVisits={} rawVisits needed={}\n", bi, rawVisits, min_visits ),debug);
  }
  return false;
}

void adjustMonomerCorrectionsBias(Energies& energies)
{
  // Reset the bias
  double bias = std::numeric_limits<double>::max();
  for ( size_t ii = 0; ii<energies._MonomerCorrectionsLength; ii++ ) {
    bias = std::min((*energies._MonomerCorrections)[ii],bias);
  }
  for ( size_t ii = 0; ii<energies._MonomerCorrectionsLength; ii++ ) {
    (*energies._MonomerCorrections)[ii] -= bias;
  }
}

double Energies::physicalEnergy(const State& state,double lambda) {
  const double* intramolecularSingleTerms = &(*this->_IntramolecularSingleTerms)[0];
  const double* intermolecularSingleTerms = &(*this->_IntermolecularSingleTerms)[0];
  const double* pairTerms = &(*this->_PairTerms)[0];
  const uint8_t* intermolecularPPairTerms = &(*this->_IntermolecularPPairTerms)[0];

  KahanAccumulator singleSum, pairSum;

  for (int ii = 0; ii < state._State.size(); ++ii) {
    singleSum.Add(intramolecularSingleTerms[state._State[ii]]);
    singleSum.Add(intermolecularSingleTerms[state._State[ii]]*lambda);
  }

  for (int xii = 0; xii < state._State.size() - 1; ++xii) {
    int xx = state._State[xii];
    for (int yii = xii + 1; yii < state._State.size(); ++yii) {
      int yy = state._State[yii];
      size_t lti = this->lowerTriangularIndex(xx,yy);
      double pairTerm = pairTerms[lti];
      int8_t intermolecular_p = intermolecularPPairTerms[lti];
      pairSum.Add(intermolecular_p ? pairTerm*lambda : pairTerm );
    }
  }

  KahanAccumulator total;
  total.Add(this->_IntramolecularBackboneEnergy);
  total.Add(this->_IntermolecularBackboneEnergy*lambda);
  total.Add(singleSum.sum);
  total.Add(pairSum.sum);
  return total.sum;
}

double Energies::reducedEnergy(State& state, double lambda ) {
  double physicalEnergy = this->physicalEnergy(state,lambda);
  KahanAccumulator correctionSum;
  correctionSum.Add(physicalEnergy);
  for ( int ii=0; ii<state._State.size(); ii++ ) {
    int index = state._State[ii];
    uint32_t lmkey_index = (*this->_MrkeyIndexToLmkeyIndex)[index];
    //    double correction = this->_MonomerCorrections[lmkey_index]/beta;
    double correction = this->_MonomerCorrections[lmkey_index];
    correctionSum.Add(correction);
  }
  return correctionSum.sum;
}


inline double acc_prob_from_log(core::T_sp debug, double beta, double testEnergy, double curEnergy ) {
  // x = log(prob) = -βΔE  (Metropolis)  or  (βi-βj)(Ej-Ei) (REX)
  double delta = beta * (testEnergy - curEnergy);
  double x = - delta;
  double p = (x >= 0.0) ? 1.0 : std::exp(x);   // never exp(large +)
  if (debug.notnilp()) core::clasp_write_string(fmt::format("acc_prob test: testE={:.6} prevE={:.6} beta={:.3} delta={:.6} p={:.3}\n", testEnergy, curEnergy, beta, delta, p ),debug);
  return p;
}

inline double acc_beta_prob_from_log(core::T_sp debug, double betai,double u_i,double betaj,double u_j) {
  // double x  = (betas[i]-betas[j])*(Ej - Ei);
  // double p  = (x >= 0.0) ? 1.0 : std::exp(-x);
  double delta = (betai-betaj)*(u_j - u_i);
  double x = - delta;
  double p = (x >= 0.0) ? 1.0 : std::exp(x);   // never exp(large +)
  if (debug.notnilp()) core::clasp_write_string(fmt::format("acc_beta_prob test: delta={:.6} p={:.3}\n", delta, p ),debug);
  return p;
}

CL_LAMBDA(energies temperature &key seed (max-iterations 1000000) flatness-callback accept-callback (flatness-steps 100) (wl-scaling 0.8) (wl-increment-stop 0.1) (wl-increment-start 1.0) (flatness-threshold 0.95) debug);
CL_DEFUN core::T_mv chem__voelz_optimize_monomer_corrections_single_temperature(core::T_sp tenergies, double temperature, core::T_sp seed, size_t max_iterations, core::T_sp flatnessCallback, core::T_sp acceptCallback, size_t flatness_steps, double wl_scaling, double wl_increment_stop, double wl_increment_start, double flatness_threshold, core::T_sp debug ) {
  double wl_increment = wl_increment_start;
  double beta = 1.0/(0.0019872*temperature);
  Energies energies(tenergies);
  if (seed.notnilp()) {
    uint32_t s = core::clasp_to_uint32_t(seed);
    energies._Rng.seed(s);
    srand(s);
  }
  core::T_sp stop_reason = nil<core::T_O>();
  // Number of slots in state
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  if (useAcceptCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  State currentState;
  size_t slotIndex = 0;
  currentState.randomState(energies);
  State testState;
  std::vector<size_t> histogram(energies._MonomerCorrectionsLength,0);
  std::vector<size_t> rawHistogram(energies._MonomerCorrectionsLength,0);
  size_t flatness_iter = 0;
  size_t iterCount = 0;
  size_t accepts = 0;
  size_t flatness_reached = 0;

  // In chem__voelz_optimize_monomer_corrections_single_temperature
  std::vector<size_t> reachable(histogram.size(), 0);
  for (size_t ii = 0; ii < energies._MrkeyIndexToLmkeyIndex->length(); ++ii) {
    uint32_t lm = (*energies._MrkeyIndexToLmkeyIndex)[ii];
    reachable[lm] += 1;
  }
  size_t reachable_count = std::count(reachable.begin(), reachable.end(), 1);

  for (size_t iter = 0; iter < max_iterations; ++iter) {
    if (debug.notnilp()) core::clasp_write_string(fmt::format("============== Top of loop at iter {}\n", iter),debug);
    size_t testSlotIndex;
    double curEnergy = energies.reducedEnergy(currentState); // calculated using new MonomerCorrections
    testState = energies.randomStep( testSlotIndex, currentState );
    if (debug.notnilp()) {
      core::clasp_write_string(fmt::format("testState = "),debug);
      for ( size_t ii=0; ii<testState._State.size(); ii++ ) {
        core::clasp_write_string(fmt::format(" {:3}", testState._State[ii]),debug);
        if (ii==testSlotIndex) core::clasp_write_string("<",debug);
      }
      core::clasp_write_string(fmt::format("\n"),debug);
    }
    double testEnergy = energies.reducedEnergy(testState);
    double p = acc_prob_from_log( debug, beta, testEnergy, curEnergy );
    //      double x = -(*betas)[betai] * (testEnergy - curEnergy);
    //      double p = acc_prob_from_log(debug,x);
    double rnd = energies.U01();
    slotIndex = testSlotIndex;
    if ( rnd <= p) {
      accepts++;
      currentState = testState;
      if (debug.notnilp()) core::clasp_write_string(fmt::format(" +++ accepted step\n"),debug);
    }
    // Update the histogram and the monomer corrections
    {
      // Increment the histogram and the MonomerCorrection for the current monomer
      uint32_t curLmkeyIndex;
#define ALL_MONOMERS_UPDATED 0
      // Update monomerCorrections for only curLmkeyIndex
      {
        curLmkeyIndex = (*energies._MrkeyIndexToLmkeyIndex)[currentState._State[slotIndex]];
        size_t numberOfMonomers = energies._LocusMonomerCount[energies._LmkeyIndexToLocus[curLmkeyIndex]];
        histogram[curLmkeyIndex] += numberOfMonomers;
        rawHistogram[curLmkeyIndex]++;
        // double beta_increment = wl_increment*(*betas)[betai];
        double beta_increment = wl_increment;
        (*energies._MonomerCorrections)[curLmkeyIndex] += beta_increment;
        if (debug.notnilp())
          core::clasp_write_string(fmt::format(" >> curLmkeyIndex {} beta={:.4} wl_increment {:.4} beta_increment {:.4}\n",
                                               curLmkeyIndex, beta, wl_increment, beta_increment ),debug);
      }
      if (debug.notnilp()) {
        core::clasp_write_string(fmt::format(" Locus/monomer-name: " ),debug);
        for ( size_t lmkey_index=0; lmkey_index<histogram.size(); lmkey_index++ ) {
          std::string label = fmt::format("{}/{}", energies._LmkeyIndexToLocus[lmkey_index], energies._LmkeyIndexToMonomerName[lmkey_index]);
          core::print(fmt::format("{:>12}", label), debug);
          if (ALL_MONOMERS_UPDATED || lmkey_index==curLmkeyIndex) core::clasp_write_string("<",debug);
        }
        core::print("\n",debug);
        core::clasp_write_string(fmt::format("          Histogram: " ),debug);
        for ( size_t ii=0; ii<histogram.size(); ii++ ) {
          core::clasp_write_string(fmt::format(" {:11}", histogram[ii]),debug);
          if (ALL_MONOMERS_UPDATED || ii==curLmkeyIndex) core::clasp_write_string("<",debug);
        }
        core::clasp_write_string(fmt::format("\n"),debug);
        core::clasp_write_string(fmt::format("       rawHistogram: " ),debug);
        for ( size_t ii=0; ii<rawHistogram.size(); ii++ ) {
          core::clasp_write_string(fmt::format(" {:11}", rawHistogram[ii]),debug);
          if (ALL_MONOMERS_UPDATED || ii==curLmkeyIndex) core::clasp_write_string("<",debug);
        }
        core::clasp_write_string(fmt::format("\n"),debug);
        core::clasp_write_string(fmt::format("Monomer corrections: "),debug);
        for ( size_t ii=0; ii<energies._MonomerCorrectionsLength; ii++ ) {
          core::clasp_write_string(fmt::format("{:12.2f}", (*energies._MonomerCorrections)[ii]),debug);
          if (ALL_MONOMERS_UPDATED || ii==curLmkeyIndex) core::clasp_write_string("<",debug);
        }
        core::clasp_write_string(fmt::format("\n"),debug);
      }
      // Reset the bias
      adjustMonomerCorrectionsBias(energies);

      if (useAcceptCallback) {
        memcpy(&(*saveState)[0],&currentState._State[0],sizeof(int32_t)*energies._NumberOfSlots);
        core::eval::funcall( acceptCallback, core::make_fixnum(iter), saveState, mk_double_float(testEnergy), energies._MonomerCorrections );
      }
    }
    {
      // Check the histogram for flatness
      flatness_iter++;
      if (flatness_iter%flatness_steps==0) {
        size_t visits = 0, hmin = SIZE_MAX, hmax = 0;
        size_t flat_count = 0;
        size_t bi = 0;
        if ( flat_enough( debug,  histogram, rawHistogram, reachable, bi, flatness_threshold ) ) flat_count++;
        bi++;
        if (debug.notnilp()) core::clasp_write_string(fmt::format(" - - - - Flat_count = {}\n", flat_count ),debug);
        if (flat_count == 1) {
          flatness_reached++;
          if (debug.notnilp()) core::clasp_write_string(fmt::format("All histograms flat enough wl_increment updated to {:.4}\n", wl_increment ),debug);
          for ( size_t ii=0; ii<histogram.size(); ii++ ) {
            histogram[ii] = 0;
            rawHistogram[ii] = 0;
          }
          wl_increment *= wl_scaling;
          if (flatnessCallback.notnilp()) {
            core::eval::funcall( flatnessCallback, core::make_fixnum(iter), mk_double_float(wl_increment), energies._MonomerCorrections );
          }
          // Write the MonomerCorrections whenever flatness is reached
          memcpy(&(*energies._MonomerCorrectionsOriginal)[0],&(*energies._MonomerCorrections)[0],energies._MonomerCorrectionsLength*sizeof(double));
        }
      }
    }
    if (wl_increment<wl_increment_stop) {
      if (debug.notnilp()) core::clasp_write_string(fmt::format("wl_increment {:.6} is < wl_increment_stop {:.6}\n",
                                                                wl_increment, wl_increment_stop ),debug);
      stop_reason = _sym_wl_increment_stop;
      goto DONE;
    }
    iterCount++;
  }
  stop_reason = _sym_max_iterations;
 DONE:
  return Values(core::make_fixnum(iterCount),
                core::make_fixnum(accepts),
                core::make_fixnum(flatness_reached),
                mk_double_float(wl_increment),
                stop_reason
                );
}
#if 0
CL_LAMBDA(energies betas &key (temperature-swap-steps 100) (max-iterations 1000000) flatness-callback accept-callback (flatness-steps 100) (wl-scaling 0.8) (wl-increment-stop 0.1) (wl-increment-start 1.0) (flatness-threshold 0.95) debug);
CL_DEFUN core::T_mv chem__voelz_optimize_monomer_corrections(core::T_sp tenergies, core::SimpleVector_double_sp betas, size_t temperature_swap_steps, size_t max_iterations, core::T_sp flatnessCallback, core::T_sp acceptCallback, size_t flatness_steps, double wl_scaling, double wl_increment_stop, double wl_increment_start, double flatness_threshold, core::T_sp debug ) {
  size_t beta_size = core::cl__length(betas);
  double wl_increment = wl_increment_start;
  Energies energies(tenergies);
  core::T_sp stop_reason = nil<core::T_O>();
  // Number of slots in state
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  if (useAcceptCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  std::vector<State> currentStates(beta_size);
  for ( size_t ii = 0; ii<currentStates.size(); ii++ ) {
    currentStates[ii].randomState(energies);
  }
  State testState;
  std::vector<std::vector<size_t>> histograms(beta_size,std::vector<size_t>(energies._MonomerCorrectionsLength,0));
  size_t flatness_iter = 0;
  size_t iterCount = 0;
  size_t accepts = 0;
  size_t flatness_reached = 0;
  for (size_t iter = 0; iter < max_iterations; ++iter) {
    if (debug.notnilp()) core::clasp_write_string(fmt::format("============== Top of loop at iter {}\n", iter),debug);
    for ( size_t betai = 0; betai<beta_size; betai++ ) {
      double beta = (*betas)[betai];
      size_t slotIndex;
      double curEnergy = energies.reducedEnergy(currentStates[betai]); // calculated using new MonomerCorrections
      testState = energies.randomStep( slotIndex, currentStates[betai] );
      if (debug.notnilp()) {
        core::clasp_write_string(fmt::format("testState = "),debug);
        for ( size_t ii=0; ii<testState._State.size(); ii++ ) {
          core::clasp_write_string(fmt::format(" {:3}", testState._State[ii]),debug);
          if (ii==slotIndex) core::clasp_write_string("<",debug);
        }
        core::clasp_write_string(fmt::format("\n"),debug);
      }
      double testEnergy = energies.reducedEnergy(testState);
      double p = acc_prob_from_log( debug, (*betas)[betai], testEnergy, curEnergy );
      //      double x = -(*betas)[betai] * (testEnergy - curEnergy);
      //      double p = acc_prob_from_log(debug,x);
      double rnd = energies.U01();
      if ( rnd <= p) {
        accepts++;
        currentStates[betai] = testState;
        if (debug.notnilp()) core::clasp_write_string(fmt::format(" +++ accepted step\n"),debug);
      }
      // Update the histogram and the monomer corrections
      {
        // Increment the histogram and the MonomerCorrection for the current monomer
        uint32_t monomerIndex;
#define ALL_MONOMERS_UPDATED 0
        // Update monomerCorrections for only monomerIndex
        {
          monomerIndex = (*energies._MrkeyIndexToLmkeyIndex)[currentStates[betai]._State[slotIndex]];
          histogram[monomerIndex]++;
          // double beta_increment = wl_increment*(*betas)[betai];
          double beta_increment = wl_increment;
          (*energies._MonomerCorrections)[monomerIndex] += beta_increment;
          if (debug.notnilp())
            core::clasp_write_string(fmt::format(" >> betai {} monomerIndex {} beta={:.4} wl_increment {:.4} beta_increment {:.4}\n",
                                                 betai, monomerIndex, beta, wl_increment, beta_increment ),debug);
        }
        if (debug.notnilp()) {
          core::clasp_write_string(fmt::format("Histogram(betai={:2}): ", betai ),debug);
          for ( size_t ii=0; ii<histograms[betai].size(); ii++ ) {
            core::clasp_write_string(fmt::format(" {:8}", histograms[betai][ii]),debug);
            if (ALL_MONOMERS_UPDATED || ii==monomerIndex) core::clasp_write_string("<",debug);
          }
          core::clasp_write_string(fmt::format("\n"),debug);
          core::clasp_write_string(fmt::format("Monomer corrections: "),debug);
          for ( size_t ii=0; ii<energies._MonomerCorrectionsLength; ii++ ) {
            core::clasp_write_string(fmt::format(" {:8.2f}", (*energies._MonomerCorrections)[ii]),debug);
            if (ALL_MONOMERS_UPDATED || ii==monomerIndex) core::clasp_write_string("<",debug);
          }
          core::clasp_write_string(fmt::format("\n"),debug);
        }
        // Reset the bias
        adjustMonomerCorrectionsBias(energies);

        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentStates[betai]._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, core::make_fixnum(iter), core::make_fixnum(betai), saveState, mk_double_float(testEnergy), energies._MonomerCorrections );
        }
      }
    }
    {
      // Check the histogram for flatness
      flatness_iter++;
      if (flatness_iter%flatness_steps==0) {
        size_t visits = 0, hmin = SIZE_MAX, hmax = 0;
        size_t flat_count = 0;
        size_t bi = 0;
        for (auto const& histogram : histograms ) {
          if ( flat_enough( debug,  histogram, bi, flatness_threshold ) ) flat_count++;
          bi++;
        }
        if (debug.notnilp()) core::clasp_write_string(fmt::format(" - - - - Flat_count = {}  histograms.size() = {}\n", flat_count, histograms.size()),debug);
        if (flat_count == histograms.size()) {
          flatness_reached++;
          if (debug.notnilp()) core::clasp_write_string(fmt::format("All histograms flat enough wl_increment updated to {:.4}\n", wl_increment ),debug);
          for ( size_t hi=0; hi<histograms.size(); hi++ ) {
            for ( size_t ii=0; ii<histograms[hi].size(); ii++ ) {
              histograms[hi][ii] = 0;
            }
          }
          wl_increment *= wl_scaling;
          if (flatnessCallback.notnilp()) {
            core::eval::funcall( flatnessCallback, core::make_fixnum(iter), mk_double_float(wl_increment), energies._MonomerCorrections );
          }
          // Write the MonomerCorrections whenever flatness is reached
          memcpy(&(*energies._MonomerCorrectionsOriginal)[0],&(*energies._MonomerCorrections)[0],energies._MonomerCorrectionsLength*sizeof(double));
        }
      }
    }
    if (wl_increment<wl_increment_stop) {
      if (debug.notnilp()) core::clasp_write_string(fmt::format("wl_increment {:.6} is < wl_increment_stop {:.6}\n",
                                                                wl_increment, wl_increment_stop ),debug);
      stop_reason = _sym_wl_increment_stop;
      goto DONE;
    }
#if 0
    if (temperature_swap_steps && (currentStates.size()>=2) && (iter % temperature_swap_steps)==0) {
      // Try the temperature swap
      size_t i = tswapPickSlot(rng), j = i+1;
#if 1
      double u_i = energies.reducedEnergy(currentStates[i]);
      double u_j = energies.reducedEnergy(currentStates[j]);
      double p = acc_beta_prob_from_log(debug,(*betas)[i],u_i,(*betas)[j],u_j);
#else
      // ChatGPT says use this
      double E_i = energies.physicalEnergy(currentStates[i]);
      double E_j = energies.physicalEnergy(currentStates[j]);
      double p = acc_beta_prob_from_log(debug,(*betas)[i],E_i,(*betas)[j],E_j);
#endif
      double rnd = U01(rng);
      if (rnd < p) {
        if (debug.notnilp()) core::clasp_write_string(fmt::format(" <><><> Swapping betai states {} and {}\n", i, j ),debug);
        std::swap(currentStates[i], currentStates[j]);
      }
    }
#endif
    iterCount++;
  }
  stop_reason = _sym_max_iterations;
 DONE:
  return Values(core::make_fixnum(iterCount),
                core::make_fixnum(accepts),
                core::make_fixnum(flatness_reached),
                mk_double_float(wl_increment),
                stop_reason
                );
}
#endif

#if 0
CL_LAMBDA(energies betas &key (temperature-swap-steps 100) (max-iterations 1000000) flatness-callback accept-callback (flatness-steps 100) (wl-scaling 0.8) (wl-increment-stop 0.1) (wl-increment 20.0) (flatness 0.95) debug);
CL_DEFUN core::T_mv chem__voelz_optimize_monomer_corrections(core::T_sp tenergies, core::SimpleVector_double_sp betas, size_t temperature_swap_steps, size_t max_iterations, core::T_sp flatnessCallback, core::T_sp acceptCallback, size_t flatness_steps, double wl_scaling, double wl_increment_stop, double wl_increment, double flatness, core::T_sp debug ) {
  size_t beta_size = core::cl__length(betas);
  Energies energies(tenergies);
  // Number of slots in state
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  if (useAcceptCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  std::vector<State> currentStates(beta_size);
  for ( size_t ii = 0; ii<currentStates.size(); ii++ ) {
    currentStates[ii].randomState(energies);
  }
  State testState;
  std::vector<std::vector<size_t>> histograms(beta_size,std::vector<size_t>(energies._MonomerCorrectionsLength,0));
  size_t flatness_iter = 0;
  size_t iterCount = 0;
  size_t accepts = 0;
  size_t flatness_reached = 0;
  for (size_t iter = 0; iter < max_iterations; ++iter) {
    if (debug.notnilp()) core::clasp_write_string(fmt::format("============== Top of loop at iter {}\n", iter),debug);
    for ( size_t betai = 0; betai<beta_size; betai++ ) {
      double beta = (*betas)[betai];
      size_t slotIndex = pickSlot(rng);
      double curEnergy = energies.reducedEnergy(currentStates[betai]); // calculated using new MonomerCorrections
      testState = energies.randomStep( slotIndex, currentStates[betai] );
      if (debug.notnilp()) {
        core::clasp_write_string(fmt::format("testState = "),debug);
        for ( size_t ii=0; ii<testState._State.size(); ii++ ) {
          core::clasp_write_string(fmt::format(" {:3}", testState._State[ii]),debug);
          if (ii==slotIndex) core::clasp_write_string("<",debug);
        }
        core::clasp_write_string(fmt::format("\n"),debug);
      }
      double testEnergy = energies.reducedEnergy(testState);
      double delta = beta * (testEnergy - curEnergy);
      double x = - delta;
      double p = (x >= 0.0) ? 1.0 : std::exp(x);   // never exp(large +)
      if (debug.notnilp()) core::clasp_write_string(fmt::format("acc_prob test: testE={:.6} prevE={:.6} beta={:.3} delta={:.6} p={:.3}\n", testEnergy, curEnergy, beta, delta, p ),debug);
      double rnd = U01(rng);
      if ( rnd <= p) {
        accepts++;
        currentStates[betai] = testState;
        if (debug.notnilp()) core::clasp_write_string(fmt::format(" +++ accepted step\n"),debug);
      }
      // Update the histogram and the monomer corrections
      {
        // Increment the histogram and the MonomerCorrection for the current monomer
        uint32_t monomerIndex;
#define ALL_MONOMERS_UPDATED 0
        // Update monomerCorrections for only monomerIndex
        {
          monomerIndex = (*energies._MrkeyIndexToLmkeyIndex)[currentStates[betai]._State[slotIndex]];
          histograms[betai][monomerIndex]++;

          (*energies._MonomerCorrections)[monomerIndex] += wl_increment;
          if (debug.notnilp())
            core::clasp_write_string(fmt::format(" >> betai {} monomerIndex {} beta={:.4} wl_increment {:.4}\n",
                                                 betai, monomerIndex, beta, wl_increment ),debug);
        }
        if (debug.notnilp()) {
          core::clasp_write_string(fmt::format("Histogram(betai={:2}): ", betai ),debug);
          for ( size_t ii=0; ii<histograms[betai].size(); ii++ ) {
            core::clasp_write_string(fmt::format(" {:8}", histograms[betai][ii]),debug);
            if (ALL_MONOMERS_UPDATED || ii==monomerIndex) core::clasp_write_string("<",debug);
          }
          core::clasp_write_string(fmt::format("\n"),debug);
          core::clasp_write_string(fmt::format("Monomer corrections: "),debug);
          for ( size_t ii=0; ii<energies._MonomerCorrectionsLength; ii++ ) {
            core::clasp_write_string(fmt::format(" {:8.2f}", (*energies._MonomerCorrections)[ii]),debug);
            if (ALL_MONOMERS_UPDATED || ii==monomerIndex) core::clasp_write_string("<",debug);
          }
          core::clasp_write_string(fmt::format("\n"),debug);
        }
        // Reset the bias
        adjustMonomerCorrectionsBias(energies);

        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentStates[betai]._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, core::make_fixnum(iter), core::make_fixnum(betai), saveState, mk_double_float(testEnergy), energies._MonomerCorrections );
        }
      }
    }
    {
      // Check the histogram for flatness
      flatness_iter++;
      if (flatness_iter%flatness_steps==0) {
        size_t visits = 0, hmin = SIZE_MAX, hmax = 0;
        size_t flat_count = 0;
        size_t bi = 0;
        for (auto const& histogram : histograms ) {
          if ( flat_enough( debug,  histogram, bi, flatness ) ) flat_count++;
          bi++;
        }
        if (debug.notnilp()) core::clasp_write_string(fmt::format(" - - - - Flat_count = {}  histograms.size() = {}\n", flat_count, histograms.size()),debug);
        if (flat_count == histograms.size()) {
          flatness_reached++;
          if (debug.notnilp()) core::clasp_write_string(fmt::format("All histograms flat enough wl_increment updated to {:.4}\n", wl_increment ),debug);
          for ( size_t hi=0; hi<histograms.size(); hi++ ) {
            for ( size_t ii=0; ii<histograms[hi].size(); ii++ ) {
              histograms[hi][ii] = 0;
            }
          }
          wl_increment *= wl_scaling;
          if (flatnessCallback.notnilp()) {
            core::eval::funcall( flatnessCallback, core::make_fixnum(iter), mk_double_float(wl_increment), energies._MonomerCorrections );
          }
        }
      }
    }
    if (wl_increment<wl_increment_stop) {
      if (debug.notnilp()) core::clasp_write_string(fmt::format("wl_increment {:.6} is < wl_increment_stop {:.6}\n",
                                                                wl_increment, wl_increment_stop ),debug);
      goto DONE;
    }
    if (temperature_swap_steps && (currentStates.size()>=2) && (iter % temperature_swap_steps)==0) {
      // Try the temperature swap
      size_t i = tswapPickSlot(rng), j = i+1;

      double E_i = energies.physicalEnergy(currentStates[i]);
      double E_j = energies.physicalEnergy(currentStates[j]);
      double betai = (*betas)[i];
      double betaj = (*betas)[j];
      double delta = (betai-betaj)*(E_j - E_i);
      double x = - delta;
      double p = (x >= 0.0) ? 1.0 : std::exp(x);   // never exp(large +)
      if (debug.notnilp()) core::clasp_write_string(fmt::format("acc_beta_prob test: delta={:.6} p={:.3}\n", delta, p ),debug);
      double rnd = U01(rng);
      if (rnd < p) {
        if (debug.notnilp()) core::clasp_write_string(fmt::format(" <><><> Swapping betai states {} and {}\n", i, j ),debug);
        std::swap(currentStates[i], currentStates[j]);
      }
    }
    iterCount++;
  }
 DONE:
  return Values(core::make_fixnum(iterCount),
                core::make_fixnum(accepts),
                core::make_fixnum(flatness_reached)
                );
}
#endif

#if 1
CL_DOCSTRING(R"doc(Perform a constant temperature Hamiltonian replica exchange monte carlo on the
ENERGIES object using the LAMBDA-WINDOWS at TEMPERATURE taking LAMBDA-STEPS with each set of lambdas
before attempting to swap windows.
Return (values lowest-energy-state accepts swap-accepts-for-each-lambda-window swap-attempts-for-each-lambda-window max-iterations lowest-energy initial-state initial-energy))doc");
CL_LAMBDA(energies lambdaWindows &key (temperature 300.0) (lambda-steps 10) (max-iterations 1000) (warm-up-iterations 100) step-callback exchange-callback debug);
CL_DEFUN core::T_mv chem__constantTemperatureHamiltonianReplicaExchangeMonteCarlo(core::T_sp tenergies,
                                                                                  core::T_sp tlambdaWindows,
                                                                                  double temperature,
                                                                                  size_t lambdaSteps,
                                                                                  size_t maxIterations,
                                                                                  size_t warmUpIterations,
                                                                                  core::T_sp stepCallback,
                                                                                  core::T_sp exchangeCallback,
                                                                                  core::T_sp debug ) {

  size_t accepts = 0;

  double beta = 1.0/(0.0019872*temperature);

  core::SimpleVector_double_sp lambdaWindows = gc::As<core::SimpleVector_double_sp>(tlambdaWindows);
  size_t numberOfLambdaWindows = core::cl__length(lambdaWindows);

  Energies energies(tenergies);
  core::SimpleVector_byte32_t_sp saveState;
  bool useStepCallback = stepCallback.notnilp();
  if (useStepCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  bool useExchangeCallback = exchangeCallback.notnilp();

  std::vector<State> currentStates(numberOfLambdaWindows);
  for ( size_t ii = 0; ii<currentStates.size(); ii++ ) currentStates[ii].randomState(energies);
  State initialState = currentStates[numberOfLambdaWindows-1];

  std::vector<byte64_t> swapAccepts(numberOfLambdaWindows,0);
  std::vector<byte64_t> swapAttempts(numberOfLambdaWindows,0);

  size_t numSlotsInState = core::cl__length(energies._MonomerLocusMaxMrkindex );

  State testState;
  State lowestState;
  double lowestEnergy = std::numeric_limits<double>::max();
  for (int ii = 0; ii < maxIterations; ++ii) {
    for (size_t lambdaWindowIdx=0; lambdaWindowIdx<numberOfLambdaWindows; lambdaWindowIdx++ ) {
      double lambda = (*lambdaWindows)[lambdaWindowIdx];
      for ( size_t lambdaStepi=0; lambdaStepi<lambdaSteps; lambdaStepi++ ) {
        double curEnergy = energies.reducedEnergy(currentStates[lambdaWindowIdx],lambda);
        size_t testSlotIndex;
        testState = energies.randomStep( testSlotIndex, currentStates[lambdaWindowIdx] );
        double testEnergy = energies.reducedEnergy(testState,lambda);
        double delta = beta * (testEnergy - curEnergy);
        double x = - delta;
        double p = (x >= 0.0) ? 1.0 : std::exp(x);   // never exp(large +)
        if (debug.notnilp()) core::clasp_write_string(fmt::format("acc_prob test: testE={:.6} prevE={:.6} beta={:.3} delta={:.6} p={:.3}\n", testEnergy, curEnergy, beta, delta, p ),debug);
        double rnd = energies.U01();
        if ( rnd <= p) {
          accepts++;
          currentStates[lambdaWindowIdx] = testState;
          curEnergy = testEnergy;
          if (debug.notnilp()) core::clasp_write_string(fmt::format(" +++ accepted step\n"),debug);
        }
        if (useStepCallback) {
          memcpy(&(*saveState)[0],&currentStates[lambdaWindowIdx]._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( stepCallback, saveState, core::make_fixnum(testSlotIndex), mk_double_float(curEnergy) );
        }
      }
      // In lambdaWindowIdx loop
    }
    // Here I would swap windows
    if (debug.notnilp()) core::clasp_write_string(fmt::format("Try to swap lambda windows now\n"),debug);

    for ( size_t swapi = ii%2; swapi<numberOfLambdaWindows; swapi += 2 ) {
      if (swapi+1<numberOfLambdaWindows) {
        State& k0 = currentStates[swapi];
        State& k1 = currentStates[swapi+1];
        double l0 = (*lambdaWindows)[swapi];
        double l1 = (*lambdaWindows)[swapi+1];
        double ul1k0 = energies.reducedEnergy(k0,l1);
        double ul0k0 = energies.reducedEnergy(k0,l0);
        double ul1k1 = energies.reducedEnergy(k1,l1);
        double ul0k1 = energies.reducedEnergy(k1,l0);
        double dul1k0_ul0k0 = (ul1k0-ul0k0);
        double dul0k1_ul1k1 = (ul0k1-ul1k1);
        double deltau = dul1k0_ul0k0 +dul0k1_ul1k1;
        double x = beta*(- deltau);
        double p = (x >= 0.0) ? 1.0 : std::exp(x);
        double rnd = energies.U01();
        if (rnd <= p) {
          swapAccepts[swapi]++;
          State tempState = currentStates[swapi];
          currentStates[swapi] = currentStates[swapi+1];
          currentStates[swapi+1] = tempState;
        }
        swapAttempts[swapi]++;
        if (useExchangeCallback && ii>warmUpIterations) {
          core::eval::funcall(exchangeCallback,
                              core::make_fixnum(ii),
                              core::make_fixnum(swapi),
                              mk_double_float(l0),
                              mk_double_float(l1),
                              mk_double_float(dul1k0_ul0k0),
                              mk_double_float(dul0k1_ul1k1));
        }
      }
    }
    double energy = energies.reducedEnergy(currentStates[numberOfLambdaWindows-1],1.0);
    if (energy < lowestEnergy) {
      lowestEnergy = energy;
      lowestState = currentStates[numberOfLambdaWindows-1];
    }
    // In maxIterations loop
  }
  double initialEnergy = energies.reducedEnergy(initialState,1.0);
  return Values(core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots,
                                                    0, false,
                                                    energies._NumberOfSlots,   // initialContentsSize
                                                    &lowestState._State[0] ),  // initialContents
                core::make_fixnum(accepts),
                core::SimpleVector_byte64_t_O::make(numberOfLambdaWindows-1,
                                                  0, false,
                                                  numberOfLambdaWindows-1,
                                                  &swapAccepts[0]),
                core::SimpleVector_byte64_t_O::make(numberOfLambdaWindows-1,
                                                  0, false,
                                                  numberOfLambdaWindows-1,
                                                  &swapAttempts[0]),
                core::make_fixnum(maxIterations),
                mk_double_float(lowestEnergy),
                core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots,
                                                    0, false,
                                                    energies._NumberOfSlots,   // initialContentsSize
                                                    &initialState._State[0] ),  // initialContents
                mk_double_float(initialEnergy)

                );
}
#endif

}; // namespace chem


