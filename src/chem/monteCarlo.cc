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


SYMBOL_EXPORT_SC_( ChemPkg, backbone_energy );
SYMBOL_EXPORT_SC_( ChemPkg, monomer_vector );
SYMBOL_EXPORT_SC_( ChemPkg, monomer_corrections );
SYMBOL_EXPORT_SC_( ChemPkg, mrkey_index_to_lmkey_index );
SYMBOL_EXPORT_SC_( ChemPkg, monomer_locus_max_mrkindex );
SYMBOL_EXPORT_SC_( ChemPkg, single_scan_energy );
SYMBOL_EXPORT_SC_( ChemPkg, pair_scan_energy_lower_triangular_matrix );

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
  size_t                         _NumberOfSlots;
  double                         _BackboneEnergy;
  core::SimpleVector_sp          _MonomerVector;
  core::SimpleVector_byte32_t_sp _MonomerLocusMaxMrkindex;
  core::SimpleVector_byte32_t_sp _MrkeyIndexToLmkeyIndex;
  core::SimpleVector_double_sp   _MonomerCorrections;
  size_t                         _MonomerCorrectionsLength;
  core::SimpleVector_double_sp   _SingleTerms;
  core::SimpleVector_double_sp   _PairTerms;
  Energies(core::T_sp energies) {
    core::T_sp val;
    val = core::eval::funcall( _sym_backbone_energy, energies );
    this->_BackboneEnergy = gc::As<core::DoubleFloat_sp>(val)->get();
    val = core::eval::funcall( _sym_monomer_vector, energies );
    this->_MonomerVector = gc::As<core::SimpleVector_sp>(val);
    val = core::eval::funcall( _sym_monomer_locus_max_mrkindex, energies );
    this->_MonomerLocusMaxMrkindex = gc::As<core::SimpleVector_byte32_t_sp>(val);
    val = core::eval::funcall( _sym_mrkey_index_to_lmkey_index, energies );
    this->_MrkeyIndexToLmkeyIndex = gc::As<core::SimpleVector_byte32_t_sp>(val);
    this->_NumberOfSlots = core::cl__length(this->_MonomerLocusMaxMrkindex);
    val = core::eval::funcall( _sym_monomer_corrections, energies );
    this->_MonomerCorrections = gc::As<core::SimpleVector_double_sp>(val);
    this->_MonomerCorrectionsLength = core::cl__length(this->_MonomerCorrections);
    val = core::eval::funcall( _sym_single_scan_energy, energies );
    this->_SingleTerms = gc::As<core::SimpleVector_double_sp>(val);
    val = core::eval::funcall( _sym_pair_scan_energy_lower_triangular_matrix, energies );
    this->_PairTerms = gc::As<core::SimpleVector_double_sp>(val);
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
    total.Add(this->_BackboneEnergy);
    total.Add(singleSum.sum);
    total.Add(pairSum.sum);
    return total.sum;
  }

  double energyFunction(State& state, double beta) {
    double* singleTerms = &(*this->_SingleTerms)[0];
    double* pairTerms = &(*this->_PairTerms)[0];
    KahanAccumulator singleSum;
    for ( int ii=0; ii<state._State.size(); ii++ ) {
      int index = state._State[ii];
      double energyTerm = singleTerms[index];
      singleSum.Add(energyTerm);
      uint32_t lmkey_index = (*this->_MrkeyIndexToLmkeyIndex)[index];
      double correction = beta*this->_MonomerCorrections[lmkey_index];
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
    total.Add(this->_BackboneEnergy);
    total.Add(singleSum.sum);
    total.Add(pairSum.sum);
    return total.sum;
  }

  State randomStep(size_t slotIndex, const State& state) {
    State tempState(state);
    int max = (*this->_MonomerLocusMaxMrkindex)[slotIndex];
    int prevMax = (slotIndex==0) ? -1 : (*this->_MonomerLocusMaxMrkindex)[slotIndex-1];
    int range = max - prevMax;
    tempState._State[slotIndex] = (rand() % range) + prevMax + 1;
    return tempState;
  }

// Temperature decay function
  double tempDecay(const double& temperature) {
      // Replace with your temperature decay logic
    return temperature * 0.99;
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
CL_DEFUN core::DoubleFloat_sp chem__mcstate_energy(core::T_sp tmcstate, double beta, core::T_sp tenergies )
{
  Energies energies(tenergies);
  State state(gc::As<core::SimpleVector_byte32_t_sp>(tmcstate));
  double testEnergy = energies.energyFunction(state,beta);
  return core::DoubleFloat_O::create(testEnergy);
}

#if 0
CL_LAMBDA(energies &key (start-temperature 100.0) (max-iterations 1000) initial-state accept-callback);
CL_DEFUN core::T_mv chem__simulatedAnnealing(core::T_sp tenergies, double startTemperature, size_t max_iterations, core::T_sp tinitial_state, core::T_sp acceptCallback ) {
  srand(time(NULL)); // Initialize random seed
  Energies energies(tenergies);
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  if (useAcceptCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
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
  double currentEnergy = energies.energyFunction(currentState);
  double temperature = startTemperature;
  State testState;
  for (int i = 0; i < max_iterations; ++i) {
    for ( size_t slotIndex = 0; slotIndex< energies._NumberOfSlots; slotIndex++ ) {
      testState = energies.randomStep( slotIndex, currentState );
      double testEnergy = energies.energyFunction(testState);
      if (testEnergy < currentEnergy) {
        currentState = testState;
        currentEnergy = testEnergy;
        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentState._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, saveState, mk_double_float(currentEnergy) );
        }
      } else {
        double acceptance_probability = exp(-(testEnergy - currentEnergy) / temperature);
        double rand_val = static_cast<double>(rand()) / RAND_MAX;
        if (rand_val < acceptance_probability) {
          currentState = testState;
          currentEnergy = testEnergy;
          if (useAcceptCallback) {
            memcpy(&(*saveState)[0],&currentState._State[0],sizeof(int32_t)*energies._NumberOfSlots);
            core::eval::funcall( acceptCallback, saveState, mk_double_float(currentEnergy) );
          }
        }
      }
    }
    temperature = energies.tempDecay(temperature);
  }
  return Values(core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots,
                                                    0, false,
                                                    energies._NumberOfSlots, // initialContentsSize
                                                    &currentState._State[0] ), // initialContents
                mk_double_float(currentEnergy),
                mk_double_float(temperature));
}
#endif

#if 0
CL_LAMBDA(energies &key (temperature 100.0) (max-iterations 1000) initial-state accept-callback);
CL_DEFUN core::T_mv chem__constantTemperatureMonteCarlo(core::T_sp tenergies, double temperature, size_t max_iterations, core::T_sp tinitial_state, core::T_sp acceptCallback ) {
  srand(time(NULL)); // Initialize random seed
  Energies energies(tenergies);
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  if (useAcceptCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
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
  double currentEnergy = energies.energyFunction(currentState);
  State testState;
  for (int i = 0; i < max_iterations; ++i) {
    for ( size_t slotIndex = 0; slotIndex< energies._NumberOfSlots; slotIndex++ ) {
      testState = energies.randomStep( slotIndex, currentState );
      double testEnergy = energies.energyFunction(testState);
      if (testEnergy < currentEnergy) {
        currentState = testState;
        currentEnergy = testEnergy;
        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentState._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, saveState, mk_double_float(currentEnergy) );
        }
      } else {
        double acceptance_probability = exp(-(testEnergy - currentEnergy) / temperature);
        double rand_val = static_cast<double>(rand()) / RAND_MAX;
        if (rand_val < acceptance_probability) {
          currentState = testState;
          currentEnergy = testEnergy;
          if (useAcceptCallback) {
            memcpy(&(*saveState)[0],&currentState._State[0],sizeof(int32_t)*energies._NumberOfSlots);
            core::eval::funcall( acceptCallback, saveState, mk_double_float(currentEnergy) );
          }
        }
      }
    }
  }
  return Values(core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots,
                                                    0, false,
                                                    energies._NumberOfSlots, // initialContentsSize
                                                    &currentState._State[0] ), // initialContents
                mk_double_float(currentEnergy) );
}
#endif


inline double acc_prob_from_log(double x) {
  // x = log(prob) = -βΔE  (Metropolis)  or  (βi-βj)(Ej-Ei) (REX)
  return (x >= 0.0) ? 1.0 : std::exp(x);   // never exp(large +)
}

CL_LAMBDA(energies betas &key (temperature-swap-steps 100) (max-iterations 1000000) accept-callback (flatness-steps 100) (wl-scaling 0.8) (wl-increment-stop 0.1) debug);
CL_DEFUN core::T_mv chem__old_optimize_monomer_corrections(core::T_sp tenergies, core::SimpleVector_double_sp betas, size_t temperature_swap_steps, size_t max_iterations, core::T_sp acceptCallback, size_t flatness_steps, double wl_scaling, double wl_increment_stop, bool debug ) {
  size_t beta_size = core::cl__length(betas);
  double wl_increment = 20.0;
  double flatness = 0.7;
  size_t flat_count = 0;
  size_t accepts = 0;
  Energies energies(tenergies);
  // Number of slots in state
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  std::mt19937 rng(123456u);
  std::uniform_real_distribution<double> U01(0.0,1.0);
  std::uniform_int_distribution<size_t>  pickSlot(0, num_slots_in_state-1 );
  std::uniform_int_distribution<size_t>  tswapPickSlot(0, beta_size-2 );
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  if (useAcceptCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  std::vector<State> currentStates(beta_size);
  for ( size_t ii = 0; ii<currentStates.size(); ii++ ) {
    currentStates[ii].randomState(energies);
  }
  State testState;
  std::vector<size_t> histogram(energies._MonomerCorrectionsLength,0);
  size_t flatness_iter = 0;
  size_t iterCount = 0;
  for (size_t iter = 0; iter < max_iterations; ++iter) {
    if (debug) core::clasp_write_string(fmt::format("Top of loop at iter {}\n", iter));
    for ( size_t betai = 0; betai<beta_size; betai++ ) {
      double beta = (*betas)[betai];
      size_t slotIndex = pickSlot(rng);
      double currentEnergy = energies.oldEnergyFunction(currentStates[betai]);
      testState = energies.randomStep( slotIndex, currentStates[betai] );
      double testEnergy = energies.oldEnergyFunction(testState);
      double x = -(*betas)[betai] * (testEnergy - currentEnergy);
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if ( rnd <= p) {
        accepts++;
        uint32_t curLmkeyIndex = (*energies._MrkeyIndexToLmkeyIndex)[testState._State[slotIndex]];
        if (debug) core::clasp_write_string(fmt::format("Accepted mc step iter {} betai {}  curLmkeyIndex {}\n", iter, betai, curLmkeyIndex ));
        {
          // Increment the histogram and the MonomerCorrection for the current monomer
          histogram[curLmkeyIndex]++;
          if (debug) {
            core::clasp_write_string(fmt::format("          Histogram: "));
            for ( size_t ii=0; ii<histogram.size(); ii++ ) {
              if (ii==curLmkeyIndex) core::clasp_write_string("*");
              core::clasp_write_string(fmt::format("{} ", histogram[ii]));
            }
            core::clasp_write_string(fmt::format("\n"));
            core::clasp_write_string(fmt::format("wl_increment {}\n", wl_increment ));
            core::clasp_write_string(fmt::format("Monomer corrections: "));
            for ( size_t ii=0; ii<energies._MonomerCorrectionsLength; ii++ ) {
              if (ii==curLmkeyIndex) core::clasp_write_string("*");
              core::clasp_write_string(fmt::format("{:.4f} ", (*energies._MonomerCorrections)[ii]));
            }
            core::clasp_write_string(fmt::format("\n"));
          }
          (*energies._MonomerCorrections)[curLmkeyIndex] += wl_increment;
          // Reset the bias
          double bias = std::numeric_limits<double>::max();
          for ( size_t ii = 0; ii<energies._MonomerCorrectionsLength; ii++ ) {
            bias = std::min((*energies._MonomerCorrections)[ii],bias);
          }
          for ( size_t ii = 0; ii<energies._MonomerCorrectionsLength; ii++ ) {
            (*energies._MonomerCorrections)[ii] -= bias;
          }
          if (debug) core::clasp_write_string(fmt::format("Updated MonomerCorrections at iter {} betai {} curLmkeyIndex {} with {}\n", iter, betai, curLmkeyIndex, wl_increment ));
        }
        currentStates[betai] = testState;
        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentStates[betai]._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, core::make_fixnum(iter), core::make_fixnum(betai), saveState, mk_double_float(testEnergy), energies._MonomerCorrections );
        }
        // Check the histogram for flatness
        flatness_iter++;
        if (flatness_iter%flatness_steps==0) {
          size_t visits = 0, hmin = SIZE_MAX;
          for (auto h : histogram) { visits += h; hmin = std::min(hmin,h); }
          if (visits >= 50* histogram.size()) {
            double avg = double(visits) / histogram.size();
            if (debug) core::clasp_write_string(fmt::format("Enough visits hmin = {}  avg = {} hmin/avg={}  flatness = {}\n", hmin, avg, hmin/avg, flatness ));
            if (double(hmin)/avg > flatness) {
              std::fill( histogram.begin(), histogram.end(), 0 );
              wl_increment *= wl_scaling;
              flat_count++;
              if (debug) {
                core::clasp_write_string(fmt::format("wl_increment updated to {}\n", wl_increment ));
              }
            }
          }
        }
        if (wl_increment<wl_increment_stop) {
          if (debug) core::clasp_write_string(fmt::format("wl_increment {} is < wl_increment_stop {}\n", wl_increment, wl_increment_stop ));
          goto DONE;
        }
      }
    }
    if (temperature_swap_steps && (currentStates.size()>=2) && (iter % temperature_swap_steps)==0) {
      // Try the temperature swap
      size_t i = tswapPickSlot(rng), j = i+1;
      double currentEnergy_i = energies.oldEnergyFunction(currentStates[i]);
      double currentEnergy_j = energies.oldEnergyFunction(currentStates[j]);
      double x = ( (*betas)[i] - (*betas)[j] ) * ( currentEnergy_j - currentEnergy_i );
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if (rnd < p) {
        if (debug) core::clasp_write_string(fmt::format("Swapping states {} and {}\n", i, j ));
        std::swap(currentStates[i], currentStates[j]);
      }
    }
    iterCount++;
  }
 DONE:
  return Values(
      core::make_fixnum(iterCount),
      core::make_fixnum(accepts),
      core::make_fixnum(flat_count));
}




CL_LAMBDA(energies betas &key (temperature-swap-steps 100) (max-iterations 1000000) accept-callback (flatness-steps 100) (wl-scaling 0.8) (wl-increment-stop 0.1) debug);
CL_DEFUN core::T_mv chem__original_optimize_monomer_corrections(core::T_sp tenergies, core::SimpleVector_double_sp betas, size_t temperature_swap_steps, size_t max_iterations, core::T_sp acceptCallback, size_t flatness_steps, double wl_scaling, double wl_increment_stop, bool debug ) {
  size_t beta_size = core::cl__length(betas);
  double wl_increment = 20.0;
  double flatness = 0.7;
  Energies energies(tenergies);
  // Number of slots in state
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  std::mt19937 rng(123456u);
  std::uniform_real_distribution<double> U01(0.0,1.0);
  std::uniform_int_distribution<size_t>  pickSlot(0, num_slots_in_state-1 );
  std::uniform_int_distribution<size_t>  tswapPickSlot(0, beta_size-2 );
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  if (useAcceptCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  std::vector<State> currentStates(beta_size);
  for ( size_t ii = 0; ii<currentStates.size(); ii++ ) {
    currentStates[ii].randomState(energies);
  }
  State testState;
  std::vector<size_t> histogram(energies._MonomerCorrectionsLength,0);
  size_t flatness_iter = 0;
  size_t accepts = 0;
  size_t flatness_counts = 0;
  size_t iterCount = 0;
  for (size_t iter = 0; iter < max_iterations; ++iter) {
    if (debug) core::clasp_write_string(fmt::format("Top of loop at iter {}\n", iter));
    for ( size_t betai = 0; betai<beta_size; betai++ ) {
      double beta = (*betas)[betai];
      size_t slotIndex = pickSlot(rng);
      double currentEnergy = energies.energyFunction(currentStates[betai],1.0);
      testState = energies.randomStep( slotIndex, currentStates[betai] );
      double testEnergy = energies.energyFunction(testState,1.0);
      double x = -(*betas)[betai] * (testEnergy - currentEnergy);
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if ( rnd <= p) {
        accepts++;
        uint32_t curLmkeyIndex = (*energies._MrkeyIndexToLmkeyIndex)[testState._State[slotIndex]];
        if (debug) core::clasp_write_string(fmt::format("Accepted mc step iter {} betai {}  curLmkeyIndex {}\n", iter, betai, curLmkeyIndex ));
        {
          // Increment the histogram and the MonomerCorrection for the current monomer
          histogram[curLmkeyIndex]++;
          if (debug) {
            core::clasp_write_string(fmt::format("          Histogram: "));
            for ( size_t ii=0; ii<histogram.size(); ii++ ) {
              if (ii==curLmkeyIndex) core::clasp_write_string("*");
              core::clasp_write_string(fmt::format("{} ", histogram[ii]));
            }
            core::clasp_write_string(fmt::format("\n"));
            core::clasp_write_string(fmt::format("wl_increment {}\n", wl_increment ));
            core::clasp_write_string(fmt::format("Monomer corrections: "));
            for ( size_t ii=0; ii<energies._MonomerCorrectionsLength; ii++ ) {
              if (ii==curLmkeyIndex) core::clasp_write_string("*");
              core::clasp_write_string(fmt::format("{:.4f} ", (*energies._MonomerCorrections)[ii]));
            }
            core::clasp_write_string(fmt::format("\n"));
          }
          (*energies._MonomerCorrections)[curLmkeyIndex] += wl_increment;
          // Reset the bias
          double bias = std::numeric_limits<double>::max();
          for ( size_t ii = 0; ii<energies._MonomerCorrectionsLength; ii++ ) {
            bias = std::min((*energies._MonomerCorrections)[ii],bias);
          }
          for ( size_t ii = 0; ii<energies._MonomerCorrectionsLength; ii++ ) {
            (*energies._MonomerCorrections)[ii] -= bias;
          }
          if (debug) core::clasp_write_string(fmt::format("Updated MonomerCorrections at iter {} betai {} curLmkeyIndex {} with {}\n", iter, betai, curLmkeyIndex, wl_increment ));
        }
        currentStates[betai] = testState;
        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentStates[betai]._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, core::make_fixnum(iter), core::make_fixnum(betai), saveState, mk_double_float(testEnergy), energies._MonomerCorrections );
        }
        // Check the histogram for flatness
        flatness_iter++;
        if (flatness_iter%flatness_steps==0) {
          size_t visits = 0, hmin = SIZE_MAX;
          for (auto h : histogram) { visits += h; hmin = std::min(hmin,h); }
          if (visits >= 50* histogram.size()) {
            double avg = double(visits) / histogram.size();
            if (debug) core::clasp_write_string(fmt::format("Enough visits hmin = {}  avg = {} hmin/avg={}  flatness = {}\n", hmin, avg, hmin/avg, flatness ));
            if (double(hmin)/avg > flatness) {
              std::fill( histogram.begin(), histogram.end(), 0 );
              wl_increment *= wl_scaling;
              flatness_counts++;
              if (debug) {
                core::clasp_write_string(fmt::format("wl_increment updated to {}\n", wl_increment ));
              }
            }
          }
        }
        if (wl_increment<wl_increment_stop) {
          if (debug) core::clasp_write_string(fmt::format("wl_increment {} is < wl_increment_stop {}\n", wl_increment, wl_increment_stop ));
          goto DONE;
        }
      }
    }
    if (temperature_swap_steps && (currentStates.size()>=2) && (iter % temperature_swap_steps)==0) {
      // Try the temperature swap
      size_t i = tswapPickSlot(rng), j = i+1;
      double currentEnergy_i = energies.energyFunction(currentStates[i],1.0);
      double currentEnergy_j = energies.energyFunction(currentStates[j],1.0);
      double x = ( (*betas)[i] - (*betas)[j] ) * ( currentEnergy_j - currentEnergy_i );
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if (rnd < p) {
        if (debug) core::clasp_write_string(fmt::format("Swapping states {} and {}\n", i, j ));
        std::swap(currentStates[i], currentStates[j]);
      }
    }
    iterCount++;
  }
 DONE:
  return Values(core::make_fixnum(iterCount),
                core::make_fixnum(accepts),
                core::make_fixnum(flatness_counts));
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

CL_LAMBDA(energies betas &key (temperature-swap-steps 100) (max-iterations 1000000) accept-callback (flatness-steps 100) (wl-scaling 0.8) (wl-increment-stop 0.1) debug);
CL_DEFUN core::T_mv chem__updated_optimize_monomer_corrections(core::T_sp tenergies, core::SimpleVector_double_sp betas, size_t temperature_swap_steps, size_t max_iterations, core::T_sp acceptCallback, size_t flatness_steps, double wl_scaling, double wl_increment_stop, bool debug ) {
  size_t beta_size = core::cl__length(betas);
  double wl_increment = 20.0;
  double flatness = 0.7;
  Energies energies(tenergies);
  // Number of slots in state
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  std::mt19937 rng(123456u);
  std::uniform_real_distribution<double> U01(0.0,1.0);
  std::uniform_int_distribution<size_t>  pickSlot(0, num_slots_in_state-1 );
  std::uniform_int_distribution<size_t>  tswapPickSlot(0, beta_size-2 );
  core::SimpleVector_byte32_t_sp saveState;
  bool useAcceptCallback = acceptCallback.notnilp();
  if (useAcceptCallback) saveState = core::SimpleVector_byte32_t_O::make(energies._NumberOfSlots);
  std::vector<State> currentStates(beta_size);
  for ( size_t ii = 0; ii<currentStates.size(); ii++ ) {
    currentStates[ii].randomState(energies);
  }
  State testState;
  std::vector<size_t> histogram(energies._MonomerCorrectionsLength,0);
  size_t flatness_iter = 0;
  size_t iterCount = 0;
  size_t accepts = 0;
  size_t flatness_reached = 0;
  for (size_t iter = 0; iter < max_iterations; ++iter) {
    if (debug) core::clasp_write_string(fmt::format("Top of loop at iter {}\n", iter));
    for ( size_t betai = 0; betai<beta_size; betai++ ) {
      double beta = (*betas)[betai];
      size_t slotIndex = pickSlot(rng);
      double currentEnergy = energies.energyFunction(currentStates[betai],1.0);
      testState = energies.randomStep( slotIndex, currentStates[betai] );
      double testEnergy = energies.energyFunction(testState,1.0);
      double x = -(*betas)[betai] * (testEnergy - currentEnergy);
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if ( rnd <= p) {
        accepts++;
        uint32_t curLmkeyIndex = (*energies._MrkeyIndexToLmkeyIndex)[testState._State[slotIndex]];
        if (debug) {
          core::clasp_write_string(fmt::format("Accepted (rnd={} <= p={}) mc step iter {} betai {}  curLmkeyIndex {}\n", rnd, p, iter, betai, curLmkeyIndex ));
        }
        // Increment the histogram and the MonomerCorrection for the current monomer
        histogram[curLmkeyIndex]++;
        if (debug) {
          core::clasp_write_string(fmt::format("          Histogram: "));
          for ( size_t ii=0; ii<histogram.size(); ii++ ) {
            if (ii==curLmkeyIndex) core::clasp_write_string("*");
            core::clasp_write_string(fmt::format("{} ", histogram[ii]));
          }
          core::clasp_write_string(fmt::format("\n"));
          core::clasp_write_string(fmt::format("Monomer corrections: "));
          for ( size_t ii=0; ii<energies._MonomerCorrectionsLength; ii++ ) {
            if (ii==curLmkeyIndex) core::clasp_write_string("*");
            core::clasp_write_string(fmt::format("{:12.2f} ", (*energies._MonomerCorrections)[ii]));
          }
          core::clasp_write_string(fmt::format("\n"));
        }
        {
          (*energies._MonomerCorrections)[curLmkeyIndex] += wl_increment;
          if (debug)
            core::clasp_write_string(fmt::format("Updating _MonomerCorrections by Updated MonomerCorrections at iter {} betai {} curLmkeyIndex {} beta={} wl_increment {} \n", iter, betai, curLmkeyIndex, beta, wl_increment ));
        }
        // Reset the bias
        adjustMonomerCorrectionsBias(energies);
        if (debug) core::clasp_write_string(fmt::format("Updated MonomerCorrections at iter {} betai {} curLmkeyIndex {} with {}\n", iter, betai, curLmkeyIndex, wl_increment ));
        currentStates[betai] = testState;
        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentStates[betai]._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, core::make_fixnum(iter), core::make_fixnum(betai), saveState, mk_double_float(testEnergy), energies._MonomerCorrections );
        }
      }
    }
    // Check the histogram for flatness
    flatness_iter++;
    if (flatness_iter%flatness_steps==0) {
      size_t visits = 0, hmin = SIZE_MAX;
      for (auto h : histogram) { visits += h; hmin = std::min(hmin,h); }
      if (visits >= 50* histogram.size()) {
        double avg = double(visits) / histogram.size();
        if (debug) core::clasp_write_string(fmt::format("Enough visits hmin = {}  avg = {} hmin/avg={}  flatness = {}\n", hmin, avg, hmin/avg, flatness ));
        if (double(hmin)/avg > flatness) {
          std::fill( histogram.begin(), histogram.end(), 0 );
          wl_increment *= wl_scaling;
          if (debug) {
            core::clasp_write_string(fmt::format("wl_increment updated to {}\n", wl_increment ));
          }
        }
      }
    }
    if (wl_increment<wl_increment_stop) {
      if (debug) core::clasp_write_string(fmt::format("wl_increment {} is < wl_increment_stop {}\n", wl_increment, wl_increment_stop ));
      goto DONE;
    }
    if (temperature_swap_steps && (currentStates.size()>=2) && (iter % temperature_swap_steps)==0) {
      // Try the temperature swap
      size_t i = tswapPickSlot(rng), j = i+1;
      double currentEnergy_i = energies.energyFunction(currentStates[i],1.0);
      double currentEnergy_j = energies.energyFunction(currentStates[j],1.0);
      double x = ( (*betas)[i] - (*betas)[j] ) * ( currentEnergy_j - currentEnergy_i );
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if (rnd < p) {
        if (debug) core::clasp_write_string(fmt::format("Swapping states {} and {}\n", i, j ));
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






#if 0
CL_LAMBDA(energies betas &key (temperature-swap-steps 100) (max-iterations 1000000) accept-callback (flatness-steps 100) (wl-scaling 0.8) (wl-increment-stop 0.1) debug);
CL_DEFUN core::T_mv chem__voelz_optimize_monomer_corrections(core::T_sp tenergies, core::SimpleVector_double_sp betas, size_t temperature_swap_steps, size_t max_iterations, core::T_sp acceptCallback, size_t flatness_steps, double wl_scaling, double wl_increment_stop, bool debug ) {
  size_t beta_size = core::cl__length(betas);
  double wl_increment = 20.0;
  double flatness = 0.7;
  Energies energies(tenergies);
  // Number of slots in state
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  std::mt19937 rng(123456u);
  std::uniform_real_distribution<double> U01(0.0,1.0);
  std::uniform_int_distribution<size_t>  pickSlot(0, num_slots_in_state-1 );
  std::uniform_int_distribution<size_t>  tswapPickSlot(0, beta_size-2 );
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
    if (debug) core::clasp_write_string(fmt::format("Top of loop at iter {}\n", iter));
    for ( size_t betai = 0; betai<beta_size; betai++ ) {
      double beta = (*betas)[betai];
      size_t slotIndex = pickSlot(rng);
      double currentEnergy = energies.energyFunction(currentStates[betai],beta);
      testState = energies.randomStep( slotIndex, currentStates[betai] );
      double testEnergy = energies.energyFunction(testState,beta);
      double x = -(*betas)[betai] * (testEnergy - currentEnergy);
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if ( rnd <= p) {
        accepts++;
        uint32_t curLmkeyIndex = (*energies._MrkeyIndexToLmkeyIndex)[testState._State[slotIndex]];
        if (debug) {
          core::clasp_write_string(fmt::format("Accepted (rnd={} <= p={}) mc step iter {} betai {}  curLmkeyIndex {}\n", rnd, p, iter, betai, curLmkeyIndex ));
        }
        {
          // Increment the histogram and the MonomerCorrection for the current monomer
          histograms[betai][curLmkeyIndex]++;
          if (debug) {
            core::clasp_write_string(fmt::format("Histogram(betai={:2}): ", betai ));
            for ( size_t ii=0; ii<histograms[betai].size(); ii++ ) {
              if (ii==curLmkeyIndex) core::clasp_write_string("*");
              core::clasp_write_string(fmt::format("{:12} ", histograms[betai][ii]));
            }
            core::clasp_write_string(fmt::format("\n"));
            core::clasp_write_string(fmt::format("Monomer corrections: "));
            for ( size_t ii=0; ii<energies._MonomerCorrectionsLength; ii++ ) {
              if (ii==curLmkeyIndex) core::clasp_write_string("*");
              core::clasp_write_string(fmt::format("{:12.2f} ", (*energies._MonomerCorrections)[ii]));
            }
            core::clasp_write_string(fmt::format("\n"));
          }
          {
            double beta_increment = wl_increment/(*betas)[betai];
            if (debug)
              core::clasp_write_string(fmt::format("Updating _MonomerCorrections by Updated MonomerCorrections at iter {} betai {} curLmkeyIndex {} beta={} wl_increment {} beta_increment {}\n", iter, betai, curLmkeyIndex, beta, wl_increment, beta_increment ));
            (*energies._MonomerCorrections)[curLmkeyIndex] += beta_increment;
          }
        }
        // Reset the bias
        double bias = std::numeric_limits<double>::max();
        for ( size_t ii = 0; ii<energies._MonomerCorrectionsLength; ii++ ) {
          bias = std::min((*energies._MonomerCorrections)[ii],bias);
        }
        for ( size_t ii = 0; ii<energies._MonomerCorrectionsLength; ii++ ) {
          (*energies._MonomerCorrections)[ii] -= bias;
        }
        if (debug) core::clasp_write_string(fmt::format("Updated MonomerCorrections at iter {} betai {} curLmkeyIndex {} with {}\n", iter, betai, curLmkeyIndex, wl_increment ));
      }
      currentStates[betai] = testState;
      if (useAcceptCallback) {
        memcpy(&(*saveState)[0],&currentStates[betai]._State[0],sizeof(int32_t)*energies._NumberOfSlots);
        core::eval::funcall( acceptCallback, core::make_fixnum(iter), core::make_fixnum(betai), saveState, mk_double_float(testEnergy), energies._MonomerCorrections );
      }
    }
    // Check the histogram for flatness
    flatness_iter++;
    if (flatness_iter%flatness_steps==0) {
      size_t visits = 0, hmin = SIZE_MAX, hmax = 0;
      size_t flat_count = 0;
      size_t bi = 0;
      for (auto histogram : histograms ) {
        for (auto h : histogram) { visits += h; hmin = std::min(hmin,h); hmax = std::max(hmax,h); }
        if (visits >= 50* histogram.size()) {
          double avg = double(visits) / histogram.size();
          if (debug) core::clasp_write_string(fmt::format("Enough visits for bi={} hmin = {}  avg = {} hmin/avg={}  flatness = {}  {}\n", bi, hmin, avg, hmin/avg, flatness, (double(hmin)/avg>flatness) ? "FLAT" : "not-flat" ));
          if (double(hmin)/avg > flatness) flat_count++;
        }
        bi++;
      }
      if (debug) core::clasp_write_string(fmt::format("Flat_count = {}  histograms.size() = {}\n", flat_count, histograms.size()));
      if (flat_count == histograms.size()) {
        flatness_reached++;
        if (debug) core::clasp_write_string(fmt::format("All histograms flat enoughwl_increment updated to {}\n", wl_increment ));
        for ( auto histogram : histograms ) {
          std::fill( histogram.begin(), histogram.end(), 0 );
        }
        wl_increment *= wl_scaling;
      }
    }
    if (wl_increment<wl_increment_stop) {
      if (debug) core::clasp_write_string(fmt::format("wl_increment {} is < wl_increment_stop {}\n", wl_increment, wl_increment_stop ));
      goto DONE;
    }
    if (temperature_swap_steps && (currentStates.size()>=2) && (iter % temperature_swap_steps)==0) {
      // Try the temperature swap
      size_t i = tswapPickSlot(rng), j = i+1;
      double currentEnergy_i = energies.energyFunction(currentStates[i],(*betas)[i]);
      double currentEnergy_j = energies.energyFunction(currentStates[j],(*betas)[j]);
      double x = ( (*betas)[i] - (*betas)[j] ) * ( currentEnergy_j - currentEnergy_i );
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if (rnd < p) {
        if (debug) core::clasp_write_string(fmt::format("Swapping states {} and {}\n", i, j ));
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


CL_LAMBDA(energies betas &key (temperature-swap-steps 100) (max-iterations 1000000) accept-callback (flatness-steps 100) (wl-scaling 0.8) (wl-increment-stop 0.1) debug);
CL_DEFUN core::T_mv chem__voelz_optimize_monomer_corrections(core::T_sp tenergies, core::SimpleVector_double_sp betas, size_t temperature_swap_steps, size_t max_iterations, core::T_sp acceptCallback, size_t flatness_steps, double wl_scaling, double wl_increment_stop, bool debug ) {
  size_t beta_size = core::cl__length(betas);
  double wl_increment = 20.0;
  double flatness = 0.7;
  Energies energies(tenergies);
  // Number of slots in state
  size_t num_slots_in_state = core::cl__length(energies._MonomerLocusMaxMrkindex );
  std::mt19937 rng(123456u);
  std::uniform_real_distribution<double> U01(0.0,1.0);
  std::uniform_int_distribution<size_t>  pickSlot(0, num_slots_in_state-1 );
  std::uniform_int_distribution<size_t>  tswapPickSlot(0, beta_size-2 );
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
    if (debug) core::clasp_write_string(fmt::format("Top of loop at iter {}\n", iter));
    for ( size_t betai = 0; betai<beta_size; betai++ ) {
      double beta = (*betas)[betai];
      size_t slotIndex = pickSlot(rng);
      double currentEnergy = energies.energyFunction(currentStates[betai],beta);
      testState = energies.randomStep( slotIndex, currentStates[betai] );
      double testEnergy = energies.energyFunction(testState,beta);
      double x = -(*betas)[betai] * (testEnergy - currentEnergy);
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if ( rnd <= p) {
        accepts++;
        uint32_t curLmkeyIndex = (*energies._MrkeyIndexToLmkeyIndex)[testState._State[slotIndex]];
        if (debug) {
          core::clasp_write_string(fmt::format("Accepted (rnd={} <= p={}) mc step iter {} betai {}  curLmkeyIndex {}\n", rnd, p, iter, betai, curLmkeyIndex ));
        }
        // Increment the histogram and the MonomerCorrection for the current monomer
        histograms[betai][curLmkeyIndex]++;
        if (debug) {
          core::clasp_write_string(fmt::format("Histogram(betai={:2}): ", betai ));
          for ( size_t ii=0; ii<histograms[betai].size(); ii++ ) {
            if (ii==curLmkeyIndex) core::clasp_write_string("*");
            core::clasp_write_string(fmt::format("{:12} ", histograms[betai][ii]));
          }
          core::clasp_write_string(fmt::format("\n"));
          core::clasp_write_string(fmt::format("Monomer corrections: "));
          for ( size_t ii=0; ii<energies._MonomerCorrectionsLength; ii++ ) {
            if (ii==curLmkeyIndex) core::clasp_write_string("*");
            core::clasp_write_string(fmt::format("{:12.2f} ", (*energies._MonomerCorrections)[ii]));
          }
          core::clasp_write_string(fmt::format("\n"));
        }
        {
          double beta_increment = wl_increment/(*betas)[betai];
          (*energies._MonomerCorrections)[curLmkeyIndex] += beta_increment;
          if (debug)
            core::clasp_write_string(fmt::format("Updating _MonomerCorrections by Updated MonomerCorrections at iter {} betai {} curLmkeyIndex {} beta={} wl_increment {} beta_increment {}\n", iter, betai, curLmkeyIndex, beta, wl_increment, beta_increment ));
        }
        // Reset the bias
        adjustMonomerCorrectionsBias(energies);
        if (debug) core::clasp_write_string(fmt::format("Updated MonomerCorrections at iter {} betai {} curLmkeyIndex {} with {}\n", iter, betai, curLmkeyIndex, wl_increment ));
        currentStates[betai] = testState;
        if (useAcceptCallback) {
          memcpy(&(*saveState)[0],&currentStates[betai]._State[0],sizeof(int32_t)*energies._NumberOfSlots);
          core::eval::funcall( acceptCallback, core::make_fixnum(iter), core::make_fixnum(betai), saveState, mk_double_float(testEnergy), energies._MonomerCorrections );
        }
      }
    }
    // Check the histogram for flatness
    flatness_iter++;
    if (flatness_iter%flatness_steps==0) {
      size_t visits = 0, hmin = SIZE_MAX, hmax = 0;
      size_t flat_count = 0;
      size_t bi = 0;
      for (auto histogram : histograms ) {
        for (auto h : histogram) { visits += h; hmin = std::min(hmin,h); hmax = std::max(hmax,h); }
        if (visits >= 50* histogram.size()) {
          double avg = double(visits) / histogram.size();
          if (debug) core::clasp_write_string(fmt::format("Enough visits for bi={} hmin = {}  avg = {} hmin/avg={}  flatness = {}  {}\n", bi, hmin, avg, hmin/avg, flatness, (double(hmin)/avg>flatness) ? "FLAT" : "not-flat" ));
          if (double(hmin)/avg > flatness) flat_count++;
        }
        bi++;
      }
      if (debug) core::clasp_write_string(fmt::format("Flat_count = {}  histograms.size() = {}\n", flat_count, histograms.size()));
      if (flat_count == histograms.size()) {
        flatness_reached++;
        if (debug) core::clasp_write_string(fmt::format("All histograms flat enoughwl_increment updated to {}\n", wl_increment ));
        for ( auto histogram : histograms ) {
          std::fill( histogram.begin(), histogram.end(), 0 );
        }
        wl_increment *= wl_scaling;
      }
    }
    if (wl_increment<wl_increment_stop) {
      if (debug) core::clasp_write_string(fmt::format("wl_increment {} is < wl_increment_stop {}\n", wl_increment, wl_increment_stop ));
      goto DONE;
    }
    if (temperature_swap_steps && (currentStates.size()>=2) && (iter % temperature_swap_steps)==0) {
      // Try the temperature swap
      size_t i = tswapPickSlot(rng), j = i+1;
      double currentEnergy_i = energies.energyFunction(currentStates[i],(*betas)[i]);
      double currentEnergy_j = energies.energyFunction(currentStates[j],(*betas)[j]);
      double x = ( (*betas)[i] - (*betas)[j] ) * ( currentEnergy_j - currentEnergy_i );
      double p = acc_prob_from_log(x);
      double rnd = U01(rng);
      if (rnd < p) {
        if (debug) core::clasp_write_string(fmt::format("Swapping states {} and {}\n", i, j ));
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




}; // namespace chem


