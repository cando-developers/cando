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
SYMBOL_EXPORT_SC_( ChemPkg, monomer_locus_max_mrkindex );
SYMBOL_EXPORT_SC_( ChemPkg, single_scan_energy );
SYMBOL_EXPORT_SC_( ChemPkg, pair_scan_energy_lower_triangular_matrix );


struct State {
  typedef std::vector<uint32_t> StateType;
  StateType _State;

  State() {};
  
  // Use a defined state
  State(core::SimpleVector_byte32_t_sp mcstate) {
    for ( size_t ii=0; ii<mcstate->length(); ii++ ) {
      this->_State[ii] = (*mcstate)[ii];
    }
  }

  void resize(size_t sz) {
    this->_State.resize(sz);
  }

  void randomState(core::SimpleVector_byte32_t_sp sv_monomerLocusMaxMrkindex ) {
    this->_State.resize(sv_monomerLocusMaxMrkindex->length());
    int prevMax = -1;
    uint32_t* monomerLocusMaxMrkindex = &(*sv_monomerLocusMaxMrkindex)[0]; 
    for ( size_t ii=0; ii<this->_State.size(); ii++ ) {
      int max = monomerLocusMaxMrkindex[ii];
      int range = max - prevMax;
      this->_State[ii] = (rand() % range) + prevMax + 1;
      prevMax = max;
    }
  }

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
    this->_NumberOfSlots = core::cl__length(this->_MonomerLocusMaxMrkindex);
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

  double energyFunction(State& state) {
    double* singleTerms = &(*this->_SingleTerms)[0];
    double* pairTerms = &(*this->_PairTerms)[0];
    KahanAccumulator sum;
    sum.Add(this->_BackboneEnergy);
    for ( int ii=0; ii<state._State.size(); ii++ ) {
      int index = state._State[ii];
      double energyTerm = singleTerms[index];
      sum.Add(energyTerm);
    }
    for ( int xii=0; xii<state._State.size()-1; xii++ ) {
      int xxindex = state._State[xii];
      for ( int yii=xii+1; yii<state._State.size(); yii++ ) {
        int yyindex = state._State[yii];
        int ltmIndex = this->lowerTriangularIndex(xxindex,yyindex);
        double energyTerm = pairTerms[ltmIndex];
        sum.Add(energyTerm);
      }
    }
    return sum.sum;
  }

  State randomStep(size_t slotIndex, const State& state) {
    State tempState(state);
    int index = slotIndex;
    int max = (*this->_MonomerLocusMaxMrkindex)[index];
    int prevMax = (index==0) ? -1 : (*this->_MonomerLocusMaxMrkindex)[index-1];
    int range = max - prevMax;
    tempState._State[index] = (rand() % range) + prevMax + 1;
    return tempState;
  }

// Temperature decay function
  double tempDecay(const double& temperature) {
      // Replace with your temperature decay logic
    return temperature * 0.99;
  }
};

CL_DEFUN core::DoubleFloat_sp chem__mcstate_energy(core::T_sp tenergies, core::T_sp tmcstate ) {
  Energies energies(tenergies);
  State state(gc::As<core::SimpleVector_byte32_t_sp>(tmcstate));
  double testEnergy = energies.energyFunction(state);
  return core::DoubleFloat_O::create(testEnergy);
}

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
    initial_state.randomState(energies._MonomerLocusMaxMrkindex);
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
    initial_state.randomState(energies._MonomerLocusMaxMrkindex);
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

}; // namespace chem


