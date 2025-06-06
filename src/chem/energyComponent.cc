/*
    File: energyComponent.cc
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

#include <clasp/core/common.h>
#include <cando/chem/scoringFunction.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/energyComponent.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#
// last include is wrappers.h
#include <clasp/core/wrappers.h>

namespace chem {








void	EnergyComponent_O::initialize()
{
  this->Base::initialize();
}


CL_LISPIFY_NAME("debugLogAsString");
CL_DEFMETHOD string EnergyComponent_O::debugLogAsString()
{
#if TURN_ENERGY_FUNCTION_DEBUG_ON == 1
  return this->_DebugLog.str();
#else
  return this->className()+"-no_log";
#endif
}

void EnergyComponent_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,enabled), this->_Enabled );
  node->field( INTERN_(kw,scale), this->_Scale);
  this->Base::fields(node);
}


string	EnergyComponent_O::enabledAsString() 
{
  stringstream	ss;
  ss.str("");
  if ( this->isEnabled() ) 
  {
    ss<<this->className();
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
    if ( this->_Debug_NumberOfTermsToCalculate < 0 ) {
      ss<<"(debug_All)";
    } else {
      ss<<"(debug_";
      ss<<this->_Debug_NumberOfTermsToCalculate<<")";
    }
#endif
    ss<<" ";
  }
  return ss.str();
};


CL_DOCSTRING(R"dx(Evaluate the energy of a component)dx");
CL_LAMBDA(energy-function component pos &key energy-scale active-atom-mask debug-interactions);
DOCGROUP(cando);
CL_DEFUN
num_real chem__energy_component_evaluate_energy(EnergyFunction_sp energy_function,
                                                EnergyComponent_sp component,
                                                chem::NVector_sp pos,
                                                core::T_sp energyScale,
                                                core::T_sp activeAtomMask,
                                                core::T_sp debugInteractions 
                                                )
{
  num_real val = component->evaluateAllComponent(energy_function,
                                                 pos,
                                                 energyScale,
                                                 nil<core::T_O>(),
                                                 false,nil<NVector_O>(),
                                                 false,false,
                                                 nil<AbstractLargeSquareMatrix_O>(),
                                                 nil<NVector_O>(),
                                                 nil<NVector_O>(),
                                                 activeAtomMask,
                                                 debugInteractions );
  return val;
};

CL_DOCSTRING(R"dx(Evaluate the energy and force of a component)dx");
CL_LAMBDA(energy-function component pos &key energy-scale force active-atom-mask debug-interactions);
DOCGROUP(cando);
CL_DEFUN
num_real chem__energy_component_evaluate_energy_force(EnergyFunction_sp energy_function,
                                                      EnergyComponent_sp component,
                                                      NVector_sp pos,
                                                      core::T_sp energyScale,
                                                      NVector_sp force,
                                                      core::T_sp activeAtomMask,
                                                      core::T_sp debugInteractions )
{
  num_real val = component->evaluateAllComponent(energy_function,
                                                 pos,
                                                 energyScale,
                                                 nil<core::T_O>(),
                                                 true,force,
                                                 false,false,
                                                 nil<AbstractLargeSquareMatrix_O>(),
                                                 nil<NVector_O>(),
                                                 nil<NVector_O>(),
                                                 activeAtomMask,
                                                 debugInteractions );
  return val;
};

EnergyComponent_sp EnergyComponent_O::filterCopyComponent(core::T_sp keepInteractionFactory) {
  IMPLEMENT_ME();
}

void copyEnergyComponent(EnergyComponent_sp newComponent, EnergyComponent_sp orig )
{
  newComponent->_Enabled = orig->_Enabled;
  newComponent->_Scale = orig->_Scale;
}

};
