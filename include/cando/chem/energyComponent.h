/*
    File: energyComponent.h
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
#define	DEBUG_LEVEL_FULL

//
// (C) 2004 Christian E. Schafmeister
//

/*
 *	energyComponent.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyComponent_H 
#define	EnergyComponent_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/quickDom.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/nVector.h>
#include <cando/chem/scoringFunction.fwd.h>
#include <cando/chem/residue.h>

		//
		// SCALE FACTOR FOR ELECTROSTIC INTERACTIONS???
		//
const double ELECTROSTATIC_MODIFIER = (93865.10/3.53*93865.10/93958.78);

#define	DECLARE_FLOAT(x) num_real x = 0.0

	// If we declare mathematica generated term variables
	// on the fly then set ON_THE_FLY_DECLARES
	// and variable declarations won't be included
	//
	// This should avoid a lot of variables being declared but
	// not used and eliminate a lot of compiler warnings
	//
#define	USE_EXPLICIT_DECLARES 1



#ifndef	delta2
#define		delta2	0.00001
#endif
#include <cando/adapt/quickDom.fwd.h>// energyComponent.h wants QDomNode needs quickDom.fwd.h

namespace       chem {


#define MAYBE_SETUP_ACTIVE_ATOM_MASK() \
  bool hasActiveAtomMask = (activeAtomMask.notnilp()); \
  auto bitvectorActiveAtomMask = gc::As_unsafe<core::SimpleBitVector_sp>(activeAtomMask);

SMART(QDomNode);
SMART(AbstractLargeSquareMatrix);
FORWARD(EnergyFunction);

#if TURN_ENERGY_FUNCTION_DEBUG_ON == 1
#define	LOG_ENERGY_CLEAR() {this->_DebugLog.str("");};
#define LOG_ENERGY(...) {fmt::print("{}:{}:{} ", __FILE__, __LINE__, __FUNCTION__ ); fmt::print(__VA_ARGS__);}
#else
#define	LOG_ENERGY_CLEAR() {}
#define LOG_ENERGY(...) {}
#endif


class EnergyTerm 
{
 private:
 public:
  virtual string	className() {THROW_HARD_ERROR("Subclass must implement");};

  EnergyTerm() {};
};


class BeyondThresholdEnergyTerm
{
 public:
};



inline string	atomLabel(Atom_sp a)
{
  stringstream	sstr;
  sstr << _rep_(a->getName());
  return sstr.str();
}

#ifndef VERYSMALL
#define VERYSMALL       0.000000000001
#endif
#define VERYSMALLSQUARED       (VERYSMALL*VERYSMALL)
#define	TENM3		0.001


#ifndef	SIGN
#define	SIGN(x)	((x)<0.0?-1.0:1.0)
#endif
#define	myPow(x,y)	(pow(x,y))
#define	ArcCos(x)	(acos(x))
#define	Cos(x)		(cos(x))
#define	Sin(x)		(sin(x))

#define	TWOPI		(3.14159265*2.0)
#define	reciprocal(x)	(1.0/(x))
#define	mysqrt(x)	sqrt(x)
#define	power1(x)	(x)
#define	power2(x)	((x)*(x))
#define	c_two		2.0
#define	c_six		6.0
#define	oneHalf		0.5
#define	c_twelve	12.0
#define	c_True		true
#define	c_False		false

	// ChiralOffset is a term added in the chiral restraint energy function
	// to try and prevent the four atoms in a chiral restraint from
	// settling into the same plane
#define	ChiralOffset	0.0



#define	TEST_FORCE( func, delta, argLow, argHigh, term, idx ) {	\
    num_real eLow = func argLow;\
    num_real eHigh = func argHigh;\
    num_real numForce = -(eHigh-eLow)/(delta);\
    LOG(("eHigh = %le eLow = %le delta = %le") , (eHigh) , (eLow) , (delta) ); \
    LOG(("numForce = %le") , numForce );			\
    if ( !_areValuesClose( numForce, term, #func, #term, idx ) ) { fails++; };\
  }

#define	TEST_DIAGONAL_HESSIAN( func, delta, argLow, argMiddle, argHigh, term, idx ) {\
    num_real eLow = func argLow;\
    num_real eMiddle = func argMiddle;\
    num_real eHigh = func argHigh;\
    num_real numHessian = (eHigh+eLow-2.0*eMiddle)/(((delta)/2.0)*((delta)/2.0));\
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
  }

#define	TEST_OFF_DIAGONAL_HESSIAN( func, delta, argMM, argPM, argMP, argPP, term, idx ) {\
    num_real eMM = func argMM;\
    num_real eMP = func argMP;\
    num_real ePM = func argPM;\
    num_real ePP = func argPP;\
    num_real numHessian = ((ePP-ePM)-(eMP-eMM))/(delta*delta);\
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
  }


#define	ForceAcc(i,o,v) {\
      if ( hasForce ) {\
        ASSERT_NOT_NAN(v); \
        force->setElement((i)+(o),(v)+force->getElement((i)+(o)));\
      }\
  }

//
// Accumulate an off diagonal Hessian element
//
#define	OffDiagHessAcc(i1,o1,i2,o2,v) {\
    if ( hasHessian ) {\
      ASSERT_NOT_NAN(v); \
      hessian->addToElement((i1)+(o1),(i2)+(o2),v);\
    }\
    if ( hasHdAndD ) {\
      auto v22 = v*dvec->element((i2)+(o2));\
      auto v11 = v*dvec->element((i1)+(o1));\
      ASSERT_NOT_NAN(v22);\
      ASSERT_NOT_NAN(v11);\
      hdvec->addToElement((i1)+(o1),v22); \
      hdvec->addToElement((i2)+(o2),v11); \
    }\
  }

//
// Accumulate a diagonal Hessian element
//
#define	DiagHessAcc(i1,o1,i2,o2,v) {\
    if ( hasHessian ) {\
      ASSERT_NOT_NAN(v); \
      hessian->addToElement((i1)+(o1),(i2)+(o2),v);\
    }\
    if ( hasHdAndD ) {\
      auto vd = v*dvec->element((i1)+(o1));\
      ASSERT_NOT_NAN(vd); \
      hdvec->addToElement((i1)+(o1),vd); \
    }\
  }



#define	ALL_ENERGY_COMPONENTS(msg) {		\
    this->_Stretch->msg; 			\
      hdvec->addToElement((i1)+(o1),v*dvec->element((i2)+(o2)));\
      hdvec->addToElement((i2)+(o2),v*dvec->element((i1)+(o1)));\
    }\
  }

#if 0
// This is a duplicate of the one above - why is it here?  Delete it
//
// Accumulate a diagonal Hessian element
//
#define	DiagHessAcc(i1,o1,i2,o2,v) {\
      if ( hasHessian ) {\
          hessian->addToElement((i1)+(o1),(i2)+(o2),v);\
      }\
      if ( hasHdAndD ) {\
          hdvec->addToElement((i1)+(o1),v*dvec->element((i1)+(o1)));\
      }\
  }
#endif


#define	ALL_ENERGY_COMPONENTS(msg) {		\
    this->_Stretch->msg; 			\
    this->_Angle->msg; 			\
    this->_Dihedral->msg; 			\
    this->_Nonbond->msg; 			\
    this->_DihedralRestraint->msg; 		\
    this->_ChiralRestraint->msg; 		\
    this->_AnchorRestraint->msg; 		\
    this->_FixedNonbondRestraint->msg; 	\
  }

template <class ComponentType, class EntryType>
  string	component_beyondThresholdInteractionsAsString(ComponentType& component) 
{
  int	bt;
  stringstream	ss;
  bt = component._BeyondThresholdTerms.end()-
    component._BeyondThresholdTerms.begin();
  ss << component.className() << "(#"<< bt <<") ";
  return ss.str();
};




SMART(EnergyComponent );
class EnergyComponent_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,EnergyComponent_O,"EnergyComponent",core::CxxObject_O);
 public:
 public: // virtual functions inherited from Object
  void initialize();
 protected: // instance variables
  bool		        _Enabled;
  double		_Scale;
  double		_ErrorThreshold;
  bool		_DebugEnergy;
  int		_Debug_NumberOfTermsToCalculate;
  core::StringOutputStream_sp 	_DebugLog;
public:
  size_t        _Evaluations;
//protected:		// Define these in subclasses
//	vector<TermClass>	_Terms;
//	vector<TermClass>	_BeyondThresholdTerms;
 public:
  CL_DEFMETHOD virtual size_t numberOfTerms() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
  void setScale(double s) {this->_Scale = s; };
  double getScale() { return this->_Scale ; };
  CL_DEFMETHOD bool isEnabled() const { return this->_Enabled; };
  CL_LISPIFY_NAME("enable");
  CL_DEFMETHOD 	void enable() {this->_Enabled = true; };
  CL_LISPIFY_NAME("disable");
  CL_DEFMETHOD 	void disable() {this->_Enabled = false; };
  CL_LISPIFY_NAME("getDebugOn");
  CL_DEFMETHOD 	bool getDebugOn() {return this->_DebugEnergy; };
  CL_LISPIFY_NAME("enableDebug");
  CL_DEFMETHOD 	void enableDebug() {this->_DebugEnergy = true; };
  CL_LISPIFY_NAME("disableDebug");
  CL_DEFMETHOD 	void disableDebug() {this->_DebugEnergy = false; };
  CL_LISPIFY_NAME("setDebug_NumberOfTermsToCalculate");
  CL_DEFMETHOD 	void setDebug_NumberOfTermsToCalculate(int i) {this->_Debug_NumberOfTermsToCalculate=i;}
  CL_LISPIFY_NAME("setErrorThreshold");
  CL_DEFMETHOD 	void	setErrorThreshold(double tr) { this->_ErrorThreshold = tr; };
  CL_LISPIFY_NAME("getErrorThreshold");
  CL_DEFMETHOD 	double	getErrorThreshold() { return this->_ErrorThreshold; };

  CL_DOCSTRING("Some energy-components are restraints - and should be disabled to calculate conformational energies");
  CL_DEFMETHOD
    virtual bool is_restraint() const { return false; };
  
  CL_LISPIFY_NAME("energy-component-evaluations");
  CL_DEFMETHOD 	size_t	evaluations() const { return this->_Evaluations; };

  CL_DEFMETHOD virtual core::List_sp extract_vectors_as_alist() const { SUBCLASS_MUST_IMPLEMENT(); };
 
  string enabledAsString();
  string debugLogAsString();

 public:	// Virtual methods

  CL_DEFMETHOD virtual void dumpTerms(core::HashTable_sp atomTypes) {_OF();SUBCLASS_MUST_IMPLEMENT();};

  virtual EnergyComponent_sp filterCopyComponent(core::T_sp keepInteraction);
  
  CL_DEFMETHOD virtual	num_real evaluateAllComponent( ScoringFunction_sp scorer,
                                                       NVector_sp 	pos,
                                                       core::T_sp componentEnergy,
                                                       bool 		calcForce,
                                                       gc::Nilable<NVector_sp> 	force,
                                                       bool		calcDiagonalHessian,
                                                       bool		calcOffDiagonalHessian,
                                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                       gc::Nilable<NVector_sp>	hdvec,
                                                       gc::Nilable<NVector_sp> dvec,
                                                       core::T_sp activeAtomMask,
                                                       core::T_sp debugInteractions ) = 0;

  virtual core::List_sp checkForBeyondThresholdInteractionsWithPosition(NVector_sp pos, double threshold ) {_OF();SUBCLASS_MUST_IMPLEMENT();};

  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos ) {_OF();SUBCLASS_MUST_IMPLEMENT();};
 public:
  EnergyComponent_O( const EnergyComponent_O& ss ); //!< Copy constructor

  EnergyComponent_O() : _Enabled(true),
                        _Scale(1.0),
                        _DebugEnergy(false), _Evaluations(0) {};
};
template <typename SP>
SP safe_alist_lookup(core::List_sp list, core::T_sp key) {
  if (list.consp()) {
    core::T_sp entry = list.unsafe_cons()->assoc(key,nil<core::T_O>(),nil<core::T_O>(),nil<core::T_O>());
    if (entry.consp()) {
      return gctools::As<SP>(CONS_CDR(entry));
    }
  }
  SIMPLE_ERROR("Could not find {}" , _rep_(key));
}
  
;
};


TRANSLATE(chem::EnergyComponent_O);
#endif
