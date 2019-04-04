/*
    File: energyNonbond.h
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
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyNonbond.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyNonbond_H  //[
#define	EnergyNonbond_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
//#include "geom/render.fwd.h"// energyNonbond.h wants DisplayList needs render.fwd.h

namespace       chem {
  FORWARD(EnergyFunction); // Declares class EnergyFunction_O {} and EnergyFunction_sp
  FORWARD(EnergyNonbond);
  FORWARD(FFNonbondDb);
  FORWARD(AtomTable);

struct TermNonBond
{
	double		dQ1Q2;
	double		dA;
	double		dC;
	int		I1; //!< i*3 index into coordinate vector, must match Mathematica code!
	int		I2; //!< i*3 index into coordinate vector, must match Mathematica code!
};


inline	string	XmlTag_Nonbond() { return "Nonbond"; };
inline	string	XmlTag_EnergyNonbond() { return "EnergyNonbond"; };

/*! A Nonbond term, this is subclassed in AmberFunction to deal with Amber specific stuff.
 * I should combine energyFunction.cc and amberFunction.cc together.
 */
class EnergyNonbond : public EnergyTerm {
public:
	string	className()	{ return "EnergyNonbond"; };
public:
                // Parameters
                // Variables
	bool		_Is14;
        Atom_sp      	_Atom1;
        Atom_sp      	_Atom2;
	double		_A;
	double		_C;
	double		_Charge1;
	double		_Charge2;
	double		_RStar;
	TermNonBond	term;
//	double		_NonBondScale; // Used in Amber
//	double		_ElectrostaticScale; // Used in Amber
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_Nonbond_debugEvalDeclares.cc>
#endif

  Atom_sp	getAtom1() { return this->_Atom1; };
  Atom_sp	getAtom2() { return this->_Atom2; };
  double	getDistance();
  bool	defineFrom(core::T_sp	forceField,
                   bool		is14,
                   EnergyAtom	*iea1,
                   EnergyAtom	*iea2,
                   EnergyNonbond_sp nb);

public:
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlRelativeToContainer( adapt::QDomNode_sp xml, Matter_sp parent );
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );

public:
        core::List_sp encode() const;
        void decode(core::List_sp alist);
};
};

namespace translate {

  template <>
struct	to_object<chem::EnergyNonbond>
  {
    typedef	core::Cons_sp ExpectedType;
    typedef	core::Cons_sp DeclareType;
    static core::T_sp convert(const chem::EnergyNonbond& nonbond)
    {
      return nonbond.encode();
    }
  };
  
  template <>
    struct	from_object<chem::EnergyNonbond>
  {
    typedef	chem::EnergyNonbond	ExpectedType;
    typedef	ExpectedType 		DeclareType;
    DeclareType _v;
    from_object(core::T_sp o)
    {
      SIMPLE_ERROR(BF("Implement me"));
    }
  };
};

namespace chem {

double	_evaluateEnergyOnly_Nonbond(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double dA, double dC, double dQ1Q2 );


class EnergyNonbond_O : public EnergyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyNonbond_O,"EnergyNonbond",EnergyComponent_O);
  
 public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
  
 public: // virtual functions inherited from Object
  void	initialize();
//	string	__repr__() const;
 public:
  typedef EnergyNonbond TermType;
 public: // instance variables
  double		_DielectricConstant;
  double		_ScaleVdw;
  double		_ScaleElectrostatic;
  double		_EnergyVdw;
  double		_EnergyElectrostatic;
  bool                _UsesExcludedAtoms;
    // Original way of defining nonbonds with list of nonbond terms
  gctools::Vec0<TermType>	_Terms;
  gctools::Vec0<TermType>	_BeyondThresholdTerms;
    // Correct way of defining nonbonds using excluded atom indices
  // FIXME:  Get rid of _FFNonbondDb - and the old code that uses it.
  core::T_sp        _FFNonbondDb;
  AtomTable_sp          _AtomTable;
  //  Tables for nonbonded calculation using excluded atoms and the Amber way
  size_t                     _ntypes;          // ntypes
  core::SimpleVector_sp      _atom_name_vector;  // atom-name-vector
  core::SimpleVector_double_sp     _charge_vector;          // charge-vector
  core::SimpleVector_double_sp     _mass_vector;            // masses
  core::SimpleVector_int32_t_sp   _atomic_number_vector;    // vec

  // Nonbond parameters generated by generate-nonbond-parameters
  core::SimpleVector_int32_t_sp   _ico_vec;             // ico-vec
  core::SimpleVector_int32_t_sp   _iac_vec;             // iac-vec
  core::SimpleVector_sp      _atom_type_vector;    // Amber atom type names
  core::SimpleVector_int32_t_sp   _local_typej_vec;      // local-typej-vec
  core::SimpleVector_double_sp     _cn1_vec;
  core::SimpleVector_double_sp     _cn2_vec;
  // Excluded atom table
  core::SimpleVector_int32_t_sp   _NumberOfExcludedAtomIndices;
  core::SimpleVector_int32_t_sp   _ExcludedAtomIndices;
 public:	
  typedef gctools::Vec0<TermType>::iterator iterator;
  iterator begin() { return this->_Terms.begin(); };
  iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
 public:
  //core::List_sp termAtIndex(size_t index) const;
  virtual size_t numberOfTerms() { return this->_Terms.size();};
 public:
  void	setDielectricConstant(double d) { this->_DielectricConstant = d; };
  double	getDielectricConstant() { return this->_DielectricConstant; };
  CL_DEFMETHOD void setVdwScale(double d) { this->_ScaleVdw = d; };
  CL_DEFMETHOD double	getVdwScale()	{return this->_ScaleVdw; };
  CL_DEFMETHOD void	setElectrostaticScale(double d) { this->_ScaleElectrostatic = d; };
  CL_DEFMETHOD double	getElectrostaticScale()	{return this->_ScaleElectrostatic; };

  double	getVdwEnergy() { return this->_EnergyVdw; };
  double	getElectrostaticEnergy() { return this->_EnergyElectrostatic; };
  CL_DEFMETHOD core::SimpleVector_int32_t_sp number_excluded_atoms() const { return this->_NumberOfExcludedAtomIndices;}
  CL_DEFMETHOD core::SimpleVector_int32_t_sp excluded_atom_list() const { return this->_ExcludedAtomIndices;}
 public:
  void addTerm(const TermType& term);
  virtual void dumpTerms();
  virtual	void	zeroEnergy();

  virtual core::List_sp extract_vectors_as_alist() const;
  
  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m );
    
  virtual double evaluateAll( NVector_sp 	pos,
                              bool 		calcForce,
                              gc::Nilable<NVector_sp> 	force,
                              bool		calcDiagonalHessian,
                              bool		calcOffDiagonalHessian,
                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                              gc::Nilable<NVector_sp>	hdvec,
                              gc::Nilable<NVector_sp> dvec);

  void evaluateUsingExcludedAtoms( NVector_sp 	pos,
                                   bool 		calcForce,
                                   gc::Nilable<NVector_sp> 	force,
                                   bool		calcDiagonalHessian,
                                   bool		calcOffDiagonalHessian,
                                   gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                   gc::Nilable<NVector_sp>	hdvec,
                                   gc::Nilable<NVector_sp> dvec);
  void evaluateTerms( NVector_sp 	pos,
                      bool 		calcForce,
                      gc::Nilable<NVector_sp> 	force,
                      bool		calcDiagonalHessian,
                      bool		calcOffDiagonalHessian,
                      gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                      gc::Nilable<NVector_sp>	hdvec,
                      gc::Nilable<NVector_sp> dvec);
  
  void expandExcludedAtomsToTerms();

  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                               NVector_sp pos );

  virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

  virtual string	beyondThresholdInteractionsAsString();

//    int countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn, core::Lisp_sp );

  virtual	double	getEnergy();
  core::T_sp getFFNonbondDb();

  void constructFromAtomTable(bool useExcludedAtoms, AtomTable_sp atomTable, core::T_sp nbforceField, bool show_progress);
  void constructNonbondTermsFromAtomTable(bool ignore14s, AtomTable_sp atomTable, core::T_sp nbforceField, bool show_progress);
  void construct14InteractionTerms(AtomTable_sp atomTable, Matter_sp matter, core::T_sp nbforceField, core::T_sp activeAtoms, bool show_progress);
  void constructExcludedAtomListFromAtomTable(AtomTable_sp atomTable, core::T_sp nbforceField, bool show_progress);

  void constructNonbondTermsFromAtomTableUsingExcludedAtoms(EnergyFunction_sp energyFunction,
                                                            core::T_sp prepareAmberEnergyNonbond );

  void constructNonbondTermsFromAList(core::List_sp values);
  void setNonbondExcludedAtomInfo(AtomTable_sp atom_table, core::SimpleVector_int32_t_sp excluded_atoms_list, core::SimpleVector_int32_t_sp number_excluded_atoms);
 public:
  EnergyNonbond_O( const EnergyNonbond_O& ss ); //!< Copy constructor

  EnergyNonbond_O() : _FFNonbondDb(_Unbound<core::T_O>()), _UsesExcludedAtoms(true) {};
};

};

TRANSLATE(chem::EnergyNonbond_O);
#endif //]
