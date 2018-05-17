/*
    File: coupling.h
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
 	coupling.h


 */

#ifndef COUPLING_H
#define COUPLING_H

#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <cando/geom/vector3.h>
#include <cando/chem/residue.fwd.h>


#include <cando/chem/chemPackage.h>


namespace chem {

  SMART(Coupling);
  FORWARD(Monomer);
  SMART(Plug);
  SMART(Oligomer);
  SMART(EntityNameSet);
/*!
	\class Coupling_O
	\brief Coupling_O stores a residue and links to Monomers
*/

  SMART(Coupling);
  class Coupling_O : public core::CxxObject_O //public core::Model_O
  {
    LISP_CLASS(chem,ChemPkg,Coupling_O,"Coupling",core::CxxObject_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    virtual core::Symbol_sp getName() const {_OF();SUBCLASS_MUST_IMPLEMENT();};
    bool		containsMonomer(Monomer_sp mon);
    virtual Monomer_sp	getMonomer1() {_OF();SUBCLASS_MUST_IMPLEMENT();};
    virtual Monomer_sp	getMonomer1_const() const {_OF();SUBCLASS_MUST_IMPLEMENT();};
    virtual Monomer_sp	getMonomer2() {_OF();SUBCLASS_MUST_IMPLEMENT();};
    virtual Monomer_sp	getMonomer2_const() const {_OF();SUBCLASS_MUST_IMPLEMENT();};

    virtual bool isRingClosing() { return false; };
    CL_LISPIFY_NAME("throwIfBadConnections");
    CL_DEFMETHOD         virtual void throwIfBadConnections() {_OF();SUBCLASS_MUST_IMPLEMENT();};
    CL_LISPIFY_NAME("Coupling-checkForErrors");
    CL_DEFMETHOD         virtual void checkForErrors(core::T_sp errorStream) {_OF();SUBCLASS_MUST_IMPLEMENT();};

    virtual Monomer_sp getOtherSideMonomer(Monomer_sp mon) {_OF();SUBCLASS_MUST_IMPLEMENT();};

    CL_LISPIFY_NAME("doCoupling");
    CL_DEFMETHOD         virtual void doCoupling(Residue_sp inResidue, Residue_sp outResidue ) {_OF();SUBCLASS_MUST_IMPLEMENT();};


    DEFAULT_CTOR_DTOR(Coupling_O);
  };


  typedef	enum	{ NoSide, InSide, OutSide } DirectionalCouplingSide;


  SMART(DirectionalCoupling);
  class DirectionalCoupling_O : public Coupling_O
  {
    LISP_CLASS(chem,ChemPkg,DirectionalCoupling_O,"DirectionalCoupling",Coupling_O);

  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    static char otherDirectionalCouplingSide( char side);
    	//! Return just the coupling name with the prefix removed
    static core::Symbol_sp couplingName(core::Symbol_sp name);
        //! Return the inPlugName
    static core::Symbol_sp inPlugName(core::Symbol_sp couplingName);
        //! Return the outPlugName
    static core::Symbol_sp outPlugName(core::Symbol_sp couplingName);
        //! Return the other PlugName
    static core::Symbol_sp otherPlugName(core::Symbol_sp plugName);
    static bool isInPlugName(core::Symbol_sp plugName);
    static bool isOutPlugName(core::Symbol_sp plugName) {
      return plugName->symbolName()->rowMajorAref(0).unsafe_character() == OUT_PLUG_PREFIX;
    };
    static bool isPlugName(core::Symbol_sp plugName) {
      return ((plugName->symbolName()->rowMajorAref(0).unsafe_character() == OUT_PLUG_PREFIX) || (plugName->symbolName()->rowMajorAref(0).unsafe_character() == IN_PLUG_PREFIX) );
    };

  private:
    core::Symbol_sp	_Name;
    Monomer_sp		_InMonomer;
    Monomer_sp		_OutMonomer;
  public:
    bool		containsMonomer();
    Monomer_sp	getInMonomer();
    Monomer_sp	getInMonomer_const() const;
    Monomer_sp	getOutMonomer();
    Monomer_sp	getOutMonomer_const() const;

    virtual Monomer_sp	getMonomer1() { return this->getInMonomer();};
    virtual Monomer_sp	getMonomer1_const() const { return this->getInMonomer_const();};
    virtual Monomer_sp	getMonomer2() { return this->getOutMonomer();};
    virtual Monomer_sp	getMonomer2_const() const { return this->getOutMonomer_const();};

    core::T_sp  getInPlugName() const;
    core::T_sp  getOutPlugName() const;
    bool	isInCouplingToMonomer(Monomer_sp mon);
    void	setInMonomer_NoSignal(Monomer_sp sin);
    void	setOutMonomer_NoSignal(Monomer_sp sout);
    void	setInMonomer(Monomer_sp sin);
    void	setOutMonomer(Monomer_sp sout);
    DirectionalCouplingSide	couplingSideOfMonomer( Monomer_sp mon);

    bool	isRingClosing() { return false; };

    string	description() const;
    void	setName( core::Symbol_sp nm);
    CL_LISPIFY_NAME("getName");
    CL_DEFMETHOD 	core::Symbol_sp	getName() const {return this->_Name;};

    void	couple( Monomer_sp sin, Monomer_sp sout );

	/*! Get the plug name of this Coupling from the point of view of the InMonomer */
    core::Symbol_sp	getInMonomerPlugName();
	/*! Get the plug name of this Coupling from the point of view of the OutMonomer */
    core::Symbol_sp	getOutMonomerPlugName();

    Plug_sp		getPlugForMonomer(Monomer_sp mon);
    Plug_sp		getPlugForOtherMonomer(Monomer_sp mon);

    bool		containsMonomer(Monomer_sp mon);
    Monomer_sp	getOtherSideMonomer(Monomer_sp mon);

    void	doCoupling(Residue_sp inResidue, Residue_sp outResidue );

    void	checkForErrors(core::T_sp errorStream);
    void	throwIfBadConnections();


  public:
    static DirectionalCoupling_sp make(core::Symbol_sp name, Monomer_sp in, Monomer_sp out) {
      GC_ALLOCATE_VARIADIC(DirectionalCoupling_O,dc,name,in,out);
      return dc;
    }
  public:
  DirectionalCoupling_O(core::Symbol_sp name, Monomer_sp in, Monomer_sp out) :
    _Name(name),
      _InMonomer(in),
      _OutMonomer(out) {};
  };

  typedef	enum	{ NoRingSide, Monomer1Side, Monomer2Side } RingCouplingSide;

  SMART(RingCoupling);
  class RingCoupling_O : public Coupling_O
  {
    LISP_CLASS(chem,ChemPkg,RingCoupling_O,"RingCoupling",Coupling_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  private:
    core::Symbol_sp		_Plug1;
    Monomer_sp	                _Monomer1;
    core::Symbol_sp		_Plug2;
    Monomer_sp	                _Monomer2;
  protected:

  public:
    core::Symbol_sp getName() const;

    virtual Monomer_sp	getMonomer1();
    virtual Monomer_sp	getMonomer1_const() const;
    virtual Monomer_sp	getMonomer2();
    virtual Monomer_sp	getMonomer2_const() const;

    void	resetMonomer1();
    void	resetMonomer2();
    void	setMonomer1_NoSignal(Monomer_sp sin);
    void	setMonomer2_NoSignal(Monomer_sp sout);
    void	setMonomer1(Monomer_sp sin);
    void	setMonomer2(Monomer_sp sout);
    void	setPlug1(core::Symbol_sp p);
    core::Symbol_sp getPlug1();
    void	setPlug2(core::Symbol_sp p);
    core::Symbol_sp getPlug2();

    RingCouplingSide	couplingSideOfMonomer( Monomer_sp mon);

    void	throwIfBadConnections();
    void	checkForErrors(core::T_sp errorStream);

    bool	isRingClosing() { return true; };

    string	description() const;

//	void	couple( Monomer_sp sin, Monomer_sp sout );


    Plug_sp		getPlugForMonomer(Monomer_sp mon);
    Plug_sp		getPlugForOtherMonomer(Monomer_sp mon);

    bool		containsMonomer(Monomer_sp mon);
    Monomer_sp	getOtherSideMonomer(Monomer_sp mon);

    void	doCoupling(Residue_sp inResidue, Residue_sp outResidue );

  public:
    static RingCoupling_sp make(Monomer_sp mon1, core::Symbol_sp plug1Name, Monomer_sp mon2, core::Symbol_sp plug2Name) {
      GC_ALLOCATE_VARIADIC(RingCoupling_O,dc,mon1,plug1Name,mon2,plug2Name);
      return dc;
    }

  RingCoupling_O(Monomer_sp mon1, core::Symbol_sp plug1, Monomer_sp mon2, core::Symbol_sp plug2 ) :
    _Plug1(plug1),
      _Monomer1(mon1),
      _Plug2(plug2),
      _Monomer2(mon2) {};
  };



};
#endif
