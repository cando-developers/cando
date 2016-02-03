       
       

//
// (C) 2004 Christian E. Schafmeister
//


/*
 	coupling.h


 */

#ifndef COUPLING_H
#define COUPLING_H

#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <cando/geom/vector3.h>
#include <cando/chem/residue.fwd.h>


#include <cando/chem/chemPackage.h>


namespace chem {

SMART(Monomer);
SMART(Coupling);
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
	friend class O_WeakOligomer;
    public:
	void initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
    public:

    private:
        gc::Nilable<Oligomer_sp>	_Oligomer;
	bool			_HasError;
    private:
	bool			_Selected;
    public:

	virtual core::Symbol_sp getName() const {_OF();SUBCLASS_MUST_IMPLEMENT();};

	bool		containsMonomer(Monomer_sp mon);

	virtual Monomer_sp	getMonomer1() {_OF();SUBCLASS_MUST_IMPLEMENT();};
	virtual Monomer_sp	getMonomer1_const() const {_OF();SUBCLASS_MUST_IMPLEMENT();};
	virtual Monomer_sp	getMonomer2() {_OF();SUBCLASS_MUST_IMPLEMENT();};
	virtual Monomer_sp	getMonomer2_const() const {_OF();SUBCLASS_MUST_IMPLEMENT();};

	void	setOligomer(Oligomer_sp o);
	void	resetOligomer();
	Oligomer_sp	getOligomer();

CL_NAME("isSelected");
CL_DEFMETHOD 	bool isSelected() { return this->_Selected; };
CL_NAME("setSelected");
CL_DEFMETHOD 	void setSelected(bool b) { this->_Selected = b; };

	void	setHasError(bool b) { this->_HasError = b;};
CL_NAME("getHasError");
CL_DEFMETHOD 	bool	getHasError() { return this->_HasError; };

        virtual bool isRingClosing() { return false; };
CL_NAME("throwIfBadConnections");
CL_DEFMETHOD         virtual void throwIfBadConnections() {_OF();SUBCLASS_MUST_IMPLEMENT();};
CL_NAME("Coupling-checkForErrors");
CL_DEFMETHOD         virtual void checkForErrors(core::T_sp errorStream) {_OF();SUBCLASS_MUST_IMPLEMENT();};

        virtual Monomer_sp getOtherSideMonomer(Monomer_sp mon) {_OF();SUBCLASS_MUST_IMPLEMENT();};

CL_NAME("doCoupling");
CL_DEFMETHOD         virtual void doCoupling(Residue_sp inResidue, Residue_sp outResidue ) {_OF();SUBCLASS_MUST_IMPLEMENT();};


        DEFAULT_CTOR_DTOR(Coupling_O);
    };





typedef	enum	{ NoSide, InSide, OutSide } DirectionalCouplingSide;


SMART(DirectionalCoupling);
class DirectionalCoupling_O : public Coupling_O
{
  LISP_CLASS(chem,ChemPkg,DirectionalCoupling_O,"DirectionalCoupling",Coupling_O);

 public:
  friend class O_WeakOligomer;
//	void	archiveBase(core::ArchiveP node);
  void initialize();


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
  static bool isOutPlugName(core::Symbol_sp plugName)
  {
    return plugName->symbolName()->schar(0) == OUT_PLUG_PREFIX;
  };
  static bool isPlugName(core::Symbol_sp plugName)
  {
    return ((plugName->symbolName()->schar(0) == OUT_PLUG_PREFIX) || (plugName->symbolName()->schar(0) == IN_PLUG_PREFIX) );
  };

 private:
  gc::Nilable<Monomer_sp>		_InMonomer;
  core::Symbol_sp	_Name;
  gc::Nilable<Monomer_sp>		_OutMonomer;

 protected:

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

  void	resetIn();
  void	resetOut();
  bool	isInCouplingToMonomer(Monomer_sp mon);
  void	setInMonomer_NoSignal(Monomer_sp sin);
  void	setOutMonomer_NoSignal(Monomer_sp sout);
  void	setInMonomer(Monomer_sp sin);
  void	setOutMonomer(Monomer_sp sout);
  DirectionalCouplingSide	couplingSideOfMonomer( Monomer_sp mon);

  bool	isRingClosing() { return false; };

  string	description() const;
  void	setName( core::Symbol_sp nm);
  CL_NAME("getName");
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


  DEFAULT_CTOR_DTOR(DirectionalCoupling_O);
};


typedef	enum	{ NoRingSide, Monomer1Side, Monomer2Side } RingCouplingSide;



SMART(RingCoupling);
class RingCoupling_O : public Coupling_O
{
    LISP_CLASS(chem,ChemPkg,RingCoupling_O,"RingCoupling",Coupling_O);

public:
	friend class O_WeakOligomer;
//	void	archiveBase(core::ArchiveP node);
	void initialize();

public:
private:

        gc::Nilable<Monomer_sp>	_Monomer1;
	core::Symbol_sp		_Plug1;
        gc::Nilable<Monomer_sp>	_Monomer2;
	core::Symbol_sp		_Plug2;
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


	DEFAULT_CTOR_DTOR(RingCoupling_O);
};



};
TRANSLATE(chem::Coupling_O);
TRANSLATE(chem::DirectionalCoupling_O);
TRANSLATE(chem::RingCoupling_O);
#endif
