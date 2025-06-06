/*
    File: restraint.h
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
 *	restraint.h
 *
 *	Maintain a Restraint object.
 *
 *
 */

#ifndef RESTRAINT_H
#define	RESTRAINT_H

#include <iostream>
#include <string>
#include <vector>
#include <clasp/core/common.h>
#include <cando/geom/ovector3.h>
#include <cando/geom/matrix.h>
#include <cando/geom/vector3.h>
#include <clasp/core/sequence.h>


#include <cando/chem/chemPackage.h>


namespace chem {


SMART(Atom);
SMART(Matter);
SMART(Aggregate);
SMART(Restraint);

SMART(Restraint );
class Restraint_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,Restraint_O,"Restraint",core::CxxObject_O);
 public:
//    void	serialize(serialize::SNode node);

 private:
  bool _Active;
 public:

  virtual Restraint_sp	copyDontRedirectAtoms() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
  virtual void redirectAtoms() {_OF(); SUBCLASS_MUST_IMPLEMENT();};

  virtual void invertStereochemistryOfRestraint() {_OF();SUBCLASS_MUST_IMPLEMENT();};

CL_LISPIFY_NAME("isActive");
CL_DEFMETHOD   bool isActive() const { return this->_Active; };
CL_LISPIFY_NAME("setActive");
CL_DEFMETHOD   void setActive(bool a) { this->_Active = a; };

 virtual void applyTransform(const Matrix& m) {};

 Restraint_O() : _Active(true) {};
  virtual ~Restraint_O() {};
};



SMART(RestraintAnchor);
class RestraintAnchor_O : public Restraint_O
{
    LISP_CLASS(chem,ChemPkg,RestraintAnchor_O,"RestraintAnchor",Restraint_O);
 public:
    static RestraintAnchor_sp make(Atom_sp atom, const Vector3& pos, double weight);

public:

	Atom_sp		_Atom;
	Vector3		_Pos;
	double		_Weight;
 public:
        bool fieldsp() const { return true; };
        void	fields(core::Record_sp node);
public:

CL_LISPIFY_NAME("setWeight");
CL_DEFMETHOD 	void setWeight(double d) { this->_Weight = d; };
CL_LISPIFY_NAME("restraint-anchor-getWeight");
CL_DEFMETHOD 	double getWeight() { return this->_Weight; };

CL_LISPIFY_NAME("setAnchorPos");
CL_DEFMETHOD 	void	setAnchorPos(const Vector3& pos )
	{
	    this->_Pos = pos;
	};
	Atom_sp	getAtom() {ASSERTNOTNULL(this->_Atom);return this->_Atom;};
CL_LISPIFY_NAME("setAtom");
CL_DEFMETHOD 	void	setAtom(Atom_sp a) {this->_Atom=a;};

	Vector3 getAnchorPos() { return this->_Pos; };

	Restraint_sp	copyDontRedirectAtoms();
	void redirectAtoms();

    	void invertStereochemistryOfRestraint();

        virtual void applyTransform(const Matrix& m);

        RestraintAnchor_O() {};
 RestraintAnchor_O(Atom_sp atom, const Vector3& pos, double weight) : _Atom(atom), _Pos(pos), _Weight(weight) {};
        virtual ~RestraintAnchor_O() {};
};



SMART(RestraintChiral);
class RestraintChiral_O : public Restraint_O
{
    LISP_CLASS(chem,ChemPkg,RestraintChiral_O,"RestraintChiral",Restraint_O);

public:
	Atom_sp		_A;
	int		_Chirality;	// S=1,R=-1

public:

	void	setParameters(int c)
	{
	    this->_Chirality = c;
	};
	Atom_sp	getAtomA() {_OF();ASSERTNOTNULL(this->_A);return this->_A;};
	void	setAtomA(Atom_sp a) {this->_A=a;};

	int	getChirality() { return this->_Chirality; };
	void	setChirality(int c) { this->_Chirality = c; };

	Restraint_sp	copyDontRedirectAtoms();
	void redirectAtoms();

    	void invertStereochemistryOfRestraint();

	RestraintChiral_O(const RestraintChiral_O& old);
	DEFAULT_CTOR_DTOR(RestraintChiral_O);
};





SMART(RestraintDistance);
class RestraintDistance_O : public Restraint_O
{
  LISP_CLASS(chem,ChemPkg,RestraintDistance_O,"RestraintDistance",Restraint_O);

 public:
  bool fieldsp() const { return true; };
  void	fields(core::Record_sp node);

 public:
  Atom_sp		_A;
  Atom_sp		_B;
  double          _K;
  double          _R0;

 public:
  Atom_sp	getAtomA() {_OF();ASSERTNOTNULL(this->_A);return this->_A;};
  void	setAtomA(Atom_sp a) {this->_A=a;};
  Atom_sp	getAtomB() {_OF();ASSERTNOTNULL(this->_B);return this->_B;};
  void	setAtomB(Atom_sp a) {this->_B=a;};

  Restraint_sp	copyDontRedirectAtoms();
  void redirectAtoms();

  void invertStereochemistryOfRestraint() { /* Nothing */ };
};


SMART(RestraintAngle);
class RestraintAngle_O : public Restraint_O
{
    LISP_CLASS(chem,ChemPkg,RestraintAngle_O,"RestraintAngle",Restraint_O);
public:
	Atom_sp		_A;
	Atom_sp		_B;
	Atom_sp		_C;
	double		_MinRad;
	double		_MaxRad;
	double		_Weight;

public:


	void	setParameters(double mn,double mx, double wt ) {
			this->_MinRad=mn;
			this->_MaxRad=mx;
			this->_Weight=wt;
	};
	Atom_sp	getAtomA() {_OF();ASSERTNOTNULL(this->_A);return this->_A;};
	void	setAtomA(Atom_sp a) {this->_A=a;};
	Atom_sp	getAtomB() {_OF();ASSERTNOTNULL(this->_B);return this->_B;};
	void	setAtomB(Atom_sp a) {this->_B=a;};
	Atom_sp	getAtomC() {_OF();ASSERTNOTNULL(this->_C);return this->_C;};
	void	setAtomC(Atom_sp a) {this->_C=a;};

	double	getMin() {return this->_MinRad;};
	double	getMax() {return this->_MaxRad;};
	double	getWeight() {return this->_Weight;};

	Restraint_sp	copyDontRedirectAtoms();
	void redirectAtoms();

    	void invertStereochemistryOfRestraint() { /* Nothing */ };


	RestraintAngle_O(const RestraintAngle_O& old);
	DEFAULT_CTOR_DTOR(RestraintAngle_O);
};




SMART(RestraintDihedral);
class RestraintDihedral_O : public Restraint_O
{
    LISP_CLASS(chem,ChemPkg,RestraintDihedral_O,"RestraintDihedral",Restraint_O);

public:
	Atom_sp		_A;
	Atom_sp		_B;
	Atom_sp		_C;
	Atom_sp		_D;
	double		_Degrees;
	double		_Weight;

public:


	void	setParameters(double deg, double wt ) {
			this->_Degrees=deg;
			this->_Weight=wt;
	};
	Atom_sp	getAtomA() {_OF();ASSERTNOTNULL(this->_A);return this->_A;};
	Atom_sp	getAtomB() {_OF();ASSERTNOTNULL(this->_B);return this->_B;};
	Atom_sp	getAtomC() {_OF();ASSERTNOTNULL(this->_C);return this->_C;};
	Atom_sp	getAtomD() {_OF();ASSERTNOTNULL(this->_D);return this->_D;};
	void	setAtomA(Atom_sp a) {this->_A=a;};
	void	setAtomB(Atom_sp a) {this->_B=a;};
	void	setAtomC(Atom_sp a) {this->_C=a;};
	void	setAtomD(Atom_sp a) {this->_D=a;};

	double	getDegrees() {return this->_Degrees;};
	double	getWeight() {return this->_Weight;};

	Restraint_sp	copyDontRedirectAtoms();
	void redirectAtoms();

    	void invertStereochemistryOfRestraint();


	RestraintDihedral_O(const RestraintDihedral_O& old);
	DEFAULT_CTOR_DTOR(RestraintDihedral_O);
};



SMART(RestraintFixedNonbond);
class RestraintFixedNonbond_O : public Restraint_O
{
    LISP_CLASS(chem,ChemPkg,RestraintFixedNonbond_O,"RestraintFixedNonbond",Restraint_O);

public:
	void	initialize();
protected:
	Matter_sp	_Matter;
public:
    static RestraintFixedNonbond_sp create(Matter_sp matter);
public:
	void setMatter(Matter_sp matter) { this->_Matter = matter;};
	Matter_sp getMatter() { return this->_Matter;};

	RestraintFixedNonbond_O(const RestraintFixedNonbond_O& old);
	DEFAULT_CTOR_DTOR(RestraintFixedNonbond_O);
};






};
#endif
