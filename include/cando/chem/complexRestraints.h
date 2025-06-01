/*
    File: complexRestraints.h
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
#ifndef	ComplexRestraints_H //[
#define ComplexRestraints_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>


namespace chem {



SMART(Residue);
SMART(Atom);
SMART(ChemInfo);

SMART(ComplexRestraint );
class ComplexRestraint_O : public core::CxxObject_O
{
  CL_DOCSTRING(R"(Subclasses represent restraints for double bonds and axial/equatorial for exocyclic atoms)");
  LISP_CLASS(chem,ChemPkg,ComplexRestraint_O,"ComplexRestraint",core::CxxObject_O);
#if INIT_TO_FACTORIES
 public:
    static ComplexRestraint_sp make();
#else
    DECLARE_INIT();
#endif
public:
public: // virtual functions inherited from Object
	void	initialize();
//	string	__repr__() const;

private: // instance variables
public:	// Creation class functions
public:

  CL_DEFMETHOD virtual size_t fillRestraints(Residue_sp residue, core::HashTable_sp cip) {_OF();SUBCLASS_MUST_IMPLEMENT();};

public:
	ComplexRestraint_O( const ComplexRestraint_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(ComplexRestraint_O);
};









SMART(RestrainedExoCyclicAtom );
class RestrainedExoCyclicAtom_O : public ComplexRestraint_O
{
    LISP_CLASS(chem,ChemPkg,RestrainedExoCyclicAtom_O,"RestrainedExoCyclicAtom",ComplexRestraint_O);
public:
    static RestrainedExoCyclicAtom_sp make();
public:
public: // virtual functions inherited from Object
    void	initialize();
//	string	__repr__() const;

private: // instance variables
    MatterName	_ExoCyclicAtomName;
    //! Its either axial or if not equatorial
    bool	_Axial;


private:
    static bool	_LazyInitializedSmarts;
//    static ChemInfo_sp	_AtomExoToSixMemberedRing;
private:
    void lazyInitializeSmarts();

public:	// Creation class functions
public:
    /*! String is case-insensitive
     * and can be "axial" or "equatorial"
     * or "A" or "E"
     */
    void setConfiguration(const string& s);
    string getConfiguration();

    void setExoCyclicAtomName(const string& name);
    string getExoCyclicAtomName();

    size_t fillRestraints(Residue_sp residue, core::HashTable_sp cip);

public:
    RestrainedExoCyclicAtom_O( const RestrainedExoCyclicAtom_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(RestrainedExoCyclicAtom_O);
};










SMART(RestrainedPiBond);
class RestrainedPiBond_O : public ComplexRestraint_O
{
    LISP_CLASS(chem,ChemPkg,RestrainedPiBond_O,"RestrainedPiBond",ComplexRestraint_O);
#if INIT_TO_FACTORIES
 public:
    static RestrainedPiBond_sp make(core::Symbol_sp config, core::List_sp piAtoms);
#else
    DECLARE_INIT();
#endif
public:
public:
    core::Symbol_sp	_Configuration;	// "E" or "Z"
    MatterName  	_PiAtomX;
    MatterName  	_PiAtomY;

public:
CL_LISPIFY_NAME("getConfiguration");
CL_DEFMETHOD     core::Symbol_sp	getConfiguration()	{return this->_Configuration;};
CL_LISPIFY_NAME("setConfiguration");
CL_DEFMETHOD     void	setConfiguration(core::Symbol_sp s)	{this->_Configuration = s;};

CL_LISPIFY_NAME("getPiAtomX");
CL_DEFMETHOD     MatterName	getPiAtomX() { return this->_PiAtomX; };
CL_LISPIFY_NAME("setPiAtomX");
CL_DEFMETHOD     void	setPiAtomX(MatterName n) { this->_PiAtomX = n; };
CL_LISPIFY_NAME("getPiAtomY");
CL_DEFMETHOD     MatterName	getPiAtomY() { return this->_PiAtomY; };
CL_LISPIFY_NAME("setPiAtomY");
CL_DEFMETHOD     void	setPiAtomY(MatterName n) { this->_PiAtomY = n; };

    void fillOneDihedralRestraint(Residue_sp residue,
				Atom_sp p1X, Atom_sp pX, Atom_sp pY, Atom_sp p1Y,
				double deg, double weight);
    virtual size_t fillRestraints(Residue_sp residue, core::HashTable_sp cip);

    DEFAULT_CTOR_DTOR(RestrainedPiBond_O);
};








};
#endif //]
