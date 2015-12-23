#ifndef	ComplexRestraints_H //[
#define ComplexRestraints_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>


namespace chem {



SMART(RestraintList);
SMART(Residue);
SMART(Atom);
SMART(ChemInfo);

SMART(ComplexRestraint );
class ComplexRestraint_O : public core::CxxObject_O
{
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
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
public:	// Creation class functions
public:

	virtual void fillRestraints(Residue_sp residue) {_OF();SUBCLASS_MUST_IMPLEMENT();};

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
//    void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
    MatterName	_ExoCyclicAtomName;
    //! Its either axial or if not equatorial
    bool	_Axial;


private:
    static bool	_LazyInitializedSmarts;
    static ChemInfo_sp	_AtomExoToSixMemberedRing;
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

    void fillRestraints(Residue_sp residue);

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
//    void	archive(core::ArchiveP node);

private:
    core::Symbol_sp	_Configuration;	// "E" or "Z"
    MatterName  	_PiAtomX;
    MatterName  	_PiAtomY;

public:
CL_NAME("getConfiguration");
CL_DEFMETHOD     core::Symbol_sp	getConfiguration()	{return this->_Configuration;};
CL_NAME("setConfiguration");
CL_DEFMETHOD     void	setConfiguration(core::Symbol_sp s)	{this->_Configuration = s;};

CL_NAME("getPiAtomX");
CL_DEFMETHOD     MatterName	getPiAtomX() { return this->_PiAtomX; };
CL_NAME("setPiAtomX");
CL_DEFMETHOD     void	setPiAtomX(MatterName n) { this->_PiAtomX = n; };
CL_NAME("getPiAtomY");
CL_DEFMETHOD     MatterName	getPiAtomY() { return this->_PiAtomY; };
CL_NAME("setPiAtomY");
CL_DEFMETHOD     void	setPiAtomY(MatterName n) { this->_PiAtomY = n; };

    void fillOneDihedralRestraint(Residue_sp residue,
				Atom_sp p1X, Atom_sp pX, Atom_sp pY, Atom_sp p1Y,
				double min, double max, double weight);
    virtual void fillRestraints(Residue_sp residue);

    DEFAULT_CTOR_DTOR(RestrainedPiBond_O);
};








};
TRANSLATE(chem::ComplexRestraint_O);
TRANSLATE(chem::RestrainedExoCyclicAtom_O);
TRANSLATE(chem::RestrainedPiBond_O);
#endif //]
