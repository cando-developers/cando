#ifndef	_kinematics_Rotamer_H
#define _kinematics_Rotamer_H

#include <clasp/core/common.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <clasp/core/holder.h>

namespace kinematics
{

    FORWARD(Rotamer);
    class Rotamer_O : public core::T_O
    {
	LISP_CLASS(kinematics,KinPkg,Rotamer_O,"Rotamer",core::T_O);
#if INIT_TO_FACTORIES
    public:
	static Rotamer_sp make(core::List_sp dihedrals, core::List_sp sigmas, core::List_sp indices, const double probability, const int count);
#else
	DECLARE_INIT();
#endif
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(Rotamer_O);
    public:
	void initialize();
    public:
    private: // instance variables here
	vector<double>	_Dihedrals;
	vector<double>	_Sigmas;
	vector<int>	_Indices;
	double		_Probability;
	int		_Count;
    public: // Functions here
	string __repr__() const;
	double probability() const { return this->_Probability;};
	
    };


    FORWARD(RotamerSetBase);
    class RotamerSetBase_O : public core::T_O
    {
	LISP_CLASS(kinematics,KinPkg,RotamerSetBase_O,"RotamerSetBase",core::T_O);
#if INIT_TO_FACTORIES
    public:
	static RotamerSetBase_sp make();
#else
	DECLARE_INIT();
#endif
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(RotamerSetBase_O);
    public:
	void initialize();
    public:
    private: // instance variables here
	gctools::Vec0<Rotamer_sp>	_Rotamers;
    public: // Functions here
	Rotamer_sp addRotamer(Rotamer_sp rotamer);
	core::List_sp asList() const;
	int size() const { return this->_Rotamers.size();};
	Rotamer_sp get(int i) const;
    };




    FORWARD(RotamerSet);
    class RotamerSet_O : public RotamerSetBase_O
    {
	LISP_CLASS(kinematics,KinPkg,RotamerSet_O,"RotamerSet",RotamerSetBase_O);
#if INIT_TO_FACTORIES
    public:
	static RotamerSet_sp make();
#else
	DECLARE_INIT();
#endif
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(RotamerSet_O);
    public:
	void initialize();
    public:
    };




    FORWARD(BackboneDependentRotamerSet);
    class BackboneDependentRotamerSet_O : public RotamerSetBase_O
    {
	friend class BackboneDependentRotamerLibrary_O;
	LISP_CLASS(kinematics,KinPkg,BackboneDependentRotamerSet_O,"BackboneDependentRotamerSet",RotamerSetBase_O);
#if INIT_TO_FACTORIES
    public:
	static BackboneDependentRotamerSet_sp make(const int phi, const int psi);
#else
	DECLARE_INIT();
#endif
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(BackboneDependentRotamerSet_O);
    public:
	void initialize();
    public:
    private: // instance variables here
	int			_Phi;
	int			_Psi;
    public: // Functions here
	int phi() const { return this->_Phi;};
	int psi() const { return this->_Psi;};
    };












    FORWARD(BackboneDependentRotamerLibrary);
    class BackboneDependentRotamerLibrary_O : public core::T_O
    {
	LISP_CLASS(kinematics,KinPkg,BackboneDependentRotamerLibrary_O,"BackboneDependentRotamerLibrary",core::T_O);
#if INIT_TO_FACTORIES
    public:
	static BackboneDependentRotamerLibrary_sp make(const int phiStep, const int phiStart, const int phiCount, const int psiStep, const int psiStart, const int psiCount);
#else
	DECLARE_INIT();
#endif
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(BackboneDependentRotamerLibrary_O);
    public:
	void initialize();
    public:
	int	_PhiStep;
	int	_PhiStart;
	int	_PhiCount;
	int	_PsiStep;
	int	_PsiStart;
	int	_PsiCount;
	gctools::Vec0<BackboneDependentRotamerSet_sp>	_RotamerSets;

	bool validPhiPsi(int phi, int psi) const;
	int _index(int phi, int psi, int& mphi, int& mpsi ) const;

	core::List_sp rotamerSetsAsList() const;
	void addRotamerSet(BackboneDependentRotamerSet_sp rotamerSet);
    };

}; /* kinematics */

TRANSLATE(kinematics::Rotamer_O);
TRANSLATE(kinematics::RotamerSetBase_O);
TRANSLATE(kinematics::RotamerSet_O);
TRANSLATE(kinematics::BackboneDependentRotamerSet_O);
TRANSLATE(kinematics::BackboneDependentRotamerLibrary_O);

#endif /* _kinematics_Rotamer_H */


