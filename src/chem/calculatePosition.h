#ifndef	CalculatePosition_H //[
#define CalculatePosition_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/candoBase/vector3.h>


#include <cando/chem/chemPackage.h>


namespace chem {


    SMART(CalculatePosition);
    class CalculatePosition_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,CalculatePosition_O,"CalculatePosition");
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables

public:	// Creation class functions

public:
	virtual Vector3 calculatePosition(core::T_sp positionLookup) {_OF();SUBCLASS_MUST_IMPLEMENT();};

	virtual MatterName closestAtomName() {_OF();SUBCLASS_MUST_IMPLEMENT();};

	CalculatePosition_O( const CalculatePosition_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CalculatePosition_O);
    };



    SMART(CalculatePositionUsingInternals );
    class CalculatePositionUsingInternals_O : public CalculatePosition_O
    {
	LISP_BASE1(CalculatePosition_O);
	LISP_CLASS(chem,ChemPkg,CalculatePositionUsingInternals_O,"CalculatePositionUsingInternals");
#if INIT_TO_FACTORIES
    public:
	static CalculatePositionUsingInternals_sp make(MatterName distanceAtomName,
						       MatterName angleAtomName,
						       MatterName dihedralAtomName,
						       double distance, double angleDegrees, double dihedralDegrees);
#else
	DECLARE_INIT();
#endif
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
	MatterName	_DistanceAtomName;
	MatterName	_AngleAtomName;
	MatterName	_DihedralAtomName;
	double	_Distance;
	double	_AngleDegrees;
	double	_DihedralDegrees;
public:	// Creation class functions

public:
	virtual Vector3 calculatePosition(core::T_sp positionLookup);

	virtual MatterName closestAtomName() {return this->_DistanceAtomName;};

	CalculatePositionUsingInternals_O( const CalculatePosition_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CalculatePositionUsingInternals_O);
    };


    SMART(CalculatePositionAlongBond );
    class CalculatePositionAlongBond_O : public CalculatePosition_O
    {
	LISP_BASE1(CalculatePosition_O);
	LISP_CLASS(chem,ChemPkg,CalculatePositionAlongBond_O,"CalculatePositionAlongBond");
#if INIT_TO_FACTORIES
    public:
	static CalculatePositionAlongBond_sp make(MatterName fromAtomName, MatterName toAtomName, double extendDistance);
#else
	DECLARE_INIT();
#endif
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);

private: // instance variables
	MatterName	_FromAtomName;
	MatterName	_ToAtomName;
	double	_ExtendDistance;
public:	// Creation class functions

public:
	virtual Vector3 calculatePosition(core::T_sp positionLookup);

	virtual MatterName closestAtomName() {return this->_ToAtomName;};

	CalculatePositionAlongBond_O( const CalculatePosition_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CalculatePositionAlongBond_O);
    };


    SMART(CalculatePositionRelativeToOrigin );
    class CalculatePositionRelativeToOrigin_O : public CalculatePosition_O
    {
	LISP_BASE1(CalculatePosition_O);
	LISP_CLASS(chem,ChemPkg,CalculatePositionRelativeToOrigin_O,"CalculatePositionRelativeToOrigin");
#if INIT_TO_FACTORIES
    public:
	static CalculatePositionRelativeToOrigin_sp make(MatterName originAtomName, MatterName toAtomName, double extendDistance);
#else
	DECLARE_INIT();
#endif
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);

private: // instance variables
	MatterName	_OriginAtomName;
	MatterName	_ToAtomName;
	double	_ExtendDistance;
public:	// Creation class functions

public:
	virtual Vector3 calculatePosition(core::T_sp positionLookup);

	virtual MatterName closestAtomName() {return this->_ToAtomName;};

	CalculatePositionRelativeToOrigin_O( const CalculatePosition_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CalculatePositionRelativeToOrigin_O);
    };





};
TRANSLATE(chem::CalculatePosition_O);
TRANSLATE(chem::CalculatePositionUsingInternals_O);
TRANSLATE(chem::CalculatePositionAlongBond_O);
TRANSLATE(chem::CalculatePositionRelativeToOrigin_O);
#endif //]
