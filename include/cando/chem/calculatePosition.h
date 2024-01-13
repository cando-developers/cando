/*
    File: calculatePosition.h
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
#ifndef	CalculatePosition_H //[
#define CalculatePosition_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>


#include <cando/chem/chemPackage.h>


namespace chem {


    SMART(CalculatePosition);
    class CalculatePosition_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,CalculatePosition_O,"CalculatePosition",core::CxxObject_O);
public: // virtual functions inherited from Object
	void	initialize();
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
	LISP_CLASS(chem,ChemPkg,CalculatePositionUsingInternals_O,"CalculatePositionUsingInternals",CalculatePosition_O);
#if INIT_TO_FACTORIES
    public:
	static CalculatePositionUsingInternals_sp make(MatterName distanceAtomName,
						       MatterName angleAtomName,
						       MatterName dihedralAtomName,
						       double distance, double angleDegrees, double dihedralDegrees);
#else
#endif
public: // virtual functions inherited from Object
	void	initialize();
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
	LISP_CLASS(chem,ChemPkg,CalculatePositionAlongBond_O,"CalculatePositionAlongBond",CalculatePosition_O);
#if INIT_TO_FACTORIES
    public:
	static CalculatePositionAlongBond_sp make(MatterName fromAtomName, MatterName toAtomName, double extendDistance);
#else
#endif
public: // virtual functions inherited from Object
	void	initialize();

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
	LISP_CLASS(chem,ChemPkg,CalculatePositionRelativeToOrigin_O,"CalculatePositionRelativeToOrigin",CalculatePosition_O);
#if INIT_TO_FACTORIES
    public:
	static CalculatePositionRelativeToOrigin_sp make(MatterName originAtomName, MatterName toAtomName, double extendDistance);
#else
#endif
public: // virtual functions inherited from Object
	void	initialize();

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
#endif //]
