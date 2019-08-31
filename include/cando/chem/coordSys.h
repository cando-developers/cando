/*
    File: coordSys.h
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


//
// coordSys.h
//


#ifndef	COORDINATESYSTEM_H
#define	COORDINATESYSTEM_H

#include <iostream>
#include <string>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/geom/matrix.h>
//#include "matter.h"
//#include "atom.h"

#include <cando/chem/chemPackage.h>

//#include "geom/render.fwd.h"// coordSys.h wants Render needs render.fwd.h

namespace chem {

SMART(Render);
SMART(Atom);


class	CoordinateSystem_O;
    typedef	gctools::smart_ptr<CoordinateSystem_O>	CoordinateSystem_sp;

FORWARD(CoordinateSystem);
class CoordinateSystem_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,CoordinateSystem_O,"CoordinateSystem",core::CxxObject_O);

public:
	void initialize();
public:
CL_LISPIFY_NAME("setOrigin");
CL_DEFMETHOD 	void	setOrigin(const Vector3& o)	{ this->origin=o; };
CL_LISPIFY_NAME("getOrigin");
CL_DEFMETHOD 	Vector3		getOrigin()	{ return (this->origin); };
CL_LISPIFY_NAME("getX");
CL_DEFMETHOD 	Vector3		getX()		{ return (this->x); };
CL_LISPIFY_NAME("getY");
CL_DEFMETHOD 	Vector3		getY()		{ return (this->y); };
CL_LISPIFY_NAME("getZ");
CL_DEFMETHOD 	Vector3		getZ()		{ return (this->z); };

	void	canonical();

	void	defineForAtoms( Atom_sp	origin,
				Atom_sp 	x,
				Atom_sp 	xy );

	void	defineForAtomVectors( const Vector3& vo,
				      const Vector3& vx,
				      const Vector3& vxy);

	void	defineForOneAtomVector( const Vector3& vo );

	void	defineForTwoAtomVectors( const Vector3& vo,
					 const Vector3& vx );

	void	defineForVectorsOriginXDirYDir( const Vector3& origin,
						const Vector3& xDir,
						const Vector3& yDir );

	void	defineForVectorsOriginXDirXYPlane( const Vector3& origin,
						   const Vector3& xDir,
						   const Vector3& xyPlaneDir );

	void	defineForVectorsOriginZDirXZPlane( const Vector3& origin,
						   const Vector3& zDir,
						   const Vector3& xzPlaneDir );

	void	defineForVectorsOriginXDirZDir( const Vector3& origin,
						const Vector3& xDir,
						const Vector3& zDir );

	void	defineFromCoordinateSystem(CoordinateSystem_sp sys);
	void	defineCanonical();

	Vector3	vectorRelativeToYou(const Vector3& pos );

	void	defineRandom();

	Matrix	matrixFromCanonical();
	Matrix	matrixToCanonical();

	Matrix	matrixForTransformTo(CoordinateSystem_sp cs );


	void	dump();
	bool		sameAs(CoordinateSystem_sp coord);

	CoordinateSystem_sp	copyCoordinateSystem();

	void	transformWithMatrixPointer( Matrix*	mp );
	void	transformWithMatrix( const Matrix&	m );

	string asString();

        virtual bool canRender() { return true; };
#ifdef RENDER
        virtual geom::Render_sp rendered(core::List_sp options);
#endif

	CoordinateSystem_O(const CoordinateSystem_O& orig);

private:
	Vector3		origin;
	Vector3		x;
	Vector3		y;
	Vector3		z;


	DEFAULT_CTOR_DTOR(CoordinateSystem_O);
};


};

TRANSLATE(chem::CoordinateSystem_O);
#endif
