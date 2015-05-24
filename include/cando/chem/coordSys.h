       
       
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
class CoordinateSystem_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,CoordinateSystem_O,"CoordinateSystem");

public:
	void initialize();
public:
//	void	archive(core::ArchiveP node);
public:
public:
	void	setOrigin(const Vector3& o)	{ this->origin=o; };
	Vector3		getOrigin()	{ return (this->origin); };
	Vector3		getX()		{ return (this->x); };
	Vector3		getY()		{ return (this->y); };
	Vector3		getZ()		{ return (this->z); };

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

	CoordinateSystem_sp	copy();

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
