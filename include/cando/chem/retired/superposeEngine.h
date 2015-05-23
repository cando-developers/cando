       
       

//
// (C) 2004 Christian E. Schafmeister
//

#ifndef	SuperposeEngine_H
#define	SuperposeEngine_H


#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/geom/matrix.h>
#include <cando/geom/coordinateArray.h>


#include <cando/chem/chemPackage.h>

#include "core/intArray.fwd.h"// superposeEngine.h wants IntArray needs intArray.fwd.h
#include "geom/coordinateArray.fwd.h"// superposeEngine.h wants CoordinateArray needs coordinateArray.fwd.h

namespace chem
{





SMART(IntArray);
SMART(CoordinateArray);
SMART(Matter);
SMART(Atom);
typedef	enum	{ init, gotFixed, gotMoveable } SuperposeEngineStates;

class	SuperposeEngine_O;
    typedef	gctools::smart_ptr<SuperposeEngine_O>	SuperposeEngine_sp;


/*! Superpose a moveable set of points onto a fixed set of points.
 * You can specify a larger list of moveable points and then index into them
 * using setMoveableIndex.
 * This will superpose the points in the following way.
 * 	- FixedIndices is an array of indeces into FixedCoordinates
 * 	- MoveableIndices is an array of indices into MoveableCoordinates
 *	Superpose Fixed points onto Moveable points.
 *	The points that are superposed are only the ones indexed by _FixedIndices and _MoveableIndices.
 */
SMART(SuperposeEngine);
class SuperposeEngine_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,SuperposeEngine_O,"SuperposeEngine");

public:
	static void lisp_initGlobals(core::Lisp_sp lisp);

	void	initialize();
//	void	archiveBase(core::ArchiveP node);
private:
	core::IntArray_sp		_FixedIndices;
	geom::CoordinateArray_sp	_FixedCoordinates;
	core::IntArray_sp		_MoveableIndices;
	geom::CoordinateArray_sp	_MoveableCoordinates;
private://do not archive
	Matrix			_Transform;
public:
	string	debugString();
	void	setFixedPoints(core::IntArray_sp indices, geom::CoordinateArray_sp fp );
	void	setMoveablePoints(core::IntArray_sp indices, geom::CoordinateArray_sp mp);

	void	setFixedAllPoints(geom::CoordinateArray_sp fp );
	void	setMoveableAllPoints(geom::CoordinateArray_sp mp);

	void	eraseMoveablePoints();
	void	appendMoveablePoint(const Vector3& pos);
	void	eraseFixedPoints();
	void	appendFixedPoint(const Vector3& pos);

//	geom::CoordinateArray_sp	getFixedCoordinates();
	geom::CoordinateArray_sp	getMoveableCoordinates();

//	void	setMoveableIndex(const vector<int> vi);
	int	getNumberOfFixedPoints();
	int	getNumberOfMoveablePoints();

//	geom::CoordinateArray_sp	getMoveablePoints() {_OF(); ASSERTNOTNULL(this->_Moveable);return this->_Moveable; };
	Matrix	superpose();

//	double	sumOfSquaresOfDifferences(ScorerState_sp scorerState);
	double	sumOfSquaresOfDifferences();
	double	rootMeanSquareDifference();




	void	doSuperpose();

	DEFAULT_CTOR_DTOR(SuperposeEngine_O);
};

inline	geom::CoordinateArray_sp	SuperposeEngine_O::getMoveableCoordinates()
{_OF();
    ASSERTNOTNULL(this->_MoveableCoordinates);
    return this->_MoveableCoordinates;
}




SMART(SuperposeSelectedAtoms );
class SuperposeSelectedAtoms_O : public SuperposeEngine_O
{
    LISP_BASE1(SuperposeEngine_O);
    LISP_CLASS(chem,ChemPkg,SuperposeSelectedAtoms_O,"SuperposeSelectedAtoms");

public:
public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
private:
	Matter_sp		_Matter;
    gctools::Vec0<Atom_sp>		_SuperposeAtoms;
private:
	geom::CoordinateArray_sp extractCoordinates(Matter_sp matter);
public:
		/*! Set the Matter to superpose
		 * Initially all _SuperposeAtoms will be empty and
		 * if you try to superpose anything an exception will be thrown
		 */
	void	setMatter(Matter_sp matter);

		/*! Get the Matter
		 */
	Matter_sp getMatter();

		/*!
		 * Regenerate _SuperposeAtoms using the selected atoms
		 * of our _Matter.
		 */
	void updateSuperposeAtoms();

		/*!
		 * Create a copy of the Matter and extract the
		 * selected atoms into _HeavyAtoms
		 */
	void	setMatterWithSelectedAtoms(Matter_sp matter);

		/*!
		 * Copy the coordinates from (matter) into our _Matter
		 * and transfer those selected atoms into the FixedCoordinates
		 * of the superposer.
		 */
	void	copyMatterCoordinatesIntoFixedCoordinates(Matter_sp matter);
		/*!
		 * Copy the coordinates from (matter) into our _Matter
		 * and transfer those selected atoms into the MoveableCoordinates
		 * of the superposer.
		 */
	void	copyMatterCoordinatesIntoMoveableCoordinates(Matter_sp matter);



	DEFAULT_CTOR_DTOR(SuperposeSelectedAtoms_O);
};

};

TRANSLATE(chem::SuperposeEngine_O);
TRANSLATE(chem::SuperposeSelectedAtoms_O);
#endif //]
