/*
    File: superposeEngine.h
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

#ifndef	SuperposeEngine_H
#define	SuperposeEngine_H


#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/geom/matrix.h>
#include <cando/geom/coordinateArray.h>


#include <cando/chem/chemPackage.h>

#include <clasp/core/intArray.fwd.h>
#include <cando/geom/coordinateArray.fwd.h>// superposeEngine.h wants SimpleVectorCoordinate needs coordinateArray.fwd.h

namespace chem
{
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
class SuperposeEngine_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,SuperposeEngine_O,"SuperposeEngine",core::CxxObject_O);

public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
  
	void	initialize();
private:
	core::ComplexVector_byte32_t_sp		_FixedIndices;
	geom::ComplexVectorCoordinate_sp	_FixedCoordinates;
	core::ComplexVector_byte32_t_sp		_MoveableIndices;
	geom::ComplexVectorCoordinate_sp	_MoveableCoordinates;
private://do not archive
	Matrix			_Transform;
public:
	string	debugString();
  void	setFixedPoints(core::ComplexVector_byte32_t_sp indices, geom::SimpleVectorCoordinate_sp fp );
	void	setMoveablePoints(core::ComplexVector_byte32_t_sp indices, geom::SimpleVectorCoordinate_sp mp);

	void	setFixedAllPoints(geom::SimpleVectorCoordinate_sp fp );
	void	setMoveableAllPoints(geom::SimpleVectorCoordinate_sp mp);

	void	eraseMoveablePoints();
	void	appendMoveablePoint(const Vector3& pos);
	void	eraseFixedPoints();
	void	appendFixedPoint(const Vector3& pos);

//	geom::SimpleVectorCoordinate_sp	getFixedCoordinates();
  inline core::Array_sp	getMoveableCoordinates() {
    return this->_MoveableCoordinates;
  }


//	void	setMoveableIndex(const vector<int> vi);
	int	getNumberOfFixedPoints();
	int	getNumberOfMoveablePoints();

//	geom::SimpleVectorCoordinate_sp	getMoveablePoints() {_OF(); ASSERTNOTNULL(this->_Moveable);return this->_Moveable; };
	Matrix	superpose();

//	double	sumOfSquaresOfDifferences(ScorerState_sp scorerState);
	double	sumOfSquaresOfDifferences();
	double	rootMeanSquareDifference();




	void	doSuperpose();

	DEFAULT_CTOR_DTOR(SuperposeEngine_O);
};




SMART(SuperposeSelectedAtoms );
class SuperposeSelectedAtoms_O : public SuperposeEngine_O
{
    LISP_CLASS(chem,ChemPkg,SuperposeSelectedAtoms_O,"SuperposeSelectedAtoms",SuperposeEngine_O);

public:
public:
	void	initialize();
private:
	Matter_sp		_Matter;
    gctools::Vec0<Atom_sp>		_SuperposeAtoms;
private:
	geom::SimpleVectorCoordinate_sp extractCoordinates(Matter_sp matter);
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
