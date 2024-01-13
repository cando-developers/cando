/*
    File: largeSquareMatrix.h
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
 *      File:   largeMatrix.h
 *
 *      Description:
 *              Routines to maintain NxN dimensional Matrices
 */

#ifndef LARGE_MATRIX_H
#define LARGE_MATRIX_H

#include <clasp/core/common.h>
#include <cando/chem/nVector.h>

#include <cando/chem/chemPackage.h>


namespace chem {

typedef	enum	{ UpperDiagonalLower, UpperDiagonal, DiagonalLower, SymmetricUpperDiagonal, SymmetricDiagonalLower } TriangleType;
typedef	enum	{ AbstractMatrix, FullMatrix, SparseMatrix } MatrixClass;

inline	string	XmlAttribute_TriangleType(TriangleType t) {
			if ( t==UpperDiagonalLower ) {
			    return "UpperDiagonalLower";
			} else if ( t==SymmetricDiagonalLower ) {
			    return "SymmetricDiagonalLower";
			} else if ( t==SymmetricUpperDiagonal ) {
			    return "SymmetricUpperDiagonal";
			} else if ( t==DiagonalLower ) {
			    return "DiagonalLower";
			} else if ( t==UpperDiagonal ) {
			    return "UpperDiagonal";
			}
			return "-undefinedTriangleType-";
		};




SMART(AbstractLargeSquareMatrix);
class AbstractLargeSquareMatrix_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,AbstractLargeSquareMatrix_O,"AbstractLargeSquareMatrix",core::CxxObject_O);

public:
  TriangleType 	_Triangle;
  uint   		_Columns;
  uint		_Rows;

protected:
  virtual	uint	indexBegin() {SUBIMP();};
  virtual	uint	indexEnd() {SUBIMP();}
  virtual	uint	indexAdvance(uint ii) {SUBIMP();};
  virtual	void	coordinatesFromIndex(uint ii,uint& row, uint& col ) {SUBIMP();};

public:
  AbstractLargeSquareMatrix_O( const AbstractLargeSquareMatrix_O& orig );


  void	setup(uint dim, TriangleType type);

//	void	describe();

  CL_DEFMETHOD uint columns()	{ return this->_Columns; };
  CL_DEFMETHOD uint rows()	{ return this->_Rows; };
  uint	dimension()	{return this->_Rows;};

  CL_DEFMETHOD virtual	vecreal	element(uint x, uint y) {SUBIMP();};
  virtual	void	setElement(uint x, uint y, vecreal d) {SUBIMP();};
  virtual	void	addToElement(uint x, uint y, vecreal d) {SUBIMP();};
  virtual	void	setAtIndex(uint ii, vecreal d) {SUBIMP();};
  virtual	vecreal	getAtIndex(uint ii) {SUBIMP();};
  virtual	void	insertionIsComplete() {};


  CL_DEFMETHOD virtual vecreal maxAbsValue() {SUBIMP(); };
  CL_DEFMETHOD virtual vecreal mostPositiveValue() {SUBIMP(); };
  CL_DEFMETHOD virtual vecreal mostNegativeValue() {SUBIMP(); };
  virtual void	fill(vecreal d) {SUBIMP();};
  void	writeMathematica(const string& fileName);
  void	writeMathematicaSymbolic(const string& fileName);
  void	dumpMatrix(uint x);


  virtual	void	insertElement(uint x, uint y) {};
  virtual	void	reset() {};
  virtual	void	zero() {this->fill(0.0);};

  virtual	MatrixClass	myClass()		{return AbstractMatrix;};
  virtual	AbstractLargeSquareMatrix_sp copy_matrix( ) {SUBIMP();};
  virtual	uint	indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y) {SUBIMP();};
  virtual	uint	col_indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y) {SUBIMP();};
  virtual	uint	col_indexFromCoordinates(uint x, uint y);
  CL_DEFMETHOD virtual	bool	hasElement(uint x, uint y) { return this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y)!=UndefinedUnsignedInt; };
  virtual	void	debug()	{};

  virtual	void	multiplyByVector(NVector_sp result, NVector_sp x);

  uint	indexFromCoordinates(uint x, uint y);

  DEFAULT_CTOR_DTOR(AbstractLargeSquareMatrix_O);
};












SMART(FullLargeSquareMatrix)
class FullLargeSquareMatrix_O : public  AbstractLargeSquareMatrix_O 
{
  LISP_CLASS(chem,ChemPkg,FullLargeSquareMatrix_O,"FullLargeSquareMatrix", AbstractLargeSquareMatrix_O );

private:
  uint		_ActiveElements;
  gctools::Vec0<vecreal>	_Values;

private:
  void	setup(uint dim, TriangleType type);

public:
  static FullLargeSquareMatrix_sp	create(uint dim, TriangleType type);

  virtual	uint	indexBegin()	{ return 0;};
  virtual	uint	indexEnd() {return this->_ActiveElements;};
  virtual	uint	indexAdvance(uint dp) {return dp+1;};
  virtual	void	coordinatesFromIndex(uint i,uint& x, uint& y );

public:

  FullLargeSquareMatrix_O( const FullLargeSquareMatrix_O& orig );

  virtual vecreal maxAbsValue();
  virtual void	fill(vecreal d);
  virtual	MatrixClass	myClass()	{return FullMatrix;};
  virtual	AbstractLargeSquareMatrix_sp copy_matrix();
  virtual	uint	indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y);
  virtual	void	insertElement(uint x, uint y) {};
  virtual	void	reset() {};
  virtual	void	setAtIndex(uint ii, vecreal d) { this->_Values[ii] = d;};
  virtual	vecreal	getAtIndex(uint ii) { return this->_Values[ii];};
  virtual	void	setElement(uint x, uint y, vecreal d) { this->_Values[this->indexFromCoordinates(x,y)]=d;};
  virtual	void	addToElement(uint x, uint y, vecreal d) { this->_Values[this->indexFromCoordinates(x,y)]+=d;};
  virtual	vecreal	element(uint x, uint y) { return this->_Values[this->indexFromCoordinates(x,y)];};
  virtual	void	debug()	{};
  DEFAULT_CTOR_DTOR(FullLargeSquareMatrix_O);
};


#if 0
inline	FullLargeSquareMatrix_sp	new_FullLargeSquareMatrix_sp(uint dim, TriangleType m) {
    return FullLargeSquareMatrix_sp(new FullLargeSquareMatrix_O(dim,m));
}

inline	FullLargeSquareMatrix_sp new_FullLargeSquareMatrix_sp(const FullLargeSquareMatrix_O& m) {
    return FullLargeSquareMatrix_sp(new FullLargeSquareMatrix_O(m));
}
#endif





SMART(SparseLargeSquareMatrix);
class SparseLargeSquareMatrix_O : public AbstractLargeSquareMatrix_O
{
  LISP_CLASS(chem,ChemPkg,SparseLargeSquareMatrix_O,"SparseLargeSquareMatrix",AbstractLargeSquareMatrix_O);

public:
  bool                          _InsertionIsComplete;
  uint		                _RowStartEntries;
  uint		                _ActiveElements;
  uint		                _ReservedElements;
  gctools::Vec0<uint> 	        _RowStarts;
  gctools::Vec0<uint> 	        _ColumnForValue;
  gctools::Vec0<vecreal> 	_Values;
  bool                          _col_OptimizationDone;
  uint		                _col_ColumnStartEntries;
  uint		                _col_ActiveElements;
  uint		                _col_ReservedElements;
  gctools::Vec0<uint> 	        _col_ColumnStarts;
  gctools::Vec0<uint> 	        _col_RowForValue;
  gctools::Vec0<vecreal> 	_col_Values;
public:
  static SparseLargeSquareMatrix_sp	create(uint dim, TriangleType type);
private:
  void	initializeStorage();
  void	expandStorage();
  void	col_expandStorage();
  void	releaseStorage();
  void	setup(uint dim, TriangleType type);

public:

  void ensureIdentical(SparseLargeSquareMatrix_sp other);

  virtual       vecreal maxAbsValue();
  virtual       vecreal mostPositiveValue();
  virtual       vecreal mostNegativeValue();
  virtual       void	fill(vecreal d);

  void	walkMatrix(core::T_sp callback);
  
  virtual	uint	indexBegin()	{ return 0; };
  virtual	uint	indexEnd() {return this->_ActiveElements;};
  virtual	uint	indexAdvance(uint i) {return i+1;};
  virtual	void	coordinatesFromIndex(uint i,uint& x, uint& y );
  virtual	void	addToElement(uint x, uint y, vecreal d);

  virtual	void	setElement(uint x, uint y, vecreal d);
  virtual	void	col_setElement(uint x, uint y, vecreal d);
  virtual	vecreal element(uint x, uint y);
  virtual	vecreal col_element(uint x, uint y);


  virtual	MatrixClass	myClass()	{return SparseMatrix;};
  virtual	AbstractLargeSquareMatrix_sp copy_matrix();
  virtual	uint	indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y);
  virtual	uint	col_indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y);
//virtual	bool	hasElement(uint x, uint y) { return this->indexFromCoordinates(x,y)>=0; };
  virtual	void	insertElement(uint x, uint y);
  virtual	void	col_insertElement(uint x, uint y);
  virtual	void	reset();
  virtual	void	insertionIsComplete();
  virtual       void    doColumnOptimization();
  virtual	void	debug();
  virtual	SparseLargeSquareMatrix_sp optimized();
  virtual	void	setAtIndex(uint ii, vecreal d) { this->_Values[ii] = d;};
  virtual	vecreal getAtIndex(uint ii) { return this->_Values[ii];};

  template <typename Func>
  void loopOverRowLowerDiagonal(uint row, uint start_col, uint end_col, Func f) {
    uint ib = this->_RowStarts[row];
    uint ie = this->_RowStarts[row+1]-1;
    for ( uint ii = ib; ii<=ie; ii++ ) {
      uint col = this->_ColumnForValue[ii];
      if (start_col<=col && col < end_col) {
        f(col,row,this->_Values[ii]);
      }
    }
  }

#if 1
  template <typename Func>
  void loopOverColumnLowerDiagonal(uint col, uint start_row, uint end_row, Func f) {
    if (!this->_col_OptimizationDone) this->doColumnOptimization();
    uint ib = this->_col_ColumnStarts[col];
    uint ie = this->_col_ColumnStarts[col+1]-1;
    for ( uint ii = ib; ii<=ie; ii++ ) {
      uint row = this->_col_RowForValue[ii];
      if ( start_row <= row && row < end_row) {
        f(col,row,this->_Values[ii]);
      }
    }
  }
#else
  template <typename Func>
  void loopOverColumnLowerDiagonal(uint column, uint row_start, uint row_end, Func f) {
    for ( int row = row_start; row < row_end; row++ ) {
      uint ib = this->_RowStarts[row];
      uint ie = this->_RowStarts[row+1];
      auto it = std::lower_bound(this->_ColumnForValue.begin()+ib,
                                 this->_ColumnForValue.begin()+ie,
                                 column );
      if (it != this->_ColumnForValue.begin()+ie && *it==column) {
        int idx = std::distance(this->_ColumnForValue.begin(),it);
        uint col = this->_ColumnForValue[idx];
        vecreal val = this->_Values[idx];
        f(col,row,val);
      }
    }
  }
#endif
  
  uint	columnForIndex(uint ii) { return this->_ColumnForValue[ii];};

		/*! Return the index of the element in row(y)
		 * that immediatly follows column (x)
		 */
  uint	indexOfFirstElementAtOrAfterX(uint x, uint y);
		/*! Return the index of the last element in row(y)
  		*/
  uint	indexOfLastElementOnRow(uint y );

public:
  SparseLargeSquareMatrix_O( const SparseLargeSquareMatrix_O& orig );
  DEFAULT_CTOR(SparseLargeSquareMatrix_O);
  virtual ~SparseLargeSquareMatrix_O(); // Non-trivial
};
};

template <>
struct gctools::GCInfo<chem::SparseLargeSquareMatrix_O> {
  static bool constexpr CanAllocateWithNoArguments = true;
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = true;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {
#if 0
inline	SparseLargeSquareMatrix_sp	new_SparseLargeSquareMatrix_sp(uint dim, TriangleType m) {
    return SparseLargeSquareMatrix_sp(new SparseLargeSquareMatrix_O(dim,m));
}

inline	SparseLargeSquareMatrix_sp new_SparseLargeSquareMatrix_sp(const SparseLargeSquareMatrix_O& m) {
    return SparseLargeSquareMatrix_sp(new SparseLargeSquareMatrix_O(m));
}
#endif






};


#endif          /* ifndef LARGE_MATRIX_H */
