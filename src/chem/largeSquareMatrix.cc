/*
    File: largeSquareMatrix.cc
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
#define		DEBUG_LEVEL_NONE
//#define	MEISTER_MALLOC_DEBUG


//#define	DEBUG_VIA_PRINTF
#include <cando/chem/largeSquareMatrix.h>
#include <iostream>
#include <iomanip>
#include <clasp/core/lispStream.h>
#include <clasp/core/ql.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/wrappers.h>


namespace chem {


#define	SWAP(x,y,t) { t=x; x=y; y=t; };






/*
 *      AbstractLargeSquareMatrix_O
 *
 *	Author:	Christian Schafmeister (1991)
 *
 */


void	AbstractLargeSquareMatrix_O::setup(uint dim, TriangleType type)
{
  this->_Triangle = type;
  ASSERT(this->_Triangle == SymmetricDiagonalLower);
  this->_Columns = dim;
  this->_Rows = dim;
};

AbstractLargeSquareMatrix_O::AbstractLargeSquareMatrix_O(const AbstractLargeSquareMatrix_O& orig) : core::CxxObject_O(orig)
{
  this->_Triangle = orig._Triangle;
  this->_Columns = orig._Columns;
  this->_Rows = orig._Rows;
}



uint	AbstractLargeSquareMatrix_O::indexFromCoordinates(uint x, uint y )
{
  uint idx = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
  if ( idx == UndefinedUnsignedInt )
  {
    SIMPLE_ERROR("Illegal coordinates for {}[{},{}]" , this->className() , x , y );
  }
  return idx;
}

uint	AbstractLargeSquareMatrix_O::col_indexFromCoordinates(uint x, uint y )
{
  uint idx = this->col_indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
  if ( idx == UndefinedUnsignedInt )
  {
    SIMPLE_ERROR("Illegal coordinates for {}[{},{}]" , this->className() , x , y );
  }
  return idx;
}




void	AbstractLargeSquareMatrix_O::writeMathematica(const string& fileName)
{
  std::ofstream	fout;
  uint		x,y;
  fout.open(fileName.c_str());
  fout << "{" << std::endl;
  fout << std::setiosflags(std::ios::fixed)
       << std::setw(16)
       << std::setprecision(8);
  for ( y = 0; y<this->dimension(); y++ ) {
    if ( y!= 0 ) {
      fout << "," << std::endl;
    }
    fout << "{";
    for ( x = 0; x<this->dimension(); x++ ) {
      if ( x!=0 ) {
        fout << ",";
      }
      fout << this->element(x,y);
    }
    fout << "}";
  }
  fout << "}" << std::endl;
  fout.close();
}

void	AbstractLargeSquareMatrix_O::writeMathematicaSymbolic(const string& fileName)
{
  std::ofstream	fout;
  uint		x,y;
  fout.open(fileName.c_str());
  fout << "{" << std::endl;
  fout << std::setiosflags(std::ios::fixed)
       << std::setw(16)
       << std::setprecision(8);
  for ( y = 0; y<this->dimension(); y++ ) {
    if ( y!= 0 ) {
      fout << "," << std::endl;
    }
    fout << "{";
    for ( x = 0; x<this->dimension(); x++ ) {
      if ( x!=0 ) {
        fout << ",";
      }
      if ( this->hasElement(x,y) ) {
        fout << "1.0";
      } else {
        fout << "0.0";
      }
    }
    fout << "}";
  }
  fout << "}" << std::endl;
  fout.close();
}

//! Dump the matrix from 0,0-upto,upto
void	AbstractLargeSquareMatrix_O::dumpMatrix( uint upto)
{
  uint		x,y;
  for ( y = 0; y<this->dimension(); y++ ) {
    for ( x = 0; x<=y; x++ ) {
      uint idx = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
      if (idx != UndefinedUnsignedInt) {
        vecreal val = this->element(x,y);
        core::clasp_write_string(fmt::format("{} {} -> {:12.7f}\n" , x, y, val ));
      }
    }
  }
}



void	AbstractLargeSquareMatrix_O::multiplyByVector(NVector_sp result, NVector_sp d)
{
  uint	x,y;
  vecreal	sum;
  if ( result->size() != this->_Rows ) {
    SIMPLE_ERROR("Result vector does not have the correct dimension");
  }
  if ( d->size() != this->_Rows ) {
    SIMPLE_ERROR("Argument vector does not have the correct dimension");
  }

  for ( y=0; y<this->_Rows; y++ ) {
    sum = 0.0;
    for ( x=0; x<this->_Rows; x++ ) {
      sum += this->element(x,y)*d->element(x);
    }
    result->setElement(y,sum);
  }
}





// ****************************************************
// ****************************************************
// ****************************************************
// ****************************************************
// ****************************************************

FullLargeSquareMatrix_sp FullLargeSquareMatrix_O::create(uint dim, TriangleType type)
{
  auto  res  = gctools::GC<FullLargeSquareMatrix_O>::allocate();
  res->setup(dim,type);
  return res;
}



void FullLargeSquareMatrix_O::setup(uint dim, TriangleType type)
{
  this->AbstractLargeSquareMatrix_O::setup(dim,type);
  if ( type == UpperDiagonalLower ) {
    this->_ActiveElements = dim*dim;
  } else if ( type == SymmetricDiagonalLower ) {
    this->_ActiveElements = ((dim+1)*dim)/2;
  } else if ( type == SymmetricUpperDiagonal ) {
    this->_ActiveElements = ((dim+1)*dim)/2;
  }
  this->_Values.resize(this->_ActiveElements,0);
}

FullLargeSquareMatrix_O::FullLargeSquareMatrix_O(const FullLargeSquareMatrix_O& orig)
    : AbstractLargeSquareMatrix_O(orig)
{
  DEPRECATED();
#if 0
  // This will fail because it calls an allocator Schafmeister Nov 2017
  this->_ActiveElements = orig._ActiveElements;
  this->_Values.assign(orig._Values.begin(),orig._Values.end());
#endif
}

void	FullLargeSquareMatrix_O::fill(vecreal d)
{
  uint	ii, iEnd;
  ii = this->indexBegin();
  iEnd = this->indexEnd();
  while (ii<iEnd) {
    this->setAtIndex(ii,d);
    ii = this->indexAdvance(ii);
  }
}

vecreal	FullLargeSquareMatrix_O::maxAbsValue()
{
  uint	ii, iEnd;
  vecreal	dMaxAbs;
  dMaxAbs = 0.0;
  ii = this->indexBegin();
  iEnd = this->indexEnd();
  while (ii<iEnd) {
    dMaxAbs = MAX(dMaxAbs,fabs(this->getAtIndex(ii)));
    ii = this->indexAdvance(ii);
  }
  return dMaxAbs;
}


uint	FullLargeSquareMatrix_O::indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y ) {
  uint		i0,i;
  if ( x >= this->_Columns || y >= this->_Rows ) {
    return UndefinedUnsignedInt;
  }
  if ( this->_Triangle == UpperDiagonalLower ) {
    i = y*this->_Columns+x;
    if ( i >= this->_ActiveElements ) {
      return UndefinedUnsignedInt;
    }
    return i;
  } else if ( this->_Triangle == SymmetricDiagonalLower ) {
    LOG("SymmetricDiagonalLower matrix" );
    if ( x > y ) {
      LOG("Swapping x,y" );
      uint t = y;
      y = x;
      x = t;
    }
    i0 = ((y+1)*y)/2;
    i = i0 + x;
    LOG("i0({}) i({})" , i0 , i  );
    if ( i>=this->_ActiveElements ) {
      return UndefinedUnsignedInt;
    }
    return i;
  } else if ( this->_Triangle == SymmetricUpperDiagonal ) {
    if ( x < y ) {
      uint t = y;
      y = x;
      x = t;
    }
    i0 = ((y+1)*y)/2;
    i = i0 + (x-y);
    if ( i>=this->_ActiveElements ) {
      return UndefinedUnsignedInt;
    }
    return i;
  }
  return UndefinedUnsignedInt;
}





AbstractLargeSquareMatrix_sp FullLargeSquareMatrix_O::copy_matrix()
{
  auto  dest  = gctools::GC<FullLargeSquareMatrix_O>::copy( *this); // = RP_Copy<FullLargeSquareMatrix_O>(this);
  return dest;
}


void	FullLargeSquareMatrix_O::coordinatesFromIndex(uint ii, uint& x, uint& y)
{
  uint		f2;
  uint		i;
  if ( ii < this->indexBegin() || ii >= this->indexEnd() ) {
    SIMPLE_ERROR("Out of bounds in indexesFromValuePtr");
  }
  if ( this->_Triangle == UpperDiagonalLower ) {
    y = ii/this->_Columns;
    x = ii%this->_Columns;
  } else if ( this->_Triangle == SymmetricDiagonalLower ) {
    f2 = ii*2;
    for (i=0; i<=this->_Rows; i++ ) {
      if ( (i*i+i)>f2 ) break;
    }
    y = i-1;
    x = ii-((y+1)*y)/2;
  } else if ( this->_Triangle == SymmetricUpperDiagonal ) {
    f2 = ii*2;
    for (i=0; i<=this->_Rows; i++ ) {
      if ( (i*i+i)>f2 ) break;
    }
    y = i-1;
    x = ii-((y+1)*y)/2+y;
		// I hacked this routine from the one above
		// I didn't test it.
    SIMPLE_ERROR("Test this function rigorously before using it");
  } else {
    SIMPLE_ERROR("Illegal triangle type");
  }
}










// ****************************************************
// ****************************************************
// ****************************************************
// ****************************************************
// ****************************************************

SYMBOL_EXPORT_SC_(ChemPkg, UpperDiagonalLower);
SYMBOL_EXPORT_SC_(ChemPkg, UpperDiagonal );
SYMBOL_EXPORT_SC_(ChemPkg, DiagonalLower );
SYMBOL_EXPORT_SC_(ChemPkg, SymmetricUpperDiagonal );
SYMBOL_EXPORT_SC_(ChemPkg, SymmetricDiagonalLower );

CL_DEFUN SparseLargeSquareMatrix_sp chem__make_sparse_large_square_matrix(size_t dim, core::T_sp type) {
  if (type == _sym_UpperDiagonalLower) {
    return SparseLargeSquareMatrix_O::create(dim,UpperDiagonalLower);
  }
  if (type == _sym_UpperDiagonal ) {
    return SparseLargeSquareMatrix_O::create(dim,UpperDiagonal);
  }
  if (type == _sym_DiagonalLower ) {
    return SparseLargeSquareMatrix_O::create(dim,DiagonalLower);
  }
  if (type == _sym_SymmetricUpperDiagonal ) {
    return SparseLargeSquareMatrix_O::create(dim,SymmetricUpperDiagonal);
  }
  if (type == _sym_SymmetricDiagonalLower ) {
    return SparseLargeSquareMatrix_O::create(dim,SymmetricDiagonalLower);
  }
  SIMPLE_ERROR("Illegal type {}", _rep_(type));
}


SparseLargeSquareMatrix_sp SparseLargeSquareMatrix_O::create(uint dim, TriangleType type)
{
  auto  res  = gctools::GC<SparseLargeSquareMatrix_O>::allocate();
  res->setup(dim,type);
  return res;
}

void SparseLargeSquareMatrix_O::fill(vecreal val)
{
  if (this->_InsertionIsComplete) {
    for ( uint idx=0; idx<this->_Values.size(); idx++ ) {
      this->_Values[idx] = 0.0;
    }
  }
}

vecreal	SparseLargeSquareMatrix_O::maxAbsValue()
{
  vecreal dMax = 0.0;
  for ( uint idx=0; idx<this->_Values.size(); idx++ ) {
    dMax = MAX(std::fabs(dMax),this->_Values[idx]);
  }
  return dMax;
}

vecreal	SparseLargeSquareMatrix_O::mostPositiveValue()
{
  vecreal dMax = 0.0;
  for ( uint idx=0; idx<this->_Values.size(); idx++ ) {
    dMax = MAX(dMax,this->_Values[idx]);
  }
  return dMax;
}


vecreal	SparseLargeSquareMatrix_O::mostNegativeValue()
{
  vecreal dMin = 0.0;
  for ( uint idx=0; idx<this->_Values.size(); idx++ ) {
    dMin = MIN(dMin,this->_Values[idx]);
  }
  return dMin;
}


SparseLargeSquareMatrix_sp SparseLargeSquareMatrix_O::optimized() 
{
  auto  res  = gctools::GC<SparseLargeSquareMatrix_O>::allocate();
  res->setup(this->_Columns,this->_Triangle);
  uint		x,y;
  for ( y = 0; y<this->dimension(); y++ ) {
    for ( x = 0; x<=y; x++ ) {
      uint idx = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
      if (idx != UndefinedUnsignedInt) {
        vecreal val = this->element(x,y);
        if (val!=0.0) {
          res->insertElement(x,y);
        }
      }
    }
  }
  for ( y = 0; y<this->dimension(); y++ ) {
    for ( x = 0; x<=y; x++ ) {
      uint idx = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
      if (idx != UndefinedUnsignedInt) {
        vecreal val = this->element(x,y);
        if (val!=0.0) {
          res->setElement(x,y,val);
        }
      }
    }
  }
  return res;
}


SparseLargeSquareMatrix_O::~SparseLargeSquareMatrix_O( )
{
#ifdef	MEISTER_MALLOC_DEBUG
  printf("Doing dtor for SparseLargeSquareMatrix_O\n");
#endif
  this->releaseStorage();
}


void	SparseLargeSquareMatrix_O::setup(uint dim, TriangleType type)
{
  this->AbstractLargeSquareMatrix_O::setup(dim,type);
  this->_InsertionIsComplete = false;
  this->_col_OptimizationDone = false;
  this->initializeStorage();
#ifdef	MEISTER_MALLOC_DEBUG
  printf("Done ctor for SparseLargeSquareMatrix_O\n");
#endif
}

bool similar(vecreal a, vecreal b) {
  if (a==0.0 && b==0.0) return true;
  if (std::fabs(a - b) < 0.00001) return true;
  return false;
}

void SparseLargeSquareMatrix_O::ensureIdentical(SparseLargeSquareMatrix_sp other)
{
  if (this->_Triangle != other->_Triangle) SIMPLE_ERROR("{}: _Triangle doesn't match", __FUNCTION__);
  if (this->_Columns != other->_Columns)
    SIMPLE_ERROR("{}: _Columns don't match", __FUNCTION__ );
  if (this->_Rows != other->_Rows)
    SIMPLE_ERROR("{}: _Rows don't match", __FUNCTION__ );
  for ( uint row=0; row<this->_Rows; row++ ) {
    for ( uint col=0; col<this->_Columns; col++ ) {
      vecreal mine = this->element(col,row);
      vecreal oth = other->element(col,row);
      if (!similar(mine,oth)) {
        SIMPLE_ERROR("{}: _Values does't match col,row={},{} this->element(col,row) = {} other->element(col,row) = {}",
                     __FUNCTION__, col, row,
                     mine, oth );
      }
    }
  }
}




SparseLargeSquareMatrix_O::SparseLargeSquareMatrix_O(const SparseLargeSquareMatrix_O& orig)
    : AbstractLargeSquareMatrix_O(orig)
{
  // THis will fail because it calls an allocator
  this->_InsertionIsComplete = orig._InsertionIsComplete;
  this->_RowStartEntries = orig._RowStartEntries;
  this->_RowStarts.assign(orig._RowStarts.begin(),orig._RowStarts.end());
  this->_ActiveElements = orig._ActiveElements;
  this->_ReservedElements = orig._ReservedElements;
  this->_ColumnForValue.assign(orig._ColumnForValue.begin(),orig._ColumnForValue.end());
  this->_Values.assign(orig._Values.begin(),orig._Values.end());
  this->_col_OptimizationDone = false;
}


void	SparseLargeSquareMatrix_O::initializeStorage()
{
  uint	rows = this->_Rows;
  this->_InsertionIsComplete = false;
  this->_RowStartEntries = rows+1;
  this->_RowStarts.assign((rows+1),0);
  this->_ActiveElements = 0;
  this->_ReservedElements = 100;
  this->_ColumnForValue.reserve(this->_ReservedElements);
  this->_Values.reserve(this->_ReservedElements);
  this->_col_ColumnStartEntries = rows+1;
  this->_col_ColumnStarts.assign((rows+1),0);
  this->_col_ActiveElements = 0;
  this->_col_ReservedElements = 100;
  this->_col_RowForValue.reserve(this->_col_ReservedElements);
  this->_col_Values.reserve(this->_col_ReservedElements);
  LOG("Initialized everything, returning" );
}


void	SparseLargeSquareMatrix_O::releaseStorage()
{
  this->_RowStarts.clear();
  this->_ColumnForValue.clear();
  this->_Values.clear();
  this->_RowStartEntries = 0;
  this->_ReservedElements = 0;
  this->_ActiveElements = 0;
  this->_col_ColumnStarts.clear();
  this->_col_RowForValue.clear();
  this->_col_Values.clear();
  this->_col_ColumnStartEntries = 0;
  this->_col_ReservedElements = 0;
  this->_col_ActiveElements = 0;
}


void	SparseLargeSquareMatrix_O::expandStorage()
{
  this->_ReservedElements *= 2;
  this->_ColumnForValue.reserve(this->_ReservedElements);
  this->_Values.reserve(this->_ReservedElements);
}

void	SparseLargeSquareMatrix_O::col_expandStorage()
{
  this->_col_ReservedElements *= 2;
  this->_col_RowForValue.reserve(this->_col_ReservedElements);
  this->_col_Values.reserve(this->_col_ReservedElements);
}


void	SparseLargeSquareMatrix_O::reset()
{
  this->releaseStorage();
  this->initializeStorage();
}


CL_DEFUN core::T_sp chem__sparse_large_square_matrix_layout(SparseLargeSquareMatrix_sp matrix)
{
  ql::list l;
  for ( int row=0; row<matrix->_Rows; row++ ) {
    uint start_row_index = matrix->_RowStarts[row];
    uint end_row_index = matrix->_RowStarts[row+1];
    ql::list lc;
    lc << core::make_fixnum(row);
    for ( int col_index = start_row_index; col_index<end_row_index; col_index++ ) {
      lc << core::make_fixnum(matrix->_ColumnForValue[col_index]);
    }
    l << lc.cons();
  }
  return l.cons();
}

void	SparseLargeSquareMatrix_O::insertElement(uint x, uint y)
{
  uint	ib,ie,i,t;
#ifdef	DEBUG_VIA_PRINTF
  printf("insert at %d,%d\n",x,y);
#endif
  if ( this->_InsertionIsComplete ) {
    SIMPLE_ERROR("SparseMatrix InsertionIsComplete so no more entries may be inserted");
  }
  if ( y>=this->_Rows || x>=this->_Columns ) {
    SIMPLE_ERROR("Overflow in matrix operation");
  }
  if ( this->_Triangle == SymmetricDiagonalLower && x>y ) {
    SWAP(x,y,t);
#ifdef	DEBUG_VIA_PRINTF
    printf("Swapping indexes to %d,%d\n",x,y);
#endif
  } else if ( this->_Triangle == SymmetricUpperDiagonal && x<y ) {
    SWAP(x,y,t);
#ifdef	DEBUG_VIA_PRINTF
    printf("Swapping indexes to %d,%d\n",x,y);
#endif
  }
  ib = this->_RowStarts[y];
  ie = this->_RowStarts[y+1];
  for (i=ib; i<ie; i++ ) {
    if ( this->_ColumnForValue[i] == x ) {
#ifdef	DEBUG_VIA_PRINTF
      printf("The element %d,%d is already in the sparse matrix\n",x,y);
#endif
	    // It's already in the sparse matrix
      return;
    }
    if ( this->_ColumnForValue[i] > x ) {
      break;
    }
  }
		// (i) contains the first element to move up
		// increase the size of the Values and ColumnForValue arrays
		// if they don't have enough space
  if ( this->_ActiveElements>=this->_ReservedElements ) {
    this->expandStorage();
  }
	// Slide the _Values and _ColumnForValue arrays up one to make
	// room for this new element
#ifdef	DEBUG_VIA_PRINTF
  printf("Sliding from %d to %d num=%d\n", i, i+1, (this->_ActiveElements-i) );
  printf("Befor slide value@%d = %d\n", i, this->_ColumnForValue[i]);
#endif
  this->_ColumnForValue.insert(this->_ColumnForValue.begin()+i,x);
  this->_Values.insert(this->_Values.begin()+i,0.0);
#ifdef	DEBUG_VIA_PRINTF
  printf("After slide value@%d = %d\n", i+1, this->_ColumnForValue[i+1]);
#endif
  this->_ActiveElements++;
		//
		// Adjust the row start entries after this one
		//
  for ( uint rr = y+1; rr<this->_RowStartEntries;rr++ ) {
    this->_RowStarts[rr]++;
  }
}


void	SparseLargeSquareMatrix_O::col_insertElement(uint x, uint y)
{
  uint	ib,ie,i,t;
#ifdef	DEBUG_VIA_PRINTF
  printf("insert at %d,%d\n",x,y);
#endif
  if ( this->_InsertionIsComplete ) {
    SIMPLE_ERROR("SparseMatrix InsertionIsComplete so no more entries may be inserted");
  }
  if ( y>=this->_Rows || x>=this->_Columns ) {
    SIMPLE_ERROR("Overflow in matrix operation");
  }
  if ( this->_Triangle == SymmetricDiagonalLower && x>y ) {
    SWAP(x,y,t);
#ifdef	DEBUG_VIA_PRINTF
    printf("Swapping indexes to %d,%d\n",x,y);
#endif
  } else if ( this->_Triangle == SymmetricUpperDiagonal && x<y ) {
    SWAP(x,y,t);
#ifdef	DEBUG_VIA_PRINTF
    printf("Swapping indexes to %d,%d\n",x,y);
#endif
  }
  ib = this->_col_ColumnStarts[y];
  ie = this->_col_ColumnStarts[y+1];
  for (i=ib; i<ie; i++ ) {
    if ( this->_col_RowForValue[i] == x ) {
#ifdef	DEBUG_VIA_PRINTF
      printf("The element %d,%d is already in the sparse matrix\n",x,y);
#endif
	    // It's already in the sparse matrix
      return;
    }
    if ( this->_ColumnForValue[i] > x ) {
      break;
    }
  }
		// (i) contains the first element to move up
		// increase the size of the Values and ColumnForValue arrays
		// if they don't have enough space
  if ( this->_col_ActiveElements>=this->_col_ReservedElements ) {
    this->col_expandStorage();
  }
	// Slide the _Values and _ColumnForValue arrays up one to make
	// room for this new element
#ifdef	DEBUG_VIA_PRINTF
  printf("Sliding from %d to %d num=%d\n", i, i+1, (this->_col_ActiveElements-i) );
  printf("Befor slide value@%d = %d\n", i, this->_col_ColumnForValue[i]);
#endif
  this->_col_RowForValue.insert(this->_col_RowForValue.begin()+i,x);
  this->_col_Values.insert(this->_col_Values.begin()+i,0.0);
#ifdef	DEBUG_VIA_PRINTF
  printf("After slide value@%d = %d\n", i+1, this->_col_RowForValue[i+1]);
#endif
  this->_col_ActiveElements++;
		//
		// Adjust the row start entries after this one
		//
  for ( uint rr = y+1; rr<this->_col_ColumnStartEntries;rr++ ) {
    this->_col_ColumnStarts[rr]++;
  }
}

//
// Return the index for the element at x,y
// If there is no element return UndefinedUnsignedInt 
uint	SparseLargeSquareMatrix_O::indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y )
{
  uint	im, t, ret;
  if ( x >= this->_Columns || y >= this->_Rows ) {
    return UndefinedUnsignedInt;
  }
  if ( this->_Triangle == SymmetricDiagonalLower && x>y ) {
    SWAP(x,y,t);
    LOG("Swapped row/column" );
  } else if ( this->_Triangle == SymmetricUpperDiagonal && x<y ) {
    SWAP(x,y,t);
    LOG("Swapped row/column" );
  }
  LOG("Looking at row: {}" , y  );
  uint ib = this->_RowStarts[y];
  uint ie = this->_RowStarts[y+1]-1;
  LOG("This row has columns in the list between indexes: {}-{}" , ib , ie );
  if ( ie<ib ) return UndefinedUnsignedInt;
  if ( ib >= this->_ColumnForValue.size() ) return UndefinedUnsignedInt;
  if ( x < this->_ColumnForValue[ib] ) return UndefinedUnsignedInt;
  if ( ie >= this->_ColumnForValue.size() ) return UndefinedUnsignedInt;
  if ( x > this->_ColumnForValue[ie] ) return UndefinedUnsignedInt;
  if ( x == this->_ColumnForValue[ib] )
  {
    ret = ib;
    goto gotOne;
  }
  if ( x == this->_ColumnForValue[ie] )
  {
    ret = ie;
    goto gotOne;
  }
  while ( ib<ie-1 ) {
    im = (ib+ie)/2;
    if ( this->_ColumnForValue[im] == x ) 
    {
      ret = im;
      goto gotOne;
    }
    if ( this->_ColumnForValue[im] < x ) {
      ib = im;
    } else {
      ie = im;
    }
  }
  LOG("Exhausted columns, returning UndefinedUnsignedInt" );
  return UndefinedUnsignedInt;
 gotOne:
  if ( ret >= this->_Values.size() )
  {
    core::clasp_write_string(fmt::format("Out of bounds at: {} line: {}\n" , __FILE__ , __LINE__ ));
    core::clasp_write_string(fmt::format("Looking for point: {}, {}\n" , x,  y ));
    core::clasp_write_string(fmt::format("Calculated: {}\n" , ret ));
    core::clasp_write_string(fmt::format("this->_Values.size() = {}\n" , this->_Values.size() ));
    this->debug();
    SIMPLE_ERROR("Calculated a position that is beyond the range of Values");
  }
  return ret;
}

uint	SparseLargeSquareMatrix_O::col_indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y )
{
  uint	im, t, ret;
  if ( x >= this->_Columns || y >= this->_Rows ) {
    return UndefinedUnsignedInt;
  }
  if ( this->_Triangle == SymmetricDiagonalLower && x>y ) {
    SWAP(x,y,t);
    LOG("Swapped row/column" );
  } else if ( this->_Triangle == SymmetricUpperDiagonal && x<y ) {
    SWAP(x,y,t);
    LOG("Swapped row/column" );
  }
  LOG("Looking at row: {}" , y  );
  uint ib = this->_col_ColumnStarts[y];
  uint ie = this->_col_ColumnStarts[y+1]-1;
  if ( ie<ib ) return UndefinedUnsignedInt;
  if ( ib >= this->_col_RowForValue.size() ) return UndefinedUnsignedInt;
  if ( x < this->_col_RowForValue[ib] ) return UndefinedUnsignedInt;
  if ( ie >= this->_col_RowForValue.size() ) return UndefinedUnsignedInt;
  if ( x > this->_col_RowForValue[ie] ) return UndefinedUnsignedInt;
  if ( x == this->_col_RowForValue[ib] )
  {
    ret = ib;
    goto gotOne;
  }
  if ( x == this->_col_RowForValue[ie] )
  {
    ret = ie;
    goto gotOne;
  }
  while ( ib<ie-1 ) {
    im = (ib+ie)/2;
    if ( this->_col_RowForValue[im] == x ) 
    {
      ret = im;
      goto gotOne;
    }
    if ( this->_col_RowForValue[im] < x ) {
      ib = im;
    } else {
      ie = im;
    }
  }
  LOG("Exhausted columns, returning UndefinedUnsignedInt" );
  return UndefinedUnsignedInt;
 gotOne:
  if ( ret >= this->_col_Values.size() )
  {
    core::clasp_write_string(fmt::format("Out of bounds at: {} line: {}\n" , __FILE__ , __LINE__ ));
    core::clasp_write_string(fmt::format("Looking for point: {}, {}\n" , x,  y ));
    core::clasp_write_string(fmt::format("Calculated: {}\n" , ret ));
    core::clasp_write_string(fmt::format("this->_Values.size() = {}\n" , this->_Values.size() ));
    this->debug();
    SIMPLE_ERROR("Calculated a position that is beyond the range of Values");
  }
  return ret;
}



CL_DEFMETHOD AbstractLargeSquareMatrix_sp	SparseLargeSquareMatrix_O::copy_matrix()
{
  auto  dest  = gctools::GC<SparseLargeSquareMatrix_O>::copy( *this); // = RP_Copy<SparseLargeSquareMatrix_O>(this->sharedThis<SparseLargeSquareMatrix_O>());
  return dest;
}

void SparseLargeSquareMatrix_O::insertionIsComplete() {
  this->_InsertionIsComplete = true;
}

void SparseLargeSquareMatrix_O::doColumnOptimization() {
  if (!this->_col_OptimizationDone) {
    if (!this->_InsertionIsComplete) {
      for ( int row = 0; row < this->_Rows; row++ ) {
        for ( int col = 0; col <= row; col++ ) {
          uint idx = this->indexFromCoordinatesOrUndefinedUnsignedInt(col,row);
          if (idx != UndefinedUnsignedInt) {
            vecreal val = this->getAtIndex(idx);
            this->col_setElement(col,row,val);
          }
        }
      }
#if 1
      for ( int row = 0; row < this->_Rows; row++ ) {
        for ( int col = 0; col <= row; col++ ) {
          uint idx = this->indexFromCoordinatesOrUndefinedUnsignedInt(col,row);
          if (idx != UndefinedUnsignedInt) {
            vecreal row_val = this->getAtIndex(idx);
            vecreal col_val = this->col_element(col,row);
            if ( row_val != col_val ) {
              SIMPLE_ERROR("There is a mismatch between the SparseLargeSquareMatrix entry %d, %d  row_val = %f, col_val = %f",
                           col, row, row_val, col_val );
            }
          }
        }
      }
#endif
    }
  }
  this->_col_OptimizationDone = true;
}


void	SparseLargeSquareMatrix_O::coordinatesFromIndex(uint is, uint& x, uint& y )
{
  uint	r;
  for ( r=0; r<this->_Rows; r++ ) {
    if ( this->_RowStarts[r] > is ) break;
  }
  if ( this->_RowStarts[r]<=is ) {
    y=UndefinedUnsignedInt;
    x=UndefinedUnsignedInt;
  }
  y = r-1;
  x = this->_ColumnForValue[is];
}


void	SparseLargeSquareMatrix_O::debug()
{
  uint	ib, ie;
  core::clasp_write_string(fmt::format("this->_Rows = {}\n" , this->_Rows ));
  core::clasp_write_string(fmt::format("this->_Columns = {}\n" , this->_Columns ));
  for ( uint r=0; r<this->_Rows; r++ ) {
    ib = this->_RowStarts[r];
    ie = this->_RowStarts[r+1];
    core::clasp_write_string(fmt::format("Row{:3} [{:3}-{:3}]: \n" , r , ib , ie ));
    for (uint i=ib;i<ie;i++ ) 
    {
      core::clasp_write_string(fmt::format(" x{}@{:<3}\n" , this->_ColumnForValue[i] , i ));
    }
    core::cl__terpri();
  }
  core::clasp_write_string("Columns: \n");
  for ( uint i=0; i<this->_ActiveElements; i++ ) {
    core::clasp_write_string(fmt::format("x{}@{} ", this->_ColumnForValue[i] , i));
  }
  core::clasp_write_string(fmt::format("Number of values: {}\n" , this->_Values.size() ));
}

CL_DEFMETHOD void	SparseLargeSquareMatrix_O::walkMatrix(core::T_sp callback)
{
  uint	ib, ie;
  for ( uint r=0; r<this->_Rows; r++ ) {
    ib = this->_RowStarts[r];
    ie = this->_RowStarts[r+1];
    for (uint i=ib;i<ie;i++ ) {
#if VECREAL==VECREAL_FLOAT
      core::T_sp val = core::make_single_float(this->_Values[i]);
#else
      core::T_sp val = mk_double_float(this->_Values[i]);
#endif
      core::eval::funcall( callback, core::make_fixnum(this->_ColumnForValue[i]), core::make_fixnum(r), val );
    }
  }
}



vecreal	SparseLargeSquareMatrix_O::element(uint x, uint y)
{
  uint ii= this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
  if ( ii == UndefinedUnsignedInt ) return 0.0;
  return this->_Values[ii];
}

vecreal	SparseLargeSquareMatrix_O::col_element(uint x, uint y)
{
  uint ii= this->col_indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
  if ( ii == UndefinedUnsignedInt ) return 0.0;
  return this->_col_Values[ii];
}

void	SparseLargeSquareMatrix_O::setElement(uint x, uint y, vecreal d)
{
  if (d!=0.0) {
    uint	dp;
    dp = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
    if ( dp == UndefinedUnsignedInt ) {
      LOG("Could not find indexFromCoordinates" );
      if ( this->_InsertionIsComplete ) {
        stringstream ss;
        ss << "Trying to write to non-existent element (x,y)=("<< x << "," << y <<") in sparse matrix";
        SIMPLE_ERROR("{}" , ss.str());
      } else {
        this->insertElement(x,y);
        dp = this->indexFromCoordinates(x,y);
      }
    } else {
      LOG("Found indexFromCoordinates dp={}" , dp );
    }
    this->_Values[dp] = d;
  }
}

void	SparseLargeSquareMatrix_O::col_setElement(uint x, uint y, vecreal d)
{
  if (d!=0.0) {
    uint	dp;
    dp = this->col_indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
    if ( dp == UndefinedUnsignedInt ) {
      LOG("Could not find indexFromCoordinates" );
      this->col_insertElement(x,y);
      dp = this->col_indexFromCoordinates(x,y);
    } else {
      LOG("Found indexFromCoordinates dp={}" , dp );
    }
    this->_col_Values[dp] = d;
  }
}

CL_DEFUN void chem__sparse_large_square_matrix_set_element(SparseLargeSquareMatrix_sp matrix, uint col, uint row, vecreal value) {
  matrix->setElement(col,row,value);
}

CL_DEFUN vecreal chem__sparse_large_square_matrix_element(SparseLargeSquareMatrix_sp matrix, uint col, uint row) {
  return matrix->element(col,row);
}

CL_DEFUN void chem__sparse_large_square_matrix_insertion_is_complete(SparseLargeSquareMatrix_sp matrix) {
  matrix->insertionIsComplete();
}

void	SparseLargeSquareMatrix_O::addToElement(uint x, uint y, vecreal d)
{
  uint	dp, idp;
  bool	inserted = false;
  dp = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
  idp = dp;
  if ( dp == UndefinedUnsignedInt ) {
    if ( this->_InsertionIsComplete ) {
      stringstream ss;
      ss << "Trying to write to non-existent element (x,y)=("<< x << "," << y <<") in sparse matrix";
      SIMPLE_ERROR("{}" , ss.str());
    } else {
      inserted = true;
      this->insertElement(x,y);
      dp = this->indexFromCoordinates(x,y);
    }
  }
  if ( this->_Values.size() == 0 )
  {
	    // This should never happen! The _Values array should have had one element inserted
    core::clasp_write_string(fmt::format("Failing in {} line: {}\n" , __FILE__ , __LINE__ ));
    core::clasp_write_string(fmt::format("Trying to add to element at: {}, {}\n" , x,  y ));
    core::clasp_write_string(fmt::format( "inserted dp = {}\n" , dp ));
    core::clasp_write_string(fmt::format( "initial dp = {}\n" , idp ));
    core::clasp_write_string(fmt::format( "inserted = {}\n" , inserted ));
    core::clasp_write_string(fmt::format( "this->_Values.size() = {}\n" , this->_Values.size() ));
    core::clasp_write_string(fmt::format( "this->_InsertionIsComplete = {}\n" , this->_InsertionIsComplete ));
  }
  ASSERTP(dp!=UndefinedUnsignedInt,"Illegal index for SparseLargeSquareMatrix");
  ASSERT_lessThan(dp,this->_Values.size());
  this->_Values[dp] += d;
}


uint	SparseLargeSquareMatrix_O::indexOfLastElementOnRow(uint y)
{
  ASSERT_lessThan(y,this->_Rows);
  return this->_RowStarts[y+1]-1;
}

// Return the index for the element at x,y
// If there is no element return UndefinedUnsignedInt
uint	SparseLargeSquareMatrix_O::indexOfFirstElementAtOrAfterX(uint x, uint y )
{
  uint	im;
  ASSERT_lessThan(x,this->_Columns);
  ASSERT_lessThan(y,this->_Rows);
  ASSERT(this->_Triangle == SymmetricDiagonalLower);
  LOG("Looking at row: {}" , y  );
  uint ib = this->_RowStarts[y];
  uint ie = this->_RowStarts[y+1]-1;
  LOG("This row has columns in the list between indexes: {}-{}" , ib , ie );
  if ( ie<ib ) return UndefinedUnsignedInt;
  if ( x < this->_ColumnForValue[ib] ) return ib;
  if ( x > this->_ColumnForValue[ie] ) return this->_RowStarts[y+1];
  if ( x == this->_ColumnForValue[ib] ) return ib;
  if ( x == this->_ColumnForValue[ie] ) return ie;
  while ( ib<ie-1 ) {
    im = (ib+ie)/2;
    if ( this->_ColumnForValue[im] == x ) return im;
    if ( this->_ColumnForValue[im] < x ) {
      ib = im;
    } else {
      ie = im;
    }
  }
  return ie;
}


};
