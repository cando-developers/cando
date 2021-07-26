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
    this->_Columns = dim;
    this->_Rows = dim;
};

    AbstractLargeSquareMatrix_O::AbstractLargeSquareMatrix_O(const AbstractLargeSquareMatrix_O& orig) : core::CxxObject_O(orig)
{
    this->_Triangle = orig._Triangle;
    this->_Columns = orig._Columns;
    this->_Rows = orig._Rows;
}


void	AbstractLargeSquareMatrix_O::fill(double d)
{
uint	ii, iEnd;
    ii = this->indexBegin();
    iEnd = this->indexEnd();
    while (ii<iEnd) {
	this->setAtIndex(ii,d);
	ii = this->indexAdvance(ii);
    }
}


double	AbstractLargeSquareMatrix_O::maxAbsValue()
{
uint	ii, iEnd;
double	dMaxAbs;
    dMaxAbs = 0.0;
    ii = this->indexBegin();
    iEnd = this->indexEnd();
    while (ii<iEnd) {
	dMaxAbs = MAX(dMaxAbs,fabs(this->getAtIndex(ii)));
	ii = this->indexAdvance(ii);
    }
    return dMaxAbs;
}


uint	AbstractLargeSquareMatrix_O::indexFromCoordinates(uint x, uint y )
{
    uint idx = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
    if ( idx == UndefinedUnsignedInt )
    {
	SIMPLE_ERROR(BF("Illegal coordinates for %s[%u,%u]") % this->className() % x % y );
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
	for ( x = 0; x<this->dimension(); x++ ) {
	    if ( x!=0 ) {
		_lisp->print(BF(" "));
	    }
	    _lisp->print(BF("%12.7lf") % this->element(x,y));
	}	
	_lisp->print(BF(""));
    }
}



void	AbstractLargeSquareMatrix_O::multiplyByVector(NVector_sp result, NVector_sp d)
{_OF();
uint	x,y;
double	sum;
    if ( result->size() != this->_Rows ) {
	SIMPLE_ERROR(BF("Result vector does not have the correct dimension"));
    }
    if ( d->size() != this->_Rows ) {
	SIMPLE_ERROR(BF("Argument vector does not have the correct dimension"));
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
  auto  res  = gctools::GC<FullLargeSquareMatrix_O>::allocate_with_default_constructor();
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
	LOG(BF("SymmetricDiagonalLower matrix") );
	if ( x > y ) {
	    LOG(BF("Swapping x,y") );
	    uint t = y;
	    y = x;
	    x = t;
	}
	i0 = ((y+1)*y)/2;
	i = i0 + x;
	LOG(BF("i0(%d) i(%d)") % i0 % i  );
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





AbstractLargeSquareMatrix_sp	FullLargeSquareMatrix_O::copy()
{
  auto  dest  = gctools::GC<FullLargeSquareMatrix_O>::copy( *this); // = RP_Copy<FullLargeSquareMatrix_O>(this);
    return dest;
}


void	FullLargeSquareMatrix_O::coordinatesFromIndex(uint ii, uint& x, uint& y)
{_OF();
uint		f2;
uint		i;
    if ( ii < this->indexBegin() || ii >= this->indexEnd() ) {
	SIMPLE_ERROR(BF("Out of bounds in indicesFromValuePtr"));
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
	SIMPLE_ERROR(BF("Test this function rigorously before using it"));
    } else {
	SIMPLE_ERROR(BF("Illegal triangle type"));
    }
}










// ****************************************************
// ****************************************************
// ****************************************************
// ****************************************************
// ****************************************************

SparseLargeSquareMatrix_sp SparseLargeSquareMatrix_O::create(uint dim, TriangleType type)
{
  auto  res  = gctools::GC<SparseLargeSquareMatrix_O>::allocate_with_default_constructor();
    res->setup(dim,type);
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
    this->initializeStorage();
#ifdef	MEISTER_MALLOC_DEBUG
    printf("Done ctor for SparseLargeSquareMatrix_O\n");
#endif
}





SparseLargeSquareMatrix_O::SparseLargeSquareMatrix_O(const SparseLargeSquareMatrix_O& orig)
	: AbstractLargeSquareMatrix_O(orig)
{
  DEPRECATED();
#if 0
  // THis will fail because it calls an allocator
    this->_InsertionIsComplete = orig._InsertionIsComplete;
    this->_RowStartEntries = orig._RowStartEntries;
    this->_RowStarts.assign(orig._RowStarts.begin(),orig._RowStarts.end());
    this->_ActiveElements = orig._ActiveElements;
    this->_ReservedElements = orig._ReservedElements;
    this->_ColumnForValue.assign(orig._ColumnForValue.begin(),orig._ColumnForValue.end());
    this->_Values.assign(orig._Values.begin(),orig._Values.end());
#endif
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
    LOG(BF("Initialized everything, returning") );
}


void	SparseLargeSquareMatrix_O::releaseStorage()
{
    this->_RowStarts.clear();
    this->_ColumnForValue.clear();
    this->_Values.clear();
    this->_RowStartEntries = 0;
    this->_ReservedElements = 0;
    this->_ActiveElements = 0;
}


void	SparseLargeSquareMatrix_O::expandStorage()
{
    this->_ReservedElements *= 2;
    this->_ColumnForValue.reserve(this->_ReservedElements);
    this->_Values.reserve(this->_ReservedElements);
}


void	SparseLargeSquareMatrix_O::reset()
{
    this->releaseStorage();
    this->initializeStorage();
}


void	SparseLargeSquareMatrix_O::insertElement(uint x, uint y)
{
uint	ib,ie,i,t;
#ifdef	DEBUG_VIA_PRINTF
    printf("insert at %d,%d\n",x,y);
#endif
    if ( this->_InsertionIsComplete ) {
	SIMPLE_ERROR(BF("SparseMatrix InsertionIsComplete so no more entries may be inserted"));
    }
    if ( y>=this->_Rows || x>=this->_Columns ) {
	SIMPLE_ERROR(BF("Overflow in matrix operation"));
    }
    if ( this->_Triangle == SymmetricDiagonalLower && x>y ) {
	SWAP(x,y,t);
#ifdef	DEBUG_VIA_PRINTF
	    printf("Swapping indices to %d,%d\n",x,y);
#endif
    } else if ( this->_Triangle == SymmetricUpperDiagonal && x<y ) {
	SWAP(x,y,t);
#ifdef	DEBUG_VIA_PRINTF
	    printf("Swapping indices to %d,%d\n",x,y);
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
	LOG(BF("Swapped row/column") );
    } else if ( this->_Triangle == SymmetricUpperDiagonal && x<y ) {
	SWAP(x,y,t);
	LOG(BF("Swapped row/column") );
    }
    LOG(BF("Looking at row: %d") % y  );
    uint ib = this->_RowStarts[y];
    uint ie = this->_RowStarts[y+1]-1;
    LOG(BF("This row has columns in the list between indices: %d-%d") % ib % ie );
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
    LOG(BF("Exhausted columns, returning UndefinedUnsignedInt") );
    return UndefinedUnsignedInt;
gotOne:
    if ( ret >= this->_Values.size() )
    {
        _lisp->print(BF("Out of bounds at: %s line: %d") % __FILE__% __LINE__ );
	_lisp->print(BF("Looking for point: %u, %u") % x% y );
	_lisp->print(BF("Calculated: %u") % ret );
	_lisp->print(BF("this->_Values.size() = %u") % this->_Values.size() );
	this->debug();
	SIMPLE_ERROR(BF("Calculated a position that is beyond the range of Values"));
    }
    return ret;
}



AbstractLargeSquareMatrix_sp	SparseLargeSquareMatrix_O::copy()
{
  auto  dest  = gctools::GC<SparseLargeSquareMatrix_O>::copy( *this); // = RP_Copy<SparseLargeSquareMatrix_O>(this->sharedThis<SparseLargeSquareMatrix_O>());
    return dest;
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
    _lisp->print(BF("this->_Rows = %d") % this->_Rows );
    _lisp->print(BF("this->_Columns = %d") % this->_Columns );
    for ( uint r=0; r<this->_Rows; r++ ) {
	ib = this->_RowStarts[r];
	ie = this->_RowStarts[r+1];
	_lisp->print(BF("Row%3d [%3d-%3d]: ") % r % ib % ie );
	for (uint i=ib;i<ie;i++ ) 
	{
	    _lisp->print(BF(" x%d@%-3d") % this->_ColumnForValue[i] % i );
	}
	_lisp->print(BF( "" ));
    }
    _lisp->print(BF("Columns: "));
    for ( uint i=0; i<this->_ActiveElements; i++ ) {
	_lisp->prin1(BF("x%d@%d ") % this->_ColumnForValue[i] % i);
    }
    _lisp->print(BF("Number of values: %u") % this->_Values.size() );
}



    double	SparseLargeSquareMatrix_O::element(uint x, uint y)
    {
	uint ii= this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
	if ( ii == UndefinedUnsignedInt ) return 0.0;
	return this->_Values[ii];
    }

    void	SparseLargeSquareMatrix_O::setElement(uint x, uint y, double d)
    {
	uint	dp;
	dp = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
	if ( dp == UndefinedUnsignedInt ) {
	    LOG(BF("Could not find indexFromCoordinates") );
	    if ( this->_InsertionIsComplete ) {
		stringstream ss;
		ss << "Trying to write to non-existent element (x,y)=("<< x << "," << y <<") in sparse matrix";
		SIMPLE_ERROR(BF("%s")%ss.str());
	    } else {
		this->insertElement(x,y);
		dp = this->indexFromCoordinates(x,y);
	    }
	} else {
	    LOG(BF("Found indexFromCoordinates dp=%d") % dp );
	}
	this->_Values[dp] = d;
    }

    void	SparseLargeSquareMatrix_O::addToElement(uint x, uint y, double d)
    {
	uint	dp, idp;
	bool	inserted = false;
	dp = this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y);
	idp = dp;
	if ( dp == UndefinedUnsignedInt ) {
	    if ( this->_InsertionIsComplete ) {
		stringstream ss;
		ss << "Trying to write to non-existent element (x,y)=("<< x << "," << y <<") in sparse matrix";
		SIMPLE_ERROR(BF("%s")%ss.str());
	    } else {
		inserted = true;
		this->insertElement(x,y);
		dp = this->indexFromCoordinates(x,y);
	    }
	}
	if ( this->_Values.size() == 0 )
	{
	    // This should never happen! The _Values array should have had one element inserted
	    _lisp->print(BF("Failing in %s line: %d") % __FILE__% __LINE__ );
	    _lisp->print(BF("Trying to add to element at: %u, %u") % x% y );
	    _lisp->print(BF( "inserted dp = %u") % dp );
	    _lisp->print(BF( "initial dp = %u") % idp );
	    _lisp->print(BF( "inserted = %d") % inserted );
	    _lisp->print(BF( "this->_Values.size() = %u") % this->_Values.size() );
	    _lisp->print(BF( "this->_InsertionIsComplete = %d") % this->_InsertionIsComplete );
	}
	ASSERTP(dp!=UndefinedUnsignedInt,"Illegal index for SparseLargeSquareMatrix");
	ASSERT_lessThan(dp,this->_Values.size());
	this->_Values[dp] += d;
    }


    uint	SparseLargeSquareMatrix_O::indexOfLastElementOnRow(uint y)
    {_OF();
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
	ASSERT(this->_Triangle == SymmetricUpperDiagonal);
	LOG(BF("Looking at row: %d") % y  );
	uint ib = this->_RowStarts[y];
	uint ie = this->_RowStarts[y+1]-1;
	LOG(BF("This row has columns in the list between indices: %d-%d") % ib % ie );
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
