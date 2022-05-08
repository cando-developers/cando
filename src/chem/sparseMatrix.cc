/*
    File: sparseMatrix.cc
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
       
#define	TURN_DEBUG_OFF

//
// (C) 2004 Christian E. Schafmeister
//


/*
 *      File:   sparseMatrix.h
 *
 *      Description:
 *              Maintain NxM dimensional Sparse Matrices.
 */

#include <clasp/core/common.h>
//#include	<malloc.h>
#include <iostream>
#include <math.h>
#include <cando/chem/sparseMatrix.h>


using namespace chem;
#ifdef	USEBOOSTPYTHON //[
void	
I nitPython_SparseMatrix()
{

    boost::python::def("XmlTag_SparseMatrix",XmlTag_SparseMatrix);

#ifdef USEBOOSTPYTHON
    PYTHON_CLASS(ChemPkg,Dumb,"","",_lisp)
	.def(boost::python::init<int,int> () )
	;
#endif
}
#endif //]

/*
 *      Dumb_SparseMatrix 
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Allocate memory for the Dumb_SparseMatrixP and initialize it
 *      to zero.
 */
Dumb_SparseMatrix::Dumb_SparseMatrix( SparseMatrixType type, int cols, int rows)
{
    this->_RowMatrix = true;
    this->_MatrixType = type;
    this->_Columns = cols;
    this->_Rows = rows;
    this->_Values.empty();
    this->_SliceStarts.resize(rows+1);
    this->_SliceStarts.fill(0);	// Fill with zero
}


/*
 *      Dumb_SparseMatrix 
 *
 *	Author:	Christian Schafmeister (1991)
 *
 */
Dumb_SparseMatrix::Dumb_SparseMatrix( )
{
    this->_Columns = 0;
    this->_Rows = 0;
    this->_Values.empty();
    this->_SliceStarts.empty(0);

}

//    THE NONZERO ENTRIES OF THE SPARSE MATRIX ARE STORED ROW-BY-ROW IN THE
//    ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO ENTRIES IN EACH ROW,
//    WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY LIES.  THESE COLUMN
//    INDICES ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN
//    JA(K) = J.  TO IDENTIFY THE INDIVIDUAL ROWS, WE NEED TO KNOW WHERE
//    EACH ROW STARTS.  THESE ROW POINTERS ARE STORED IN THE ARRAY IA;
//    I.E., IF M(I,J) IS THE FIRST NONZERO ENTRY (STORED) IN THE I-TH ROW
//    AND  A(K) = M(I,J),  THEN  IA(I) = K.  MOREOVER, IA(N+1) POINTS TO
//    THE FIRST LOCATION FOLLOWING THE LAST ELEMENT IN THE LAST ROW.
//    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS  IA(I+1) - IA(I),
//    THE NONZERO ENTRIES IN THE I-TH ROW ARE STORED CONSECUTIVELY IN
//
//            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),
//
//    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN
//
//            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).
//
//    SINCE THE COEFFICIENT MATRIX IS SYMMETRIC, ONLY THE NONZERO ENTRIES
//    IN THE UPPER TRIANGLE NEED BE STORED, FOR EXAMPLE, THE MATRIX
//
//             ( 1  0  2  3  0 )
//             ( 0  4  0  0  0 )
//         M = ( 2  0  5  6  0 )
//             ( 3  0  6  7  8 )
//             ( 0  0  0  8  9 )
//
//    COULD BE STORED AS
//
//            \ 1  2  3  4  5  6  7  8  9 10 11 12 13
//         ---+--------------------------------------
//         IA \ 1  4  5  8 12 14
//         JA \ 1  3  4  2  1  3  4  1  3  4  5  4  5
//          A \ 1  2  3  4  2  5  6  3  6  7  8  8  9
//
//    OR (SYMMETRICALLY) AS
//
//            \ 1  2  3  4  5  6  7  8  9
//         ---+--------------------------
//         IA \ 1  4  5  7  9 10
//         JA \ 1  3  4  2  3  4  4  5  5
//          A \ 1  2  3  4  5  6  7  8  9          .
//
//


double*	Dumb_SparseMatrix::elementPointer(int c, int r )
{
int	sliceStart, sliceEnd, ss, se, sm;
    this->_orderIndices(c,r);

    if ( c >= this->_Columns || r >= this->_Rows ) {
	SIMPLE_ERROR(("Overflow in matrix operation"));
    }
    ssstart = &(this->_Values[this->_SliceStarts[r]];
    send = &(this->_Values[this->_SliceStarts[r+1]-1];
    if ( c < sstart->_Offset ) return NULL;
    if ( c > send->_Offset ) return NULL;
    if ( c == ss->this->_Values[sliceStart]._Offset ) return &this->_Values[sliceStart]._Value;
    if ( c == this->_Values[sliceEnd]._Offset ) return &this->_Values[sliceEnd]._Value;

    ss = sliceStart;
    se = sliceEnd;
    if ( // Check if rs is the value
    while ( se!=ss ) {
	sm = (se+ss)/2;
	if ( c == this->_Values
	if ( c < rm ) {
	    re = rm;
	} else {
	    rs = rm;
	}
    }


    return this->_dPValues+(r)*this->_Columns+(c);
}


bool	Dumb_SparseMatrix::
