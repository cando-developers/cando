/*
    File: sparseMatrix.h
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
 *      File:   sparseSparseMatrix.h
 *
 *      Description:
 *              Routines to maintain sparse NxM dimensional Matrices
 */

#ifndef SparseSparseMatrix_H
#define	SparseSparseMatrix_H

#include <clasp/core/common.h>

#include <cando/chem/chemPackage.h>


namespace chem {

typedef	enum	{ FullMatrx, UpperTriangular, LowerTriangular } SparseMatrixType;

typedef	struct	{
	int	_Offset;
	double	_Value;
} SparseElement;

typedef	vector<SparseElement>	SparseSlice;

inline string	XmlTag_SparseMatrix() { return "SparseMatrix";};
class	Dumb_SparseMatrix;
    typedef	gctools::smart_ptr<Dumb_SparseMatrix>	RPSparseMatrix;

class	Dumb_SparseMatrix {
public:
	bool			_RowMatrix;
//	SparseMatrixType	_MatrixType;
	int     		_Columns;
	int			_Rows;
	vector<int>		_SliceStarts;
	vector<SparseElement>	_Values;

		Dumb_SparseMatrix( int columns, int rows );
		Dumb_SparseMatrix( );
	virtual	~Dumb_SparseMatrix();

	void	describe();

	void	copy( RPSparseMatrix matrix);


	void	zero();		// SparseMatrixZero();
	void	fill(double d);

//	void	add(RPSparseMatrix a, RPSparseMatrix b);
//	void	sub(RPSparseMatrix a, RPSparseMatrix b);

	double	maxAbsValue();
	double	minAbsValue();
	int	columns()	{ return this->_Columns; };
	int	rows()	{ return this->_Rows; };

	void	expand(int c, int r );

	double	element(int c, int r);
	void	setElement(int c, int r, double d);
	void	addToElement(int c, int r, double d);


};


inline	RPSparseMatrix	new_RPSparseMatrix(int cols, int rows) {
    return RPSparseMatrix(new Dumb_SparseMatrix(cols,rows));
}

};

#endif          /* ifndef SparseMatrix_H */
