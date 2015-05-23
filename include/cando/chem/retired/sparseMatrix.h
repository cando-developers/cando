       
       
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
