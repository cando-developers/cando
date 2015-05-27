       
       
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
class AbstractLargeSquareMatrix_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,AbstractLargeSquareMatrix_O,"AbstractLargeSquareMatrix");

protected:
    TriangleType 	_Triangle;
    uint   		_Columns;
    uint		_Rows;

public:
protected:
    virtual	uint	indexBegin() {SUBIMP();};
    virtual	uint	indexEnd() {SUBIMP();}
    virtual	uint	indexAdvance(uint ii) {SUBIMP();};
    virtual	void	coordinatesFromIndex(uint ii,uint& row, uint& col ) {SUBIMP();};

public:
    AbstractLargeSquareMatrix_O( const AbstractLargeSquareMatrix_O& orig );


    void	setup(uint dim, TriangleType type);

//	void	describe();

    void	zero() {this->fill(0.0);};


    uint	columns()	{ return this->_Columns; };
    uint	rows()	{ return this->_Rows; };
    uint	dimension()	{return this->_Rows;};

    virtual	double	element(uint x, uint y) {SUBIMP();};
    virtual	void	setElement(uint x, uint y, double d) {SUBIMP();};
    virtual	void	addToElement(uint x, uint y, double d) {SUBIMP();};
    virtual	void	setAtIndex(uint ii, double d) {SUBIMP();};
    virtual	double	getAtIndex(uint ii) {SUBIMP();};
    virtual	void	insertionIsComplete() {};


    double	maxAbsValue();
    void	fill(double d);
    void	writeMathematica(const string& fileName);
    void	writeMathematicaSymbolic(const string& fileName);
    void	dumpMatrix(uint x);


    virtual	void	insertElement(uint x, uint y) {};
    virtual	void	reset() {};

    virtual	MatrixClass	myClass()		{return AbstractMatrix;};
    virtual	AbstractLargeSquareMatrix_sp copy( ) {SUBIMP();};
    virtual	uint	indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y) {SUBIMP();};
    virtual	bool	hasElement(uint x, uint y) { return this->indexFromCoordinatesOrUndefinedUnsignedInt(x,y)!=UndefinedUnsignedInt; };
    virtual	void	debug()	{};

    virtual	void	multiplyByVector(NVector_sp result, NVector_sp x);

    uint	indexFromCoordinates(uint x, uint y);

    DEFAULT_CTOR_DTOR(AbstractLargeSquareMatrix_O);
};












SMART(FullLargeSquareMatrix)
    class FullLargeSquareMatrix_O : public  AbstractLargeSquareMatrix_O 
{
    LISP_BASE1( AbstractLargeSquareMatrix_O );
    LISP_CLASS(chem,ChemPkg,FullLargeSquareMatrix_O,"FullLargeSquareMatrix");

private:
	uint		_ActiveElements;
	vector<double>	_Values;

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

virtual	MatrixClass	myClass()	{return FullMatrix;};
virtual	AbstractLargeSquareMatrix_sp copy();
virtual	uint	indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y);
virtual	void	insertElement(uint x, uint y) {};
virtual	void	reset() {};
virtual	void	setAtIndex(uint ii, double d) { this->_Values[ii] = d;};
virtual	double	getAtIndex(uint ii) { return this->_Values[ii];};
virtual	void	setElement(uint x, uint y, double d) { this->_Values[this->indexFromCoordinates(x,y)]=d;};
virtual	void	addToElement(uint x, uint y, double d) { this->_Values[this->indexFromCoordinates(x,y)]+=d;};
virtual	double	element(uint x, uint y) { return this->_Values[this->indexFromCoordinates(x,y)];};
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
    LISP_BASE1(AbstractLargeSquareMatrix_O);
    LISP_CLASS(chem,ChemPkg,SparseLargeSquareMatrix_O,"SparseLargeSquareMatrix");

private:
	bool		_InsertionIsComplete;
	uint		_RowStartEntries;
	vector<uint> 	_RowStarts;
	uint		_ActiveElements;
	uint		_ReservedElements;
	vector<uint> 	_ColumnForValue;
	vector<double> 	_Values;
public:
    static SparseLargeSquareMatrix_sp	create(uint dim, TriangleType type);
private:
	void	initializeStorage();
	void	expandStorage();
	void	releaseStorage();
	void	setup(uint dim, TriangleType type);

public:
virtual	uint	indexBegin()	{ return 0; };
virtual	uint	indexEnd() {return this->_ActiveElements;};
virtual	uint	indexAdvance(uint i) {return i+1;};
virtual	void	coordinatesFromIndex(uint i,uint& x, uint& y );
virtual	void	addToElement(uint x, uint y, double d);

virtual	void	setElement(uint x, uint y, double d);
virtual	double	element(uint x, uint y);


virtual	MatrixClass	myClass()	{return SparseMatrix;};
virtual	AbstractLargeSquareMatrix_sp copy();
virtual	uint	indexFromCoordinatesOrUndefinedUnsignedInt(uint x, uint y);
//virtual	bool	hasElement(uint x, uint y) { return this->indexFromCoordinates(x,y)>=0; };
virtual	void	insertElement(uint x, uint y);
virtual	void	reset();
virtual	void	insertionIsComplete() {this->_InsertionIsComplete = true;};
virtual	void	debug();
virtual	void	setAtIndex(uint ii, double d) { this->_Values[ii] = d;};
virtual	double	getAtIndex(uint ii) { return this->_Values[ii];};

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
	virtual ~SparseLargeSquareMatrix_O();

};

#if 0
inline	SparseLargeSquareMatrix_sp	new_SparseLargeSquareMatrix_sp(uint dim, TriangleType m) {
    return SparseLargeSquareMatrix_sp(new SparseLargeSquareMatrix_O(dim,m));
}

inline	SparseLargeSquareMatrix_sp new_SparseLargeSquareMatrix_sp(const SparseLargeSquareMatrix_O& m) {
    return SparseLargeSquareMatrix_sp(new SparseLargeSquareMatrix_O(m));
}
#endif






};


TRANSLATE(chem::AbstractLargeSquareMatrix_O);
TRANSLATE(chem::FullLargeSquareMatrix_O);
TRANSLATE(chem::SparseLargeSquareMatrix_O);
#endif          /* ifndef LARGE_MATRIX_H */
