//
// (C) 2004 Christian E. Schafmeister
//


/*
 *      File:   coordinateArray.h
 *
 *      Description:
 *              Routines to maintain a vector of Vector3's
 */

#ifndef ColoredVertexList_H
#define ColoredVertexList_H

#include <vector>
#include "core/foundation.h"
#include "core/object.h"
#include "candoBasePackage.fwd.h"
#include "vector3.h"
#include "ovector3.fwd.h"
#include "matrix.h"


namespace candoBase 
{

    SMART(Color);

    SMART(DisplayList);

    SMART(ColoredVertex);
    class ColoredVertex_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(candoBase,CandoBasePkg,ColoredVertex_O,"ColoredVertex");
	DECLARE_INIT();
    public:
	void initialize();
//	void archiveBase(core::ArchiveP node);
    private:
	Vector3	_Pos;
	Color_sp	_Color;
    public:
	static ColoredVertex_sp make(OVector3_sp pos, Color_sp color);
    public:
	void setPosition(const Vector3& pos) { this->_Pos = pos;};
	Vector3	getPosition_value() {return this->_Pos;};
	Vector3& getPosition() { return this->_Pos;};
	void setColor(Color_sp c) { this->_Color = c;};
	Color_sp getColor() { return this->_Color;};
	void transform(const Matrix& m);
	DEFAULT_CTOR_DTOR(ColoredVertex_O);
    };




    SMART(ColoredVertexList);
    class ColoredVertexList_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(candoBase,CandoBasePkg,ColoredVertexList_O,"ColoredVertexList");
    public:
	void initialize();
    public:
//	void	archiveBase(core::ArchiveP node);
    private:
        gctools::Vec0<ColoredVertex_sp>	_Vertices;
    public:
	static ColoredVertexList_sp create(uint size)
	{
	    GC_ALLOCATE(ColoredVertexList_O,a);
	    a->resize(size);
	    return a;
	}
	static ColoredVertexList_sp create(core::Cons_sp vertexCons, core::Lisp_sp e);
    public:
	typedef gctools::Vec0<ColoredVertex_sp>::iterator	iterator;
    public:

	ColoredVertexList_O( const ColoredVertexList_O& orig );


	iterator	begin();
	iterator	end();
	uint	size() const { return this->_Vertices.size(); };
	ColoredVertex_sp getElement(uint i);
	void	setElement(uint i, ColoredVertex_sp pos );
	void	append(ColoredVertex_sp pos) { this->_Vertices.push_back(pos); };
	void	resize(int i);
	void	clear() { this->_Vertices.clear();};
	void	transform(const Matrix& transMatrix);
	DEFAULT_CTOR_DTOR(ColoredVertexList_O);
    };


};
TRANSLATE(candoBase::ColoredVertex_O);
TRANSLATE(candoBase::ColoredVertexList_O);
#endif          /* ifndef ColoredVertexList_H */
