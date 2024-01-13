/*
    File: coloredVertexList.h
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
 *      File:   coordinateArray.h
 *
 *      Description:
 *              Routines to maintain a vector of Vector3's
 */

#ifndef ColoredVertexList_H
#define ColoredVertexList_H

#include <vector>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/vector3.h>
#include <cando/geom/ovector3.fwd.h>
#include <cando/geom/matrix.h>


namespace geom 
{

    SMART(Color);

    SMART(DisplayList);

    SMART(ColoredVertex);
    class ColoredVertex_O : public core::CxxObject_O
    {
	LISP_CLASS(geom,GeomPkg,ColoredVertex_O,"ColoredVertex",core::CxxObject_O);
	DECLARE_INIT();
    public:
	void initialize();
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
    class ColoredVertexList_O : public core::CxxObject_O
    {
	LISP_CLASS(geom,GeomPkg,ColoredVertexList_O,"ColoredVertexList",core::CxxObject_O);
    public:
	void initialize();
    public:
    private:
        gctools::Vec0<ColoredVertex_sp>	_Vertices;
    public:
	static ColoredVertexList_sp create(uint size)
	{
          auto a = gctools::GC<ColoredVertexList_O>::allocate_with_default_constructor();
	    a->resize(size);
	    return a;
	}
	static ColoredVertexList_sp create(core::Cons_sp vertexCons, core::LispPtr e);
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
#endif          /* ifndef ColoredVertexList_H */
