/*
    File: coloredVertexList.cc
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
#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/geom/coloredVertexList.h>
//#include "core/archiveNode.h"
#include <cando/geom/ovector3.h>
#include <clasp/core/environment.h>
#include <cando/geom/color.h>
#ifdef RENDER
#include <geom/render.h>
#endif
#include <clasp/core/wrappers.h>



namespace geom 
{



#define ARGS_ColoredVertex_O_make "(pos color)"
#define DECL_ColoredVertex_O_make ""
#define DOCS_ColoredVertex_O_make "Create a colored vertex"
CL_LISPIFY_NAME(make-colored-vertex);
CL_DEFUN ColoredVertex_sp ColoredVertex_O::make(OVector3_sp opos, Color_sp color)
    {
	ColoredVertex_sp me = ColoredVertex_O::create();
	Vector3 pos = opos->get();
	me->setPosition(pos);
	me->setColor(color);
	return me;
    }






    
    
    void ColoredVertex_O::initialize()
    {
	this->Base::initialize();
	this->_Pos.set(0.0,0.0,0.0);
	this->_Color = _Nil<Color_O>();
    }


#ifdef XML_ARCHIVE
    void ColoredVertex_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->archivePlainObject<Vector3>("v","VPos",this->_Pos);
	node->archiveObject("c",this->_Color);
    }
#endif


    void ColoredVertex_O::transform(const Matrix& m)
    {
	this->_Pos = m.multiplyByVector3(this->_Pos);
    }






void	ColoredVertexList_O::initialize()
{
    this->Base::initialize();
    this->_Vertices.clear();
}


    ColoredVertexList_sp ColoredVertexList_O::create(core::Cons_sp elements, core::LispPtr env)
{
    uint sz = elements->length();
    ColoredVertexList_sp res = ColoredVertexList_O::create(sz);
    uint idx = 0;
    for ( ; elements.notnilp(); elements = elements->cdr() )
    {
	res->setElement(idx,core::oCar(elements).as<ColoredVertex_O>());
	idx++;
    }
    return res;
}

void ColoredVertexList_O::resize(int iElement)
{
    ColoredVertex_sp cv = _Nil<ColoredVertex_O>();
    this->_Vertices.resize(iElement,cv);
}

ColoredVertexList_O::iterator	ColoredVertexList_O::begin()
{
    return this->_Vertices.begin();
}

ColoredVertexList_O::iterator	ColoredVertexList_O::end()
{
    return this->_Vertices.end();
}

ColoredVertex_sp	ColoredVertexList_O::getElement(uint i)
{
    ASSERT_lessThan(i,this->_Vertices.size() );
    ColoredVertex_sp res = this->_Vertices[i];
    return res;
}

void	ColoredVertexList_O::setElement(uint i, ColoredVertex_sp pos)
{
    ASSERT_lessThan(i,this->_Vertices.size() );
    this->_Vertices[i] = pos;
}


/*! Copy constructor
 */



void	ColoredVertexList_O::transform( const Matrix& tm )
{
int  sz;
    sz = this->_Vertices.size();
    for ( int i=0; i<sz; i++ ) 
    {
	ColoredVertex_sp cv = this->_Vertices[i];
	cv->transform(tm);
    }
}



#ifdef XML_ARCHIVE
    void	ColoredVertexList_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->archiveVector0("vertices",this->_Vertices);
}
#endif







};
