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


    EXPOSE_CLASS(geom,ColoredVertex_O);

#define ARGS_ColoredVertex_O_make "(pos color)"
#define DECL_ColoredVertex_O_make ""
#define DOCS_ColoredVertex_O_make "Create a colored vertex"
CL_LISPIFY_NAME(make-colored-vertex);
CL_DEFUN ColoredVertex_sp ColoredVertex_O::make(OVector3_sp opos, Color_sp color)
    {_G();
	ColoredVertex_sp me = ColoredVertex_O::create();
	Vector3 pos = opos->get();
	me->setPosition(pos);
	me->setColor(color);
	return me;
    }


    void ColoredVertex_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ColoredVertex_O>()
	    .def("setPosition",&ColoredVertex_O::setPosition)
	    .def("getPosition",&ColoredVertex_O::getPosition_value)
	    .def("setColor",&ColoredVertex_O::setColor)
	    .def("getColor",&ColoredVertex_O::getColor)
	    ;
//	af_def(CorePkg,"make-colored-vertex",&ColoredVertex_O::make,ARGS_ColoredVertex_O_make,DECL_ColoredVertex_O_make,DOCS_ColoredVertex_O_make);
    }

    void ColoredVertex_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(CorePkg,ColoredVertex,"","",_lisp)
	    .def("setPosition",&ColoredVertex_O::setPosition)
	    .def("getPosition",&ColoredVertex_O::getPosition_value)
	    .def("setColor",&ColoredVertex_O::setColor)
	    .def("getColor",&ColoredVertex_O::getColor)
	    ;
#endif
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




    EXPOSE_CLASS(geom,ColoredVertexList_O);
    void ColoredVertexList_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ColoredVertexList_O>()
	.def("getElement",&ColoredVertexList_O::getElement)
	.def("setElement",&ColoredVertexList_O::setElement)
	;
    }

    void ColoredVertexList_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(CorePkg,ColoredVertexList,"","",_lisp)
	.def("getElement",&ColoredVertexList_O::getElement)
	.def("setElement",&ColoredVertexList_O::setElement)
	;
#endif
    }

void	ColoredVertexList_O::initialize()
{
    this->Base::initialize();
    this->_Vertices.clear();
}


    ColoredVertexList_sp ColoredVertexList_O::create(core::Cons_sp elements, core::Lisp_sp env)
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
{_G();
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
{_G();
    ASSERT_lessThan(i,this->_Vertices.size() );
    ColoredVertex_sp res = this->_Vertices[i];
    return res;
}

void	ColoredVertexList_O::setElement(uint i, ColoredVertex_sp pos)
{_G();
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
{_G();
    this->Base::archiveBase(node);
    node->archiveVector0("vertices",this->_Vertices);
}
#endif







};



