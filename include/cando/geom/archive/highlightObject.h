/*
    File: highlightObject.h
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
#ifndef	HighlightObject_H //[
#define HighlightObject_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/highlightObject.fwd.h>

namespace geom
{

    SMART(PickableObject);
    class PickableObject_O : public core::CxxObject_O
    {
	LISP_CLASS(geom,GeomPkg,PickableObject_O,"PickableObject",core::CxxObject_O);
    public:
	void	initialize();
    public:
	virtual bool isPicked(uint id) {return false;};


	DEFAULT_CTOR_DTOR(PickableObject_O);
    };



    SMART(ControlObject);
    class ControlObject_O : public PickableObject_O
    {
	LISP_CLASS(geom,GeomPkg,ControlObject_O,"ControlObject",PickableObject_O);
    private:
    public:
	void	initialize();
    public:
	DEFAULT_CTOR_DTOR(ControlObject_O);
    };



    SMART(DragObject);
    class DragObject_O : public ControlObject_O
    {
	LISP_CLASS(geom,GeomPkg,DragObject_O,"DragObject",ControlObject_O);
    protected:
	core::Function_sp	_StartCode;
	core::Function_sp	_DragCode;
	core::Function_sp	_ReleaseCode;
    public:
	static DragObject_sp make(core::Function_sp startCode, core::Function_sp dragCode, core::Function_sp releaseCode); // 
    public:
	void	initialize();
    public:
	void setStartCode(core::Function_sp proc) { this->_StartCode = proc;};
	core::Function_sp getStartCode() { return this->_StartCode;};
	void setDragCode(core::Function_sp proc) { this->_DragCode = proc;};
	core::Function_sp getDragCode() { return this->_DragCode;};
	void setReleaseCode(core::Function_sp proc) { this->_ReleaseCode = proc;};
	core::Function_sp getReleaseCode() { return this->_ReleaseCode;};

	DEFAULT_CTOR_DTOR(DragObject_O);
    };




    FORWARD(HighlightedObject);
class HighlightedObject_O : public PickableObject_O
{
    LISP_CLASS(geom,GeomPkg,HighlightedObject_O,"HighlightedObject",PickableObject_O);
public: // virtual functions inherited from Object
    void	initialize();
//	string	__repr__() const;

private: // instance variables - not serialized
    uint     	_HighlightedId;

public:
    HighlightedObject_O( const HighlightedObject_O& ss ); //!< Copy constructor

    uint getHighlightedId() { return this->_HighlightedId;};
    void setHighlightedId(uint id);
    virtual bool isPicked(uint id);

    explicit HighlightedObject_O();
    virtual ~HighlightedObject_O();

};



class HighlightedHolder_O : public HighlightedObject_O
{
    LISP_CLASS(geom,GeomPkg,HighlightedHolder_O,"HighlightedHolder",HighlightedObject_O);
    DECLARE_MAKE_INIT();
public: // virtual functions inherited from Object
    void	initialize();
private:
    core::T_sp	_HeldObject;
public:
    string __repr__() const;
    DEFAULT_CTOR_DTOR(HighlightedHolder_O);
};




};



#endif //]
