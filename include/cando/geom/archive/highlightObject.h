#ifndef	HighlightObject_H //[
#define HighlightObject_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/executables.fwd.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/highlightObject.fwd.h>

namespace geom
{

    SMART(PickableObject);
    class PickableObject_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(geom,GeomPkg,PickableObject_O,"PickableObject");
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
	virtual bool isPicked(uint id) {return false;};


	DEFAULT_CTOR_DTOR(PickableObject_O);
    };



    SMART(ControlObject);
    class ControlObject_O : public PickableObject_O
    {
	LISP_BASE1(PickableObject_O);
	LISP_CLASS(geom,GeomPkg,ControlObject_O,"ControlObject");
    private:
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
	DEFAULT_CTOR_DTOR(ControlObject_O);
    };



    SMART(DragObject);
    class DragObject_O : public ControlObject_O
    {
	LISP_BASE1(ControlObject_O);
	LISP_CLASS(geom,GeomPkg,DragObject_O,"DragObject");
    protected:
	core::Function_sp	_StartCode;
	core::Function_sp	_DragCode;
	core::Function_sp	_ReleaseCode;
    public:
	static DragObject_sp make(core::Function_sp startCode, core::Function_sp dragCode, core::Function_sp releaseCode); // 
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
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
    LISP_BASE1(PickableObject_O);
    LISP_CLASS(geom,GeomPkg,HighlightedObject_O,"HighlightedObject");
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);
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
    LISP_BASE1(HighlightedObject_O);
    LISP_CLASS(geom,GeomPkg,HighlightedHolder_O,"HighlightedHolder");
    DECLARE_MAKE_INIT();
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);
private:
    core::T_sp	_HeldObject;
public:
    string __repr__() const;
    DEFAULT_CTOR_DTOR(HighlightedHolder_O);
};




};


TRANSLATE(geom::PickableObject_O);
TRANSLATE(geom::ControlObject_O);
TRANSLATE(geom::DragObject_O);
TRANSLATE(geom::HighlightedObject_O);
TRANSLATE(geom::HighlightedHolder_O);

#endif //]