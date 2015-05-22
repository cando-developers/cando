#define	DEBUG_LEVEL_NONE


#include <clasp/core/common.h>
#include <cando/candoBase/omatrix.h>
#include <cando/candoBase/vector3.h>
#include <cando/candoBase/ovector3.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>


namespace candoBase {


#ifdef XML_ARCHIVE
    void OMatrix_O::archiveBase(core::ArchiveP node)
    {
	node->archivePlainObject<Matrix>( "matrix",
					  "Matrix", this->_Value);
    }
#endif


    string	OMatrix_O::__repr__() const
    {
	stringstream ss;
	ss << this->_Value.asString();
	return ss.str();
    }
    string	OMatrix_O::__str__()
    {
	return this->__repr__();
    }


    void	OMatrix_O::setFromString(const string& s)
    {_G();
	this->_Value.setFromString(s,_lisp);
    }

#if 0
    void	OMatrix_O::setFromStringFast(const string& s)
    {_G();
	this->_Value.setFromStringFast(s,_lisp);
    }
#endif

    void	OMatrix_O::setAll(const Matrix& m)
    {_G();
	this->_Value.setAll(m);
    }

//string OMatrix_O::asStringFast() { return this->_Value.asStringFast(); }


    string OMatrix_O::asString() { return this->_Value.asString(); }
    void OMatrix_O::setToIdentity() { this->_Value.setToIdentity(); }
    /*! Set the matrix element index as the value for testing */
    void OMatrix_O::setToIndexAsValue() { this->_Value.setToIndexAsValue();};
    void OMatrix_O::translate(Vector3 pos) { this->_Value.translate(&pos); }
    void OMatrix_O::rotationX(double radians) { this->_Value.rotationX(radians); }
    void OMatrix_O::rotationY(double radians) { this->_Value.rotationY(radians); }
    void OMatrix_O::rotationZ(double radians) { this->_Value.rotationZ(radians); }


    OMatrix_sp OMatrix_O::clone() const
    {_G();
	GC_COPY(OMatrix_O,clone,*this);
	return clone;
    }

    Matrix OMatrix_O::mul(const Matrix& other) const
    {_G();
	Matrix res = this->_Value*other;
	return res;
    };

    Vector3 OMatrix_O::mul(const Vector3& vec) const
    {_G();
	Vector3 res = this->_Value*vec;
	return res;
    }




    void OMatrix_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<OMatrix_O>()
	    .def("setToIdentity",&OMatrix_O::setToIdentity)
	    .def("setToIndexAsValue",&OMatrix_O::setToIndexAsValue)
	    .def("asString",&OMatrix_O::asString)
//	.def("asStringFast",&OMatrix_O::asStringFast)
	    .def("setFromString",&OMatrix_O::setFromString)
//	.def("setFromStringFast",&OMatrix_O::setFromStringFast)
	    .def("translate",&OMatrix_O::translate)
	    .def("rotationX",&OMatrix_O::rotationX)
	    .def("rotationY",&OMatrix_O::rotationY)
	    .def("rotationZ",&OMatrix_O::rotationZ)
	    .def("atRowColPut",&OMatrix_O::atRowColPut)
	    .def("atRowColGet",&OMatrix_O::atRowColGet)
	    .def("extract3x3",&OMatrix_O::extract3x3)
	    .def("transposed3x3",&OMatrix_O::transposed3x3)
	    .def("is3x3Orthogonal",&OMatrix_O::is3x3Orthogonal)
	    ;
    }

    void OMatrix_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(CorePkg,OMatrix,"","",_lisp)
	    .def("setToIdentity",&OMatrix_O::setToIdentity)
	    .def("setToIndexAsValue",&OMatrix_O::setToIndexAsValue)
	    .def("asString",&OMatrix_O::asString)
//	.def("asStringFast",&OMatrix_O::asStringFast)
	    .def("setFromString",&OMatrix_O::setFromString)
//	.def("setFromStringFast",&OMatrix_O::setFromStringFast)
	    .def("translate",&OMatrix_O::translate)
	    .def("rotationX",&OMatrix_O::rotationX)
	    .def("rotationY",&OMatrix_O::rotationY)
	    .def("rotationZ",&OMatrix_O::rotationZ)
	    .def("atRowColPut",&OMatrix_O::atRowColPut)
	    .def("atRowColGet",&OMatrix_O::atRowColGet)
	    .def("extract3x3",&OMatrix_O::extract3x3)
	    .def("transposed3x3",&OMatrix_O::transposed3x3)
	    .def("is3x3Orthogonal",&OMatrix_O::is3x3Orthogonal)
	    ;
#endif
    }

    EXPOSE_CLASS(candoBase,OMatrix_O);

};

