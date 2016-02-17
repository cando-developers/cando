#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//
#include <clasp/core/common.h>
#include <cando/geom/coordinateArray.h>
//#include "core/archiveNode.h"
//#include "core/serializerNode.h"
#include <cando/geom/ovector3.h>
#include <cando/geom/omatrix.h>
#include <clasp/core/lispStream.h>
#ifdef RENDER
#include <cando/geom/color.h>
#include <geom/render.h>
#endif
#include <clasp/core/wrappers.h>



namespace geom 
{


#define ARGS_CoordinateArray_O_make "(fnsize &optional vals)"
#define DECL_CoordinateArray_O_make ""
#define DOCS_CoordinateArray_O_make "fnsize vals"
CL_LAMBDA(fnsize &optional vals);
CL_LISPIFY_NAME(make-coordinate-array);
CL_DEFUN CoordinateArray_sp CoordinateArray_O::make(gc::Nilable<core::Fixnum_sp> fnsize, core::List_sp vals)
    {
	if ( fnsize.notnilp() && vals.notnilp() )
	{
	    SIMPLE_ERROR(BF("You can only pass size or vals"));
	}
	GC_ALLOCATE(CoordinateArray_O,me);
	if ( fnsize.notnilp() )
	{
          core::Fixnum size = core::clasp_to_fixnum(fnsize);
	    me->resize(size);
	} else if ( vals.notnilp() ) {
	    me->fillFromList(vals);
	}
	return me;
    }





void CoordinateArray_O::fillFromList(core::List_sp vals)
{
  if ( this->size() != 0 ) {
    SIMPLE_ERROR(BF("Coordinate array must be empty to fillFromCons"));
  }
  if ( vals.notnilp() ) {
    for ( auto cur : vals ) {
      core::T_sp obj = core::oCar(cur);
      if ( core::cl__consp(obj) )
      {
        core::Cons_sp consVec = obj.as<core::Cons_O>();
        Vector3 vec;
        vec.fillFromCons(consVec);
        this->appendElement(vec);
      } else if ( OVector3_sp ovec = obj.asOrNull<OVector3_O>()) {
        this->appendElement(ovec->value());
      } else {
        SIMPLE_ERROR(BF("Cannot convert value[%s] to Vector3") % obj->__repr__() );
      }
    }
  }
}

    void CoordinateArray_O::initialize()
    {
	this->Base::initialize();
	this->_Points.clear();
    }


CoordinateArray_sp CoordinateArray_O::create(core::List_sp elements)
{
  size_t sz = core::cl__length(elements);
  CoordinateArray_sp res = CoordinateArray_O::create(sz);
  uint idx = 0;
  for ( auto cur : elements ) {
    res->setElement(idx,core::oCar(cur).as<OVector3_O>()->get());
    idx++;
  }
  return res;
}

void CoordinateArray_O::resize(int iElement)
{
  Vector3	zero;
  zero.set(0.0,0.0,0.0);
  this->_Points.resize(iElement,zero);
}

    CoordinateArray_O::iterator CoordinateArray_O::begin()
    {
	return this->_Points.begin();
    }

    CoordinateArray_O::iterator CoordinateArray_O::end()
    {
	return this->_Points.end();
    }

    Vector3& CoordinateArray_O::getElement(uint i)
    {
	ASSERT_lessThan(i,this->_Points.size() );
	Vector3 &res = this->_Points[i];
	return res;
    }

    core::T_sp CoordinateArray_O::getElementObject(uint i) const
    {
	ASSERT_lessThan(i,this->_Points.size() );
	Vector3 res = this->_Points[i];
	OVector3_sp ov = OVector3_O::create(res);
	return ov;
    }

    core::T_sp CoordinateArray_O::setElementObject(uint i,core::T_sp obj)
    {
	ASSERT_lessThan(i,this->_Points.size() );
	if ( obj->isAssignableTo<OVector3_O>() )
	{
	    OVector3_sp ov = obj.as<OVector3_O>();
	    this->_Points[i] = ov->get();
	} else
	{
	    SIMPLE_ERROR(BF("setElementObject for CoordinateArray expects an OVector3 - got a %s") % obj->__class()->classNameAsString() );
	}
	return obj;
    }


CL_LISPIFY_NAME("getElement");
CL_DEFMETHOD     Vector3 CoordinateArray_O::value_getElement(uint i)
    {
	ASSERT_lessThan(i,this->_Points.size() );
	Vector3 res = this->_Points[i];
	return res;
    }
CL_LISPIFY_NAME("setElement");
CL_DEFMETHOD     void CoordinateArray_O::setElement(uint i, Vector3 pos)
    {
	ASSERT_lessThan(i,this->_Points.size() );
	this->_Points[i] = pos;
    }


/*! Copy constructor
 */
    CoordinateArray_O::CoordinateArray_O( const CoordinateArray_O& c ) : Base(c)
    {
	int  sz;
	sz = c._Points.size();
	this->resize(sz);
	for ( int i=0; i<sz; i++ ) {
	    this->_Points[i] = c._Points[i];
	}
    }

    void CoordinateArray_O::zeroAllCoordinates()
    {
	int  sz;
	sz = this->_Points.size();
	for ( int i=0; i<sz; i++ ) {
	    this->_Points[i].set(0.0,0.0,0.0);
	}
    }

    void CoordinateArray_O::add(CoordinateArray_sp y)
    {_OF();
	int  sz;
	ASSERT_eq(this->_Points.size(),y->size());
	sz = this->_Points.size();
	for ( int i=0; i<sz; i++ ) {
	    this->_Points[i] = this->_Points[i].add(y->getElement(i));
	}
    }


    double CoordinateArray_O::distanceTo(CoordinateArray_sp y)
    {_OF();
	int  sz;
	double dist = 0.0;
	sz = this->_Points.size();
	ASSERT_eq(this->_Points.size(),y->size());
	ASSERT(this->_Points.size()>0);
	for ( int i=0; i<sz; i++ ) {
	    Vector3 v = this->_Points[i].sub(y->getElement(i));
	    double d = v.dotProduct(v);
	    dist += d;
	}
	dist /= (double)(sz);
	return sqrt(dist);
    }


CL_LISPIFY_NAME("multiplyByScalar");
CL_DEFMETHOD     void CoordinateArray_O::multiplyByScalar(double d)
    {
	int  sz;
	sz = this->_Points.size();
	for ( int i=0; i<sz; i++ ) 
	{
	    this->_Points[i] = this->_Points[i].multiplyByScalar(d);
	}
    }

#ifdef RENDER
    Render_sp CoordinateArray_O::asDisplayListOfPoints()
    {
	int  sz;
	DisplayList_sp dl = DisplayList_O::create();
	GrPoints_sp pl = GrPoints_O::create();
	sz = this->_Points.size();
	for ( int i=0; i<sz; i++ ) 
	{
	    pl->appendPoint(this->_Points[i],Color_O::white(_lisp));
	}
	dl->append(pl);
	return dl;
    }
#endif
CL_LISPIFY_NAME("transform");
CL_DEFMETHOD     void CoordinateArray_O::transform( const Matrix& tm )
    {
	int  sz;
	sz = this->_Points.size();
	for ( int i=0; i<sz; i++ ) {
	    this->_Points[i] = tm.multiplyByVector3(this->_Points[i]);
	}
    }

    void CoordinateArray_O::multiplyBy( double mul )
    {
	int  sz;
	sz = this->_Points.size();
	for ( int i=0; i<sz; i++ ) {
	    this->_Points[i] = this->_Points[i].multiplyByScalar(mul);
	}
    }

    string CoordinateArray_O::__repr__() const
    {
	return this->asString();
    }

CL_LISPIFY_NAME("asString");
CL_DEFMETHOD     string CoordinateArray_O::asString( ) const
    {
	int  sz;
	stringstream	ss;
	sz = this->_Points.size();
	for ( int i=0; i<sz; i++ ) {
	    ss << this->_Points[i].asString() << std::endl;
	}
	return ss.str();
    }


#ifdef OLD_SERIALIZE
    void CoordinateArray_O::serialize(serialize::SNode node)
    {_OF();
	node->orderedPOD(this->_Points);
    }
#endif

#if 0
void CoordinateArray_O::fields(core::Record_sp node)
{
  node->field(
	if ( node->saving() ) {
	    uint sz = this->size();
	    node->attribute("size",sz);
	    int i;
	    core::ArchiveP vnode;
	    vector<Vector3>::iterator vi;
	    for ( i=0,vi=this->_Points.begin(); vi!=this->_Points.end(); i++,vi++ ) {
		vnode = node->createChildNode("V");
		vnode->attribute("i",i);
		(*vi).archive(vnode);
	    }
	} else {
	    uint sz=0, i=0;
	    node->attribute("size",sz);
	    this->resize(sz);
	    vector<core::ArchiveP>::iterator	vi;
	    i = 0;
	    for ( vi=node->begin_Children(); vi!=node->end_Children(); vi++ ) {
		(*vi)->setRecognized(true);
		(*vi)->attribute("i",i);
		this->_Points[i].archive(*vi);
		i++;
	    }
	}
    }
#endif


    CoordinateArray_sp CoordinateArray_O::copy() const
    {
	GC_COPY(CoordinateArray_O,na,*this);
	return na;
    }


CL_LISPIFY_NAME("write-to-stream");
CL_DEFMETHOD     void CoordinateArray_O::writeToStream(string const& info, core::T_sp sout)
    {_OF();
        core::clasp_write_string((BF("+++CoordinateArray %d\n") % this->_Points.size() ).str(),sout);
        core::clasp_write_string((BF("+ %s\n") % info).str(),sout);
	for (int i=0; i<(int)this->_Points.size(); i++ )
	{
	    this->_Points[i].write(sout);
            core::clasp_terpri(sout);
	}
    }

#if 0
    string CoordinateArray_O::parseFromStream(core::Stream_sp sin)
    {_OF();
	streampos spos;
	string line, info;
	int numLines;
	bool hitEof;
	spos = sin->tell();
	sin->readLine(line,hitEof);
	if ( !hitEof )
	{
	    SIMPLE_ERROR(BF("Could not read header for CoordinateArray - got[%s] eof[%d] ")
			       % line % hitEof );
	}
	if ( line.substr(0,18)=="+++CoordinateArray")
	{
	    stringstream snum;
	    snum.str(line.substr(18,999));
	    snum >> numLines;
	    ASSERTF(numLines>=0,BF("Illegal number of entries[%d]") % numLines);
	    this->resize(numLines);
	    sin->readLine(line,hitEof);
//	    getline(sin->stream(),line);
	    if ( !hitEof )
	    {
		SIMPLE_ERROR(BF("Could not read info - got[%s]") % line );
	    }
	    info = line.substr(2,9999);
	    for ( int i=0; i<numLines; i++ )
	    {
		Vector3 pos;
		sin->readLine(line,hitEof);
//		getline(sin->stream(),line);
		if ( !hitEof )
		{
		    SIMPLE_ERROR(BF("Error reading relative line[%d] contents[%s] for CoordinateArray") % i % line );
		}
		pos.parseFromString(line);
		this->getElement(i)=pos;
	    }
	} else
	{
	    sin->seek(spos);
	    SIMPLE_ERROR(BF("Could not read CoordinateArray - header was[%s]") % line );
	}
	return info;
    }
#endif






    CoordinateArrayWithHash_sp CoordinateArrayWithHash_O::create(uint sz)
    {
	GC_ALLOCATE(CoordinateArrayWithHash_O,h);
	h->resize(sz);
	return h;
    }







    void CoordinateArrayWithHash_O::initialize()
    {
	this->Base::initialize();
	this->_Hash = UndefinedUnsignedInt;
    }

#ifdef XML_ARCHIVE
    void CoordinateArrayWithHash_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute("hash",this->_Hash);
    }
#endif


#ifdef OLD_SERIALIZE
    void CoordinateArrayWithHash_O::serialize(serialize::SNode node)
    {
	this->Base::serialize(node);
	node->attribute("hash",this->_Hash);
    }
#endif


core::List_sp CoordinateArray_O::encode() {
  core::Vector_sp v = core::core__make_vector(cl::_sym_DoubleFloat_O,3*this->_Points.size());
  size_t cur(0);
  for ( size_t i(0); i<this->_Points.size(); ++i ) {
    (*v)[cur++] = core::clasp_make_double_float(this->_Points[i].getX());
    (*v)[cur++] = core::clasp_make_double_float(this->_Points[i].getY());
    (*v)[cur++] = core::clasp_make_double_float(this->_Points[i].getZ());
  }
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,pnts),core::clasp_make_fixnum(this->_Points.size())),
                                  core::Cons_O::create(INTERN_(kw,pnts),v));
}

void CoordinateArray_O::decode(core::List_sp c) {
  core::Fixnum_sp fn_size = gc::As<core::Fixnum_sp>(oCdr(oFirst(c)));
  core::Vector_sp v = gc::As<core::Vector_sp>(oCdr(oSecond(c)));
  this->resize(core::clasp_to_fixnum(fn_size));
  size_t cur(0);
  for ( size_t i(0); i<this->_Points.size(); ++i ) {
    // TODO: Optimize this if v is a specialized vector of doubles
    double x = core::clasp_to_double((*v)[cur++]);
    double y = core::clasp_to_double((*v)[cur++]);
    double z = core::clasp_to_double((*v)[cur++]);
    Vector3 v3(x,y,z);
    this->_Points[i] = v3;
  }
}



};



