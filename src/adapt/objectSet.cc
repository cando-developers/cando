/*
    File: objectSet.cc
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

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/environment.h>
#include <clasp/core/multipleValues.h>
#include <cando/adapt/objectSet.h>
#include <clasp/core/wrappers.h>

namespace adapt {

void ObjectSet_O::initialize() {
  this->Base::initialize();
  this->_Set = HashTableEq_O::create_default();
}

CL_LISPIFY_NAME("addObjects");
CL_DEFMETHOD void ObjectSet_O::addObjects(ObjectSet_sp other) {
  
  other->map([this](T_sp obj) {this->insert(obj); });
}

CL_LISPIFY_NAME("asCons");
CL_DEFMETHOD List_sp ObjectSet_O::asCons() {
  List_sp res = nil<T_O>();
  this->map([&res](T_sp o) {
	    res = Cons_O::create(o,res);
  });
  return res;
}

CL_LISPIFY_NAME("adapt:objectSetUnion");
CL_DEFMETHOD ObjectSet_sp ObjectSet_O::setUnion(ObjectSet_sp other) {
  ObjectSet_sp os = ObjectSet_O::create();
  os->addObjects(this->sharedThis<ObjectSet_O>());
  os->addObjects(other);
  return os;
}

CL_LISPIFY_NAME("adapt:objectSetIntersection");
CL_DEFMETHOD ObjectSet_sp ObjectSet_O::intersection(ObjectSet_sp b) {
  
  ObjectSet_sp nset;
  nset = ObjectSet_O::create();

  this->map([&b, &nset](T_sp o) {
            if ( b->contains(o) ) {
                LOG("Found it!!!" );
                nset->insert(o);
            } else {
                LOG("Not found" );
            }
  });
  return nset;
}

CL_LISPIFY_NAME("relativeComplement");
CL_DEFMETHOD ObjectSet_sp ObjectSet_O::relativeComplement(ObjectSet_sp b) {
  ObjectSet_sp nset;
  nset = ObjectSet_O::create();
  this->map([&b, &nset](T_sp o) {
	if ( !b->contains(o) ) {
	    nset->insert(o);
	}
  });
  return nset;
}

string ObjectSet_O::asString() const {
  stringstream ss;
  this->map([&ss](T_sp si) {
                ss << _rep_(si) << " ";
  });
  return ss.str();
}

/*! Return a new set that takes every element of (this) in combination
	with every element in b separated by a comma
*/
ObjectSet_sp ObjectSet_O::cartesianProduct(ObjectSet_sp b) {
  
  ObjectSet_sp nset;
  stringstream sstr;
  nset = ObjectSet_O::create();
  this->map([&b, &nset, this](T_sp si) {
            this->map( [&si,&nset] (T_sp bi) {
                    Cons_sp op = _lisp->create<Cons_O>(si,bi);
                    nset->insert(op);
                });
  });
  return nset;
}

/*! Return a new set that takes every element of (this) in combination
	with every element in b separated by a comma
*/
ObjectSet_sp ObjectSet_O::cartesianProductWrapped(ObjectSet_sp b, const ObjectSetCartesianProductWrapper &wrapper) {
  
  ObjectSet_sp nset;
  stringstream sstr;
  nset = ObjectSet_O::create();

  this->map([&b, &nset, this, &wrapper](T_sp si) {
            this->map( [&si,&nset,&wrapper] (T_sp bi) {
                    T_sp op = wrapper(si,bi);
                    nset->insert(op);
                });
  });
  return nset;
}

void ObjectSet_O::map(std::function<void(T_sp)> const &fn) {
  this->_Set->mapHash([&fn](T_sp key, T_sp val) {
                fn(key);
  });
}

void ObjectSet_O::map(std::function<void(T_sp)> const &fn) const {
  this->_Set->mapHash([&fn](T_sp key, T_sp val) {
                fn(key);
  });
}

CL_LISPIFY_NAME("addObjectsInCons");
CL_DEFMETHOD void ObjectSet_O::addObjectsInCons(List_sp c) {
  while (c.consp()) {
    this->insert(oCar(c));
    c = oCdr(c);
  }
}

#if defined(XML_ARCHIVE)
void ObjectSet_O::archive(ArchiveP node) {
  
  stringstream suid;
  if (node->saving()) {
    
    if (this->_Set.size() != 0) {
      int i = 0;
      stringstream suid;
      T_sp obj;
      set<gctools::smart_ptr<T_O>>::iterator oi;
      for (oi = this->_Set.begin(); oi != this->_Set.end(); i++, oi++) {
        obj = (*oi);
        suid.str("");
        suid << i;
        node->archiveObject(suid.str(), obj);
      }
    }
  } else {
    
    VectorNodes::iterator ci;
    T_sp object;
    this->_Set.clear();
    for (ci = node->begin_Children(); ci != node->end_Children(); ci++) {
      object = node->getArchive()->loadObjectDirectly((*ci));
      ASSERTNOTNULL(object);
      this->_Set.insert(object);
    }
  }
}
#endif // defined(XML_ARCHIVE)



};
