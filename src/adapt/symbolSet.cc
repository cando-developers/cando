/*
    File: symbolSet.cc
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
#define DEBUG_LEVEL_FULL

//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/adapt/symbolSet.h>
#include <clasp/core/multipleValues.h>
#include <cando/adapt/stringList.h>
#include <cando/adapt/objectSet.h>
#include <clasp/core/environment.h>
#include <clasp/core/cons.h>
#include <cando/adapt/symbolList.h>
#include <clasp/core/numbers.h>
#include <clasp/core/wrappers.h>

namespace adapt {

extern Symbol_sp _sym_entries;

/*
  __BEGIN_DOC(classes.SymbolSet.!class.SymbolSet)
  \requiredKeyed{entries:}{Cons::listOfStrings}

  Create a SymbolSet containing the strings in \sa{listOfStrings}.
  __END_DOC
*/
SymbolSet_sp SymbolSet_O::make(List_sp entries) {
  
  SymbolSet_sp me = SymbolSet_O::create();
  me->insertConsSymbols(entries);
  return me;
}

CL_LISPIFY_NAME("insertConsSymbols");
CL_DEFMETHOD SymbolSet_sp SymbolSet_O::insertConsSymbols(List_sp vals) {
  for (auto p : vals) {
    Symbol_sp t = gc::As<Symbol_sp>(oCar(p));
    this->insert(t);
  }
  return this->sharedThis<SymbolSet_O>();
}

CL_LISPIFY_NAME("insertSymbolSet");
CL_DEFMETHOD void SymbolSet_O::insertSymbolSet(SymbolSet_sp s) {
  s->map([this](Symbol_sp k) {this->insert(k); });
}

void SymbolSet_O::insertSymbolList(SymbolList_sp slist) {
  slist->map([this](Symbol_sp e) {this->insert(e); });
}

SymbolSet_sp SymbolSet_O::copy() {
  
  GC_COPY(SymbolSet_O, snew, *this);
  snew->clear();
  this->map([&snew](Symbol_sp k) {
                snew->insert(k);
  });
  return snew;
}

string SymbolSet_O::__repr__() const {
  stringstream ss;
  string nm;
  ss.str("");
  ss << "(" << this->_instanceClass()->classNameAsString() << " ";
  this->map([&ss](Symbol_sp s) {
                ss << _rep_(s) << " ";
  });
  ss << " )";
  return ss.str();
}
CL_LISPIFY_NAME("asString");
CL_DEFMETHOD string SymbolSet_O::asString() {
  stringstream ss;
  string nm;
  ss.str("");
  this->map([&ss](Symbol_sp s) {
                ss << _rep_(s) << " ";
  });
  return ss.str();
}

void SymbolSet_O::fields(core::Record_sp node) {
  // this->Base::fields(node);
  switch (node->stage()) {
  case Record_O::initializing: 
  case Record_O::loading: {
    Vector_sp valueVec;
    node->field(INTERN_(kw,data), valueVec );
    this->_Symbols->clrhash();
    for ( size_t i(0), iEnd(cl__length(valueVec)); i<iEnd; ++i ) {
      T_sp key = valueVec->elt(i);
      this->_Symbols->hash_table_setf_gethash(key,_lisp->_true());
    }
  }
      break;
  case Record_O::saving: {
    Vector_sp valueVec = core__make_vector(cl::_sym_T_O, this->_Symbols->hashTableCount() );
    size_t idx = 0;
    this->_Symbols->mapHash([&idx,&valueVec](T_sp key, T_sp val) {
        valueVec->setf_elt(idx++,key);
      });
    node->field(INTERN_(kw,data), valueVec);
  }
      break;
  case Record_O::patching: {
    IMPLEMENT_MEF(BF("Do patching"));
  }
      break;
  }
}


void SymbolSet_O::insertVectorStrings(const VectorStrings &symbolsAsStrings) {
  for (VectorStrings::const_iterator it = symbolsAsStrings.begin(); it != symbolsAsStrings.end(); it++) {
    Symbol_sp sym = _lisp->intern(*it);
    this->insert(sym);
  }
}

void SymbolSet_O::insertVector(Vector_sp vec) {
  for (int i = 0, iEnd(vec->length()); i < iEnd; ++i) {
    Symbol_sp sym = gc::As<Symbol_sp>(vec->elt(i));
    this->insert(sym);
  }
}

CL_LISPIFY_NAME("contains");
CL_DEFMETHOD bool SymbolSet_O::contains(Symbol_sp s) {
  
  return this->_Symbols->contains(s);
}

CL_LISPIFY_NAME("containsSubset");
CL_DEFMETHOD bool SymbolSet_O::containsSubset(SymbolSet_sp sub) {
  bool missed(false);
  this->map([&missed, this](Symbol_sp s) {
                if ( !this->contains(s) ) {
                    missed = true;
                    return;
                }
  });
  return !missed;
}

CL_LISPIFY_NAME("clear");
CL_DEFMETHOD void SymbolSet_O::clear() {
  _OF();
  this->_Symbols->clrhash();
}

/*! Check if the SymbolSet contains the same strings as this
 */
bool SymbolSet_O::equal(T_sp obj) const {
  if (this->eq(obj))
    return true;
  if (gc::IsA<SymbolSet_sp>(obj)) {
    SymbolSet_sp ss = gc::As<SymbolSet_sp>(obj);
    if (this->size() != ss->size())
      return false;
    bool missed(false);
    this->map([&missed, &ss](Symbol_sp s) {
                    if ( !ss->contains(s) ) 
                    {
                        missed = true;
                        return;
                    }
    });
    return !missed;
  }
  return false;
}

void SymbolSet_O::remove(Symbol_sp s) {
  _OF();
#ifdef DEBUG_ON
  if (this->strs.count(s) != 1) {
    SIMPLE_ERROR(BF("The string: %s was not found in SymbolSet") % s);
  }
#endif
  this->_Symbols->remhash(s);
}

void SymbolSet_O::initialize() {
  this->Base::initialize();
  this->_Symbols = HashTableEq_O::create_default();
}

SymbolSet_O::SymbolSet_O(const SymbolSet_O &ss) : Base(ss) {
  this->_Symbols = ss._Symbols;
}

CL_LISPIFY_NAME("adapt:SymbolSetUnion");
CL_DEFMETHOD SymbolSet_sp SymbolSet_O::setUnion(SymbolSet_sp b) {
  
  SymbolSet_sp nset;
  nset = SymbolSet_O::create();
  this->map([&nset](Symbol_sp s) { nset->insert(s); });
  b->map([&nset](Symbol_sp s) { nset->insert(s); });
  return nset;
}

CL_LISPIFY_NAME("adapt:SymbolSetIntersection");
CL_DEFMETHOD SymbolSet_sp SymbolSet_O::intersection(SymbolSet_sp b) {
  
  SymbolSet_sp nset = SymbolSet_O::create();
  this->map([&b, &nset](Symbol_sp s) {
            LOG(BF("Looking for(%s)") % _rep_((*si)) );
            if ( b->contains(s) )
            {
                LOG(BF("Found it!!!") );
                nset->insert(s);
            } else {
                LOG(BF("Not found") );
            }
  });
  return nset;
}

CL_LISPIFY_NAME("relativeComplement");
CL_DEFMETHOD SymbolSet_sp SymbolSet_O::relativeComplement(SymbolSet_sp b) {
  SymbolSet_sp nset;
  nset = SymbolSet_O::create();
  this->map([&nset, &b](Symbol_sp si) {
                if ( !b->contains(si) )
                {
                    nset->insert(si);
                }
  });
  return nset;
}

/*! Return a new set that takes every element of (this) in combination
  with every element in b separated by a comma
*/
CL_LISPIFY_NAME("cartesianProduct");
CL_DEFMETHOD ObjectSet_sp SymbolSet_O::cartesianProduct(SymbolSet_sp b) {
  
  ObjectSet_sp nset = ObjectSet_O::create();
  this->map([&b, &nset](Symbol_sp si) {
                b->map( [&si,&nset] (Symbol_sp bi) {
                        Cons_sp p = Cons_O::create(si,bi);
                        nset->insert(p);
                    } );
  });
  return nset;
}

/*! Return a new set that takes every element of (this) in combination
  with every element in b separated by a comma
*/
ObjectSet_sp SymbolSet_O::cartesianProductWrapped(SymbolSet_sp b, const SymbolSetCartesianProductWrapper &wrapper) {
  
  ObjectSet_sp nset;
  nset = ObjectSet_O::create();
  this->map([&b, &nset, &wrapper](Symbol_sp si) {
                b->map( [&si,&nset,&wrapper] (Symbol_sp bi) {
                        T_sp p = wrapper(si,bi);
                        nset->insert(p);
                    });
  });
  return nset;
}

Symbol_mv SymbolSet_O::first() {
  T_sp found = _Nil<T_O>();
  Symbol_sp result = _Nil<Symbol_O>();
  this->_Symbols->map_while_true([&found, &result](T_sp key, T_sp val) -> bool {
                result = key;
                found = _lisp->_true();
                return false;
  });
  return Values(result, found);
}

CL_LISPIFY_NAME("asCons");
CL_DEFMETHOD List_sp SymbolSet_O::asCons() {
  
  List_sp cur = _Nil<List_V>();
  this->map([&cur](Symbol_sp si) {
                cur = Cons_O::create(si,cur);
  });
  return cur;
}

void SymbolSet_O::map(std::function<void(Symbol_sp)> const &fn) {
  this->_Symbols->mapHash([&fn](T_sp key, T_sp val) {
                fn(gc::As<Symbol_sp>(key));
  });
}

void SymbolSet_O::map(std::function<void(Symbol_sp)> const &fn) const {
  this->_Symbols->mapHash([&fn](T_sp key, T_sp val) {
                fn(gc::As<Symbol_sp>(key));
  });
}


void SymbolSet_O::map_while_true(std::function<bool(Symbol_sp)> const &fn) {
  this->_Symbols->map_while_true([&fn](T_sp key, T_sp val) -> bool {
      return fn(gc::As<Symbol_sp>(key));
  });
}

void SymbolSet_O::map_while_true(std::function<bool(Symbol_sp)> const &fn) const {
  this->_Symbols->map_while_true([&fn](T_sp key, T_sp val) -> bool {
      return fn(gc::As<Symbol_sp>(key));
  });
}




};
