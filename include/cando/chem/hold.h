
#ifndef	Hold_H
#define	Hold_H

#include <cando/adapt/indexedObjectBag.h>


namespace chem {
  extern core::Symbol_sp& _sym_STARparserNodeHolderSTAR;

  template <class oClass>
    class Hold {
  public:
    core::Fixnum nodeHolderIndex;
    Hold( ) {
      gc::smart_ptr<oClass> obj = oClass::create();
      adapt::IndexedObjectBag_sp bag = gc::As<adapt::IndexedObjectBag_sp>(chem::_sym_STARparserNodeHolderSTAR->symbolValue());
      core::Fixnum index = bag->set_next_available_entry(obj);
      this->nodeHolderIndex = index;
    };
    Hold( gctools::smart_ptr<oClass> obj ) {
      adapt::IndexedObjectBag_sp bag = gc::As<adapt::IndexedObjectBag_sp>(chem::_sym_STARparserNodeHolderSTAR->symbolValue());
      core::Fixnum index = bag->set_next_available_entry(obj);
      this->nodeHolderIndex = index;
    };
    gc::smart_ptr<oClass> value() {
      adapt::IndexedObjectBag_sp bag = gc::As<adapt::IndexedObjectBag_sp>(chem::_sym_STARparserNodeHolderSTAR->symbolValue());
      core::T_sp obj = bag->entry(this->nodeHolderIndex);
      return obj;
    }
    virtual ~Hold() {
      adapt::IndexedObjectBag_sp bag = gc::As<adapt::IndexedObjectBag_sp>(chem::_sym_STARparserNodeHolderSTAR->symbolValue());
      bag->remove_entry(this->nodeHolderIndex);
    };
  };

 /*!Use a GC compatible allocator here */
  template <class T>
    Hold<T>* createHolder(gc::smart_ptr<T> obj) {
    return new Hold<T>(obj);
  }

  template <class T>
    Hold<T>* createHolder() {
    return new Hold<T>();
  }

 /*! Use a corresponding GC compatible deallocator here */
  template <class T>
    void deleteHolder(Hold<T>* objP) {
    delete objP;
  };

};
#endif
