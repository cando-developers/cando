
#ifndef	Hold_H
#define	Hold_H


namespace chem {


template <class oClass>
class Hold {
public:
  gctools::smart_ptr<oClass>	_o;
  Hold( ) { this->_o = oClass::create(); };
  Hold( gctools::smart_ptr<oClass> o ) { this->_o = o; };
  gc::smart_ptr<oClass> value() { return this->_o; };
  virtual ~Hold() {};
};

 /*!Use a GC compatible allocator here */
 template <class T>
   Hold<T>* createHolder(gc::smart_ptr<T> obj) {
   return gc::RootClassAllocator<Hold<T>>::allocate(obj);// return new Hold<T>(obj);
 }

 template <class T>
   Hold<T>* createHolder() {
   return gc::RootClassAllocator<Hold<T>>::allocate();
//   return new Hold<T>();
 }

 /*! Use a corresponding GC compatible deallocator here */
template <class T>
  void deleteHolder(Hold<T>* objP) {
  gc::RootClassAllocator<Hold<T>>::deallocate(objP);
//  delete objP;
};

};
#endif
