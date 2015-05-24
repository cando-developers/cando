
#ifndef	Hold_H
#define	Hold_H


namespace chem {


template <class oClass>
class Hold {
public:
    gctools::smart_ptr<oClass>	_obj;
    Hold( ) { this->_obj = _Nil<oClass>(); };
    Hold( gctools::smart_ptr<oClass> o ) { this->_obj = o; };
	virtual ~Hold() { this->_obj.reset(); };
};


};
#endif
