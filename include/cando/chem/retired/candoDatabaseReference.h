#ifndef	CandoDatabaseReference_H //[
#define CandoDatabaseReference_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/entity.h>

#include <cando/chem/chemPackage.h>

//#include "core/objRef.fwd.h"// candoDatabaseReference.h wants ObjRef needs objRef.fwd.h

namespace chem {

SMART(ObjRef);

SMART(CandoDatabaseReference);
class CandoDatabaseReference_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,CandoDatabaseReference_O,"CandoDatabaseReference",core::CxxObject_O);
public:
//    void	archive(core::ArchiveP node);
	void	initialize();
private:
	core::ObjRef_sp	_Reference;
public:
    static CandoDatabaseReference_sp create(core::Lisp_sp e,CandoDatabase_sp bdb, const string& link);
public:


		/*!Lookup and return the object in the database referenced by this one
		 * and cast it to the desired type
		 */
	template <class otype>
          gctools::smart_ptr<otype> get() { return this->oget().as<otype>();};

		//!Lookup and return the object in the database referenced by this one
	core::T_sp oget();


	CandoDatabaseReference_O( const CandoDatabaseReference_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CandoDatabaseReference_O);
};


// Constructor function

};
TRANSLATE(chem::CandoDatabaseReference_O);
#endif //]
