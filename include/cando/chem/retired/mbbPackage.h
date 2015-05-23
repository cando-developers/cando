#ifndef	mbbPackage_H
#define mbbPackage_H

#include <cando/chem/chemPackage.h>

#include "core/lisp.fwd.h"// mbbPackage.h wants Lisp needs lisp.fwd.h
#include "core/package.fwd.h"// mbbPackage.h wants Package needs package.fwd.h

namespace chem
{
SMART(Lisp);
SMART(Package);
//
// Non-class functions that belong in the mbbPackage are initialized here
//
inline void initialize_mbbPackage( core::Lisp_sp lisp, core::Package_sp mbbPackage )
{
//    initialize_tableOperations(lisp);
}
};
#endif
