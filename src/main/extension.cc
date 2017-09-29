/*
    File: main.cc
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
#define DEBUG_LEVEL_NONE

#ifdef _TARGET_OS_LINUX
#include <signal.h>
#include <sys/resource.h>
#endif

#ifdef USE_MPI
#include <boost/mpi.hpp>
#endif
#include <string>
#include <clasp/core/foundation.h>
#include <clasp/core/bundle.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/compiler.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/array.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/candoOpenMp.h>
#include <clasp/core/cons.h>
#include <clasp/core/commandLineOptions.h>
#include <clasp/llvmo/llvmoPackage.h>
#include <clasp/gctools/gctoolsPackage.h>
#include <clasp/clbind/clbindPackage.h>
#include <clasp/sockets/socketsPackage.h>
#include <clasp/serveEvent/serveEventPackage.h>
#include <clasp/asttooling/asttoolingPackage.h>
#include <clasp/core/pathname.h>
#include <cando/geom/geomPackage.h>
#include <cando/units/unitsPackage.h>
#include <cando/adapt/adaptPackage.h>
#include <cando/chem/readAmberParameters.h>
#include <cando/chem/chemPackage.h>
//#include <cando/kinematics/kinematicsPackage.h>

#ifdef USE_MPI
#include <clasp/mpip/mpiPackage.h>
#include <clasp/mpip/claspMpi.h>
#endif





namespace cando {

CL_INITIALIZER void cando_initializer()
{
//  printf("%s:%d In cando_initializer\n", __FILE__, __LINE__);
  geom::GeomExposer_O* GeomPkg = new geom::GeomExposer_O(_lisp);
  units::UnitsExposer_O* UnitsPkg = new units::UnitsExposer_O(_lisp);
  adapt::AdaptExposer_O* AdaptPkg = new adapt::AdaptExposer_O(_lisp);
  chem::ChemExposer_O* ChemPkg = new chem::ChemExposer_O(_lisp);
  _lisp->installPackage(GeomPkg);
  _lisp->installPackage(UnitsPkg);
  _lisp->installPackage(AdaptPkg);
  _lisp->installPackage(ChemPkg);
  core::Cons_sp features = core::Cons_O::create(_lisp->intern("CANDO",KeywordPkg),cl::_sym_STARfeaturesSTAR->symbolValue());
  core::List_sp loads = core::_sym_STARextension_startup_loadsSTAR->symbolValue();
  loads = core::Cons_O::create(core::Cons_O::createList(cl::_sym_load,core::cl__pathname(core::Str_O::create("source-dir:extensions;cando;src;lisp;start-cando.lisp"))),loads);
  core::_sym_STARextension_startup_loadsSTAR->defparameter(loads);
}


};
//core::Initializer global_cando_initializer(cando_initializer);

