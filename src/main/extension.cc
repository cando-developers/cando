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
#include <clasp/core/hwinfo.h>
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

CL_PRE_GC_STARTUP void cando_initialize_gc()
{
  printf("%s:%d In cando_initialize_gc() do something to initialize the GC\n", __FILE__, __LINE__ );
}


CL_DEFUN core::T_sp core__leap_command_line_scripts()
{
  ql::list scripts;
  int iarg = 1;
  while (iarg<core::global_options->_RawArguments.size() && core::global_options->_RawArguments[iarg] != "--") {
    if (core::global_options->_RawArguments[iarg] == "-f") {
      if (iarg+1 < core::global_options->_RawArguments.size()) {
        scripts << core::SimpleBaseString_O::make(core::global_options->_RawArguments[iarg+1]);
        iarg++;
      } else {
        SIMPLE_ERROR(BF("Illegal -f argument - missing argument"));
      }
    }
    iarg++;
  }
  return scripts.cons();
}

CL_DEFUN core::T_sp core__leap_command_line_includes()
{
  ql::list includes;
  int iarg = 1;
  while (iarg<core::global_options->_RawArguments.size() && core::global_options->_RawArguments[iarg] != "--") {
    if (core::global_options->_RawArguments[iarg] == "-I") {
      if (iarg+1 < core::global_options->_RawArguments.size()) {
        includes << core::SimpleBaseString_O::make(core::global_options->_RawArguments[iarg+1]);
        iarg++;
      } else {
        SIMPLE_ERROR(BF("Illegal -f argument - missing argument"));
      }
    }
    iarg++;
  }
  return includes.cons();
}

CL_DEFUN core::T_sp core__leap_command_line_ignore_leaprc()
{
  int iarg = 1;
  while (iarg<core::global_options->_RawArguments.size() && core::global_options->_RawArguments[iarg] != "--") {
    if (core::global_options->_RawArguments[iarg] == "-s") {
      return _lisp->_true();
    }
  }
  return nil<core::T_O>();
}


void leap_process_arguments(core::CommandLineOptions* options) {
  for ( int iarg = 1; iarg<options->_RawArguments.size(); ++iarg) {
    std::string arg = options->_RawArguments[iarg];
    if (arg == "-h" || arg == "--help") {
      printf("cleap options\n"
             "-I {dir}  - Add {dir} to the leap search path.\n"
             "-f {file} - Source the {file} script.\n"
             "-s        - Ignore the leaprc.\n");
    }
  }
}



CL_INITIALIZER void cando_initializer()
{
  geom::GeomExposer_O* GeomPkg = new geom::GeomExposer_O(_lisp);
  units::UnitsExposer_O* UnitsPkg = new units::UnitsExposer_O(_lisp);
  adapt::AdaptExposer_O* AdaptPkg = new adapt::AdaptExposer_O(_lisp);
  chem::ChemExposer_O* ChemPkg = new chem::ChemExposer_O(_lisp);
  _lisp->installPackage(GeomPkg);
  _lisp->installPackage(UnitsPkg);
  _lisp->installPackage(AdaptPkg);
  _lisp->installPackage(ChemPkg);
  core::Cons_sp features = core::Cons_O::create(_lisp->intern("CANDO",KeywordPkg),cl::_sym_STARfeaturesSTAR->symbolValue());
  std::string executable_path = core::global_options->_ExecutableName;
  std::vector<string> parts = core::split(executable_path,"/");
  std::string exec_name = parts[parts.size()-1];
  if (exec_name.find("cando")!=std::string::npos) {
    std::pair<core::LoadEvalEnum,std::string> cmd1(core::cloEval,"(in-package :cando-user))");
    core::global_options->_LoadEvalList.insert(core::global_options->_LoadEvalList.begin(),cmd1);
    std::pair<core::LoadEvalEnum,std::string> cmd0(core::cloLoad,"source-dir:extensions;cando;src;lisp;start-cando.lisp");
    core::global_options->_LoadEvalList.insert(core::global_options->_LoadEvalList.begin(),cmd0);
  }
  if (exec_name.find("leap")!=std::string::npos) {
    core::global_options->_ProcessArguments = leap_process_arguments;
    std::pair<core::LoadEvalEnum,std::string> cmdNm1(core::cloEval,"(leap:setup-default-paths)");
    core::global_options->_LoadEvalList.push_back(cmdNm1);
    std::pair<core::LoadEvalEnum,std::string> cmdN(core::cloEval,"(leap.commands:leap-repl-then-exit)");
    core::global_options->_LoadEvalList.push_back(cmdN);
    std::pair<core::LoadEvalEnum,std::string> cmd1(core::cloEval,"(in-package :cando-user)");
    core::global_options->_LoadEvalList.insert(core::global_options->_LoadEvalList.begin(),cmd1);
    std::pair<core::LoadEvalEnum,std::string> cmd0(core::cloLoad,"source-dir:extensions;cando;src;lisp;start-cando.lisp");
    core::global_options->_LoadEvalList.insert(core::global_options->_LoadEvalList.begin(),cmd0);
    core::global_options->_Features.push_back("leap-syntax");
  }
}

SYMBOL_EXPORT_SC_(KeywordPkg,wait);
SYMBOL_EXPORT_SC_(KeywordPkg,no_auto_lparallel);

CL_TERMINATOR void cando_terminator()
{
  if (gctools::As<core::Cons_sp>(cl::_sym_STARfeaturesSTAR->symbolValue())->memberEq(kw::_sym_no_auto_lparallel).nilp()) {
    core::T_sp lparallel = _lisp->findPackage("LPARALLEL");
    if (lparallel.notnilp()) {
      core::Symbol_sp end_kernel = _lisp->findSymbol("END-KERNEL",lparallel);
      core::eval::funcall(end_kernel,kw::_sym_wait,_lisp->_true());
    }
  }
}


};
//core::Initializer global_cando_initializer(cando_initializer);

