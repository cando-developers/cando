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

DOCGROUP(cando)
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
        SIMPLE_ERROR(("Illegal -f argument - missing argument"));
      }
    }
    iarg++;
  }
  return scripts.cons();
}

DOCGROUP(cando)
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
        SIMPLE_ERROR(("Illegal -f argument - missing argument"));
      }
    }
    iarg++;
  }
  return includes.cons();
}

DOCGROUP(cando)
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
    if (arg == "-h") {
      printf("Usage: cleap [options]\n"
             "Options:\n"
             " -h         Generate this message.\n"
             " -s         Ignore leaprc startup file.\n"
             " -I {dir}   Add {dir} to search path.\n"
             " -f {file}  Source {file}.\n"
             " -t c       Startup using cclasp and quicklisp for all systems\n");
      exit(0);
    } else if (arg == "-t") {
      options->_StartupStage = options->_RawArguments[++iarg][0];
      printf("Loading stage %c\n", options->_StartupStage );
    }
  }
  options->_LoadEvalList.push_back(std::pair<core::LoadEvalEnum,std::string>(core::cloEval,"(leap.commands:leap-repl-then-exit)"));
  options->_Features.insert("LEAP-SYNTAX");
}


CL_PRE_GC_STARTUP void cando_initialize_gc()
{
  std::vector<string> parts = core::split(core::global_options->_ExecutableName, "/");
  std::string exec_name = parts[parts.size()-1];
  if (exec_name.find("leap") != std::string::npos ||
      core::global_options->_Features.contains("LEAP-SYNTAX")) {
    core::global_options->_ProcessArguments = leap_process_arguments;
  }
}


CL_INITIALIZER void cando_initializer()
{
  if (!core::global_options->_Features.contains("IGNORE-EXTENSIONS")) {
    geom::GeomExposer_O* GeomPkg = new geom::GeomExposer_O(_lisp);
    units::UnitsExposer_O* UnitsPkg = new units::UnitsExposer_O(_lisp);
    adapt::AdaptExposer_O* AdaptPkg = new adapt::AdaptExposer_O(_lisp);
    chem::ChemExposer_O* ChemPkg = new chem::ChemExposer_O(_lisp);
    _lisp->installPackage(GeomPkg);
    _lisp->installPackage(UnitsPkg);
    _lisp->installPackage(AdaptPkg);
    _lisp->installPackage(ChemPkg);
    core::global_options->_Features.insert("CANDO");
  }
}

};

