import os
from subprocess import call

# configure can provide...
# cfg.extensions_includes  - Append directory names relative to top that will
def configure(cfg):
    cfg.check_cxx(stlib='expat', cflags='-Wall', uselib_store='EXPAT')
    cfg.extensions_stlib += cfg.env.STLIB_EXPAT
    cfg.extensions_stlib += cfg.env.STLIB_Z
    cfg.define("BUILD_EXTENSION",1)   # add this whenever building an extension
#    cfg.define("DEBUG_ENERGY_FUNCTION",1)


def build(bld):
    print("In Extensions build")
    cando_project_headers_name = "include/cando/main/project_headers.h"
    cando_project_headers = bld.path.find_node(cando_project_headers_name)
    assert cando_project_headers!=None, "Could not find %s in %s" % (cando_project_headers_name,bld.path.abspath())
    bld.extensions_gcinterface_include_files.append(cando_project_headers)
    bld.extensions_include_dirs.append("extensions/cando/include")
    bld.extensions_include_dirs.append("extensions/cando/include/cando/main")
    bld.recurse('include')
    bld.recurse('src')
    print("Leaving Extensions build")


def build_extension(cfg):
    print("In cando build_extension")
    print("     clasp -> %s" % cfg.path.find_node("src/corePackage.cc"))
    print("     cfg.path.abspath() -> %s" % cfg.path.abspath())
    root = cfg.path.search_node("../../")
    clasp_executable = "%s/build/clasp" % root
    command = [ clasp_executable, \
                "-l", "source-dir:extensions;cando;src;lisp;build-cando.lisp",
                "-e", "(asdf:make :build-cando :print t)",
                "-l", "source-dir:extensions;cando;src;lisp;build-cando-jupyter.lisp",
                "-e", "(asdf:make :build-cando-jupyter :print t)",
                "-e", "(core:quit)"
    ]
    call(command)

