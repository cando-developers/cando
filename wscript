import sys
from subprocess import call
import waflib.Task
import os

# configure can provide...
# cfg.extensions_includes  - Append directory names relative to top that will
def configure(cfg):
    cfg.check_cxx(stlib='expat', cflags='-Wall', uselib_store='EXPAT')
    cfg.extensions_stlib += cfg.env.STLIB_EXPAT
    cfg.extensions_stlib += cfg.env.STLIB_Z
    cfg.check_cxx(stlib='boost_graph', cflags='-Wall', uselib_store='BOOST-boost_graph')
    cfg.define("BUILD_EXTENSION",1)   # add this whenever building an extension
#    cfg.define("DEBUG_ENERGY_FUNCTION",1)


class duplicate_executable(waflib.Task.Task):
    def run(self):
        return self.exec_command('cp %s %s' % ( self.inputs[0].abspath(),self.outputs[0].abspath()))
    
def build(bld):
    # The following will copy iclasp-<gc> to icando-<gc>
    cp_1 = duplicate_executable(env=bld.env)
    cp_1.set_inputs(bld.iclasp_executable)
    cp_1.set_outputs(bld.path.find_or_declare(bld.iclasp_executable.abspath().replace("clasp","cando")))
    bld.add_to_group(cp_1)
    if (bld.stage_val >= 3):
        # The following will copy cclasp-<gc> to ccando-<gc>
        cp_2 = duplicate_executable(env=bld.env)
        cp_2.set_inputs(bld.cclasp_executable)
        cp_2.set_outputs(bld.path.find_or_declare(bld.cclasp_executable.abspath().replace("clasp","cando")))
        bld.add_to_group(cp_2)
    # The following will copy iclasp-<gc> to ileap-<gc>
    cp_3 = duplicate_executable(env=bld.env)
    cp_3.set_inputs(bld.iclasp_executable)
    cp_3.set_outputs(bld.path.find_or_declare(bld.iclasp_executable.abspath().replace("clasp","leap")))
    bld.add_to_group(cp_3)
    if (bld.stage_val >= 3):
        # The following will copy cclasp-<gc> to cleap-<gc>
        cp_4 = duplicate_executable(env=bld.env)
        cp_4.set_inputs(bld.cclasp_executable)
        cp_4.set_outputs(bld.path.find_or_declare(bld.cclasp_executable.abspath().replace("clasp","leap")))
        bld.add_to_group(cp_4)
    print("In extensions build bld.cclasp_executable = %s" % bld.cclasp_executable)
    print("      bld.stage_val = %s" % bld.stage_val)
    cando_project_headers_name = "include/cando/main/project_headers.h"
    cando_project_headers = bld.path.find_node(cando_project_headers_name)
    assert cando_project_headers!=None, "Could not find %s in %s" % (cando_project_headers_name,bld.path.abspath())
    bld.extensions_gcinterface_include_files.append(cando_project_headers)
    bld.extensions_include_dirs.append("extensions/cando/include")
    bld.extensions_include_dirs.append("extensions/cando/include/cando/main")
    bld.extensions_lisp_files = []
    bld.recurse('include')
    bld.recurse('src')
    if (bld.stage_val>=5):
        bld_extensions = build_extension(env=bld.env)
        bld_extensions.set_inputs([bld.cclasp_executable,bld.cclasp_fasl,bld.asdf_fasl_cclasp]+bld.extensions_lisp_files)
        bld_extensions.set_outputs([bld.path.parent.parent.find_or_declare("fasl/cando.fasb"),
                                    bld.path.parent.parent.find_or_declare("fasl/cando-jupyter.fasb")])
        bld.install_files('${PREFIX}/lib/clasp/', bld_extensions.outputs, relative_trick = True, cwd = bld.path.parent.parent)   #source
        
        bld.add_to_group(bld_extensions)
    print("Leaving extensions build without cclasp_executable")

class build_extension(waflib.Task.Task):
    def run(self):
        cmd = [ self.inputs[0].abspath(),
                "-N",
                "-l", "source-dir:extensions;cando;src;lisp;load-cando.lisp",
                "-l", "source-dir:extensions;cando;src;lisp;load-cando-jupyter.lisp",
                "-e", "(core:quit)" ]
        print("build_extension cmd -> %s" % cmd)
        print("build_extension outputs -> %s" % self.outputs)
        return self.exec_command(cmd)
    def exec_command(self, cmd, **kw):
        kw['stdout'] = sys.stdout
        return super(build_extension, self).exec_command(cmd, **kw)
    def keyword(self):
        return 'build extensions using... '
