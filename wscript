import sys
from subprocess import call
import waflib.Task
from waflib import Utils, Logs, Task, TaskGen
import os

def options(ctx):
    ctx.add_option('--enable-jupyter', action = 'store_true', dest = 'enable_jupyter',
                   help = 'Build jupyter-lab version of cando')

# configure can provide...
# cfg.extensions_includes  - Append directory names relative to top that will
def configure(cfg):
    cfg.check_cxx(stlib='expat', cflags='-Wall', uselib_store='EXPAT')
    cfg.extensions_stlib += cfg.env.STLIB_EXPAT
    cfg.extensions_stlib += cfg.env.STLIB_Z
    cfg.check_cxx(stlib='boost_graph', cflags='-Wall', uselib_store='BOOST-boost_graph')
    cfg.define("BUILD_EXTENSION",1)   # add this whenever building an extension
#    cfg.define("DEBUG_ENERGY_FUNCTION",1)
    fetch_git_revision("src/lisp/modules/quicklisp", "https://github.com/quicklisp/quicklisp-client.git", label="master")
    fetch_git_revision("src/lisp/modules/quicklisp/local-projects/cl-netcdf","https://github.com/clasp-developers/cl-netcdf.git",label="master")
    fetch_git_revision("src/lisp/modules/quicklisp/local-projects/static-vectors","https://github.com/clasp-developers/static-vectors.git",label="master")
    if (cfg.options.enable_jupyter):
        print("Pulling in jupyter-lab code")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/trivial-garbage", "https://github.com/clasp-developers/trivial-garbage.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/bordeaux-threads", "https://github.com/clasp-developers/bordeaux-threads.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/cffi", "https://github.com/clasp-developers/cffi.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/usocket", "https://github.com/clasp-developers/usocket.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/uuid", "https://github.com/clasp-developers/uuid.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/cl-jupyter", "https://github.com/drmeister/cl-jupyter.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/cl-ipykernel", "https://github.com/clasp-developers/cl-ipykernel.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/cl-ipywidgets", "https://github.com/clasp-developers/cl-ipywidgets.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/cl-nglview", "https://github.com/clasp-developers/cl-nglview.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/cl-bqplot", "https://github.com/clasp-developers/cl-bqplot.git", label="master")
	fetch_git_revision("src/lisp/modules/quicklisp/local-projects/trivial-backtrace", "https://github.com/clasp-developers/trivial-backtrace", label="master")


class duplicate_executable(waflib.Task.Task):
    def run(self):
        print("     creating command for %s %s" % (self.inputs[0].abspath(),self.outputs[0].abspath()))
        return self.exec_command('cp %s %s' % ( self.inputs[0].abspath(),self.outputs[0].abspath()))

def rename_executable(p,frm,to):
    parts = os.path.split(p);
    return "%s/%s" % (parts[0],parts[1].replace(frm,to))
    
def build(bld):
    print("bld.options.enable_jupyter -> %s" % bld.options.enable_jupyter)
    bld.env.options = bld.options
    print("bld.env.options.enable_jupyter -> %s" % bld.env.options.enable_jupyter)
    # The following will copy iclasp-<gc> to icando-<gc>
    cp_1 = duplicate_executable(env=bld.env)
    cp_1.set_inputs(bld.iclasp_executable)
    bld.icando_executable = bld.path.find_or_declare(rename_executable(bld.iclasp_executable.abspath(),"clasp","cando"))
    bld.ccando_executable = bld.path.find_or_declare(rename_executable(bld.cclasp_executable.abspath(),"clasp","cando"))
    cp_1.set_outputs(bld.icando_executable)
    bld.add_to_group(cp_1)
    bld.install_as('${PREFIX}/bin/%s' % bld.icando_executable.name, bld.icando_executable, chmod = Utils.O755)
    if (bld.stage_val >= 3):
        # The following will copy cclasp-<gc> to ccando-<gc>
        cp_2 = duplicate_executable(env=bld.env)
        cp_2.set_inputs(bld.cclasp_executable)
        cp_2.set_outputs(bld.ccando_executable)
        bld.add_to_group(cp_2)
        bld.install_as('${PREFIX}/bin/%s' % bld.ccando_executable.name, bld.ccando_executable, chmod = Utils.O755)
        print("Going to build cando")        
    # The following will copy iclasp-<gc> to ileap-<gc>
    bld.ileap_executable = bld.path.find_or_declare(rename_executable(bld.iclasp_executable.abspath(),"clasp","leap"))
    bld.cleap_executable = bld.path.find_or_declare(rename_executable(bld.cclasp_executable.abspath(),"clasp","leap"))
    cp_3 = duplicate_executable(env=bld.env)
    cp_3.set_inputs(bld.iclasp_executable)
    cp_3.set_outputs(bld.path.find_or_declare(rename_executable(bld.iclasp_executable.abspath(),"clasp","leap")))
    bld.add_to_group(cp_3)
    bld.install_as('${PREFIX}/bin/%s' % bld.ileap_executable.name, bld.ileap_executable, chmod = Utils.O755)
    if (bld.stage_val >= 3):
        # The following will copy cclasp-<gc> to cleap-<gc>
        cp_4 = duplicate_executable(env=bld.env)
        cp_4.set_inputs(bld.cclasp_executable)
        cp_4.set_outputs(bld.path.find_or_declare(rename_executable(bld.cclasp_executable.abspath(),"clasp","leap")))
        bld.add_to_group(cp_4)
        bld.install_as('${PREFIX}/bin/%s' % bld.cleap_executable.name, bld.cleap_executable, chmod = Utils.O755)
        print("Going to build leap")
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
        print("bld.path.abspath() -> %s" % bld.path.abspath())
        extensions_result = bld.path.parent.parent.find_or_declare("%s/src/lisp/cando-user/packages.fasl" % (bld.path.abspath()))
        print("extensions_result -> %s" % extensions_result.abspath())
        bld_extensions = build_extension(env=bld.env)
        bld_extensions.set_inputs([bld.ccando_executable,bld.cclasp_fasl,bld.asdf_fasl_cclasp]+bld.extensions_lisp_files)
#        bld_extensions.set_outputs([extensions_result])
        bld.add_to_group(bld_extensions)
    print("Leaving extensions build without cclasp_executable")

class build_extension(waflib.Task.Task):
    def run(self):
        print("In run self.env.options.enable_jupyter -> %s" % self.env.options.enable_jupyter)
        cmd = [ self.inputs[0].abspath(), "-N"]
        if (self.env.options.enable_jupyter):
            cmd = cmd + [ "-e", "(ql:quickload :cando-jupyter)" ]
        cmd = cmd + [ "-e", "(core:quit)" ]
        print("build_extension cmd -> %s" % cmd)
        print("build_extension outputs -> %s" % self.outputs)
        return self.exec_command(cmd)
    def exec_command(self, cmd, **kw):
        kw['stdout'] = sys.stdout
        return super(build_extension, self).exec_command(cmd, **kw)
    def keyword(self):
        return 'build extensions using... '

def fetch_git_revision(path, url, revision = "", label = "master"):
    print("Git repository %s  url: %s\n     revision: %s  label: %s\n" % (path, url, revision, label))
    ret = os.system("./tools-for-build/fetch-git-revision.sh '%s' '%s' '%s' '%s'" % (path, url, revision, label))
    if ( ret != 0 ):
        raise Exception("Failed to fetch git url %s" % url)
    
