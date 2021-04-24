from waflib.Tools import c
from waflib.Tools import cxx
import sys
import time
import subprocess
import waflib.Task
from waflib import Utils, Logs, Task, TaskGen
import os
from wscript_utils import *

def options(ctx):
    ctx.add_option('--enable-jupyter', action = 'store_true', dest = 'enable_jupyter',
                   help = 'Build jupyter-lab version of cando')

# configure can provide...
# cfg.extensions_includes  - Append directory names relative to top that will

def analyze_clasp(cfg):
    cfg.extensions_clasp_gc_names += [ "cando" ]

def configure(cfg):
    print("----- Configuring cando extension")
    cfg.check_cxx(stlib='expat', cflags='-Wall', uselib_store='EXPAT')
    cfg.extensions_stlib += cfg.env.STLIB_EXPAT
    cfg.extensions_stlib += cfg.env.STLIB_Z
    cfg.extensions_clasp_gc_names += "cando"
    cfg.check_cxx(stlib='boost_graph', cflags='-Wall', uselib_store='BOOST-boost_graph')
    cfg.check_cxx(lib='netcdf', cflags='-Wall', uselib_store='BOOST-netcdf',msg="Checking for netcdf [netcdf]")
    if (cfg.options.enable_jupyter):
        cfg.env.enable_jupyter = True
        cfg.check_cxx(stlib='zmq', cflags='-Wall', uselib_store='ZMQ')
    cfg.define("BUILD_EXTENSION",1)   # add this whenever building an extension
#    cfg.define("DEBUG_ENERGY_FUNCTION",1)

class duplicate_executable(waflib.Task.Task):
    def run(self):
        print("     creating command for %s %s" % (self.inputs[0].abspath(),self.outputs[0].abspath()))
        return self.exec_command('cp %s %s' % ( self.inputs[0].abspath(),self.outputs[0].abspath()))

def rename_executable(p,frm,to):
    parts = os.path.split(p);
    return "%s/%s" % (parts[0],parts[1].replace(frm,to))
    
def build(bld):
    print("Recursed into bld bld.stage_val == %d" % bld.stage_val )
    print("bld.env.enable_jupyter -> %s" % bld.env.enable_jupyter)
    bld.env.options = bld.options
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

def build3(bld):
    print("Recursed into extensions/cando/wscript:build3")
    bld.ccando_executable = bld.path.find_or_declare(rename_executable(bld.cclasp_executable.abspath(),"clasp","cando"))
    task = symlink_executable(env=bld.env)
    task.set_inputs(bld.iclasp_executable)
    task.set_outputs(bld.ccando_executable)
    bld.add_to_group(task)
    bld.symlink_as('${PREFIX}/bin/icando',bld.iclasp_executable.abspath().split("/")[-1])

def build4(bld):
    log.info("bld.stage_val cando = %d" % bld.stage_val )
    bld.add_group()
    bld_extensions = build_cando_extension(env=bld.env)
    cando_system_files = bld.path.ant_glob("src/lisp/**/*.lisp")
    log.info("cando_system_files = %s" % cando_system_files )
    start_cando_file = bld.path.find_node("src/lisp/start-cando.lisp")
    if (start_cando_file == None):
       raise Exception("Could not fild start_cando_file")       
    log.info("start_cando_file = %s" % start_cando_file )
    snapshot_file = bld.path.find_or_declare("generated/%s.snapshot" % "cando")
    log.info("snapshot_file -> %s" % snapshot_file.abspath())
    bld_extensions.set_inputs([bld.ccando_executable,
                               bld.cclasp_link_product,
                               start_cando_file] + cando_system_files )
    bld_extensions.set_outputs([snapshot_file])
    bld.add_to_group(bld_extensions)
    bld.add_group()
    bld.dcando_executable = bld.path.find_or_declare("cando")
    embed_snapshot(bld,snapshot_file,bld.ccando_executable,bld.dcando_executable,"cando")

def post_install(ctx):
    print('In post_install(ctx)')
    if (ctx.cmd[:len('install_c')] == 'install_c'):
        prefix = ctx.env.PREFIX
        cando_parts = os.path.split(ctx.ccando_executable.abspath())
        installed_cando = "%s/bin/%s" % (prefix, cando_parts[1])
        # this trick with test -e /dev/fd/3 doesn't work on FreeBSD - cracauer fixme
        # cmd = '%s -e "(sys:quit)" 2>&1 | (ls -l /dev/fd > cboehm_install.log 2>&1 ; tee -a cboehm_install.log) | (echo foo 1>&3 && tee 1>&3 || cat)' % installed_cando
        cmd = '%s -N -D -e "(sys:quit)"' % installed_cando
        print("Executing post-install command %s" % cmd)
        print("environment variable  ASDF_OUTPUT_TRANSLATIONS=%s" % os.getenv("ASDF_OUTPUT_TRANSLATIONS")) 
        print("environment variable            XDG_CACHE_HOME=%s" % os.getenv("XDG_CACHE_HOME")) 
        print("environment variable CLASP_QUICKLISP_DIRECTORY=%s" % os.getenv("CLASP_QUICKLISP_DIRECTORY")) 
        sys.stdout.flush()
        print("Start time: %s" % time.asctime())
        sys.stdout.flush()
        # Try subprocess.call instead of exec_command that eats output
        subprocess.call(["/bin/sh","-c",cmd])
        #ctx.exec_command(cmd)
        cando_symlink = "%s/bin/cando" % prefix
        leap_parts = os.path.split(ctx.cleap_executable.abspath())
        installed_cleap = "%s/bin/%s" % (prefix, leap_parts[1])
        cleap_symlink = "%s/bin/cleap" % prefix
        try:
            if (os.stat(cando_symlink)):
                os.unlink(cando_symlink)
        except OSError as e:
            pass
        os.symlink(installed_cando,cando_symlink)
        try:
            if (os.stat(cleap_symlink)):
                os.unlink(cleap_symlink)
        except OSError as e:
            pass
        os.symlink(installed_cleap,cleap_symlink)

#            ctx.exec_command('${PREFIX}/bin/icando-boehm')


class build_cando_extension(waflib.Task.Task):
    def run(self):
        print("In build_extension self.env.enable_jupyter -> %s" % self.env.enable_jupyter)
        cmd = [ self.inputs[0].abspath(), "-N"]
        cmd = cmd + [ "-e", "(ql:update-all-dists :prompt nil)" ]
        if (self.env.enable_jupyter):
            cmd = cmd + [ "-e", "(ql:quickload :cando-jupyter)" ]
        saveFile = self.outputs[0].abspath()
        cmd = cmd + [ "-e", "(setf ext:*image-save-load-startup* 'cl-user:start-cando-user-from-snapshot)"]
        cmd = cmd + [ "-e", "(gctools:save-lisp-and-die \"%s\")" % saveFile ]
        cmd = cmd + [ "-e", "(core:quit)" ]
        print("build_extension cmd -> %s" % cmd)
        print("build_extension outputs -> %s" % self.outputs)
        return self.exec_command(cmd)
    def exec_command(self, cmd, **kw):
        kw['stdout'] = sys.stdout
        return super(build_cando_extension, self).exec_command(cmd, **kw)
    def keyword(self):
        return 'build extensions using... '

