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
    
def build1(bld):
    print("Recursed into bld bld.stage_val == %d" % bld.stage_val )
    print("bld.options.enable_jupyter -> %s" % bld.options.enable_jupyter)
    bld.env.options = bld.options
    print("bld.env.options.enable_jupyter -> %s" % bld.env.options.enable_jupyter)
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
        
def build5(bld):
    if (bld.stage_val >= 5):
        print("bld.stage_val cando = %d" % bld.stage_val )
        # The following will copy cclasp-<gc> to ccando-<gc>
        bld.ccando_executable = bld.path.find_or_declare(rename_executable(bld.cclasp_executable.abspath(),"clasp","cando"))
        cp_2 = duplicate_executable(env=bld.env)
        cp_2.set_inputs(bld.cclasp_executable)
        cp_2.set_outputs(bld.ccando_executable)
        bld.add_to_group(cp_2)
        bld_extensions = build_extension(env=bld.env)
        snapshot_file = bld.path.find_or_declare("generated/cando.snapshot")
        print("snapshot_file -> %s" % snapshot_file.abspath())
        bld_extensions.set_inputs([bld.ccando_executable,
                                   bld.cclasp_link_product])
        bld_extensions.set_outputs([snapshot_file])
        bld.add_to_group(bld_extensions)
        bld.add_group()
        bld.dcando_executable = bld.path.find_or_declare("cando")
        print("dcando target -> %s" % bld.dcando_executable )
        print("bld.env[DEST_OS] = %s" % bld.env["DEST_OS"] )
        if (bld.env["DEST_OS"] == DARWIN_OS):
            print("dcando target -> %s" % bld.dcando_executable )
            env2 = bld.env.derive()
            env2.append_value("LINKFLAGS",["-Wl,-exported_symbols_list",bld.exported_symbols_file.abspath()])
            env2.append_value("LINKFLAGS",["-sectcreate", "__CLASP", "__clasp", snapshot_file.abspath()])
            link2 = cxx.cxxprogram(env=env2)
            link2.name = "final_build"
            link2.set_inputs( bld.iclasp_link_task.inputs) # snapshot_file
            link2.set_outputs( [ bld.dcando_executable ] )
            bld.add_to_group(link2)

#         cando_snapshot_product = bld.variant_obj.snapshot_node(bld,appname="cando",stage='c')
#         print("cando_snapshot_product = %s" % cando_snapshot_product)
#         leap_snapshot_product = bld.variant_obj.snapshot_node(bld,appname="leap",stage='c')
#         print("leap_snapshot_product = %s" % leap_snapshot_product)
#         task_cando_snapshot = link_snapshot(env=bld.env)
#         task_cando_snapshot.set_inputs([bld.iclasp_executable,
#                                         bld.cclasp_link_product ])
#         task_cando_snapshot.set_outputs([cando_snapshot_product])
#         bld.add_to_group(task_cando_snapshot)
#         task_leap_snapshot = link_snapshot(env=bld.env)
#         task_leap_snapshot.set_inputs([bld.ileap_executable,
#                                         bld.cclasp_link_product ])
#         task_leap_snapshot.set_outputs([leap_snapshot_product])
#         bld.add_to_group(task_leap_snapshot)
        
    #     bld.install_as('${PREFIX}/bin/%s' % bld.ccando_executable.name, bld.ccando_executable, chmod = Utils.O755)
    #     print("Going to build cando")
    # # The following will copy iclasp-<gc> to ileap-<gc>
    # bld.ileap_executable = bld.path.find_or_declare(rename_executable(bld.iclasp_executable.abspath(),"clasp","leap"))
    # bld.cleap_executable = bld.path.find_or_declare(rename_executable(bld.cclasp_executable.abspath(),"clasp","leap"))
    # cp_3 = duplicate_executable(env=bld.env)
    # cp_3.set_inputs(bld.iclasp_executable)
    # cp_3.set_outputs(bld.path.find_or_declare(rename_executable(bld.iclasp_executable.abspath(),"clasp","leap")))
    # bld.add_to_group(cp_3)
    # bld.install_as('${PREFIX}/bin/%s' % bld.ileap_executable.name, bld.ileap_executable, chmod = Utils.O755)
    # if (bld.stage_val >= 3):
    #     # The following will copy cclasp-<gc> to cleap-<gc>
    #     cp_4 = duplicate_executable(env=bld.env)
    #     cp_4.set_inputs(bld.cclasp_executable)
    #     cp_4.set_outputs(bld.path.find_or_declare(rename_executable(bld.cclasp_executable.abspath(),"clasp","leap")))
    #     bld.add_to_group(cp_4)
    #     bld.install_as('${PREFIX}/bin/%s' % bld.cleap_executable.name, bld.cleap_executable, chmod = Utils.O755)
    #     print("Going to build leap")
#     if (bld.stage_val>=5):
#         print("bld.path.abspath() -> %s" % bld.path.abspath())
#         extensions_result = bld.path.parent.parent.find_or_declare("%s/src/lisp/cando-user/packages.fasl" % (bld.path.abspath()))
#         print("extensions_result -> %s" % extensions_result.abspath())
#         bld_extensions = build_extension(env=bld.env)
#         bld_extensions.set_inputs([bld.icando_executable,bld.cclasp_link_product,bld.cclasp_asdf_fasl]+bld.extensions_lisp_files)
# #        bld_extensions.set_outputs([extensions_result])
#         bld.add_to_group(bld_extensions)
#         cando_snapshot_product = bld.variant_obj.snapshot_node(bld,appname="cando",stage='c')
#         print("cando_snapshot_product = %s" % cando_snapshot_product)
#         leap_snapshot_product = bld.variant_obj.snapshot_node(bld,appname="leap",stage='c')
#         print("leap_snapshot_product = %s" % leap_snapshot_product)
#         task_cando_snapshot = link_snapshot(env=bld.env)
#         task_cando_snapshot.set_inputs([bld.iclasp_executable,
#                                         bld.cclasp_link_product ])
#         task_cando_snapshot.set_outputs([cando_snapshot_product])
#         bld.add_to_group(task_cando_snapshot)
#         task_leap_snapshot = link_snapshot(env=bld.env)
#         task_leap_snapshot.set_inputs([bld.ileap_executable,
#                                         bld.cclasp_link_product ])
#         task_leap_snapshot.set_outputs([leap_snapshot_product])
#         bld.add_to_group(task_leap_snapshot)
#     bld.add_post_fun(post_install)
#     if (bld.stage_val >= 4):
#         # The following will copy iclasp-<gc> to icando-<gc>
#         cp_1 = duplicate_executable(env=bld.env)
#         cp_1.set_inputs(bld.cclasp_executable)
#         bld.ccando_executable = bld.path.find_or_declare(rename_executable(bld.cclasp_executable.abspath(),"clasp","cando"))
#         cp_1.set_outputs(bld.ccando_executable)
#         bld.add_to_group(cp_1)
#         bld.install_as('${PREFIX}/bin/%s' % bld.ccando_executable.name, bld.ccando_executable, chmod = Utils.O755)

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


class build_extension(waflib.Task.Task):
    def run(self):
        print("In run self.env.options.enable_jupyter -> %s" % self.env.options.enable_jupyter)
        cmd = [ self.inputs[0].abspath(), "-N"]
        if (self.env.options.enable_jupyter):
            cmd = cmd + [ "-e", "(ql:quickload :cando-jupyter)" ]
        saveFile = self.outputs[0].abspath()
        cmd = cmd + [ "-e", "(gctools:save-lisp-and-die \"%s\")" % saveFile ]
        cmd = cmd + [ "-e", "(core:quit)" ]
        print("build_extension cmd -> %s" % cmd)
        print("build_extension outputs -> %s" % self.outputs)
        return self.exec_command(cmd)
    def exec_command(self, cmd, **kw):
        kw['stdout'] = sys.stdout
        return super(build_extension, self).exec_command(cmd, **kw)
    def keyword(self):
        return 'build extensions using... '
