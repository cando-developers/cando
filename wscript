#-*- mode: python; coding: utf-8-unix -*-

import sys, os
from subprocess import call
from waflib import Task

sys.path.append('../../tools-for-build/')
from wscript_utils import log, clasp_task

def configure(cfg):
    cfg.check_cxx(stlib='expat', cflags='-Wall', uselib_store='EXPAT')
    cfg.extensions_stlib += cfg.env.STLIB_EXPAT
    cfg.extensions_stlib += cfg.env.STLIB_Z
    cfg.define("BUILD_EXTENSION",1)   # add this whenever building an extension
#    cfg.define("DEBUG_ENERGY_FUNCTION",1)

def build(bld):
    log.pprint('BLUE', 'cando build(), %s, stage %s=%s' %
                (bld.variant,
                 bld.stage,
                 bld.stage_val))
    cando_project_headers_name = "include/cando/main/project_headers.h"
    cando_project_headers = bld.path.find_node(cando_project_headers_name)
    assert cando_project_headers != None, "Could not find %s in %s" % (cando_project_headers_name, bld.path.abspath())
    bld.extensions_gcinterface_include_files.append(cando_project_headers)
    bld.extensions_include_dirs.append("extensions/cando/include")
    bld.extensions_include_dirs.append("extensions/cando/include/cando/main")
    bld.extensions_lisp_files = []
    bld.recurse('include')
    bld.recurse('src')
    if (bld.stage_val >= 5):
        task = build_cando(env=bld.env)
        task.set_inputs([bld.iclasp_executable,
                         bld.cclasp_fasl,
                         bld.asdf_fasl_cclasp] +
                        bld.extensions_lisp_files)
        task.set_outputs([bld.path.parent.parent.find_or_declare("fasl/cando.fasb"),
                          bld.path.parent.parent.find_or_declare("fasl/cando-jupyter.fasb")])
        bld.install_files('${PREFIX}/lib/clasp/', task.outputs, relative_trick = True, cwd = bld.path.parent.parent)   #source
        bld.add_to_group(task)
    log.pprint('BLUE', 'cando build() has finished')

class build_cando(clasp_task):
    #ext_out = ['.fasl']
    def run(self):
        cmd = self.clasp_command_line(self.inputs[0].abspath(),
                                      image = False,
                                      forms = ['(load "source-dir:extensions;cando;src;lisp;build-cando.lisp")',
                                               '(asdf:make :build-cando :print t :verbose t)',
                                               '(load "source-dir:extensions;cando;src;lisp;build-cando-jupyter.lisp")',
                                               '(asdf:make :build-cando-jupyter :print t :verbose t)',
                                               '(core:quit)'],
                                      *self.bld.clasp_aclasp)
        return self.exec_command(cmd)
