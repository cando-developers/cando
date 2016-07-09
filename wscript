import os

# configure can provide...
# cfg.plugins_includes  - Append directory names relative to top that will
def configure(cfg):
    cfg.plugins_include_dirs.append("plugins/cando/include")
    cfg.plugins_gcinterface_include_files.append("cando/main/project_headers.h")
    cfg.check_cxx(stlib='expat', cflags='-Wall', uselib_store='EXPAT')
    cfg.plugins_stlib += cfg.env.STLIB_EXPAT
    cfg.plugins_stlib += cfg.env.STLIB_Z
    
def build(bld):
    bld.recurse('include')
    bld.recurse('src')

