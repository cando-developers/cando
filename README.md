# CANDO
## Computer Aided Nanostructure Design and Optimization

Cando is an open-source, computational chemistry environment for
building and designing functional molecules and materials.

# Building Cando 

Cando needs [Clasp](https://github.com/clasp-developers/clasp.git) and Clasp has [dependencies](https://github.com/clasp-developers/clasp/wiki/Build-Instructions) that need to be installed first.  Cando has a few more [dependencies](https://github.com/drmeister/cando/wiki/Cando-dependencies).  Advanced use of Cando needs [AmberTools18](http://ambermd.org/GetAmber.php#ambertools) - I assume that anyone downloading Cando is familiar with Amber and Ambertools18.

1. Clone Clasp using: 

```
git clone https://github.com/clasp-developers/clasp.git
```

2. Create a wscript.conf file in the top level clasp directory as described below. (See Note 1)

An example clasp/wscript.config looks like:
```
LLVM_CONFIG_BINARY="/usr/local/opt/llvm@6/bin/llvm-config"   # path to find llvm6 llvm-config executable
PREFIX="/opt/clasp"                                          # installs everything in /opt/clasp
```

3. Build and install Cando.

```
./build-more-systems.sh configure build_cboehm install_cboehm
```

3. Run Cando or run the leap emulator. (We will simplify these instructions)
```
export AMBERHOME=/opt/amber # adjust as needed
export ASDF_EXPORT_TRANSLATIONS=/:
export CLASP_QUICKLISP_DIRECTORY=/opt/clasp/lib/clasp/src/lisp/modules/quicklisp/

/opt/clasp/bin/cleap

/opt/clasp/bin/cando

```

## Free energy perturbation demo (You need AmberTools and Amber pmemd.cuda)

1. Create an empty directory and ```cd``` into it.

2. Set up the following environment variables.
```
export AMBERHOME=/opt/amber # adjust as needed
export ASDF_EXPORT_TRANSLATIONS=/:
export CLASP_QUICKLISP_DIRECTORY=/opt/clasp/lib/clasp/src/lisp/modules/quicklisp/
export PATH=$AMBERHOME/bin:$PATH
```

2. Evaluate:

```
/opt/clasp/bin/cleap -f /opt/clasp/lib/clasp/extensions/cando/src/examples/lysozyme-fep/fep.leaprc

cd jobs

make -j<number-of-cores>    
```

This will take many hours. We have a 'distributor' that distributes the work to GPU accelerated AWS spot instances or a heterogeneous cluster of GPU machines and can carry out the calculation in about 15 min.

## Updating Cando and Clasp

Both Clasp and Cando are evolving - if you want the latest versions you need to update them both.

1. Change directory into the top Clasp directory and pull Clasp updates.  I'm assuming you are using the master branch.

```
git pull origin master
```

2. Update Cando.

```
(cd extensions/cando && git pull origin master)
```

3. Rebuild and reinstall as above.

Note 1. Clasp can be configured by copying the wscript.config.template to wscript.config and editing it. The most important option is the LLVM_CONFIG_BINARY - it is the path to the ```llvm-config``` executable for llvm6.  Once clasp has that it can find the rest of llvm.  *On macOS and Linux you must download llvm6 using homebrew or some other package manager. Clasp uses a lot of advanced features of the llvm library.*

Example wscript.config
```
LLVM_CONFIG_BINARY="/usr/local/opt/llvm/bin/llvm-config"
PREFIX="/opt/clasp"                                          # installs everything in /opt/clasp
```
