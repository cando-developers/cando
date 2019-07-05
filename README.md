# CANDO
## Computer Aided Nanostructure Design and Optimization

Cando is a computational chemistry environment for
building and designing functional molecules and materials.

# Building Cando 

Cando needs [Clasp](https://github.com/clasp-developers/clasp.git) and Clasp has [dependencies](https://github.com/clasp-developers/clasp/wiki/Build-Instructions) that need to be installed first.  Cando has a few more [dependencies](https://github.com/drmeister/cando/wiki/Cando-dependencies).  Advanced use of Cando needs [AmberTools18](http://ambermd.org/GetAmber.php#ambertools) - I assume that anyone downloading Cando is familiar with Amber and Ambertools18.

1. Clone Clasp using: 

```
git clone https://github.com/clasp-developers/clasp.git
```

2. Come up with a wscript.config as described below. (See Note 1)

3. Maybe read 
```
cat ./build-more-systems.README
```
which allows you to control some build parameters.  The default is probably what you want, which is to build Cando and install everything into /opt/clasp, in which case you do nothing.

4. Build and install Cando.

```
./build-more-systems.sh configure build_cboehm install_cboehm
```

5. Run Cando or run the leap emulator.
```
setenv AMBERHOME=/opt/amber # adjust as neede

/opt/clasp/bin/cando

/opt/clasp/bin/cleap
```

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
