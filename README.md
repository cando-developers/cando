CANDO - Computer Aided Nanostructure Design and Optimization
============

Cando is a computational chemistry environment for
building and designing functional molecules and materials.

# Building Cando 

Cando needs [Clasp](https://github.com/clasp-developers/clasp.git) and Clasp has [dependencies](https://github.com/clasp-developers/clasp/wiki/Build-Instructions) that need to be installed first.  Advanced use of Cando needs [AmberTools18](http://ambermd.org/GetAmber.php#ambertools) - I assume that anyone downloading Cando is familiar with Amber and Ambertools18.

1. Clone Clasp using: 

```
git clone https://github.com/clasp-developers/clasp.git
```

2. Add Cando to Clasp:

```
./waf add_cando
```

3. Configure Cando with the installation directory. We use ```/opt/clasp``` for example. 

```
./waf configure --prefix=/opt/clasp
```

4. Build and install Cando.

```
./waf install_cboehm
```

5. Run Cando or run the leap emulator.
```
setenv AMBERHOME=/path/to/amber

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

3. Build and install everything.

```
./waf build_iboehm
```

