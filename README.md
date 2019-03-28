CANDO - Computer Aided Nanostructure Design and Optimization
============

Cando is a computational chemistry environment for
building and designing functional molecules and materials.

Cando is an extension of the Clasp Common Lisp programming environment (github.com/clasp-developers/clasp).

# Building

Cando needs [clasp](https://github.com/clasp-developers/clasp.git) and clasp has [dependencies](https://github.com/clasp-developers/clasp/wiki/Build-Instructions) that need to be installed first.

1. Clone clasp using: 

```
git clone https://github.com/clasp-developers/clasp.git
```

2. Add cando to clasp:

```
./waf add_cando
```

3. Configure cando with the installation directory. We use ```/opt/clasp``` for example. 

```
./waf configure --prefix=/opt/clasp
```

4. Build and install cando.

```
./waf install_cboehm
```

5. Run cando or run the leap emulator.
```
setenv AMBERHOME=/path/to/amber

/opt/clasp/bin/cando

/opt/clasp/bin/cleap
```

## Updating clasp and cando

Both clasp and cando are evolving - if you want the latest versions you need to update them both.

1. Update clasp

```
git pull origin master
```

2. Update cando

```
(cd extensions/cando && git pull origin master)
```

3. Rebuild everything

```
./waf build_iboehm
```

