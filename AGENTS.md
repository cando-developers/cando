#
## Cando is built within Clasp

Cando is within a directory tree cando/extensions/cando/...

## Adding new classes

Add a header file under the cando/include directory tree and a source (.cc) file in the cando/src/ directory tree.
The .cc file needs to also be added to the cscript.lisp file for the build system.
Then the clasp/build directory tree needs to be cleared out and run ./koga --deep-clean to reset the build system.
