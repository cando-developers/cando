

all:
	make -C src/main boehm-release-cxx
	make cando-lisp
	make cando-data

cando-lisp:
	(cd src/lisp; bjam bundle )

cando-data:
	(cd src/data; bjam data )
