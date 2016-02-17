

all:
	make -C src/main boehm-release-cxx
	make cando-lisp
	make cando-data

cando-lisp:
	(cd src/lisp; bjam bundle )

cando-data:
	(cd src/data; bjam data )

cloc-files:
	find src/adapt -name '*.cc' -print >/tmp/cando.cloc
	find src/chem -name '*.cc' -print >>/tmp/cando.cloc
	find src/geom -name '*.cc' -print >>/tmp/cando.cloc
	find src/main -name '*.cc' -print >>/tmp/cando.cloc
	find src/units -name '*.cc' -print >>/tmp/cando.cloc
	find include/cando/adapt -name '*.h' -print >>/tmp/cando.cloc
	find include/cando/chem -name '*.h' -print >>/tmp/cando.cloc
	find include/cando/geom -name '*.h' -print >>/tmp/cando.cloc
	find include/cando/main -name '*.h' -print >>/tmp/cando.cloc
	find include/cando/units -name '*.h' -print >>/tmp/cando.cloc
