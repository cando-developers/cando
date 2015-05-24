


include	$(CLASP_HOME)/makefile.inc

NEEDS=initClasses.inc symbols.inc

sa subAll: $(NEEDS)
	make scrape
	$(BJAM) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single install_dbg install_opt -j $(COMPILE_PROCESSORS)
	make subBundle


sad subAllDebug: $(NEEDS)
	make scrape
	$(BJAM) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single install_dbg -j $(COMPILE_PROCESSORS)
	make subBundle


san: $(NEEDS)
	$(BJAM) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single install_dbg install_opt -j $(COMPILE_PROCESSORS) -a -n >_bjam.txt

sb subBundle:
	cp libs/[a-z]* $(CANDO_APP_LIB_DIR)
bjd:
	$(BJAM) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single debug -j $(COMPILE_PROCESSORS)

initClasses.inc: *.h
	make mi




subMeta sm:
	make cscrape
#	make m4
	make mi


symbols.inc m4: symbols.m4
	m4 symbols.m4 >symbols.inc

relink:
	rm libs/[a-z]*
	make sa

subClean sc:
	rm -rf bin/*

totalClean:
	rm -rf bin install_release
