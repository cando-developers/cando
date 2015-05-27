


include	$(CLASP_HOME)/makefile.inc


M4_NEEDS = openmmExpose.scrape.inc
NEEDS = initClasses.inc $(M4_NEEDS)



sad subAllDebug: $(NEEDS)
	make scrape
	$(BJAM) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single  -j $(COMPILE_PROCESSORS) install-lib_d
	make subBundle

sa subAll: $(NEEDS)
	make scrape
	$(BJAM) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single  -j $(COMPILE_PROCESSORS) install-lib_d install-lib
	make subBundle

san: $(NEEDS)
	$(BJAM) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single release install -j $(COMPILE_PROCESSORS) -a -n >_bjam.txt

sb subBundle:
	cp libs/[a-z]* $(CANDO_APP_LIB_DIR)
bjd:
	$(BJAM) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single debug -j $(COMPILE_PROCESSORS)

subMeta sm:
	make m4
	make mi
	make cscrape

initClasses.inc : *.h openmmExpose.scrape.inc
	make mi

openmmExpose.scrape.inc m4: openmmExpose.m4 exposeMacros.m4
	m4 -s openmmExpose.m4 >openmmExpose.scrape.inc



relink:
	rm libs/[a-z]*
	make sa

subClean sc:
	rm -rf bin/*

totalClean:
	rm -rf bin install_release



sub-prebuild:
	make scrape-symbols
	make scrape-classes

