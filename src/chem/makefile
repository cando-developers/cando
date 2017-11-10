
#include	$(CLASP_HOME)/makefile.inc

PARSER_NEEDS = _msmarts_Parser.cc _gaff_Parser.cc 
M4_NEEDS = initClasses.inc
NEEDS = $(M4_NEEDS) $(PARSER_NEEDS)

msmarts_Parser.cc : msmarts_Parser.yy
	bison -p msmarts_ -v -d -o msmarts_Parser.cc msmarts_Parser.yy

gaff_Parser.cc : gaff_Parser.yy
	bison -p gaff_ -v -d -o gaff_Parser.cc gaff_Parser.yy

parsers: msmarts_Parser.cc gaff_Parser.cc



FLEX=flex




subAllDebug sad: $(NEEDS)
	make scrape
	$(BJAM) -j $(COMPILE_PROCESSORS) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single install_dbg 
	make sb


subAll sa: $(NEEDS)
	make scrape
	$(BJAM) -j $(COMPILE_PROCESSORS) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single install_dbg install_opt
	make sb

subBundle sb:
	install -c libs/lib* $(CANDO_APP_LIB_DIR)

#
# run bjam and print all the commands it would execute but don't evaluate them
#
san: $(NEEDS)
	$(BJAM) -j $(COMPILE_PROCESSORS) toolset=$(BUILD_TOOLSET) link=$(BUILD_LINK) threading=single install_dbg install_opt -a -n >_bjam.txt


subMeta sm:
	touch initClasses.inc
	make m4
	make mi
	make cscrape

commit:
	svn commit -m "commit"


subDox:
	make msmarts_Parser.dox


subClean:
	rm -rf bin/*
	rm -f *.tab.*
	rm -f lex.*

export CT=/usr/local/bin/ctags
tags: 	$(SOURCES) $(HEADERS)
	echo $(CT)
	echo *.cc *.cpp *.h
	$(CT) --language-force=c++ *.cc *.cpp *.h

initClasses.inc: *.h
	make mi

m4: $(M4_NEEDS)
	echo Made $(M4_NEEDS)

m4Templates.h: m4Templates.m4
	m4 --synclines m4Templates.m4 >m4Templates.h
predefinedSymbols.inc pm4: predefinedSymbols.m4
	m4 predefinedSymbols.m4 >predefinedSymbols.inc
svmInstructions.inc svmm4: svmInstructions.m4
	m4 svmInstructions.m4 >svmInstructions.inc


external_m4Templates.h em4 : external_m4Templates.m4
	m4 --synclines external_m4Templates.m4 >external_m4Templates.h

userDocs: userDocs.pdf
	install -c userDocs.pdf $(CANDO_HOME)/docs

editDocs: editDocs.pdf

docs:	userDocs.pdf editDocs.pdf

userDocsView: 
	make userDocs
	open -n userDocs.pdf


editDocs.dvi:	scrapedEditDocs.tex editDocs.tex docCommands.tex
	latex editDocs.tex
	latex editDocs.tex
	latex editDocs.tex
	makeindex editDocs
	latex editDocs.tex

editDocs.pdf: editDocs.dvi
	dvipdf editDocs.dvi


userDocs.dvi:	scrapedDocs.tex userDocs.tex docCommands.tex
	latex userDocs.tex
	latex userDocs.tex
	latex userDocs.tex
	makeindex userDocs
	latex userDocs.tex

userDocs.pdf: userDocs.dvi
	dvipdf userDocs.dvi

scrapedDocs.tex: $(SOURCES) msmarts_Parser.yy $(COMMON_DIR)/scrapeDocs.py *.cc
	python $(COMMON_DIR)/scrapeDocs.py scrapedDocs.tex scrapedEditDocs.tex *.cc msmarts_Parser.yy

scrapedEditDocs.tex: $(SOURCES) msmarts_Parser.yy $(COMMON_DIR)/scrapeDocs.py *.cc
	python $(COMMON_DIR)/scrapeDocs.py scrapedDocs.tex scrapedEditDocs.tex *.cc msmarts_Parser.yy

lex:
	make lex.lisp.cc

msmarts_Parser.dox:

bisonHelp:
	"$(BISON)" -h

winParsers:
	(cmd.exe /C bison -d -v -p lisp lispParse.yy)


moeForceField.tab.cc: moeForceField.yy
	$(BISON) -v `cygpath -a moeForceField.yy`

moeForceField.tab.hh: moeForceField.tab.cc



#
# Common parser generators
#
#


_lisp_PureLexer.cc : lisp_PureLexer.ll
	$(FLEX) -Plisp --nounistd --header-file=_lisp_PureLexer.h  -o_lisp_PureLexer.cc lisp_PureLexer.ll    

_lisp_Parser.cc : lisp_Parser.yy
	$(BISON) -p lisp -v --defines=_lisp_Parser.h -d -o _lisp_Parser.cc lisp_Parser.yy

_msmarts_Parser.cc : msmarts_Parser.yy
	echo MAKE SURE BISON IS UP TO DATE VERSION MUST BE higher than 2.5
	$(BISON) --version
	$(BISON) -p msmarts -v -d -o _msmarts_Parser.cc msmarts_Parser.yy

_gaff_Parser.cc : gaff_Parser.yy
	$(BISON) -p gaff -v -t -d -o _gaff_Parser.cc gaff_Parser.yy



sub-prebuild:
	make scrape-symbols
	make scrape-classes

