SHELL = /bin/bash

#-----------------------------------------------------------------------------
ARGTABLEDIR = argtable2-6
FFTW3DIR = fftw-3.1.2
QWTDIR = qwt-5.0.1
BIN = ${PWD}/../bin
#-----------------------------------------------------------------------------

.SUFFIXES: .c .cpp .o

all: shell qwtplot

shell: make_argtable2 make_fftw3 make_shell

qwtplot: make_qwt make_qwtplot

make_argtable2: .argtableisinstalled

.argtableisinstalled:
	cd ${ARGTABLEDIR}; ./configure --prefix ${PWD};\
	make; make install;\
	rm -rf ${PWD}/share ${PWD}/man;\
	touch ../.argtableisinstalled;
	
make_fftw3: .fftw3installed

.fftw3installed:
	cd ${FFTW3DIR}; ./configure --prefix $(PWD);\
	make; make install;\
	rm -rf .libs ${PWD}/share ${PWD}/bin;\
	touch ../.fftw3installed;

make_qwt: .qwtinstalled

#setenv QTDIR towhereitis; setenv PATH ${PATH}:${QTDIR}/bin

.qwtinstalled:
	OLDPATH=$(PATH);\
	export PATH=$(QTDIR)/bin:$(PATH);\
	cd ${QWTDIR};	qmake; make;\
	cp src/*.h $(PWD)/include;\
	cp lib/lib*.* $(PWD)/lib;\
	touch ../.qwtinstalled;\
	export PATH=$(OLDPATH);\

make_shell:
	cd commshell; make -f Makefile shell; make -f Makefile installshell;

make_qwtplot:
	cd commshell; make -f Makefile qwtplot; make -f Makefile installqwtplot;

dep:
	cd commshell; make -f Makefile dep;
	
clean:
	rm -f .argtableisinstalled .qwtinstalled .mapminstalled .fftw3installed .utilitiesinstalled;
	cd commshell; make -f Makefile clean;
	cd ${ARGTABLEDIR}; make clean;
	cd ${FFTW3DIR}; make clean;
#	cd ${QWTDIR}; make clean;
