
#-----------------------------------------------------------------------------
UP_CC = g++

BIN = $(PWD)/../../bin

COPT =-g

CINC = -I. -I $(PWD)/../PPP -I $(PWD)/../include

CLIB = -L $(PWD)/../lib -largtable2 -lfftw3

OBJS = gwlCreateAxis.o gwlSignalGen.o gwlCft.o gwlCwt.o gwlConvert.o\
      gwlIwt.o gwlDiffeoLin.o gwlSignalSum.o gwlSignalRead.o gwlCwtMaxLine.o\
      gwlET2D.o gwlET2DFilter.o gwlET3D.o gwlET3DFilter.o gwlDispModel.o\
      gwlDiffeoDisp.o gwlAutoCorr.o gwlTransFK.o gwlOptiSP.o gwlOptiSI.o\
      gwlSignalFilter.o gwlNNpred.o gwlWavelets.o

EXE = $(OBJS:.o=)

#-----------------------------------------------------------------------------
# clear suffix list and set new one
.SUFFIXES:
.SUFFIXES: .cpp .o

#-----------------------------------------------------------------------------
all : shell qwtplot installshell installqwtplot

shell : $(OBJS) $(EXE)

installshell :
	cp -f $(EXE) $(BIN)

qwtplot : gwlPlot.o gwlPlot

installqwtplot :
	cp -f gwlPlot $(BIN)

gwlPlot.o : gwlPlot.cpp
	$(QTDIR)/bin/moc gwlPlot.cpp -o gwlPlot_moc.cpp;
	$(UP_CC) -c -g -Wno-deprecated gwlPlot.cpp $(CINC) -I $(QTDIR)/include

gwlPlot : gwlPlot.o
	$(UP_CC) -g  $@.o $(COPT) -L $(PWD)/../lib -L $(QTDIR)/lib -largtable2 -lfftw3 -lqwt -o $@

.cpp.o :
	$(UP_CC) -c -Wno-deprecated $< $(CINC)

$(EXE) : $(OBJS)
	$(UP_CC)  $@.o $(COPT) $(CLIB) -o $@

clear :
	rm -f $(OBJS) gwlPlot.o

clean : clear
	rm -f $(EXE) gwlPlot

dependencies.make:
	touch dependencies.make

dep:
	$(UP_CC) $(COPT) $(CINC) -MM *.cpp > dependencies.make
