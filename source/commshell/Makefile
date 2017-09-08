#-------------------------------------------------------------------------------
# GWL - Geophysical Wavelet Library
#-------------------------------------------------------------------------------
# Copyright (C) 2002-2017 Mikhail Kulesh, Matthias Holschneider
#  
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

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
	$(UP_CC) -g  $@.o $(COPT) -L $(PWD)/../lib -L $(QTDIR)/lib -largtable2 -lfftw3 -lqwt -lqt-mt -o $@

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
