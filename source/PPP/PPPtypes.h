/* 
 * This file is a part of GWL - Geophysical Wavelet Library
 * Copyright (C) 2007 Mikhail Kulesh and Matthias Holschneider
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * For more information please visit: http://users.math.uni-potsdam.de/~gwl
 * Email: mkulesh@math.uni-potsdam.de
 * ICQ: 103-405-403
 */

#ifndef _PPPTYPES
#define _PPPTYPES

/************************************************************************
 * PPP library settings. If necessary, these variables should be defined
 * before #include "PPPtypes.h"
 ************************************************************************/
// #define PPPCONF_CHECKINDEX true      // index control in all container objects
// #define PPPCONF_USEPARSING true      // using command line parsing (argtable2 external library)
// #define PPPCONF_USEFFTW3 true        // using external FFTW3 library for FFT calculation
// #define PPPCONF_USEMAPM true         // using external MAPM library
// #define PPPCONF_USEOPTIMIZATION true // include optimization objects
// #define PPPCONF_USEWAVESSOL true     // include waves solutions objects
// #define PPPCONF_USEQWT true          // include QWT and related graphical objects

/************************************************************************
 * PPP library implementation
 ************************************************************************/
// External parsing library
#ifdef PPPCONF_USEPARSING
  #include "UTparsing.h"
#endif

// PPP library version
#ifdef PPPCONF_CHECKINDEX
  #define PPPCONF_VERSION "ver. 1.5 (check index activated)"
#else
  #define PPPCONF_VERSION "ver. 1.5"
#endif

// Basic errors for functions arguments
#define ARG_VALUE       "invalid value of input parameter in procedure: "
#define ARG_POW2        "argument is not power of two in procedure: "

// Basic errors for virtual methods
#define VIRT_NOTDEF     "virtual method is not defined for this object: "
#define VIRT_NOTDEF32   "function not defined in Win32 application: "

// Basic file input-output messages
#define FILE_READASC    "reading ASC file: "
#define FILE_READSTA    "reading STA file: "
#define FILE_READWAV    "reading WAV file: "
#define FILE_WRITEASC   "writing ASC file: "
#define FILE_WRITEBIN   "writing binary file: "
#define FILE_ERROPEN    "cannot open file: "
#define FILE_ERRFORM    "version or format of input file is incorrect in procedure: "
#define FILE_ERREXP     "The procedure expects file virsion: "

// Basic memory errors
#define MEM_ERRINDEX    "out of index"
#define MEM_ERRALLOC    "allocation memory problem in procedure: "

// This is the section with incudes of standart objects library
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <complex>
#include <functional>
#include <stdlib.h>
#include <time.h>
#include <strstream>
using namespace std;

// Some small helpfull objects
template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

template < class A, class Comp > class PPPCompareLess {
  private:
    Comp _comp;
  public:
    PPPCompareLess ( Comp const & comp ): _comp ( comp ) {};
    bool inline operator ()(A const & a, A const & b ) const {
      return ( _comp(a)  < _comp(b) );
      };
  };

// External MAPM library
#ifdef PPPCONF_USEMAPM
  #include "m_apm.h"
  #include "PPPmapm.h"
  #include "PPPmcomplex.h"
#endif
  
// PPP objects
#include "PPPappliccon.h"
PPPApplicationCon ConApplication;

#include "PPPbaseobj.h"
#include "PPPcomplex.h"
class PPPellipse2D;
class PPPellipse3D;
template<class T> class PPPNullTransform;
#include "PPPbasetempl.h"
#include "PPPnulltrans.h"

template<class AType> class PPPMatrixContainer;
#include "PPPmfunc.h"
#include "PPPvector.h"
#include "PPPaxis.h"
#include "PPPmatrix.h"
#include "PPPsignal.h"
#include "PPPspectr.h"

#include "PPPwavelets.h"
#include "wavelets/PPPmorlet.h"
#include "wavelets/PPPmorletre.h"
#include "wavelets/PPPcauchy.h"
#include "wavelets/PPPshanon.h"
#include "wavelets/PPPhaar.h"
#include "wavelets/PPPdelta.h"
#include "PPPspectrpar.h"

#include "PPPlinalg.h"
#include "PPPapprox.h"
#include "approx/PPPapproxbispl.h"
#include "approx/PPPapproxcole.h"
#include "approx/PPPapproxgauss.h"
#include "approx/PPPapproxpol.h"
#include "approx/PPPapproxvel.h"
#include "approx/PPPapproxtwog.h"
#include "PPPdispModel.h"

#include "polarization/PPPellipse.h"
#include <polarization/PPPellipse2D.h>
#include "polarization/PPPellipse2Dfilter.h"
#include <polarization/PPPellipse3D.h>
#include "polarization/PPPellipse3Dfilter.h"

#include "PPPObjectIO.h"

#include "UTFFT.h"
#ifdef PPPCONF_USEFFTW3
  #include "UTFFTW3.cpp"
#else
  #include "UTFFTnr.cpp"
#endif

#include "PPPfft.h"
#include "PPPconvolution.h"
#include "PPPtransfour.h"
#include "PPPtranswav.h"
#include "PPPtranselli.h"
#include "PPPtransFK.h"

#include "PPPpropagDiss.h"
#include "PPPpropagLin.h"

#include "PPPwav.h"
#include "PPPsignalgen.h"
#include "PPPsignalread.h"
#include "PPPObjectIO.h"

#ifdef PPPCONF_USEOPTIMIZATION
  #include "PPPopti.h"
  #include "PPPoptiLM.h"
  #include "PPPOptiSP.h"
  #include "PPPOptiSI.h"
  #include "PPPOptiMult.h"
#endif

#ifdef PPPCONF_USEWAVESSOL
  #include "PPPinifile.h"
  #include "PPPwaves.h"
  #include "PPPwavesSpec.h"
  #include "PPPwavesClass.h"
  #include "PPPwavesCoss.h"
  #include "PPPwavesRedc.h"
#endif

#ifdef PPPCONF_USEQWT
  #include <qapplication.h>
  #include <qmainwindow.h>
  #include <qpopupmenu.h>
  #include <qstatusbar.h>
  #include <qmenubar.h>
  #include <qmessagebox.h>
  #include "qwt_plot_curve.h"
  #include "qwt_plot_spectrogram.h"
  #include "qwt_plot_picker.h"
  #include "qwt_scale_engine.h"
  #include "UTPlotRasterData.h"
  #include "UTMatrixOfQwtPlot.h"
  #include "PPPsignalplot.h"
#endif

#endif
