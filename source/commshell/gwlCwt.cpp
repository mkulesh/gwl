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

#include "gwlMainObject.h"

typedef gwlMainObject< PPPSignalContainer<double>, PPPSpectrContainer<PPPcomplex> > gwlMainDouble;
typedef gwlMainObject< PPPSignalContainer<PPPcomplex>, PPPSpectrContainer<PPPcomplex> > gwlMainCmpl;

template<class AType, class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_int        o_wttype;
    UTOption_file       o_freq;
    UTOption_dbl        o_cutoff;
    UTOption_file       o_trace;
    UTOption_int        o_points;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n", "name",    "<str>", "name of the spectrum (by default 'wavelet spectrum')", "wavelet spectrum"),
      o_wttype ("y", "wttype",  "<unsigned>", "type of the continuous wavelet transform: 0-time integral, 1-fft based, 2-multi-convolution method (by default 1)", 1),
      o_freq   ("f", "freq",    "<file>", "input file with frequency axis (by default 'frequency.dat')", "frequency.dat"),
      o_cutoff ("u", "cutoff",  "<real>", "wavelet cutoff for slow wavelet transform (by default 0.01)", 0.01),
      o_trace  ("a", "trace",   "<file>", "output file with calculation time trace (by default '')", ""),
      o_points ("s", "points",  "<unsigned>", "if given, reduce wavelet spectrum to 'points' length", 0)
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_progr);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(ATypeMain::o_iscmpl);
        ATypeMain::o_parser.add(o_wttype);
        ATypeMain::o_parser.add(o_freq);
        ATypeMain::o_parser.add(ATypeMain::o_wavelet);
        ATypeMain::o_parser.add(ATypeMain::o_wavpar);
        ATypeMain::o_parser.add(o_cutoff);
        ATypeMain::o_parser.add(o_trace);
        ATypeMain::o_parser.add(o_points);
        };

    void calc(void) {
      // prepare transform parameters
      PPPAxis aFreq;
      ATypeMain::read_bindata(aFreq, o_freq.getValue());
      ATypeMain::aSpectrumPar.initialize(
        (unsigned)o_wttype.getValue(),
        aFreq,
        ATypeMain::o_wavelet.getValue(),
        ATypeMain::o_wavpar.getValue(),
        "delta",
        1.0,
        o_cutoff.getValue());
      // calculation of wavelet spectrum
      PPPTransWavelet<AType> trans;
      trans.setObjectName(o_name.getValue());
      trans.setShowProgress(ATypeMain::o_progr.getValue());
      if(o_trace.isOptionGiven())
        trans.setTraceName(o_trace.getValue());
      trans.WT(ATypeMain::aDest, ATypeMain::aSource, ATypeMain::aSpectrumPar);
      ATypeMain::aDest.setObjectName(o_name.getValue());
      if(o_points.isOptionGiven() && o_points.getValue() != 0)
        {
        PPPSpectrContainer<PPPcomplex> aNewSpec;
        unsigned aNewOffSet = ATypeMain::aDest.points()/o_points.getValue();
        aNewSpec.assign(ATypeMain::aDest, 0, ATypeMain::aDest.voices(), 1, 0, o_points.getValue(), aNewOffSet);
        ATypeMain::aDest.assign(aNewSpec);
        }
      strstream str2;
      str2 << ATypeMain::aDest.getInfo() << ends;
      PPPBaseObject :: onMessage(str2.str());
      };

  };  // end of object

  
main(int  argc, char **argv)
  {
  gwlMain<double,gwlMainDouble> WT1("Continuous wavelet transform","gwlCwt","signal.dat","cwt.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  if(WT1.isComplex())
    {
    gwlMain<PPPcomplex,gwlMainCmpl> WT2("Continuous wavelet transform","gwlCwt","signal.dat","cwt.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  else WT1.evaluate();
  return 0;
  }




