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

typedef gwlMainObject< PPPAxis, PPPSignalContainer<double> > gwlMainDouble;

template<class AType, class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_dbl        o_freq;
    UTOption_dbl        o_time;
    UTOption_lit        o_four;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_freq   ("f", "freq",  "<real>", "frequency of the wavelet (by default 1.0)", 1.0),
      o_time   ("e", "time",  "<real>", "time position of the wavelet (by default 0.0)", 0.0),
      o_four   ("u", "four",  "if set, calculate and add Fourier spectrum to the wavelet representation", false)
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(ATypeMain::o_iscmpl);
        ATypeMain::o_parser.add(ATypeMain::o_wavelet);
        ATypeMain::o_parser.add(ATypeMain::o_wavpar);
        ATypeMain::o_parser.add(o_freq);
        ATypeMain::o_parser.add(o_time);
        ATypeMain::o_parser.add(o_four);
        };

    void calc(void) {
      // Wavelet
      PPPWavelet tmpwav;
      PPPWavelet *wavelet;
      PPPWavelet::WaveletType waveletType = tmpwav.strToWavelet(ATypeMain::o_wavelet.getValue());
      PPPSpectrParams sp;
      wavelet = sp.createWavelet(waveletType,ATypeMain::o_wavpar.getValue());
      wavelet->setPosition(o_time.getValue());
      wavelet->setFrequency(o_freq.getValue());
      // ceate of seismogramm
      PPPSignalContainer<PPPcomplex> sig, four;
      wavelet->evalTimeRepresentation(sig, ATypeMain::aSource);
      unsigned chan = (ATypeMain::isComplex())? 2:1;
      if(o_four.isOptionGiven()) chan+=3;
      ATypeMain::aDest.prepare(sig.points(), chan, sig.getAxis(), sig.getObjectName());
      for(unsigned i=0; i<ATypeMain::aSource.size(); i++)
        {
        ATypeMain::aDest(i,0) = real(sig(i));
        if(ATypeMain::isComplex()) ATypeMain::aDest(i,1) = imag(sig(i));
        }
      if(o_four.isOptionGiven())
        {
        double smpl = ATypeMain::aSource.getSamplingFreq();
        PPPAxis freq(ATypeMain::aSource.size(), -smpl/2.0, smpl/2.0, PPPAxis::ATlin, PPPSPECTRCONTAINER_FREQ);
        wavelet->evalFreqRepresentation(four, freq);
        for(unsigned i=0; i<ATypeMain::aSource.size(); i++)
          {
          ATypeMain::aDest(i,chan-3) = four.getAxis(i);
          ATypeMain::aDest(i,chan-2) = real(four(i));
          ATypeMain::aDest(i,chan-1) = imag(four(i));
          }
        }
      // Information
      strstream str;
      str << wavelet->getInfo() << endl << ATypeMain::aDest.getInfo() << ends;
      PPPBaseObject :: onMessage(str.str());
      delete wavelet;
      };

  };  // end of object


main(int  argc, char **argv)
  {
  gwlMain<double,gwlMainDouble> Wav("Generation of a wavelet representation","gwlWavelets","time.dat","wavelet.dat");
  Wav.parse(argc,argv);
  ConApplication.onMessage(ConApplication.getAppName());
  Wav.evaluate();
  return 0;
  }

