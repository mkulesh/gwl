/*******************************************************************************
 * GWL - Geophysical Wavelet Library
 * *****************************************************************************
 * Copyright (C) 2002-2017 Mikhail Kulesh, Matthias Holschneider
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *  
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#define PPPCONF_USEOPTIMIZATION      // include optimization objects

#include "gwlMainObject.h"

typedef gwlMainObject< PPPDispersionModel, PPPDispersionModel > gwlMainDouble;

template<class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_file       o_sig;
    UTOption_file       o_osig;
    UTOption_dbl        o_dist;
    UTOption_lit        o_isacorr;
    UTOption_dbl        o_eps;
    UTOption_str        o_chan;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n",  "name",    "<str>", "name of the inverse signal (by default 'optimized dispersion model')", "optimized dispersion model"),
      o_sig    ("s",  "sig",     "<file>", "input signal's file name (by default 'signal.dat')", "signal.dat"),
      o_osig   ("z",  "osig",    "<file>", "optimized signal's file name (by default 'signalopt.dat')", "signalopt.dat"),
      o_dist   ("d",  "dist",    "<real>", "propagation distance (by default 0)", 0.0),
      o_isacorr("a",  "acorr",   "if set, the input signal/spectrum is an autocorrelation", false),
      o_eps    ("e",  "eps",     "<real>", "maximal precision (by default 1e-08)", 1E-08),
      o_chan   ("h",  "chan",    "<str>",  "channels of input spectrum which will be optimized (by default '0,1')", "0,1")
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_progr);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_sig);
        ATypeMain::o_parser.add(o_osig);
        ATypeMain::o_parser.add(o_dist);
        ATypeMain::o_parser.add(o_isacorr);
        ATypeMain::o_parser.add(o_eps);
        ATypeMain::o_parser.add(o_chan);
        };

    inline PPPSpectrParams &getSpectrParams(void) {
      return ATypeMain::aSpectrumPar;
      };

    void calc(void) {
      // channels
      PPPVectorContainer<unsigned> chan;
      chan.strToVector("{"+string(o_chan.getValue())+"}");
      // source signal
      PPPSignalContainer<double> aSignal;
      ATypeMain::read_bindata(aSignal, o_sig.getValue(), true);
      // information
      strstream str1;
      if(o_chan.isOptionGiven())
        str1 << aSignal.getInfo() << endl << "  channels: " << chan.vectorToStr() << ends;
      else
        str1 << aSignal.getInfo() << ends;
      ATypeMain::onMessage(str1.str());
      // optimization parameters
      PPPPropagatorDiss::PropType aPropType = (o_isacorr.isOptionGiven())? PPPPropagatorDiss::STfourcross : PPPPropagatorDiss::STfour;
      PPPVectorContainer<double> aParDouble(1);
      aParDouble[0] = o_eps.getValue();
      // optimization
      strstream str2;
      if(o_chan.isOptionGiven())
        {
        PPPOptiSignal<PPPOptiLevenbergMarq> opti;
        opti.prepare(aSignal,ATypeMain::aSource,o_dist.getValue(),aPropType,chan[0],chan[1]);
        opti.setParams(aParDouble);
        opti.setShowProgress(ATypeMain::o_progr.getValue());
        opti.optimization();
        str2 << opti.getInfo();
        ATypeMain::aDest.assign(opti.getModel());
        if(o_osig.isOptionGiven())
          ATypeMain::write_bindata(opti.getResult(), o_osig.getValue());
        }
      else
        {
        PPPOptiMult< PPPOptiSignal<PPPOptiLevenbergMarq> > opti;
        opti.prepare(aSignal,ATypeMain::aSource,o_dist.getValue(),aPropType);
        opti.setParams(aParDouble);
        opti.setShowProgress(ATypeMain::o_progr.getValue());
        opti.optimization();
        str2 << opti.getInfo();
        ATypeMain::aDest.assign(opti.getModel());
        if(o_osig.isOptionGiven())
          ATypeMain::write_bindata(opti.getResultSig(), o_osig.getValue());
        }
      ATypeMain::aDest.evalModel(ATypeMain::aSource.getFreq());
      ATypeMain::aDest.setObjectName(o_name.getValue());
      str2 << endl << ATypeMain::aDest.getInfo() << ends;
      ATypeMain::onMessage(str2.str());
      };

  };  // end of object


main(int  argc, char **argv)
  {
  gwlMain<gwlMainDouble> WT1("Signal optimization","gwlOptiSI","dispmodel.dat","dispopti.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  WT1.evaluate();
  return 0;
  }




