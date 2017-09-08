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
    UTOption_file       o_spec;
    UTOption_dbl        o_dist;
    UTOption_lit        o_isacorr;
    UTOption_int        o_cmpl;
    UTOption_dbl        o_eps;
    UTOption_file       o_ospec;
    UTOption_str        o_chan;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n",  "name",    "<str>", "name of the inverse signal (by default 'optimized dispersion model')", "optimized dispersion model"),
      o_spec   ("s",  "spec",    "<file>", "input spectrum's file name (by default 'cwt.dat')", "cwt.dat"),
      o_ospec  ("z",  "ospec",   "<file>", "optimized spectrum's file name (by default 'cwtopt.dat')", "cwtopt.dat"),
      o_dist   ("d",  "dist",    "<real>", "propagation distance (by default 0)", 0.0),
      o_isacorr("a",  "acorr",   "if set, the input signal/spectrum is an autocorrelation", false),
      o_cmpl   ("p",  "cmpl",    "<unsigned>", "type of optimized component: 3-abs, 4-arg (by default 3)", 3),
      o_eps    ("e",  "eps",     "<real>", "maximal precision (by default 1e-08)", 1E-08),
      o_chan   ("h",  "chan",    "<str>",  "channels of input spectrum which will be optimized (by default '0,1')", "0,1")
        {
        ATypeMain :: o_parser.add(ATypeMain :: o_nomess);
        ATypeMain :: o_parser.add(ATypeMain :: o_progr);
        ATypeMain :: o_parser.add(ATypeMain :: o_infile);
        ATypeMain :: o_parser.add(ATypeMain :: o_outfile);
        ATypeMain :: o_parser.add(ATypeMain :: o_outtype);
        ATypeMain :: o_parser.add(o_name);
        ATypeMain :: o_parser.add(o_spec);
        ATypeMain :: o_parser.add(o_ospec);
        ATypeMain :: o_parser.add(o_dist);
        ATypeMain :: o_parser.add(o_isacorr);
        ATypeMain :: o_parser.add(o_cmpl);
        ATypeMain :: o_parser.add(o_eps);
        ATypeMain :: o_parser.add(o_chan);
        };

    inline PPPSpectrParams &getSpectrParams(void) {
      return ATypeMain :: aSpectrumPar;
      };

    void calc(void) {
      // channels
      PPPVectorContainer<unsigned> chan;
      chan.strToVector("{"+string(o_chan.getValue())+"}");
      // wavelet spectrum
      PPPSpectrContainer<PPPcomplex> aSpectr;
      ATypeMain :: read_bindata(aSpectr, o_spec.getValue(), true);
      // information
      strstream str1;
      if(o_chan.isOptionGiven())
        str1 << aSpectr.getInfo() << endl << "  channels: " << chan.vectorToStr() << ends;
      else
        str1 << aSpectr.getInfo() << ends;
      ATypeMain :: onMessage(str1.str());
      // optimization parameters
      PPPPropagatorDiss::PropType aPropType = (o_isacorr.isOptionGiven())? PPPPropagatorDiss::STwavcross : PPPPropagatorDiss::STwav;
      PPPVectorContainer<double> aParDouble(1);
      aParDouble[0] = o_eps.getValue();
      // optimization
      strstream str2;
      if(o_chan.isOptionGiven())
        {
        PPPOptiSpectrum<PPPOptiLevenbergMarq> opti;
        opti.prepare(aSpectr,ATypeMain :: aSource,o_dist.getValue(),aPropType,&getSpectrParams(),o_cmpl.getValue(),chan[0],chan[1]);
        opti.setParams(aParDouble);
        opti.setShowProgress(ATypeMain :: o_progr.getValue());
        opti.optimization();
        str2 << opti.getInfo();
        ATypeMain :: aDest.assign(opti.getModel());
        if(o_ospec.isOptionGiven())
          ATypeMain :: write_bindata(opti.getResult(), o_ospec.getValue());
        }
      else
        {
        PPPOptiMult< PPPOptiSpectrum<PPPOptiLevenbergMarq> > opti;
        opti.prepare(aSpectr,ATypeMain :: aSource,o_dist.getValue(),aPropType,&getSpectrParams(),o_cmpl.getValue());
        opti.setParams(aParDouble);
        opti.setShowProgress(ATypeMain :: o_progr.getValue());
        opti.optimization();
        str2 << opti.getInfo();
        ATypeMain :: aDest.assign(opti.getModel());
        if(o_ospec.isOptionGiven())
          ATypeMain :: write_bindata(opti.getResultSpec(), o_ospec.getValue());
        }
      ATypeMain :: aDest.evalModel(ATypeMain :: aSource.getFreq());
      ATypeMain :: aDest.setObjectName(o_name.getValue());
      str2 << endl << ATypeMain :: aDest.getInfo() << ends;
      ATypeMain :: onMessage(str2.str());
      };

  };  // end of object


main(int  argc, char **argv)
  {
  gwlMain<gwlMainDouble> WT1("Spectrum optimization","gwlOptiSP","dispmodel.dat","dispopti.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  WT1.evaluate();
  return 0;
  }




