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
#include "gwlMainObject.h"

typedef gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSignalContainer<double> > gwlMainDouble;

template<class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_type;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      gwlMainDouble(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n",  "name",    "<str>",  "name of the maximul lines (by default 'maximum lines')", "maximum lines"),
      o_type   ("y",  "type",    "<str>",  "type of maximum line calculation: amaxt, amaxf (by default amaxt)", "amaxt")
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_type);
        };

    void calc(void) {
      PPPTransWavelet<PPPcomplex> trans;
      // prepare calculation parameters
      map < string, PPPTransWavelet<PPPcomplex>::TimeFreqType > typesMap;
      typesMap["amaxt"] = PPPTransWavelet<PPPcomplex>::TFTtime;
      typesMap["amaxf"] = PPPTransWavelet<PPPcomplex>::TFTfreq;
      if(typesMap.find(o_type.getValue()) == typesMap.end())
        PPPBaseObject :: onError("type of maximum line calculation is incorrect");
      PPPAxis::AxisSign axsign = ATypeMain::aSpectrumPar.getFreq().getSign();
      // full spectrum
      if(axsign == PPPAxis::ASfull)
        {
        PPPSpectrContainer<PPPcomplex> aWgPlus,aWgMinus;
        PPPSignalContainer<double> aRplus,aRminus;
        trans.WaveletSeparate(aWgPlus,aWgMinus,ATypeMain::aSource);
        // positive ridge
        trans.Ridge(aRplus,aWgPlus,typesMap[o_type.getValue()],PPPAxis::ASplus);
        // negative ridge
        trans.Ridge(aRminus,aWgMinus,typesMap[o_type.getValue()],PPPAxis::ASminus);
        // combine positive and negative ridge
        ATypeMain::aDest.prepare(aRplus.points(),2,aRplus.getAxis(),aRplus.getObjectName());
        ATypeMain::aDest.getChannel(0).assign(aRplus.getChannel(0));
        ATypeMain::aDest.getChannel(1).assign(aRminus.getChannel(0));
        }
      else
        {
        trans.Ridge(ATypeMain::aDest,ATypeMain::aSource,typesMap[o_type.getValue()],axsign);
        }
      // information  
      ATypeMain::aDest.setObjectName(o_name.getValue());
      strstream str;
      str << ATypeMain::aDest.getInfo() << ends;
      PPPBaseObject :: onMessage(str.str());
      };

  };  // end of object


main(int  argc, char **argv)
  {
  gwlMain<gwlMainDouble> WT1("Maximum lines calculation","gwlCwtMaxLine","spectr.dat","maxline.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  WT1.evaluate();
  return 0;
  }




