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

typedef gwlMainObject< PPPAxis, PPPSignalContainer<double> > gwlMainDouble;
typedef gwlMainObject< PPPAxis, PPPSignalContainer<PPPcomplex> > gwlMainCmpl;

template<class AType, class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_par;
    UTOption_str        o_type;
    UTOption_str        o_modul;
    UTOption_dbl        o_noise;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n", "name",   "<str>",  "name of the signal (by default 'signal')", "signal"),
      o_type   ("y", "type",   "<str>",  "type of synthetic signal (by default 'harmon')", "harmon"),
      o_par    ("p", "par",    "<str>",  "parameters for signal generation (by default '')", ""),
      o_modul  ("d", "modul",  "<str>",  "parameters for windowed modulation of signal (by default '')", ""),
      o_noise  ("s", "noise",  "<real>", "value of random noise (by default 0)", 0.0)
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(ATypeMain::o_iscmpl);
        ATypeMain::o_parser.add(o_type);
        ATypeMain::o_parser.add(o_par);
        ATypeMain::o_parser.add(o_modul);
        ATypeMain::o_parser.add(o_noise);
        };

    void calc(void) {
      // prepare generation parameters
      map < string, unsigned > typesMap;
      typesMap["zero"] = 0;
      typesMap["delta"] = 1;
      typesMap["harmon"] = 2;
      typesMap["harmrot"] = 3;
      typesMap["harmphase"] = 4;
      typesMap["rickdiss"] = 5;
      if(typesMap.find(o_type.getValue()) == typesMap.end())
        PPPBaseObject :: onError(PPPSIGNALGEN_ERRTYPE+string("calc()"));

      PPPSignalGen<AType> aPar;
      aPar.getParams().strToVector("{"+string(o_par.getValue())+"}");
      // Generation of signal
      ATypeMain::aDest.prepare(ATypeMain::aSource.size(), 1, ATypeMain::aSource, o_name.getValue());
      aPar.setShowNotation(PPPBaseObject::NMappend);
      switch(typesMap[o_type.getValue()])
        {
        case 0:  aPar.genZero(ATypeMain::aDest);          break;
        case 1:  aPar.genDelta(ATypeMain::aDest);         break;
        case 2:  aPar.genHarmon(ATypeMain::aDest);        break;
        case 3:  aPar.genHarmonRot(ATypeMain::aDest);     break;
        case 4:  aPar.genHarmonPhase(ATypeMain::aDest);   break;
        case 5:  aPar.genRickerDissip(ATypeMain::aDest);  break;
        };
      // random noise
      if(o_noise.isOptionGiven())
        {
        aPar.getParams().resize(1);
        aPar.getParams()[0] = o_noise.getValue();
        aPar.onNotation("\n");
        aPar.addRandomNoise(ATypeMain::aDest);
        };
      // Modulation of signal
      if(o_modul.isOptionGiven())
        {
        aPar.getParams().strToVector("{"+string(o_modul.getValue())+"}");
        aPar.onNotation("\n");
        aPar.evalWindowedModulation(ATypeMain::aDest);
        };
      // Information
      ATypeMain::aDest.setObjectName(o_name.getValue());
      strstream str;
      str << ATypeMain::aDest.getInfo() << endl << "  " << aPar.getNotation() << ends;
      PPPBaseObject :: onMessage(str.str());
      };

  };  // end of object

  
main(int  argc, char **argv)
  {
  gwlMain<double,gwlMainDouble> WT1("Generation of a synthetic signal","gwlSignalGen","time.dat","signal.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  if(WT1.isComplex())
    {
    gwlMain<PPPcomplex,gwlMainCmpl> WT2("Generation of a synthetic signal","gwlSignalGen","time.dat","signal.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  else WT1.evaluate();
  return 0;
  }

