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

class gwlMain : public gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<PPPcomplex> >
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_params;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<PPPcomplex> >(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n",  "name",    "<str>", "name of deformed spectrum (by default 'wavelet spectrum after diffeomorphism')", "wavelet spectrum after diffeomorphism"),
      o_params ("d",  "diffpar", "<str>", "parameters of diffeomorphism (by default '')", "")
        {
        o_parser.add(o_nomess);
        o_parser.add(o_infile);
        o_parser.add(o_outfile);
        o_parser.add(o_outtype);
        o_parser.add(o_name);
        o_parser.add(o_params);
        o_parser.add(o_wavelet);
        o_parser.add(o_wavpar);
        };

    void calc(void) {
      // prepare transform parameters
      PPPVectorContainer<double> aPar;
      aPar.strToVector("{"+string(o_params.getValue())+"}");
      aSpectrumPar.setInverseParams(o_wavelet.getValue(), o_wavpar.getValue());
      // calculation of linear diffeomorphism
      PPPPropagatorLin DL;
      DL.getParams().assign(aPar);
      DL.evalWaveletPropag(aDest,aSource,aSpectrumPar);
      aDest.setObjectName(o_name.getValue());
      strstream str;
      str << DL.getInfo() << ends;
      onMessage(str.str());
      };

  };  // end of object


main(int  argc, char **argv)
  {
  gwlMain WT1("Linear diffeomorphism of the wavelet spectrum","gwlDiffeoLin","spectr.dat","diffeolin.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  WT1.evaluate();
  return 0;
  }




