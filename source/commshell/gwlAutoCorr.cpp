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

typedef gwlMainObject< PPPSignalContainer<double>, PPPSignalContainer<double> > gwlMainDouble;
typedef gwlMainObject< PPPSignalContainer<PPPcomplex>, PPPSignalContainer<PPPcomplex> > gwlMainCmpl;

template<class AType,class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_chan1;
    UTOption_str        o_chan2;

  public:
    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n", "name",    "<str>", "name of the output signal (by default 'signal')", "signal"),
      o_chan1  ("1", "chan1",   "<str>",  "first set of channels of input spectrum (by default '')", ""),
      o_chan2  ("2", "chan2",   "<str>",  "second set of channels of input spectrum (by default '')", "")
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_chan1);
        ATypeMain::o_parser.add(o_chan2);
        };

    void calc(void) {
      PPPTransFour<AType> aTrans;
      PPPSignalContainer<PPPcomplex> aFour0;
      if(o_chan1.isOptionGiven() || o_chan2.isOptionGiven())
        {
        PPPVectorContainer<unsigned> chan1,chan2;
        chan1.strToVector("{"+string(o_chan1.getValue())+"}");
        chan2.strToVector("{"+string(o_chan2.getValue())+"}");
        if(chan1.size() != chan2.size())
          PPPBaseObject :: onError("first and second set of channels must have the same length");
        PPPMatrixContainer<unsigned> chan(chan1.size(),2);
        for(unsigned i=0; i<chan1.size(); i++)
          {
          chan(i,0) = chan1[i];
          chan(i,1) = chan2[i];
          }
        aTrans.CCS(aFour0, ATypeMain::aSource, chan);
        }
      else
        {
        aTrans.CCF(aFour0, ATypeMain::aSource, ATypeMain::aSource);
        }
      double aMin = ATypeMain::aSource.getAxis().getMin();
      double aMax = ATypeMain::aSource.getAxis().getMax();
      double aDelta = ATypeMain::aSource.getAxis().getSamplingPeriod();
      aTrans.ICCF(ATypeMain::aDest,aFour0,aMin,2.0*(aMax-aMin)+aDelta);
      // information
      strstream str;
      str << ATypeMain::aDest.getInfo() << ends;
      PPPBaseObject :: onMessage(str.str());
      };

  };  // end of object

  
main(int  argc, char **argv)
  {
  gwlMain<double,gwlMainDouble> WT1("Autocorrelaton of seismogram","gwlAutoCorr","signal.dat","autocorr.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  PPPObjectIO::ObjectType head = WT1.read_binheader();
  if(head == PPPObjectIO::SIGD)
    {
    WT1.evaluate();
    }
  else
    {
    gwlMain<PPPcomplex,gwlMainCmpl> WT2("Autocorrelaton of seismogram","gwlAutoCorr","signal.dat","autocorr.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  return 0;
  }




