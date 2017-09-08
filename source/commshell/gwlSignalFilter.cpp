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

template<class AType, class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_dbl        o_fmin;
    UTOption_dbl        o_fmax;

  public:
    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n",  "name",  "<str>",  "name of the output signal (by default 'filtered signal')", "filtered signal"),
      o_fmin   ("1",  "fmin",  "<real>", "minimal frequency of band-pass filter (by default 0)", 0.0),
      o_fmax   ("2",  "fmax",  "<real>", "maxima frequency of band-pass filter (by default 100)", 100.0)
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_fmin);
        ATypeMain::o_parser.add(o_fmax);
        };

    void calc(void) {
      PPPTransFour<AType> trans;
      trans.BandPassFilter(ATypeMain::aDest,ATypeMain::aSource,o_fmin.getValue(),o_fmax.getValue());
      ATypeMain::aDest.setObjectName(o_name.getValue());
      };

  };  // end of object

  
main(int  argc, char **argv)
  {
  gwlMain<double,gwlMainDouble> WT1("Filtering of the signal","gwlSignalFilter","signal.dat","filtered.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  PPPObjectIO::ObjectType head = WT1.read_binheader();
  if(head == PPPObjectIO::SIGD)
    {
    WT1.evaluate();
    }
  else
    {
    gwlMain<PPPcomplex,gwlMainCmpl> WT2("Filtering of the signal","gwlSignalFilter","signal.dat","filtered.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  return 0;
  }




