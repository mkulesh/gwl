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

typedef gwlMainObject< PPPSignalContainer<double>, PPPSignalContainer<double> > gwlMainDouble;
typedef gwlMainObject< PPPSignalContainer<PPPcomplex>, PPPSignalContainer<PPPcomplex> > gwlMainCmpl;

template<class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;

  public:
    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n", "name",    "<str>", "name of the output signal (by default 'signal')", "signal")
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        };

    void read(void) {
      ATypeMain::parseFileNames();
      // Generation of signal
      bool aSaveMessage = ATypeMain::getShowMessage();
      for(unsigned k=0; k<ATypeMain::aFileNames.size(); k++)
        {
        ATypeMain::o_infile.setValue(ATypeMain::aFileNames[k].c_str());
        ATypeMain::setShowMessage(false);
        ATypeMain::read();
        ATypeMain::setShowMessage(aSaveMessage);
        ATypeMain::onMessage("addition of file: "+ATypeMain::aFileNames[k]);
        if(k == 0)
          {
          ATypeMain::aDest.assign(ATypeMain::aSource);
          continue;
          }
        if(ATypeMain::aSource.channels() != ATypeMain::aDest.channels() || ATypeMain::aSource.points() != ATypeMain::aDest.points())
          ATypeMain::onError("invalid size of the source signal "+ATypeMain::aFileNames[k]);
        for(unsigned i=0; i<ATypeMain::aSource.points(); i++)
          for(unsigned j=0; j<ATypeMain::aSource.channels(); j++)
            ATypeMain::aDest(i,j) = ATypeMain::aDest(i,j)+ATypeMain::aSource(i,j);
        }
      // Information
      ATypeMain::aDest.setObjectName(o_name.getValue());
      strstream str;
      str << ATypeMain::aDest.getInfo() << ends;
      PPPBaseObject :: onMessage(str.str());
      };

    void calc(void) {
      };

  };  // end of object

  
main(int  argc, char **argv)
  {
  gwlMain<gwlMainDouble> WT1("Add the signals","gwlSignalSum","input.dat","signal.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  WT1.parseFileNames();
  PPPObjectIO::ObjectType head = WT1.read_binheader();
  if(head == PPPObjectIO::SIGD)
    {
    WT1.evaluate();
    }
  else
    {
    gwlMain<gwlMainCmpl> WT2("Add the signals","gwlSignalSum","input.dat","signal.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  return 0;
  }




