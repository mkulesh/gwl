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

#define gwlET2D_ERR03  "type of polarization filter is undefined: "
#define gwlET2D_ERR04  "number of parameter of polarization filter is invalid"

typedef gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<PPPcomplex> > gwlMainSpecToSpec;
typedef gwlMainObject< PPPSignalContainer<double>, PPPSignalContainer<double> > gwlMainSigToSig;

template<class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_filter;
    UTOption_file       o_elli;
    PPPTransElli        aTrans;
    vector<PPPEllipse3Dfilter> aFilters;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n", "name",    "<str>",  "name of the elliptic parameters (by default 'filtered spectrum')", "filtered spectrum"),
      o_filter ("f", "filter",  "<str>",  "polarization filter parameters (by default '')", ""),
      o_elli   ("e", "elli",    "<file>", "input file elliptical parameters (by default 'elli.dat')", "elli.dat")
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_filter);
        ATypeMain::o_parser.add(o_elli);
        };

    void parsePolarizationFilters() {
      map < string, PPPEllipse3Dfilter::FilterType > typesMap;
      typesMap["norm"] = PPPEllipse3Dfilter::NRFNorm;
      typesMap["lin"] = PPPEllipse3Dfilter::NRFlin;
      typesMap["linxy"] = PPPEllipse3Dfilter::NRFLinXY;
      typesMap["linxz"] = PPPEllipse3Dfilter::NRFLinXZ;
      typesMap["linyz"] = PPPEllipse3Dfilter::NRFLinYZ;
      typesMap["elli"] = PPPEllipse3Dfilter::NRFelli;
      typesMap["ellixy"] = PPPEllipse3Dfilter::NRFElliXY;
      typesMap["ellixz"] = PPPEllipse3Dfilter::NRFElliXZ;
      typesMap["elliyz"] = PPPEllipse3Dfilter::NRFElliYZ;
      typesMap["prog"] = PPPEllipse3Dfilter::NRFprog;
      typesMap["retro"] = PPPEllipse3Dfilter::NRFretro;
      typesMap["azimuth"] = PPPEllipse3Dfilter::NRFazimuth;
      typesMap["dip"] = PPPEllipse3Dfilter::NRFdip;
      string aStr = o_filter.getValue();
      PPPVectorContainer<double> aPars;
      while(1)
        {
        unsigned pos = aStr.find(",");
        if(pos >= aStr.size()) { break; }
        string strname = aStr.substr(0,pos);
        if(typesMap.find(strname) == typesMap.end())
          PPPBaseObject :: onError(gwlET2D_ERR03+strname);
        aStr.erase(0,pos+1);
        if(typesMap[strname] != PPPEllipse3Dfilter::NRFNorm)
          {
          string strpar = "";
          for(unsigned i=0;i<2; i++)
            {
            pos = aStr.find(",");
            if(pos >= aStr.size()) { strpar.append(aStr); break; }
            strpar.append(aStr.substr(0,pos));
            if(i==0) strpar.append(",");
            aStr.erase(0,pos+1);
            }
          aPars.strToVector("{"+strpar+"}");
          }
        else
          aPars.strToVector("{0,0}");
        if(aPars.size() != 2)
          ATypeMain::onError(gwlET2D_ERR04);
        aFilters.push_back(PPPEllipse3Dfilter(typesMap[strname],aPars[0],aPars[1],aPars[0],aPars[1]));
        }
      };

    void calc(void) {
      // prepare transform parameters
      if(o_filter.isOptionGiven())
        parsePolarizationFilters();
      // calculation of elliptic properties
      calc(ATypeMain::aDest, ATypeMain::aSource);
      };

    void calc(PPPSignalContainer<double> &aLocDest, PPPSignalContainer<double> &aLocSource) {
      PPPSignalContainer<PPPellipse3D> aElli;
      ATypeMain::read_bindata(aElli, o_elli.getValue());
      aTrans.eval3DFilter(aLocDest, aLocSource, aElli);
      aLocDest.setObjectName(o_name.getValue());
      };

    void calc(PPPSpectrContainer<PPPcomplex> &aLocDest, PPPSpectrContainer<PPPcomplex> &aLocSource) {
      if(aFilters.size()==0) return;
      aLocDest.prepare(
        aLocSource.voices(),
        aLocSource.points(),
        3*aFilters.size(),
        aLocSource.getTime(),
        aLocSource.getFreq(),
        "filtered spectrum");
      // calculation
      PPPSpectrContainer<PPPellipse3D> aElli;
      ATypeMain::read_bindata(aElli, o_elli.getValue());
      for(unsigned k=0; k<aFilters.size(); k++)
        aTrans.eval3DFilter(aLocDest, aLocSource, aElli, aFilters[k], k);
      aLocDest.setObjectName(o_name.getValue());
      // information
      strstream str;
      str << aLocDest.getInfo() << ends;
      ATypeMain::onMessage(str.str());
      };

  };  // end of object


main(int  argc, char **argv)
  {
  string appname = "3D polarization filter";
  string appcode = "gwlET3DFilter";
  gwlMain<gwlMainSpecToSpec> WT1(appname.c_str(),appcode.c_str(),"source.dat","filtered.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  PPPObjectIO::ObjectType head = WT1.read_binheader();
  if(head == PPPObjectIO::SPECC)
    {
    WT1.evaluate();
    }
  else if(head == PPPObjectIO::SIGD)
    {
    gwlMain<gwlMainSigToSig> WT2(appname.c_str(),appcode.c_str(),"source.dat","filtered.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  else
    ConApplication.onError(gwlMain_ERR02);
  return 0;
  }




