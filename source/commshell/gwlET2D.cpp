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

#define gwlET2D_ERR01  "type of the elliptic mashine does not correspont to type of input object"
#define gwlET2D_ERR02  "frequency interval for averaging's method is invalid"
#define gwlET2D_ERR03  "type of polarization filter is undefined: "
#define gwlET2D_ERR04  "number of parameter of polarization filter is invalid"

typedef gwlMainObject< PPPSignalContainer<double>, PPPSignalContainer<PPPellipse2D> > gwlMainSigToSigDbl;
typedef gwlMainObject< PPPSignalContainer<PPPcomplex>, PPPSignalContainer<PPPellipse2D> > gwlMainSigToSigCmpl;
typedef gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<PPPellipse2D> > gwlMainSpecToSpec;
typedef gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSignalContainer<PPPellipse2D> > gwlMainSpecToSig;
typedef gwlMainObject< PPPSpectrContainer<PPPellipse2D>, PPPSignalContainer<PPPellipse2D> > gwlMainElliToSig;

template<class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_type;
    UTOption_int        o_tw;
    UTOption_dbl        o_filter;
    UTOption_str        o_freq;
    UTOption_file       o_mline;
    UTOption_int        o_channel;
    PPPTransElli        aTrans;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n", "name",    "<str>", "name of the elliptic parameters (by default 'elliptic parameters')", "elliptic parameters"),
      o_type   ("y", "type",    "<str>", "type of the elliptic mashine: rene, morozov, scovar, acovar, complex, mlinet, mlinef (by default 'complex')", "complex"),
      o_tw     ("w", "tw",      "<unsigned>", "time window length for scovar and acovar methods (by default 2)", 2),
      o_filter ("f", "filter",  "<real>", "energy filter by spectral transform (by default 0)", 0.0),
      o_freq   ("q", "freq",    "<str>", "frequency interval for averaging's method (by default '')", ""),
      o_mline  ("d", "mline",   "<file>", "input file for maximum line's method (by default '')", ""),
      o_channel("c", "chan",    "<int>",  "channel of input spectrum (by default 0)", 0)
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_type);
        ATypeMain::o_parser.add(o_tw);
        ATypeMain::o_parser.add(o_filter);
        ATypeMain::o_parser.add(o_freq);
        ATypeMain::o_parser.add(o_mline);
        ATypeMain::o_parser.add(o_channel);
        };

    PPPConstETmachine parseMachine() {
      map < string, PPPConstETmachine > typesMap;
      typesMap["rene"] = WTERene;
      typesMap["morozov"] = WTEMorozov;
      typesMap["scovar"] = WTESumCovar;
      typesMap["acovar"] = WTECovar;
      typesMap["complex"] = WTEComplex;
      typesMap["mlinet"] = WTEMlinet;
      typesMap["mlinef"] = WTEMlinef;
      if(typesMap.find(o_type.getValue()) == typesMap.end())
        PPPBaseObject :: onError(PPPTRANSELLI_ERRMASH+string("calc()"));
      return typesMap[o_type.getValue()];
      };

    void calc(void) {
      // prepare transform parameters
      aTrans.setMachine(parseMachine());
      // calculation of elliptic properties
      calc(ATypeMain::aDest, ATypeMain::aSource);
      ATypeMain::aDest.setObjectName(o_name.getValue());
      // information
      strstream str;
      str << ATypeMain::aDest.getInfo() << ends;
      ATypeMain::onMessage(str.str());
      };

    void calc(PPPSignalContainer<PPPellipse2D> &aLocDest, PPPSignalContainer<double> &aLocSource) {
      if(aLocSource.channels() != 2)
        PPPBaseObject :: onError(gwlET2D_ERR01);
      if(aTrans.getMachine() == WTECovar || aTrans.getMachine() == WTESumCovar)
        aTrans.eval2DSignalPar(aLocDest,aLocSource,o_tw.getValue());
      else
        aTrans.eval2DSignalPar(aLocDest,aLocSource);
      };

    void calc(PPPSignalContainer<PPPellipse2D> &aLocDest, PPPSignalContainer<PPPcomplex> &aLocSource) {
      PPPSignalContainer<double> aSigDbl;
      PPPSignalContainer<PPPellipse2D> aElli;
      aLocDest.prepare(aLocSource);
      aLocDest.setObjectName(PPPBASETEMPLATE_EPAR2D);
      aSigDbl.prepare(aLocSource.points(),2,aLocSource.getAxis(),aLocSource.getObjectName());
      for(unsigned j=0; j<aLocSource.channels(); j++)
        {
        for(unsigned i=0; i<aLocSource.points(); i++)
          {
          aSigDbl(i,0) = real(aLocSource(i,j));
          aSigDbl(i,1) = imag(aLocSource(i,j));
          }
        if(aTrans.getMachine() == WTECovar || aTrans.getMachine() == WTESumCovar)
          aTrans.eval2DSignalPar(aElli,aSigDbl,o_tw.getValue());
        else
          aTrans.eval2DSignalPar(aElli,aSigDbl);
        aLocDest.getChannel(j).assign(aElli.getChannel(0));
        aTrans.setShowMessage(false);
        }
      };

    void calc(PPPSpectrContainer<PPPellipse2D> &aLocDest, PPPSpectrContainer<PPPcomplex> &aLocSource) {
      if(aTrans.getMachine() == WTERene && ATypeMain::aSpectrumPar.getFreq().getSign() != PPPAxis::ASplus && aLocSource.channels() != 2)
        PPPBaseObject :: onError(gwlET2D_ERR01);
      if(aTrans.getMachine() == WTEComplex && ATypeMain::aSpectrumPar.getFreq().getSign() != PPPAxis::ASfull)
        PPPBaseObject :: onError(gwlET2D_ERR01);
      if(aTrans.getMachine() == WTECovar && ATypeMain::aSpectrumPar.getFreq().getSign() != PPPAxis::ASplus && aLocSource.channels() != 2)
        PPPBaseObject :: onError(gwlET2D_ERR01);
      aTrans.eval2DSpectrPar(aLocDest,aLocSource,o_filter.getValue(),o_tw.getValue(),o_channel.getValue());
      };

    void calc(PPPSignalContainer<PPPellipse2D> &aLocDest, PPPSpectrContainer<PPPcomplex> &aLocSource) {
      PPPSignalContainer<double> aRidge;
      ATypeMain::read_bindata(aRidge,o_mline.getValue());
      aTrans.eval2DSignalPar(aLocDest,aLocSource,aRidge,o_filter.getValue());
      };

    void calc(PPPSignalContainer<PPPellipse2D> &aLocDest, PPPSpectrContainer<PPPellipse2D> &aLocSource) {
      PPPVectorContainer<double> freq;
      freq.strToVector("{"+string(o_freq.getValue())+"}");
      if(freq.size() != 2)
        PPPBaseObject :: onError(gwlET2D_ERR02);
      aTrans.eval2DSignalPar(aLocDest,aLocSource,PPPTransElli::TFTtime,freq[0],freq[1]);
      };

  };  // end of object

  
main(int  argc, char **argv)
  {
  string appname = "2D elliptic transform";
  string appcode = "gwlET2D";
  gwlMain<gwlMainSigToSigDbl> WT1(appname.c_str(),appcode.c_str(),"source.dat","ellipar2D.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  PPPObjectIO::ObjectType head = WT1.read_binheader();
  if(head == PPPObjectIO::SIGD)
    {
    WT1.evaluate();
    }
  else if(head == PPPObjectIO::SIGC)
    {
    gwlMain<gwlMainSigToSigCmpl> WT2(appname.c_str(),appcode.c_str(),"source.dat","ellipar2D.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  else if(head == PPPObjectIO::SPECC && (WT1.parseMachine()==WTEComplex || WT1.parseMachine()==WTECovar || WT1.parseMachine()==WTERene))
    {
    gwlMain<gwlMainSpecToSpec> WT3(appname.c_str(),appcode.c_str(),"source.dat","ellipar2D.dat");
    WT3.parse(argc, argv);
    WT3.evaluate();
    }
  else if(head == PPPObjectIO::SPECC && (WT1.parseMachine()==WTEMlinet || WT1.parseMachine()==WTEMlinef))
    {
    gwlMain<gwlMainSpecToSig> WT4(appname.c_str(),appcode.c_str(),"source.dat","ellipar2D.dat");
    WT4.parse(argc, argv);
    WT4.evaluate();
    }
  else if(head == PPPObjectIO::SPECE2D)
    {
    gwlMain<gwlMainElliToSig> WT5(appname.c_str(),appcode.c_str(),"source.dat","ellipar2D.dat");
    WT5.parse(argc, argv);
    WT5.evaluate();
    }
  else
    ConApplication.onError(gwlMain_ERR02);

  return 0;
  }




