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

typedef gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<PPPcomplex> > gwlMainSpecToSpec;
typedef gwlMainObject< PPPSignalContainer<PPPcomplex>, PPPSignalContainer<PPPcomplex> > gwlMainSigToSig;

template<class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_type;
    // filter
    UTOption_str        o_filter;
    UTOption_file       o_elli;
    vector<PPPEllipse2Dfilter> aFilters;
    PPPTransElli        aTrans;
    PPPTransWavelet<PPPcomplex> aTransCmpl;
    // wavelet transform
    UTOption_file       o_freq;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n", "name",    "<str>",  "name of the elliptic parameters (by default 'elliptic parameters')", "elliptic parameters"),
      o_type   ("y", "type",    "<str>",  "type of the elliptic machine: acovar, complex (by default 'complex')", "complex"),
      o_filter ("f", "filter",  "<str>",  "polarization filter parameters (by default '')", ""),
      o_elli   ("e", "elli",    "<file>", "input file elliptical parameters (by default 'elli.dat')", "elli.dat"),
      o_freq   ("q", "freq",    "<file>", "input file with frequency axis (by default 'frequency.dat')", "frequency.dat")
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(ATypeMain::o_progr);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_type);
        ATypeMain::o_parser.add(o_filter);
        ATypeMain::o_parser.add(o_elli);
        // wavelet transform
        ATypeMain::o_parser.add(o_freq);
        ATypeMain::o_parser.add(ATypeMain::o_wavelet);
        ATypeMain::o_parser.add(ATypeMain::o_wavpar);
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

    void parsePolarizationFilters() {
      map < string, PPPEllipse2Dfilter::FilterType > typesMap;
      typesMap["linhor"] = PPPEllipse2Dfilter::TRFLinHor;
      typesMap["linvert"] = PPPEllipse2Dfilter::TRFLinVert;
      typesMap["ellihor"] = PPPEllipse2Dfilter::TRFElliHor;
      typesMap["ellivert"] = PPPEllipse2Dfilter::TRFElliVert;
      typesMap["linhors"] = PPPEllipse2Dfilter::TRFLinHorS;
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
        if(aPars.size() != 2)
          ATypeMain::onError(gwlET2D_ERR04);
        aFilters.push_back(PPPEllipse2Dfilter(typesMap[strname],aPars[0],aPars[1]));
        }
      };

    void calc(void) {
      // prepare transform parameters
      if(o_filter.isOptionGiven())
        parsePolarizationFilters();
      aTrans.setMachine(parseMachine());
      aTrans.setShowProgress(ATypeMain::o_progr.getValue());
      aTransCmpl.setShowProgress(ATypeMain::o_progr.getValue());
      // calculation of elliptic properties
      if(aFilters.size()>0)
        calc(ATypeMain::aDest, ATypeMain::aSource);
      };

    void calc(PPPSpectrContainer<PPPcomplex> &aLocDest, PPPSpectrContainer<PPPcomplex> &aLocSource) {
      unsigned newchann = (aTrans.getMachine() == WTEComplex)? aFilters.size() : 2*aFilters.size();
      aLocDest.prepare(
        aLocSource.voices(),
        aLocSource.points(),
        newchann,
        aLocSource.getTime(),
        aLocSource.getFreq(),
        "filtered spectrum");
      // calculation
      PPPSpectrContainer<PPPellipse2D> aElli;
      ATypeMain::read_bindata(aElli, o_elli.getValue());
      for(unsigned k=0; k<aFilters.size(); k++)
        aTrans.eval2DFilter(aLocDest, aLocSource, aElli, aFilters[k], k);
      aLocDest.setObjectName(o_name.getValue());
      // information
      strstream str;
      str << aLocDest.getInfo() << ends;
      ATypeMain::onMessage(str.str());
      };

    void calc(PPPSignalContainer<PPPcomplex> &aLocDest, PPPSignalContainer<PPPcomplex> &aLocSource) {
      // prepare transform parameters
      PPPAxis aFreq;
      ATypeMain::read_bindata(aFreq, o_freq.getValue());
      ATypeMain::aSpectrumPar.initialize(1,
        aFreq,
        ATypeMain::o_wavelet.getValue(),
        ATypeMain::o_wavpar.getValue(),
        "delta",
        1.0);
      // prepare inverse seismogramm
      vector< PPPSignalContainer<PPPcomplex> > aRecSeis;
      aRecSeis.resize(aFilters.size());
      for(unsigned k=0; k<aFilters.size(); k++)
        {
        aRecSeis[k].prepare(aLocSource);
        aRecSeis[k].setObjectName(o_name.getValue());
        }
      // calculation
      if(parseMachine() == WTEComplex)
        calcFullFilterComplex(aRecSeis, aLocSource);
      else
        calcFullFilterACovar(aRecSeis, aLocSource);
      // write the filtered seismogramms
      for(unsigned k=0; k<aRecSeis.size(); k++)
        {
        string newName(ATypeMain::o_outfile.getValue());
        if(aRecSeis.size() > 1)
          {
          char num[10];
          sprintf(num,"(%d)",k+1);
          newName.replace(newName.find("."),1,string(num)+".");
          }
        ATypeMain::write_bindata(aRecSeis[k], newName.c_str());
        }
      ATypeMain::o_outfile.setValue("");
      };

    void calcFullFilterComplex(vector< PPPSignalContainer<PPPcomplex> > &aRecSeis, PPPSignalContainer<PPPcomplex> &aLocSource) {
      PPPSignalContainer<PPPcomplex>    aInvSig;
      PPPSpectrContainer<PPPcomplex>    aCWT,aCWTfil;
      PPPSpectrContainer<PPPellipse2D>  aElli;
      aTransCmpl.WT(aCWT, aLocSource, ATypeMain::aSpectrumPar);
      aCWTfil.prepare(aCWT.voices(),aCWT.points(),aFilters.size(),aCWT.getTime(),aCWT.getFreq(),"filtered spectrum");
      char numchan[100];
      for(unsigned j=0; j<aLocSource.channels(); j++)
        {
        if(j>0)
          {
          sprintf(numchan, "%s%3d","current channel: ",j+1);
          aTrans.onProgress(100*j/(aLocSource.channels()-1),numchan);
          }
        aTrans.eval2DSpectrPar(aElli, aCWT, 0.0, 1, j, (j==0));
        for(unsigned k=0; k<aFilters.size(); k++)
          aTrans.eval2DFilter(aCWTfil, aCWT, aElli, aFilters[k], k, j);
        aTransCmpl.IWT(aInvSig, aCWTfil, ATypeMain::aSpectrumPar);
        for(unsigned k=0; k<aFilters.size(); k++)
          aRecSeis[k].getChannel(j).assign(aInvSig.getChannel(k));
        aTrans.setShowMessage(false);
        aTransCmpl.setShowMessage(false);
        }
      aTrans.onProgress(-1,"");
      };

    void calcFullFilterACovar(vector< PPPSignalContainer<PPPcomplex> > &aRecSeis, PPPSignalContainer<PPPcomplex> &aLocSource) {
      PPPSignalContainer<double>        aSig,aInvSig;
      PPPSpectrContainer<PPPcomplex>    aCWT,aCWTfil;
      PPPSpectrContainer<PPPellipse2D>  aElli;
      aSig.prepare(aLocSource.points(),2,aLocSource.getAxis(),aLocSource.getObjectName());
      char numchan[100];
      for(unsigned j=0; j<aLocSource.channels(); j++)
        {
        if(j>0)
          {
          sprintf(numchan, "%s%3d","current channel: ",j+1);
          aTrans.onProgress(100*j/(aLocSource.channels()-1),numchan);
          }
        for(unsigned i=0; i<aLocSource.points(); i++)
          {
          aSig(i,0) = real(aLocSource(i,j));
          aSig(i,1) = imag(aLocSource(i,j));
          }
        aTrans.WT(aCWT, aSig, ATypeMain::aSpectrumPar);
        if(j==0)
          aCWTfil.prepare(aCWT.voices(),aCWT.points(),2*aFilters.size(),aCWT.getTime(),aCWT.getFreq(),"filtered spectrum");
        aTrans.eval2DSpectrPar(aElli,aCWT,0.0,1);
        for(unsigned k=0; k<aFilters.size(); k++)
          aTrans.eval2DFilter(aCWTfil, aCWT, aElli, aFilters[k], k);
        aTrans.IWT(aInvSig, aCWTfil, ATypeMain::aSpectrumPar);
        for(unsigned k=0; k<aFilters.size(); k++)
          for(unsigned i=0; i<aLocSource.points(); i++)
            aRecSeis[k](i,j) = PPPcomplex(aInvSig(i,2*k),aInvSig(i,2*k+1));
        aTrans.setShowMessage(false);
        }
      aTrans.onProgress(-1,"");
      };

  };  // end of object


main(int  argc, char **argv)
  {
  string appname = "2D polarization filter";
  string appcode = "gwlET2DFilter";
  gwlMain<gwlMainSpecToSpec> WT1(appname.c_str(),appcode.c_str(),"source.dat","filtered.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  PPPObjectIO::ObjectType head = WT1.read_binheader();
  if(head == PPPObjectIO::SPECC)
    {
    WT1.evaluate();
    }
  else if(head == PPPObjectIO::SIGC)
    {
    gwlMain<gwlMainSigToSig> WT2(appname.c_str(),appcode.c_str(),"source.dat","filtered.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  else
    ConApplication.onError(gwlMain_ERR02);

  return 0;
  }




