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

typedef gwlMainObject< PPPAxis, PPPAxis > gwlMainAxis;
typedef gwlMainObject< PPPSignalContainer<double>, PPPSignalContainer<double> > gwlMainSignalDouble;
typedef gwlMainObject< PPPSignalContainer<PPPcomplex>, PPPSignalContainer<double> > gwlMainSignalCmpl;
typedef gwlMainObject< PPPSignalContainer<PPPellipse2D>, PPPSignalContainer<double> > gwlMainSignalElli2D;
typedef gwlMainObject< PPPSignalContainer<PPPellipse3D>, PPPSignalContainer<double> > gwlMainSignalElli3D;
typedef gwlMainObject< PPPSpectrContainer<double>, PPPSpectrContainer<double> > gwlMainSpectrDouble;
typedef gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<double> > gwlMainSpectrCmpl;
typedef gwlMainObject< PPPSpectrContainer<PPPellipse2D>, PPPSpectrContainer<double> > gwlMainSpectrElli2D;
typedef gwlMainObject< PPPSpectrContainer<PPPellipse3D>, PPPSpectrContainer<double> > gwlMainSpectrElli3D;

#define gwlConvert_ERR01 "Code of real component is unknown in procedure: calc()"

template<class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_comp;
    UTOption_dbl        o_filter;
    UTOption_int        o_channel;
    UTOption_int        o_voices;
    UTOption_int        o_points;
    UTOption_lit        o_degree;
    UTOption_lit        o_modpi;
    unsigned            aChann;
    double              aFilter;
    strstream           aInfo;
    PPPVectorContainer<unsigned> aComp;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n",  "name",   "<str>", "name of the axis (by default the name of source object)", ""),
      o_comp   ("p",  "comp",   "<str>",  "transformation codes of the source complex object (by default 3,5)", "3,5"),
      o_filter ("f",  "filter", "<real>", "phase filter limit in percent to absolute value (by default 1)", 1.0),
      o_channel("c",  "chan",   "<int>",  "channel of input signal/spectrum (by default 0)", 0),
      o_voices ("v",  "voices", "<unsigned>", "ASCII mode only: reduce spectrum's matrix to 'voices' rows (by default 256)", 256),
      o_points ("s",  "points", "<unsigned>", "ASCII mode only: reduce spectrum's matrix to 'points' colums (by default 256)", 256),
      o_degree ("d",  "degree", "if set, the angular values will be written in degrees", false),
      o_modpi  ("1",  "modpi",  "if set, the modulo pi will be applied to angular values (only valid for: 3D azimuth)", false)
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_comp);
        ATypeMain::o_parser.add(o_filter);
        ATypeMain::o_parser.add(o_channel);
        ATypeMain::o_parser.add(o_voices);
        ATypeMain::o_parser.add(o_points);
        ATypeMain::o_parser.add(o_degree);
        ATypeMain::o_parser.add(o_modpi);
        };

    void calc(void) {
      if(o_comp.isOptionGiven())
        aComp.strToVector("{"+string(o_comp.getValue())+"}");
      if(aComp.size() > 0)
        {
        aChann = o_channel.getValue();
        aFilter = o_filter.getValue();
        calc(ATypeMain::aSource);
        if(o_name.isOptionGiven())
          ATypeMain::aDest.setObjectName(o_name.getValue());
        aInfo << ATypeMain::aDest.getInfo() << ends;
        ATypeMain::onMessage(aInfo.str());
        }
      };

    void calc(PPPAxis &aData) {
      ATypeMain::aDest.assign(ATypeMain::aSource);
      };

    void calc(PPPSignalContainer<double> &aData) {
      ATypeMain::aDest.assign(ATypeMain::aSource);
      };

    void calc(PPPSignalContainer<PPPcomplex> &aData) {
      bool atype = (aComp.size()==1);
      unsigned newsize = (atype)? aData.channels() : aComp.size();
      ATypeMain::aDest.prepare(aData.points(),newsize,aData.getAxis(),aData.getObjectName());
      double maxmod = 0.0;
      for(unsigned i=0; i<newsize; i++)
        {
        unsigned ind1 = (atype)? 0 : i;
        unsigned ind2 = (atype)? i : aChann;
        if(aComp[ind1] == TCargf || aComp[ind1] == TCabsm || aComp[ind1] == TCargm)
          {
          PPPcomplexAbs wand;
          maxmod = wand(aData.getChannel(ind2).getMaxValue(wand));
          }
        switch(aComp[ind1])
          {
          case TCre:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexRe()); break;
          case TCim:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexIm()); break;
          case TCabs:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexAbs()); break;
          case TCarg:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexArg(o_degree.getValue())); break;
          case TCargf: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexArgFilter(maxmod, aFilter, -M_PI, o_degree.getValue())); break;
          case TCabsm: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexAbsM(maxmod)); break;
          case TCargm: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexArgM(maxmod)); break;
          default: ATypeMain::onError(gwlConvert_ERR01);
          }
        if(!atype)
          aInfo << "Transform complex signal to: " << ATypeMain::aDest.getChannel(i).getObjectName() << endl;
        else if(atype && i==0)
          aInfo << "Transform complex signal to: " << ATypeMain::aDest.getChannel(i).getObjectName() << "[" << ATypeMain::aDest.channels() << "]" << endl;
        }
      };

    void calc(PPPSignalContainer<PPPellipse2D> &aData) {
      bool atype = (aComp.size()==1);
      unsigned newsize = (atype)? aData.channels() : aComp.size();
      ATypeMain::aDest.prepare(aData.points(),newsize,aData.getAxis(),aData.getObjectName());
      for(unsigned i=0; i<newsize; i++)
        {
        unsigned ind1 = (atype)? 0 : i;
        unsigned ind2 = (atype)? i : aChann;
        switch(aComp[ind1])
          {
          case TErmin:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Drmin()); break;
          case TErmax:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Drmax()); break;
          case TEratio:     aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dratio()); break;
          case TEtilt:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dtilt(o_degree.getValue())); break;
          case TEphasediff: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dphasediff(o_degree.getValue())); break;
          case TEwx:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dwx()); break;
          case TEwy:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dwy()); break;
          case TEphase1:    aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dphase1()); break;
          case TEphase2:    aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dphase2()); break;
          default: ATypeMain::onError(gwlConvert_ERR01);
          }
        if(!atype)
          aInfo << "Transform 2D elliptic signal to: " << ATypeMain::aDest.getChannel(i).getObjectName() << endl;
        else if(atype && i==0)
          aInfo << "Transform 2D elliptic signal to: " << ATypeMain::aDest.getChannel(i).getObjectName() << "[" << ATypeMain::aDest.channels() << "]" << endl;
        }
      };

    void calc(PPPSignalContainer<PPPellipse3D> &aData) {
      bool atype = (aComp.size()==1);
      unsigned newsize = (atype)? aData.channels() : aComp.size();
      ATypeMain::aDest.prepare(aData.points(),newsize,aData.getAxis(),aData.getObjectName());
      for(unsigned i=0; i<newsize; i++)
        {
        unsigned ind1 = (atype)? 0 : i;
        unsigned ind2 = (atype)? i : aChann;
        switch(aComp[ind1])
          {
          case TErmin:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Drmin()); break;
          case TErmax:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Drmax()); break;
          case TErmidd:     aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Drmidd()); break;
          case TEratio:     aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dratio()); break;
          case TEratio1:    aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dratio1()); break;
          case TEratio2:    aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dratio2()); break;
          case TEwx:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dwx()); break;
          case TEwy:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dwy()); break;
          case TEwz:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dwz()); break;
          case TEplanarx:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplanarx()); break;
          case TEplanary:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplanary()); break;
          case TEplanarz:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplanarz()); break;
          case TEplancosx:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplancosx()); break;
          case TEplancosy:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplancosy()); break;
          case TEplancosz:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplancosz()); break;
          case TEsignratio: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dsignratio()); break;
          case TEdip:       aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Ddip(o_degree.getValue())); break;
          case TEazimuth:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dazimuth(o_degree.getValue(),o_modpi.getValue())); break;
          default: ATypeMain::onError(gwlConvert_ERR01);
          }
        if(!atype)
          aInfo << "Transform 3D elliptic signal to: " << ATypeMain::aDest.getChannel(i).getObjectName() << endl;
        else if(atype && i==0)
          aInfo << "Transform 3D elliptic signal to: " << ATypeMain::aDest.getChannel(i).getObjectName() << "[" << ATypeMain::aDest.channels() << "]" << endl;
        }
      };

    void calc(PPPSpectrContainer<double> &aData) {
      ATypeMain::aDest.assign(ATypeMain::aSource);
      };

    void calc(PPPSpectrContainer<PPPcomplex> &aData) {
      bool atype = (aComp.size()==1);
      unsigned newsize = (atype)? aData.channels() : aComp.size();
      ATypeMain::aDest.prepare(aData.voices(), aData.points(), newsize, aData.getTime(), aData.getFreq(), aData.getObjectName());
      double maxmod = 0.0;
      for(unsigned i=0; i<newsize; i++)
        {
        unsigned ind1 = (atype)? 0 : i;
        unsigned ind2 = (atype)? i : aChann;
        if(aComp[ind1] == TCargf || aComp[ind1] == TCabsm || aComp[ind1] == TCargm)
          {
          PPPcomplexAbs wand;
          maxmod = wand(aData.getChannel(ind2).getMaxValue(wand));
          }
        switch(aComp[ind1])
          {
          case TCre:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexRe()); break;
          case TCim:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexIm()); break;
          case TCabs:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexAbs()); break;
          case TCarg:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexArg(o_degree.getValue())); break;
          case TCargf: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexArgFilter(maxmod, aFilter, -M_PI, o_degree.getValue())); break;
          case TCabsm: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexAbsM(maxmod)); break;
          case TCargm: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPcomplexArgM(maxmod)); break;
          default: ATypeMain::onError(gwlConvert_ERR01);
          }
        if(!atype)
          aInfo << "Transform complex spectrum to: " << ATypeMain::aDest.getChannel(i).getObjectName() << endl;
        else if(atype && i==0)
          aInfo << "Transform complex spectrum to: " << ATypeMain::aDest.getChannel(i).getObjectName() << "[" << ATypeMain::aDest.channels() << "]" << endl;
        }
      };

    void calc(PPPSpectrContainer<PPPellipse2D> &aData) {
      bool atype = (aComp.size()==1);
      unsigned newsize = (atype)? aData.channels() : aComp.size();
      ATypeMain::aDest.prepare(aData.voices(), aData.points(), newsize, aData.getTime(), aData.getFreq(), aData.getObjectName());
      for(unsigned i=0; i<newsize; i++)
        {
        unsigned ind1 = (atype)? 0 : i;
        unsigned ind2 = (atype)? i : aChann;
        switch(aComp[ind1])
          {
          case TErmin:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Drmin()); break;
          case TErmax:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Drmax()); break;
          case TEratio:     aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dratio()); break;
          case TEtilt:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dtilt(o_degree.getValue())); break;
          case TEphasediff: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dphasediff(o_degree.getValue())); break;
          case TEwx:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dwx()); break;
          case TEwy:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dwy()); break;
          case TEphase1:    aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dphase1()); break;
          case TEphase2:    aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse2Dphase2()); break;
          default: ATypeMain::onError(gwlConvert_ERR01);
          }
        if(!atype)
          aInfo << "Transform 2D elliptic spectrum to: " << ATypeMain::aDest.getChannel(i).getObjectName() << endl;
        else if(atype && i==0)
          aInfo << "Transform 2D elliptic spectrum to: " << ATypeMain::aDest.getChannel(i).getObjectName() << "[" << ATypeMain::aDest.channels() << "]" << endl;
        }
      };

    void calc(PPPSpectrContainer<PPPellipse3D> &aData) {
      bool atype = (aComp.size()==1);
      unsigned newsize = (atype)? aData.channels() : aComp.size();
      ATypeMain::aDest.prepare(aData.voices(), aData.points(), newsize, aData.getTime(), aData.getFreq(), aData.getObjectName());
      for(unsigned i=0; i<newsize; i++)
        {
        unsigned ind1 = (atype)? 0 : i;
        unsigned ind2 = (atype)? i : aChann;
        switch(aComp[ind1])
          {
          case TErmin:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Drmin()); break;
          case TErmax:      aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Drmax()); break;
          case TErmidd:     aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Drmidd()); break;
          case TEratio:     aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dratio()); break;
          case TEratio1:    aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dratio1()); break;
          case TEratio2:    aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dratio2()); break;
          case TEwx:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dwx()); break;
          case TEwy:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dwy()); break;
          case TEwz:        aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dwz()); break;
          case TEplanarx:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplanarx()); break;
          case TEplanary:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplanary()); break;
          case TEplanarz:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplanarz()); break;
          case TEplancosx:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplancosx()); break;
          case TEplancosy:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplancosy()); break;
          case TEplancosz:  aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dplancosz()); break;
          case TEsignratio: aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dsignratio()); break;
          case TEdip:       aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Ddip(o_degree.getValue())); break;
          case TEazimuth:   aData.getChannel(ind2).compTransform(ATypeMain::aDest.getChannel(i), PPPellipse3Dazimuth(o_degree.getValue(),o_modpi.getValue())); break;
          default: ATypeMain::onError(gwlConvert_ERR01);
          }
        if(!atype)
          aInfo << "Transform 3D elliptic spectrum to: " << ATypeMain::aDest.getChannel(i).getObjectName() << endl;
        else if(atype && i==0)
          aInfo << "Transform 3D elliptic spectrum to: " << ATypeMain::aDest.getChannel(i).getObjectName() << "[" << ATypeMain::aDest.channels() << "]" << endl;
        }
      };

    virtual void write() {
      if(aComp.size() > 0)
        write(ATypeMain::aDest, ATypeMain::o_outfile.getValue());
      else
        write(ATypeMain::aSource, ATypeMain::o_outfile.getValue());
      };

    void write(PPPAxis &aData, const string &aName) {
      if(ATypeMain::o_outtype.getValue() == 1)
        {
        aData.write(aName);
        }
      else if(ATypeMain::o_outtype.getValue() == 2)
        {
        PPPBaseObject :: onMessage(FILE_WRITEBIN+string(aName));
        FILE *outfile = fopen(aName.c_str(), "wb");
        aData.fwrite(outfile);
        fclose(outfile);
        }
      else PPPBaseObject :: onError(gwlMain_ERR01);
      };

    template <class DataType>
    void write(PPPSignalContainer<DataType> &aData, const string &aName) {
      if(ATypeMain::o_outtype.getValue() == 1)
        {
        aData.write(aName);
        }
      else if(ATypeMain::o_outtype.getValue() == 2)
        {
        PPPBaseObject :: onMessage(FILE_WRITEBIN+string(aName));
        FILE *outfile = fopen(aName.c_str(), "wb");
        aData.fwrite(outfile);
        fclose(outfile);
        }
      else PPPBaseObject :: onError(gwlMain_ERR01);
      };

    template <class DataType>
    void write(PPPSpectrContainer<DataType> &aData, const string &aName) {
      if(ATypeMain::o_outtype.getValue() == 1)
        {
        aData.write(aName,o_voices.getValue(),o_points.getValue());
        }
      else if(ATypeMain::o_outtype.getValue() == 2)
        {
        PPPBaseObject :: onMessage(FILE_WRITEBIN+string(aName));
        FILE *outfile = fopen(aName.c_str(), "wb");
        aData.fwrite(outfile);
        ATypeMain::aSpectrumPar.fwrite(outfile);
        fclose(outfile);
        }
      else if(ATypeMain::o_outtype.getValue() == 3)
        {
        aData.writegnuplot(aName);
        }
      else PPPBaseObject :: onError(gwlMain_ERR01);
      };

  };  // end of object


main(int  argc, char **argv)
  {
  string appname = "Convert of complex binary objects to real/ASCII";
  string appcode = "gwlConvert";
  gwlMain<gwlMainAxis> WT1(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  PPPObjectIO::ObjectType head = WT1.read_binheader();
  if(head == PPPObjectIO::AXIS)
    {
    WT1.evaluate();
    }
  else if(head == PPPObjectIO::SIGD)
    {
    gwlMain<gwlMainSignalDouble> WT2(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  else if(head == PPPObjectIO::SIGC)
    {
    gwlMain<gwlMainSignalCmpl> WT3(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
    WT3.parse(argc, argv);
    WT3.evaluate();
    }
  else if(head == PPPObjectIO::SIGE2D)
    {
    gwlMain<gwlMainSignalElli2D> WT4(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
    WT4.parse(argc, argv);
    WT4.evaluate();
    }
  else if(head == PPPObjectIO::SIGE3D)
    {
    gwlMain<gwlMainSignalElli3D> WT5(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
    WT5.parse(argc, argv);
    WT5.evaluate();
    }
  else if(head == PPPObjectIO::SPECD)
    {
    gwlMain<gwlMainSpectrDouble> WT6(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
    WT6.parse(argc, argv);
    WT6.evaluate();
    }
  else if(head == PPPObjectIO::SPECC)
    {
    gwlMain<gwlMainSpectrCmpl> WT7(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
    WT7.parse(argc, argv);
    WT7.evaluate();
    }
  else if(head == PPPObjectIO::SPECE2D)
    {
    gwlMain<gwlMainSpectrElli2D> WT8(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
    WT8.parse(argc, argv);
    WT8.evaluate();
    }
  else if(head == PPPObjectIO::SPECE3D)
    {
    gwlMain<gwlMainSpectrElli3D> WT9(appname.c_str(),appcode.c_str(),"source.dat","converted.dat");
    WT9.parse(argc, argv);
    WT9.evaluate();
    }
  else
    ConApplication.onError(gwlMain_ERR02);
  return 0;
  }




