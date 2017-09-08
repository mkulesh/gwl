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

#define gwlSignalRead_ERR01  "invalid signal type in routine: "
#define gwlSignalRead_ERR02  "invalid format of input file in routine: "

typedef gwlMainObject< PPPSignalContainer<double>, PPPSignalContainer<double> > gwlMainDouble;
typedef gwlMainObject< PPPSignalContainer<PPPcomplex>, PPPSignalContainer<PPPcomplex> > gwlMainCmpl;

template<class AType, class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_str        o_type;
    UTOption_str        o_format;
    UTOption_lit        o_istime;
    UTOption_dbl        o_smplfreq;
    UTOption_dbl        o_tmin;
    UTOption_dbl        o_tmax;
    UTOption_lit        o_to2p;
    UTOption_str        o_chan;
    UTOption_int        o_geo;
    UTOption_str        o_rot;
    UTOption_int        o_resample;
    UTOption_dbl        o_mult;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name     ("n", "name",     "<str>", "name of the output signal (by default 'signal')", "signal"),
      o_format   ("f", "format",   "<str>", "format of input data (ASCII, STA. By default 'ASCII')", "ASCII"),
      o_type     ("y", "type",     "<str>", "structure of output signal (func, seis, seis2D. By default 'func')", "func"),
      o_istime   ("e", "istime",   "if set, the time column is present in the input file", false),
      o_smplfreq ("s", "smplfreq", "<real>", "sampling frequency of source signal (by default 100)", 100.0),
      o_tmin     ("1", "tmin",     "<real>", "minimum value of time of source signal (by default 0)", 0.0),
      o_tmax     ("2", "tmax",     "<real>", "maximum value of time frequency of source signal (by default 0)", 0.0),
      o_to2p     ("p", "to2p",     "if set, the signal will be expand to power of two", false),
      o_chan     ("a", "chan",     "<str>",  "channel(s) of source file which will be readed \"chan1,chan2,...chanN\" (by default '')", ""),
      o_geo      ("g", "geo",      "<unsigned>", "index of geophone for SAT format (by default 0)", 0),
      o_rot      ("r", "rot",      "<str>",  "parameters for rotation procedure of two channels \"chan1,chan2,angle\" (by default '')", ""),
      o_resample ("l", "resample", "<unsigned>",  "resampling of seismograms (by default '0')", 0),
      o_mult     ("u", "mult",     "<real>",  "multiplication of all seismograms with a constant (by default '1')", 1)
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(ATypeMain::o_iscmpl);
        ATypeMain::o_parser.add(o_format);
        ATypeMain::o_parser.add(o_type);
        ATypeMain::o_parser.add(o_istime);
        ATypeMain::o_parser.add(o_smplfreq);
        ATypeMain::o_parser.add(o_tmin);
        ATypeMain::o_parser.add(o_tmax);
        ATypeMain::o_parser.add(o_to2p);
        ATypeMain::o_parser.add(o_chan);
        ATypeMain::o_parser.add(o_geo);
        ATypeMain::o_parser.add(o_rot);
        ATypeMain::o_parser.add(o_resample);
        ATypeMain::o_parser.add(o_mult);
        };

    void read(void) {
      // prepare reading parameters
      map < string, unsigned > inputTypes;
      inputTypes["ASCII"] = 1;
      inputTypes["STA"] = 2;
      inputTypes["WAV"] = 3;
      if(inputTypes.find(o_format.getValue()) == inputTypes.end())
        PPPBaseObject :: onError(gwlSignalRead_ERR01+string("read()"));

      map < string, unsigned > outputTypes;
      outputTypes["func"] = 1;
      outputTypes["seis"] = 2;
      outputTypes["seis2D"] = 3;
      if(outputTypes.find(o_type.getValue()) == outputTypes.end())
        PPPBaseObject :: onError(gwlSignalRead_ERR01+string("read()"));

      PPPSignalRead<AType> aPar;
      aPar.setSamplingFreq(o_smplfreq.getValue());
      aPar.setTmin(o_tmin.getValue());
      aPar.setTmax(o_tmax.getValue());
      aPar.setToPowerOfTwo(o_to2p.isOptionGiven());
      aPar.setReample(o_resample.getValue());
      aPar.getChannels().strToVector("{"+string(o_chan.getValue())+"}");
      if(o_mult.isOptionGiven()) aPar.setMult(o_mult.getValue());
      bool isTime = o_istime.isOptionGiven();

      // check: is the input signal already binary?
      unsigned currType = outputTypes[o_type.getValue()];
      if(currType < 3)
        {
        aPar.setFileName(ATypeMain::o_infile.getValue());
        PPPObjectIO::ObjectType head = ATypeMain::read_binheader();
        if(head == PPPObjectIO::SIGD || head == PPPObjectIO::SIGC)
          {
          ATypeMain::read_bindata(ATypeMain::aDest,aPar.getFileName().c_str());
          currType = 0;
          }
        else if (head != PPPObjectIO::UNDEF)
          {
          PPPBaseObject :: onError(gwlSignalRead_ERR02+string("read()"));
          }
        }

      // Generation of signal
      aPar.setShowNotation(PPPBaseObject::NMappend);
      switch(currType)
        {
        case 1:
          if(aPar.getChannels().size()==0)
            aPar.getChannels().strToVector("{0,1}");
          if(inputTypes[o_format.getValue()] == 1)
             aPar.readFunctionASCII(ATypeMain::aDest,isTime);
          else if(inputTypes[o_format.getValue()] == 2)
             {
             aPar.setGeophone(o_geo.getValue());
             aPar.readFunctionGolm(ATypeMain::aDest);
             }
          else
            aPar.readFunctionWav(ATypeMain::aDest);
          break;
        case 2:
          if(o_rot.isOptionGiven())
            aPar.getRotation().strToVector("{"+string(o_rot.getValue())+"}");
          if(inputTypes[o_format.getValue()] == 1)
             aPar.readSeisASCII(ATypeMain::aDest,isTime);
          else
             {
             aPar.setGeophone(o_geo.getValue());
             aPar.readSeisGolm(ATypeMain::aDest);
             }
          break;
        case 3:
          ATypeMain::parseFileNames();
          if(ATypeMain::aFileNames.size() != 2)
            PPPBaseObject :: onError("too few or too many input signals are given");
          aPar.setFileName(ATypeMain::aFileNames[0],ATypeMain::aFileNames[1]);
          aPar.readSeis2Files(ATypeMain::aDest,isTime);
          break;
        }
      // post processing
      aPar.postProcessing(ATypeMain::aDest);
      ATypeMain::aDest.setObjectName(o_name.getValue());
      // Information
      strstream str;
      str << ATypeMain::aDest.getInfo() << endl << "  " << aPar.getNotation() << ends;
      PPPBaseObject :: onMessage(str.str());
      };

    void calc(void) {
      };

  };  // end of object

  
main(int  argc, char **argv)
  {
  gwlMain<double,gwlMainDouble> WT1("Import of a signal","gwlSignalRead","source.asc","signal.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  if(WT1.isComplex())
    {
    gwlMain<PPPcomplex,gwlMainCmpl> WT2("Import of a signal","gwlSignalRead","source.asc","signal.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  else WT1.evaluate();
  return 0;
  }

