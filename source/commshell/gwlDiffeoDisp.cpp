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

#define gwlDiffeoDisp_ERR01  "type of the propagator does not correspont to type of input object"

typedef gwlMainObject< PPPSignalContainer<double>, PPPSignalContainer<double> > gwlMainSigReal;
typedef gwlMainObject< PPPSignalContainer<PPPcomplex>, PPPSignalContainer<PPPcomplex> > gwlMainSigCmpl;
typedef gwlMainObject< PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<PPPcomplex> > gwlMainSpecCmpl;

template<class AType, class ATypeMain> class gwlMain : public ATypeMain
  {
  private:
    UTOption_str        o_name;
    UTOption_file       o_model;
    UTOption_dbl        o_dist;
    UTOption_int        o_step;
    UTOption_lit        o_isacorr;
    UTOption_dbl        o_decrs;
    UTOption_int        o_prop;
    PPPPropagatorDiss   aPropag;
    unsigned aChannels;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      ATypeMain(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n", "name",    "<str>", "name of the propagated function (by default 'propagated signal')", "propagated signal"),
      o_model  ("l", "model",   "<file>", "input file with dispersion model (by default 'dispmodel.dat')", "dispmodel.dat"),
      o_dist   ("d", "dist",    "<real>", "propagation distance (by default 0)", 0.0),
      o_step   ("s", "step",    "<unsigned>", "count of iterations (by default 1)", 1),
      o_isacorr("a", "acorr",   "if set, the input signal/spectrum is an autocorrelation", false),
      o_decrs  ("e", "decrs",   "<real>", "decrement coefficient for source signal (by default 1)", 1.0),
      o_prop   ("p", "prop",    "<unsigned>", "type of the wavelet propagator: 1-standart, 2-modulus-argument, 3-causal, 4,5-elliptic (by default 1)", 1)
        {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_model);
        ATypeMain::o_parser.add(o_dist);
        ATypeMain::o_parser.add(o_step);
        ATypeMain::o_parser.add(o_isacorr);
        ATypeMain::o_parser.add(o_decrs);
        ATypeMain::o_parser.add(o_prop);
        };

    void calc(void) {
      // prepare transform parameters
      ATypeMain::read_bindata(aPropag.getModel(), o_model.getValue());
      aPropag.setDistance(o_dist.getValue());
      aPropag.setShowMessage(false);
      strstream str1;
      str1 << aPropag.getModel().getInfo() << ends;
      PPPBaseObject :: onMessage(str1.str());
      // calculation
      aChannels = o_step.getValue()+1;
      calc(ATypeMain::aSource);
      ATypeMain::aDest.setObjectName(o_name.getValue());
      // information
      strstream str2;
      str2 << ATypeMain::aDest.getInfo() << ends;
      PPPBaseObject :: onMessage(str2.str());
      };

    void calc(PPPSignalContainer<AType> &aData) {
      PPPTransWavelet<AType> aTrans;
      PPPSignalContainer<PPPcomplex> aFour0, aFour1, aFour;
      aTrans.FT(aFour0, aData);
      aFour.prepare(aFour0.points(), aChannels, aFour0.getAxis(), aFour0.getObjectName());
      aFour.getChannel(0).assign(aFour0.getChannel(0));
      if(o_decrs.isOptionGiven())
        for(unsigned i=0; i<aFour.points(); i++)
          aFour(i,0) = aFour(i,0)/o_decrs.getValue();
      double DeltaDis = aPropag.getDistance();
      PPPPropagatorDiss::PropType aPropType = (o_isacorr.isOptionGiven())? PPPPropagatorDiss::STfourcross : PPPPropagatorDiss::STfour;
      for(unsigned i=1; i<aChannels; i++)
        {
        if(aPropag.getModel().isAnalytical() && aPropag.getModel().getWn().getType() == PPPApproximate::APTtwogauss)
          aPropag.getModel().getWn()[6] = aPropag.getDistance();
        aPropag.evalFouriePropag(aPropType,aFour1,aFour0);
        aFour.getChannel(i).assign(aFour1.getChannel(0));
        if(i == 1) PPPBaseObject :: onMessage(aPropag.getNotation());
        aPropag.setDistance(DeltaDis+aPropag.getDistance());
        }
      double aMin = aData.getAxis().getMin();
      double aMax = aData.getAxis().getMax();
      aTrans.IFT(ATypeMain::aDest,aFour,aMin,aMax);
      };

    void calc(PPPSpectrContainer<PPPcomplex> &aData) {
      PPPSpectrContainer<PPPcomplex> aCWT;
      ATypeMain::aDest.prepare(aData.voices(), aData.points(), aChannels, aData.getTime(), aData.getFreq(), aData.getObjectName());
      ATypeMain::aDest.getChannel(0).assign(aData.getChannel(0));
      if(o_decrs.isParsed())
        for(unsigned i=0; i<ATypeMain::aDest.voices(); i++)
          for(unsigned j=0; j<ATypeMain::aDest.points(); j++)
            ATypeMain::aDest(i,j,0) = ATypeMain::aDest(i,j,0)/o_decrs.getValue();
      double DeltaDis = aPropag.getDistance();
      PPPPropagatorDiss::PropType aPropType = (o_isacorr.isOptionGiven())? PPPPropagatorDiss::STwavcross : PPPPropagatorDiss::STwav;
      for(unsigned i=1; i<aChannels; i++)
        {
        aPropag.evalWaveletPropag(aPropType,aCWT,aData,ATypeMain::aSpectrumPar,o_prop.getValue());
        ATypeMain::aDest.getChannel(i).assign(aCWT.getChannel(0));
        if(i == 1) PPPBaseObject :: onMessage(aPropag.getNotation());
        aPropag.setDistance(DeltaDis+aPropag.getDistance());
        }
      };

  };  // end of object

  
main(int  argc, char **argv)
  {
  string appname = "Dispersive propagator";
  string appcode = "gwlDiffeoDisp";

  gwlMain<double,gwlMainSigReal> WT1(appname.c_str(),appcode.c_str(),"source.dat","diffeodisp.dat");
  WT1.parse(argc, argv);
  ConApplication.onMessage(ConApplication.getAppName());
  PPPObjectIO::ObjectType head = WT1.read_binheader();
  if(head == PPPObjectIO::SIGD)
    {
    WT1.evaluate();
    }
  else if(head == PPPObjectIO::SIGC)
    {
    gwlMain<PPPcomplex,gwlMainSigCmpl> WT2(appname.c_str(),appcode.c_str(),"source.dat","diffeodisp.dat");
    WT2.parse(argc, argv);
    WT2.evaluate();
    }
  else if(head == PPPObjectIO::SPECC)
    {
    gwlMain<PPPcomplex,gwlMainSpecCmpl> WT3(appname.c_str(),appcode.c_str(),"source.dat","diffeodisp.dat");
    WT3.parse(argc, argv);
    WT3.evaluate();
    }
  else
    ConApplication.onError(gwlMain_ERR02);
  return 0;
  }




