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
#ifndef _PPPPROPAGATORDISS
#define _PPPPROPAGATORDISS

#define PPPPROPAGATORDISS_NAME       "dispersive propagator"
#define PPPPROPAGATORDISS_VEL        "phase and group velocity"
#define PPPPROPAGATORDISS_SYNT       "diffeomorphism of Ricker wavelet in Fourier space"
#define PPPPROPAGATORDISS_FOUR       "calculation diffeomorphism in Fourier space"
#define PPPPROPAGATORDISS_FOURCR     "calculation diffeomorphism in cross-correlation Fourier space"
#define PPPPROPAGATORDISS_FOURCRS    "calculation diffeomorphism in cross-correlation Fourier space for seismogram"
#define PPPPROPAGATORDISS_CRMATR     "cross-correlations matrix"
#define PPPPROPAGATORDISS_WAV        "calculation diffeomorphism in wavelet space"
#define PPPPROPAGATORDISS_WAVCR      "calculation diffeomorphism in cross-correlation wavelet space"
#define PPPPROPAGATORDISS_ERRNOTDEF  "propagation type is not defined in procedure: "

/************************************************************************
 * PPPPropagatorDiss
 * @article{Kulesh2005PaAG,
 *   author = {M. Kulesh and M. Holschneider and M. S. Diallo and Q. Xie and F. Scherbaum},
 *   title = {Modeling of wave dispersion using continuous wavelet transforms},
 *   journal = {Pure and Applied Geophysics},
 *   year = {2005},
 *   volume = {162},
 *   number = {5},
 *   pages = {843-855},
 *   urllink = {http://www.springerlink.com/content/g23g1t8436541341}
 * }
 ************************************************************************/
class PPPPropagatorDiss : public PPPBaseObject
  {
  public:

    typedef enum {STnone, STfour, STfourcross, STwav, STwavcross} PropType;

  private:
    PPPcomplex                  _j;
    PPPDispersionModel          _model;
    double                      _distance;
    PPPTransWavelet<PPPcomplex> _trans;

  public:

    PPPPropagatorDiss(void) : _j(0.0,1.0), _distance(0.0) {
      setObjectName(PPPPROPAGATORDISS_NAME);
      };

    friend ostream& operator << (ostream& aDest, PPPPropagatorDiss &aSour) {
      aDest << aSour.getObjectName() << ": " << aSour.getNotation() << endl;
      aDest << "  " << aSour.getModel().getInfo() << ends;
      return aDest;
      };

    inline PPPDispersionModel &getModel() { return _model; };
    inline double getDistance(void) const { return _distance; };
    inline void setDistance(double adX) { _distance=adX; };

    void evalFouriePropag(PropType aPropType, PPPSignalContainer<PPPcomplex> &aDest,
      PPPSignalContainer<PPPcomplex> &aSource, double dX1=0, double dX2=0) {
      switch(aPropType)
        {
        case STfour: onMessage(PPPPROPAGATORDISS_FOUR);
             dX1 = getDistance(); dX2 = 0.0;
             break;
        case STfourcross: onMessage(PPPPROPAGATORDISS_FOURCR);
             if(dX1 == 0.0 && dX2 == 0.0) dX2 = getDistance();
             break;
        default: onError(PPPPROPAGATORDISS_ERRNOTDEF+string("evalFouriePropag"));
             break;
        }
      double f,atnval;
      PPPSignalContainer<PPPcomplex> Splus,Sminus;
      _trans.FTSeparate(Splus,Sminus,aSource);
      for(unsigned i=0; i<aSource.points()/2; i++)
        {
        f = aSource.getAxis(i);
        if(f == 0.0 || f<getModel().getFrMin() || f>getModel().getFrMax())
          {
          Splus(i) = 0.0;
          Sminus(i) = 0.0;
          continue;
          }
        atnval = exp(-(dX1+dX2)*getModel().Atn(f));
        Splus(i) = atnval*exp(-2.0*_j*M_PI*(dX1-dX2)*getModel().Phi(f))*Splus(i);
        Sminus(i) = atnval*exp(+2.0*_j*M_PI*(dX1-dX2)*getModel().Phi(f))*Sminus(i);
        }
      aDest.resize(aSource.points());
      aDest.setAxis(aSource.getAxis());
      _trans.FTCut(aDest,Splus,Sminus);
      // Notation
      strstream str;
      if(aPropType == STfour)
        {
        str << PPPPROPAGATORDISS_FOUR << endl;
        str << "  YS2(f) = exp(-dX*Atn(f))*exp(-2.0*Pi*i*dX*Phi(f))*YS1(f), dX=" << dX1 << ends;
        }
      if(aPropType == STfourcross)
        {
        str << PPPPROPAGATORDISS_FOURCR << endl;
        str << "  YS2(f) = exp(-(dX1+dX2)*Atn(f))*exp(-2.0*Pi*i*(dX1-dX2)*Phi(f))*YS1(f)" << ends;
        }
      onNotation(str.str());  
      return;
      };

    void evalFourieCrossPropag(PPPSignalContainer<PPPcomplex> &aDest, PPPSignalContainer<PPPcomplex> &aSource,
      PPPMatrixContainer<unsigned> &aNumb) {
      onMessage(PPPPROPAGATORDISS_FOURCRS);
      if(aNumb.cols() != 2)
        onError(ARG_VALUE+string("evalFourieCrossPropag"));
      aDest.resize(aSource.points(), aNumb.rows());
      aDest.setAxis(aSource.getAxis());
      aDest.getChannel(0).assign(aSource.getChannel(0));
      PPPSignalContainer<PPPcomplex> aYRes;
      double dX1, dX2;
      for(unsigned i=1; i<aNumb.rows(); i++)
        {
        dX1 = (double)(aNumb(i,0)-aNumb(0,0))*getDistance();
        dX2 = (double)(aNumb(i,1)-aNumb(0,1))*getDistance();
        setShowMessage(false);
        evalFouriePropag(STfourcross,aYRes,aSource,dX1,dX2);
        setShowMessage(true);
        aDest.getChannel(i).assign(aYRes.getChannel(0));
        }
      // Notation
      strstream str;
      str << PPPPROPAGATORDISS_FOURCRS << endl;
      str << "  YS2(f) = exp(-(dX1+dX2)*Atn(f))*exp(-2.0*Pi*i*(dX1-dX2)*Phi(f))*YS1(f)" << endl;
      str << "  " << PPPPROPAGATORDISS_CRMATR << ": " << aNumb.matrixToStr() << ends;
      onNotation(str.str());
      };

    void evalWaveletPropag(PropType aPropType, PPPSpectrContainer<PPPcomplex> &aDest,
      PPPSpectrContainer<PPPcomplex> &aSour, PPPSpectrParams &aTransPar,
      unsigned aAlg, double dX1=0, double dX2=0) {
      // Distance modifacation
      switch(aPropType)
        {
        case STwav: onMessage(PPPPROPAGATORDISS_WAV);
             dX1 = getDistance(); dX2 = 0.0;
             break;
        case STwavcross: onMessage(PPPPROPAGATORDISS_WAVCR);
             if(dX1 == 0.0 && dX2 == 0.0) dX2 = getDistance();
             break;
        default: onError(PPPPROPAGATORDISS_ERRNOTDEF+string("evalWaveletPropag"));
             break;
        }
      // Separate into positive and negative spectrum
      PPPSpectrContainer<PPPcomplex> APlus, AMinus;
      bool aCatFlag = (aTransPar.getFreq().getSign() == PPPAxis::ASfull);
      if(aCatFlag)
        {
        _trans.WaveletSeparate(APlus, AMinus, aSour);
        for(unsigned i=0; i<AMinus.voices(); i++) AMinus.getFreq()[i] = fabs(AMinus.getFreq(i));
        AMinus.getFreq().assign(APlus.getFreq());
        }
      else APlus.assign(aSour);
      // dispersion model
      PPPSignalContainer<double> PropPar;
      getModel().evalContainer(PropPar, APlus.getFreq());
      // Propagation
      PPPSpectrContainer<PPPcomplex> ANewPlus(APlus), ANewMinus(AMinus);
      strstream str;
      if(aAlg == 1) // Standart propagator
        {
        for(unsigned i=0; i<ANewPlus.voices(); i++)
          for(unsigned j=0; j<ANewPlus.points(); j++)
          {
          double f = APlus.getFreq(i);
          double b = APlus.getTime(j);
          double atnval = exp(-(dX1+dX2)*PropPar(i,PPPDispersionModel::iAtn));
          ANewPlus(i,j) = atnval*exp(-2.0*_j*M_PI*(dX1-dX2)*
             (PropPar(i,PPPDispersionModel::iWn) - f*PropPar(i,PPPDispersionModel::iWndiv)))*
             APlus.Get(f,b-(dX1-dX2)*PropPar(i,PPPDispersionModel::iWndiv));
          if(aTransPar.getFreq().getSign() == PPPAxis::ASfull)
            ANewMinus(i,j) = atnval*exp(+2.0*_j*M_PI*(dX1-dX2)*
            (PropPar(i,PPPDispersionModel::iWn) - f*PropPar(i,PPPDispersionModel::iWndiv)))*
            AMinus.Get(f,b-(dX1-dX2)*PropPar(i,PPPDispersionModel::iWndiv));
          }
        if(aPropType == STwav)
          {
          str << PPPPROPAGATORDISS_WAV << endl;
          str << "  WgS2(t,f) = exp(-dX*Atn(f))*exp(-2.0*Pi*i*dX*(Phi(f)-f*Phidiv(f)))*WgS1(t-dX*Phidiv(f),f)" << ends;
          }
        if(aPropType == STwavcross)
          {
          str << PPPPROPAGATORDISS_WAVCR << endl;
          str << "  WgS2(t,f) = exp(-(dX+dX2)*Atn(f))*exp(-2.0*Pi*i*(dX1-dX2)*(Phi(f)-f*Phidiv(f)))*WgS1(t-(dX1-dX2)*Phidiv(f),f)" << ends;
          }
        }
      else if(aAlg == 2) // Modulus-Argument propagator
        {
        for(unsigned i=0; i<ANewPlus.voices(); i++)
          for(unsigned j=0; j<ANewPlus.points(); j++)
          {
          double f = APlus.getFreq(i);
          double b = APlus.getTime(j);
          double atnval = exp(-(dX1+dX2)*PropPar(i,PPPDispersionModel::iAtn));
          double zp1 = abs(APlus.Get(f,b-(dX1-dX2)/PropPar(i,PPPDispersionModel::iCg)));
          double zp2 = arg(APlus.Get(f,b-(dX1-dX2)/PropPar(i,PPPDispersionModel::iCp))); // in general -1.0/f
          double zm1 = abs(AMinus.Get(f,b-(dX1-dX2)/PropPar(i,PPPDispersionModel::iCg)));
          double zm2 = arg(AMinus.Get(f,b-(dX1-dX2)/PropPar(i,PPPDispersionModel::iCp))); // in general -1.0/f
          ANewPlus(i,j) = atnval*zp1*exp(_j*zp2);
          if(aTransPar.getFreq().getSign() == PPPAxis::ASfull)
            ANewMinus(i,j) = atnval*zm1*exp(_j*zm2);
          }
        if(aPropType == STwav)
          {
          str << PPPPROPAGATORDISS_WAV << endl;
          str << "  WgS2(t,f) = exp(-dX*Atn(f))*|WgS1(t-dX/VSg(f),f)|*exp[i arg WgS1(t-dX/VSp(f),f)]" << ends;
          }
        if(aPropType == STwavcross)
          {
          str << PPPPROPAGATORDISS_WAVCR << endl;
          str << "  WgS2(t,f) = exp(-(dX1+dX2)*Atn(f))*|WgS1(t-(dX1-dX2)/VSg(f),f)|*exp[i arg WgS1(t-(dX1-dX2)/VSp(f),f)]" << ends;
          }
        }
      else if(aAlg == 3) // Causal propagator
        {
        PPPcomplex phasep, phasem;
        for(unsigned i=0; i<ANewPlus.voices(); i++)
          {
          double f = APlus.getFreq(i);
          phasep = -2.0*M_PI*_j*(PropPar(i,PPPDispersionModel::iWn)-f*PropPar(i,PPPDispersionModel::iWndiv))*(dX1-dX2) - (PropPar(i,PPPDispersionModel::iAtn)-f*PropPar(i,PPPDispersionModel::iAtndiv))*(dX1+dX2);
          phasem = +2.0*M_PI*_j*(PropPar(i,PPPDispersionModel::iWn)-f*PropPar(i,PPPDispersionModel::iWndiv))*(dX1-dX2) - (PropPar(i,PPPDispersionModel::iAtn)-f*PropPar(i,PPPDispersionModel::iAtndiv))*(dX1+dX2);
          double fs = 1.0 + (dX1+dX2)*f*PropPar(i,PPPDispersionModel::iAtndiv)/(aTransPar.getWaveletPar(PPPSpectrParams::WSTdirect));
          double atnval = 1.0/pow(fs, aTransPar.getWaveletPar(PPPSpectrParams::WSTdirect));
          for(unsigned j=0; j<ANewPlus.points(); j++)
            {
            double b = APlus.getTime(j);
            ANewPlus(i,j) = atnval*exp(phasep)*APlus.Get(f*fs, b-(dX1-dX2)*PropPar(i,PPPDispersionModel::iWndiv));
            if(aTransPar.getFreq().getSign() == PPPAxis::ASfull)
              ANewMinus(i,j) = atnval*exp(phasem)*AMinus.Get(f*fs, b-(dX1-dX2)*PropPar(i,PPPDispersionModel::iWndiv));
            }
          }
        if(aPropType == STwav)
          {
          str << PPPPROPAGATORDISS_WAV << endl;
          str << "  WgS2(t,f) = atnval*exp(phase)*WgS1(t-dX*ReKprim(f),f*fs)" << ends;
          }
        if(aPropType == STwavcross)
          {
          str << PPPPROPAGATORDISS_WAVCR << endl;
          str << "  WgS2(t,f) = atnval*exp(phase)*WgS1(t-(dX1-dX2)*ReKprim(f),f*fs)" << ends;
          }
        }
      else if(aAlg == 4) // Classic elliptic propagator
        {
        PPPellipse2D ew;
        PPPcomplex aC1,aC2;
        for(unsigned i=0; i<ANewPlus.voices(); i++)
          for(unsigned j=0; j<ANewPlus.points(); j++)
          {
          double f = APlus.getFreq(i);
          double b = APlus.getTime(j);
          double dX = getDistance();
          aC1 = APlus.Get(f, b-dX*PropPar(i,PPPDispersionModel::iWndiv));
          aC2 = AMinus.Get(f, b-dX*PropPar(i,PPPDispersionModel::iWndiv));
          ew.CWTParams(aC1, aC2);
          ew.CWTParamsDefomDiss(exp(-dX*PropPar(i,PPPDispersionModel::iAtn)), 2.0*M_PI*dX*(PropPar(i,PPPDispersionModel::iWn) - f*PropPar(i,PPPDispersionModel::iWndiv)));
          ew.CWTParamsInverse(aC1,aC2);
          ANewPlus(i,j) = 2.0*aC1;
          ANewMinus(i,j) = 2.0*aC2;
          }
        str << PPPPROPAGATORDISS_WAV << endl;
        str << "  WgS2(f,t) = ElliParInverse(Etrans1(WgS1(f,t)))" << ends;
        }
      else if(aAlg == 5) // Morlet elliptic propagator
        {
        PPPellipse2D ew1, ew2;
        PPPcomplex aC1,aC2;
        for(unsigned i=0; i<ANewPlus.voices(); i++)
          for(unsigned j=0; j<ANewPlus.points(); j++)
          {
          double f = APlus.getFreq(i);
          double b = APlus.getTime(j);
          double dX = getDistance();
          aC1 = APlus.Get(f,b-dX/PropPar(i,PPPDispersionModel::iCg));
          aC2 = AMinus.Get(f,b-dX/PropPar(i,PPPDispersionModel::iCg));
          ew1.CWTParams(aC1, aC2);
          aC1 = APlus.Get(f,b-dX/PropPar(i,PPPDispersionModel::iCp));
          aC2 = AMinus.Get(f,b-dX/PropPar(i,PPPDispersionModel::iCp));
          ew2.CWTParams(aC1, aC2);
          ew1.CWTParamsDefomDiss(exp(-dX*PropPar(i,PPPDispersionModel::iAtn)), ew2);
          ew1.CWTParamsInverse(aC1,aC2);
          ANewPlus(i,j) = 2.0*aC1;
          ANewMinus(i,j) = 2.0*aC2;
          }
        str << PPPPROPAGATORDISS_WAV << endl;
        str << "  WgS2(f,t) = ElliParInverse(Etrans2(WgS1(f,t)))" << ends;
        }
      else onError(PPPPROPAGATORDISS_ERRNOTDEF+string("evalWaveletPropag"));
      // Cat positive and negative components
      if(aCatFlag) _trans.WaveletCat(aDest, ANewPlus, ANewMinus);
      else aDest.assign(ANewPlus);
      onNotation(str.str());
      return;
      };

    void assign(PPPPropagatorDiss &aSource) {
      getModel().assign(aSource.getModel());
      setDistance(aSource.getDistance());
      onNotation(aSource.getNotation());
      };

  }; // end of object



#endif
