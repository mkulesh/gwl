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
#ifndef _PPPTRANSELLI
#define _PPPTRANSELLI

#define PPPTRANSELLI_NAME       "elliptic transform"
#define PPPTRANSELLI_CALCDIR    "calculation of polarization attributes"
#define PPPTRANSELLI_CALCINV    "inversion of polarization attributes"
#define PPPTRANSELLI_INV        "inverse spectrum"
#define PPPTRANSELLI_FILT       "coefficient for energy filter"
#define PPPTRANSELLI_TWIND      "time window"
#define PPPTRANSELLI_AVERTM     "average for time window"
#define PPPTRANSELLI_AVERFR     "average for frequency window"
#define PPPTRANSELLI_ERRMASH    "elliptic machine is undefined in procedure: "
#define PPPTRANSELLI_ERRSOURS   "invalid size of source seismogramm in procedure: "
#define PPPTRANSELLI_SIGLINF    "lineary filter for signal"
#define PPPTRANSELLI_DEL        "number of erased point"
#define PPPTRANSELLI_CHANNAME   "channel"

/************************************************************************
 * PPPTransElli
 ************************************************************************/
class PPPTransElli : public PPPTransWavelet<double>
  {
  private:
    PPPConstETmachine _machine;

  public:
    PPPTransElli(void) {
      setObjectName(PPPTRANSELLI_NAME);
      };

  inline PPPConstETmachine getMachine(void) {  return _machine;  };
  void setMachine(PPPConstETmachine amashine) { _machine = amashine;  };

  /** Direct Elliptic Transform for signal ************************************/
  void eval2DSignalPar(PPPSignalContainer<PPPellipse2D> &aElli, PPPSignalContainer<double> &aVoice, unsigned aN=1) {
    if(aVoice.channels() != 2)
      onError(PPPTRANSELLI_ERRSOURS+string("eval3DSignalPar()"));
    strstream                           str;
    PPPSignalContainer<double>          wx,wy;
    PPPSignalContainer<PPPcomplex>      Hx,Hy,aSigCmpl;
    PPPVectorContainer<double>          Zx,Zy;
    PPPTransWavelet<PPPcomplex>         aTransCmpl;
    // Prepare of parameter variable
    aElli.prepare(aVoice.points(),1,aVoice.getAxis(),PPPBASETEMPLATE_EPAR2D);
    aElli.getChannel(0).assign(PPPellipse2D(0.0));
    bool aSaveStatus = getShowMessage();
    setShowMessage(false);
    switch(getMachine())
      {
      case WTERene:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSIGNALCONTAINER_NAME << ends;
        onMessage(str.str());
        HT(Hx,aVoice,PPPAxis::ASplus);
        for(unsigned i=0; i<aElli.points(); i++)
          aElli(i).ComplexTraceParams(Hx(i,0),Hx(i,1));
        break;
      case WTEMorozov:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSIGNALCONTAINER_NAME << ends;
        onMessage(str.str());
        HT(Hx,aVoice,PPPAxis::ASplus);
        for(unsigned i=0; i<aElli.points(); i++)
          aElli(i).MorozovParams(Hx(i,0),Hx(i,1));
        break;
      case WTESumCovar:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSIGNALCONTAINER_NAME << endl << "  "
            << PPPTRANSELLI_TWIND << "=" << aN << ends;
        onMessage(str.str());
        Zx.resize(aN);
        Zy.resize(aN);
        for(unsigned i=0; i<aElli.points(); i++)
          {
          for(unsigned j=0; j<aN; j++)
            {
            if((i+j)>=aN/2 && (i+j)<(aElli.points()+aN/2))
              {
              Zx[j] = aVoice(i+j-aN/2,0);
              Zy[j] = aVoice(i+j-aN/2,1);
              }
            else
              Zx[j] = Zy[j] = 0.0;
            }
          aElli(i).CovarianceSumParams(Zx,Zy);
          }
        break;
      case WTECovar:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSIGNALCONTAINER_NAME << endl << "  "
            << PPPTRANSELLI_TWIND << "=" << aN << ends;
        onMessage(str.str());
        HT(Hx,aVoice,PPPAxis::ASplus);
        Hx.PhaseDifferentiation(wx);
        for(unsigned i=0; i<aElli.points(); i++)
          aElli(i).CovarianceParams(Hx(i,0),Hx(i,1),wx(i,0),wx(i,1),aN);
        break;
      case WTEComplex:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSIGNALCONTAINER_NAME << ends;
        onMessage(str.str());
        aSigCmpl.prepare(aVoice.points(),1,aVoice.getAxis(),aVoice.getObjectName());
        for(unsigned i=0; i<aVoice.points(); i++)
          aSigCmpl(i) = PPPcomplex(aVoice(i,0), aVoice(i,1));
        aTransCmpl.setShowMessage(false);
        aTransCmpl.HT(Hx,aSigCmpl,PPPAxis::ASminus);
        aTransCmpl.HT(Hy,aSigCmpl,PPPAxis::ASplus);
        Hx.PhaseDifferentiation(wx);
        Hy.PhaseDifferentiation(wy);
        for(unsigned i=0; i<aElli.points(); i++)
          if(aVoice(i) != 0.0)
            aElli(i).CWTParams(Hx(i),Hy(i),wx(i),wy(i));
        break;
      default:
        onError(PPPTRANSELLI_ERRMASH+string("eval2DSignalPar()"));
        break;
      }
    setShowMessage(aSaveStatus);
    };

  /** Average of 2D elliptic spectrum *****************************************/
  void eval2DSignalPar(PPPSignalContainer<PPPellipse2D> &aFunc, PPPSpectrContainer<PPPellipse2D> &aElli,
    PPPTransWavelet<double>::TimeFreqType aType, double aMin, double aMax) {
    strstream str;
    unsigned ind1,ind2,i,j,k;
    switch(aType)
      {
      case TFTtime:
        str << PPPTRANSELLI_AVERFR << ": " << aMin << ".." << aMax << ends;
        onMessage(str.str());
        aFunc.prepare(aElli.points(),1,aElli.getTime(),PPPBASETEMPLATE_EPAR2D);
        ind1 = aElli.getFreq().locateFloor(aMin);
        ind2 = aElli.getFreq().locateFloor(aMax);
        for(i=0; i<aElli.points(); i++)
          {
          for(k=0,j=ind1; j<=ind2; j++)
            {
            aFunc(i) = aFunc(i) + aElli(j,i);
            if(aElli(j,i).absmax()!=0.0 || aElli(j,i).absmin()!=0.0) k++;
            }
          if(k>0) aFunc(i) = (1.0/(double)k)*aFunc(i);
          }
        break;
      case TFTfreq:
        str << PPPTRANSELLI_AVERTM << ": " << aMin << ".." << aMax << ends;
        onMessage(str.str());
        aFunc.prepare(aElli.voices(),1,aElli.getFreq(),PPPBASETEMPLATE_EPAR2D);
        ind1 = aElli.getTime().locateFloor(aMin);
        ind2 = aElli.getTime().locateFloor(aMax);
        for(j=0; j<aElli.voices(); j++)
          {
          for(k=0,i=ind1; i<=ind2; i++)
            {
            aFunc(j) = aFunc(j) + aElli(j,i);
            if(aElli(j,i).absmax()!=0.0 || aElli(j,i).absmin()!=0.0) k++;
            }
          if(k>0) aFunc(j) = (1.0/(double)k)*aFunc(j);
          }
      break;
      }
    };

  /** Direct Elliptic Transform (CWT method for Ridge) ************************/
  void eval2DSignalPar(PPPSignalContainer<PPPellipse2D> &aElli, PPPSpectrContainer<PPPcomplex> &aWgFull,
    PPPSignalContainer<double> &aRidge, double aFilter) {
    strstream str;
    PPPSignalContainer<double>::InterpType aInterp = PPPSignalContainer<double>::ITsplin;
    double curreng = 0.0, WRplus = 0.0, WRminus = 0.0;
    PPPcomplexAbs trans;
    unsigned i;
    PPPcomplex c1, c2;
    PPPSpectrContainer<PPPcomplex> aWg1,aWg2;
    PPPTransWavelet<PPPcomplex> _trans;
    _trans.WaveletSeparate(aWg1,aWg2,aWgFull);
    double maxeng = trans(aWg1.getChannel().getMaxValue(trans))+trans(aWg2.getChannel().getMaxValue(trans));
    switch(getMachine())
      {
      case WTEMlinet:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSPECTRCONTAINER_NAME << endl << "  "
            << PPPTRANSELLI_FILT << "=" << aFilter << "%" << ends;
        onMessage(str.str());
        aElli.prepare(aWg1.points(), 1, aWg1.getTime(), PPPBASETEMPLATE_EPAR2D);
        for(i=0; i<aElli.points(); i++)
          {
          WRplus = aRidge.Get(aElli.getAxis(i),0,aInterp);
          WRminus = aRidge.Get(aElli.getAxis(i),1,aInterp);
          if(WRminus == 0.0 || WRplus == 0.0)
            {
            aElli(i) = 0.0;
            continue;
            }
          c1 = aWg2(aWg2.getFreq().locateFloor(WRminus),i);
          c2 = aWg1(aWg1.getFreq().locateFloor(WRplus),i);
          curreng = abs(c1) + abs(c2);
          if(curreng >= aFilter*maxeng/100.0)
            aElli(i).CWTParams(c1,c2);
          else
            aElli(i) = 0.0;
          }
        break;
      case WTEMlinef:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSPECTRCONTAINER_NAME << endl << "  "
            << PPPTRANSELLI_FILT << "=" << aFilter << "%" << ends;
        onMessage(str.str());
        aElli.prepare(aWg1.voices(), 1, aWg1.getFreq(), PPPBASETEMPLATE_EPAR2D);
        for(i=0; i<aElli.points(); i++)
          {
          WRplus = aRidge.Get(aElli.getAxis(i),0,aInterp);
          WRminus = aRidge.Get(aElli.getAxis(i),1,aInterp);
          if(WRminus == 0.0 || WRplus == 0.0)
            {
            aElli(i) = 0.0;
            continue;
            }
          c1 = aWg2(i,aWg2.getTime().locateFloor((WRminus+WRplus)/2.0));
          c2 = aWg1(i,aWg1.getTime().locateFloor((WRminus+WRplus)/2.0));
          curreng = abs(c1) + abs(c2);
          if(curreng >= aFilter*maxeng/100.0)
            aElli(i).CWTParams(c1,c2);
          else
            aElli(i) = 0.0;
          }
        break;
      }
    };


  /** Inverse elliptic Transform for signal ***********************************/
  void eval2DSignalParInverse(PPPSignalContainer<PPPcomplex> &aVoice, PPPSignalContainer<PPPellipse2D> &aElli) {
    strstream str;
    aVoice.prepare(aElli.points(),1,aElli.getAxis(),PPPTRANSFOUR_INV);
    PPPcomplex z1,z2;
    switch(getMachine())
      {
      case WTEComplex:
        str << PPPTRANSELLI_CALCINV << endl << "  " << getMachine() << ", "
            << PPPSIGNALCONTAINER_NAME << ends;
        onMessage(str.str());
        for(unsigned i=0; i<aElli.points(); i++)
          {
          aElli(i).CWTParamsInverse(z1,z2);
          aVoice(i) = z1+z2;
          }
        break;
      default:
        onError(PPPTRANSELLI_ERRMASH+string("eval2DSignalParInverse()"));
        break;
      }
    };

  /** Direct Elliptic Transform for spectrum **********************************/
  inline PPPcomplex& getCwtLink1(PPPSpectrContainer<PPPcomplex> &aWg, unsigned i, unsigned j, unsigned k=0) {
    return (getMachine() == WTEComplex)? aWg(aWg.voices()/2-1-i,j,k) : aWg(i,j,2*k);
    };

  inline PPPcomplex& getCwtLink2(PPPSpectrContainer<PPPcomplex> &aWg, unsigned i, unsigned j, unsigned k=0) {
    return (getMachine() == WTEComplex)? aWg(aWg.voices()/2+i,j,k) : aWg(i,j,2*k+1);
    };

  void eval2DSpectrPar(PPPSpectrContainer<PPPellipse2D> &aElli, PPPSpectrContainer<PPPcomplex> &aWg, double aFilter = 0.0, unsigned aN=1, unsigned aChann=0, bool aDoPrepare=true) {
    strstream str;
    register unsigned i,j;
    PPPcomplex c1, c2;
    double maxeng = 0.0, curreng = 0.0;
    PPPcomplexAbs trans;
    PPPSignalContainer<double> w1,w2;
    PPPSignalContainer<PPPcomplex> H1(aWg.points()),H2(aWg.points());
    switch(getMachine())
      {
      /* @article{Diallo2006GeophA,
          author = {M. S. Diallo and M. Kulesh and M. Holschneider and F. Scherbaum and F. Adler},
          title = {Characterization of polarization attributes of seismic waves using continuous wavelet transforms},
          journal = {Geophysics},
          year = {2006},
          volume = {71},
          number = {3},
          pages = {V67-V77},
          keywords = {seismic waves; wavelet transforms; geophysical signal processing},
          urllink = {http://link.aip.org/link/?GPY/71/V67/1}
        } */
      case WTEComplex:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSPECTRCONTAINER_NAME << ",  "
            << PPPTRANSELLI_CHANNAME << "=" << aChann << endl << "  "
            << PPPTRANSELLI_FILT << "=" << aFilter << "%" << ends;
        onMessage(str.str());
        if(aDoPrepare)
          {
          aElli.realloc(aWg.voices()/2, aWg.points());
          aElli.setNames(PPPBASETEMPLATE_EPAR2D,PPPSPECTRCONTAINER_TIME,PPPSPECTRCONTAINER_FREQ);
          aElli.setTime(aWg.getTime());
          aElli.getFreq().assign(aWg.getFreq(),PPPAxis::ASplus);
          }
        maxeng = trans(aWg.getChannel().getMaxValue(trans));
        for (i=0; i<aElli.voices(); i++)
          {
          for (j=0; j<aElli.points(); j++)
            {
            c1 = 2.0*getCwtLink1(aWg,i,j,aChann);
            c2 = 2.0*getCwtLink2(aWg,i,j,aChann);
            if(abs(c1) == 0.0 && abs(c2) == 0.0)
              {
              aElli(i,j) = 0.0;
              continue;
              }
            curreng = abs(c1) + abs(c2);
            if(curreng >= aFilter*maxeng/100.0)
              aElli(i,j).CWTParams(c1,c2);
            else
              aElli(i,j) = 0.0;
            }
          }
        break;
      case WTECovar:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSPECTRCONTAINER_NAME << ",  "
            << PPPTRANSELLI_CHANNAME << "=" << aChann << endl << "  "
            << PPPTRANSELLI_FILT << "=" << aFilter << "%" << endl << "  "
            << PPPTRANSELLI_TWIND << "=" << aN << ends;
        onMessage(str.str());
        if(aDoPrepare)
          aElli.prepare(aWg.voices(), aWg.points(), 1, aWg.getTime(), aWg.getFreq(), PPPBASETEMPLATE_EPAR2D);
        maxeng = trans(aWg.getChannel(0).getMaxValue(trans))+trans(aWg.getChannel(1).getMaxValue(trans));
        H1.setAxis(aWg.getTime());
        H2.setAxis(aWg.getTime());
        for (i=0; i<aElli.voices(); i++)
          {
          for (j=0; j<aElli.points(); j++)
            {
            H1(j) = 2.0*getCwtLink1(aWg,i,j,aChann);
            H2(j) = 2.0*getCwtLink2(aWg,i,j,aChann);
            }
          H1.PhaseDifferentiation(w1);
          H2.PhaseDifferentiation(w2);
          for (j=0; j<aElli.points(); j++)
            {
            c1 = H1(j);
            c2 = H2(j);
            curreng = abs(c1) + abs(c2);
            if(curreng >= aFilter*maxeng/100.0)
              aElli(i,j).CovarianceParams(c1,c2,w1(j),w2(j),aN);
            else
              aElli(i,j) = 0.0;
            }
          }
        break;
      case WTERene:
        str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
            << PPPSPECTRCONTAINER_NAME << endl << "  "
            << PPPTRANSELLI_FILT << "=" << aFilter << "%" << ends;
        onMessage(str.str());
        if(aDoPrepare)
          aElli.prepare(aWg.voices(), aWg.points(), 1, aWg.getTime(), aWg.getFreq(), PPPBASETEMPLATE_EPAR2D);
        maxeng = trans(aWg.getChannel(0).getMaxValue(trans))+trans(aWg.getChannel(1).getMaxValue(trans));
        for (i=0; i<aElli.voices(); i++)
          for (j=0; j<aElli.points(); j++)
            {
            c1 = 2.0*getCwtLink1(aWg,i,j,aChann);
            c2 = 2.0*getCwtLink2(aWg,i,j,aChann);
            curreng = abs(c1) + abs(c2);
            if(curreng >= aFilter*maxeng/100.0)
              aElli(i,j).ComplexTraceParams(c1,c2);
            else
              aElli(i,j) = 0.0;
            }
        break;
      default:
        onError(PPPTRANSELLI_ERRMASH+string("eval2DSpectrPar()"));
        break;
      }
    };

  /** Tilt-Ratio Filter (2D method for spectr) ********************************/
  /* @article{Diallo2006GeophA,
      author = {M. S. Diallo and M. Kulesh and M. Holschneider and F. Scherbaum and F. Adler},
      title = {Characterization of polarization attributes of seismic waves using continuous wavelet transforms},
      journal = {Geophysics},
      year = {2006},
      volume = {71},
      number = {3},
      pages = {V67-V77},
      keywords = {seismic waves; wavelet transforms; geophysical signal processing},
      urllink = {http://link.aip.org/link/?GPY/71/V67/1}
    } */
  void eval2DFilter(PPPSpectrContainer<PPPcomplex> &aWgFiltered, PPPSpectrContainer<PPPcomplex> &aWg,
    PPPSpectrContainer<PPPellipse2D> &aElli, PPPEllipse2Dfilter &aFilter, unsigned aChann, unsigned aSChann=0) {
    register unsigned i,j,numb=0;
    for(i=0; i<aElli.voices(); i++)
      for(j=0; j<aElli.points(); j++)
        if(aFilter.isContent(aElli(i,j)))
           {
           getCwtLink1(aWgFiltered,i,j,aChann) = getCwtLink1(aWg,i,j,aSChann);
           getCwtLink2(aWgFiltered,i,j,aChann) = getCwtLink2(aWg,i,j,aSChann);
           }
        else
           {
           getCwtLink1(aWgFiltered,i,j,aChann) = 0.0;
           getCwtLink2(aWgFiltered,i,j,aChann) = 0.0;
           numb++;
           }
    strstream str;
    str << aFilter.getInfo() << endl << "  " << PPPTRANSELLI_DEL << ": " << numb <<
      " (" << 100.0*numb/(aElli.voices()*aElli.points()) << "%)" << ends;
    onMessage(str.str());
    }

  /** 3D Elliptic Transform for signal ****************************************/
  void eval3DSignalPar(PPPSignalContainer<PPPellipse3D> &aElli, PPPSignalContainer<double> &aVoice, unsigned aN=1) {
    if(aVoice.channels() != 3)
      onError(PPPTRANSELLI_ERRSOURS+string("eval3DSignalPar"));
    PPPSignalContainer<double>     aW;
    PPPSignalContainer<PPPcomplex> aH;
    PPPVectorContainer<double>     Zx,Zy,Zz;
    PPPTransFour<double> Trans;
    // Prepare of parameter variable
    aElli.prepare(aVoice.points(), 1, aVoice.getAxis(), PPPBASETEMPLATE_EPAR3D);
    for(unsigned i=0; i<aElli.points(); i++) aElli(i) = 0;
    Trans.setShowMessage(false);
    // Calculaton of parameter variable
    strstream str;
    str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
        << PPPSIGNALCONTAINER_NAME << "[" << aVoice.channels() << "]";
    switch(getMachine())
      {
      case WTEMorozov:
        str << ends;
        onMessage(str.str());
        Trans.HT(aH,aVoice,PPPAxis::ASplus);
        for(unsigned i=0; i<aElli.points(); i++)
          aElli(i).MorozovParams(aH(i,0),aH(i,1),aH(i,2));
        break;
      case WTESumCovar:
        str << endl << "  " << PPPTRANSELLI_TWIND << "=" << aN << ends;
        onMessage(str.str());
        Zx.resize(aN);
        Zy.resize(aN);
        Zz.resize(aN);
        for(unsigned i=0; i<aElli.points(); i++)
          {
          for(unsigned j=0; j<aN; j++)
            {
            if((i+j)>=aN/2 && (i+j)<(aElli.points()+aN/2))
              {
              Zx[j] = aVoice(i+j-aN/2,0);
              Zy[j] = aVoice(i+j-aN/2,1);
              Zz[j] = aVoice(i+j-aN/2,2);
              }
            else
              Zx[j] = Zy[j] = Zz[j] = 0.0;
            }
          aElli(i).CovarianceSumParams(Zx,Zy,Zz);
          }
        break;
      case WTECovar:
        str << endl << "  " << PPPTRANSELLI_TWIND << "=" << aN << ends;
        onMessage(str.str());
        Trans.HT(aH,aVoice,PPPAxis::ASplus);  aH.PhaseDifferentiation(aW);
        for(unsigned i=0; i<aElli.points(); i++)
          aElli(i).CovarianceParams(aH(i,0),aH(i,1),aH(i,2),aW(i,0),aW(i,1),aW(i,2),aN);
        break;
      default:
        onError(PPPTRANSELLI_ERRMASH+string("eval3DSignalPar"));
        break;
      }
    };

  /** Inverse 3D Elliptic Transform for signal ********************************/
  void eval3DSignalParInverse(PPPSignalContainer<double> &aVoice, PPPSignalContainer<PPPellipse3D> &aElli) {
    aVoice.prepare(aElli.points(), 3, aElli.getAxis(), PPPTRANSFOUR_INV);
    PPPcomplex p1,p2,p3;
    strstream str;
    str << PPPTRANSELLI_CALCINV << endl << "  " << getMachine() << ", "
        << PPPSIGNALCONTAINER_NAME << "[" << aVoice.channels() << "]";
    switch(getMachine())
      {
      case WTEMorozov:
        str << ends;
        onMessage(str.str());
        for(unsigned i=0; i<aElli.points(); i++)
          {
          aElli(i).MorozovParamsInverse(p1,p2,p3);
          aVoice(i,0) = p1.real();
          aVoice(i,1) = p2.real();
          aVoice(i,2) = p3.real();
          }
        break;
      default:
        onError(PPPTRANSELLI_ERRMASH+string("TimeElliParInverse"));
        break;
      }
    };

  /** 3D Elliptic Transform for spectr ****************************************/
  void eval3DSpectrPar(PPPSpectrContainer<PPPellipse3D> &aElli, PPPSpectrContainer<PPPcomplex> &aWg,
    double aFilter = 0.0, unsigned aN=1) {
    if(aWg.channels() != 3)
      onError(PPPTRANSELLI_ERRSOURS+string("eval3DSpectrPar"));
    aElli.prepare(aWg.voices(), aWg.points(), 1, aWg.getTime(), aWg.getFreq(), PPPBASETEMPLATE_EPAR3D);
    unsigned i,j,k;
    PPPVectorContainer<PPPcomplex> aC(3);
    double maxeng = 0.0, curreng = 0.0;
    PPPcomplexAbs trans;
    maxeng = abs(aWg.getChannel(0).getMaxValue(trans))+
             abs(aWg.getChannel(1).getMaxValue(trans))+
             abs(aWg.getChannel(2).getMaxValue(trans));
    PPPSignalContainer<double>   aW;
    PPPSignalContainer<PPPcomplex> aH(aWg.points(),3);
    char StrStat[100];
    strstream str;
    str << PPPTRANSELLI_CALCDIR << endl << "  " << getMachine() << ", "
        << PPPSPECTRCONTAINER_NAME << "[" << aWg.channels() << "]";
    switch(getMachine())
      {
      /* @article{Diallo2005GP,
          author = {M. S. Diallo and M. Kulesh and M. Holschneider and F. Scherbaum},
          title = {Instantaneous polarization attributes in the time-frequency domain and wavefield separation},
          journal = {Geophysical Prospecting},
          volume = {53},
          number = {5},
          pages = {723-731},
          year = {2005},
          doi = {10.1111/j.1365-2478.2005.00500.x},
          urllink = {http://www.blackwell-synergy.com/doi/abs/10.1111/j.1365-2478.2005.00500.x},
          eprint = {http://www.blackwell-synergy.com/doi/pdf/10.1111/j.1365-2478.2005.00500.x}
        } */
      case WTEMorozov:
        str << endl << "  " << PPPTRANSELLI_FILT << "=" << aFilter << "%" << ends;
        onMessage(str.str());
        for (i=0; i<aWg.voices(); i++) for (j=0; j<aWg.points(); j++)
          {
          for (k=0; k<3; k++) aC[k] = aWg(i,j,k);
          curreng = abs(aC[0]) + abs(aC[1]) + abs(aC[2]);
          if(curreng >= aFilter*maxeng/100.0)
            aElli(i,j).MorozovParams(aC[0],aC[1],aC[2]);
          else
            aElli(i,j) = 0.0;
          }
        break;
      /* @article{Kulesh2007GJI,
          author = {M. Kulesh and M. S. Diallo and M. Holschneider and K. Kurennaya and F. Kruger and M. Ohrnberger and F. Scherbaum},
          title = {Polarization analysis in the wavelet domain based on the adaptive covariance method},
          journal = {Geophysical Journal International},
          year = {2007},
          volume = {170},
          number = {2},
          pages = {667-678},
          doi = {10.1111/j.1365-246X.2007.03417.x},
          urllink = {http://www.blackwell-synergy.com/doi/abs/10.1111/j.1365-246X.2007.03417.x},
          eprint = {http://www.blackwell-synergy.com/doi/pdf/10.1111/j.1365-246X.2007.03417.x}
        } */
      case WTECovar:
        str << endl << "  " << PPPTRANSELLI_FILT << "=" << aFilter << "%" << endl << "  " << PPPTRANSELLI_TWIND << "=" << aN << ends;
        onMessage(str.str());
        aH.setAxis(aWg.getTime());
        for (i=0; i<aWg.voices(); i++)
          {
          sprintf(StrStat,"%s: %3d",PPPTRANSWAVELET_CURR,i+1);
          onProgress(100*i/(aWg.voices()-1),StrStat);
          for (j=0; j<aWg.points(); j++)
            for (k=0; k<3; k++) aH(j,k) = aWg(i,j,k);
          aH.PhaseDifferentiation(aW);
          for (j=0; j<aWg.points(); j++)
            {
            for (k=0; k<3; k++) aC[k] = aH(j,k);
            curreng = abs(aC[0]) + abs(aC[1]) + abs(aC[2]);
            if(curreng >= aFilter*maxeng/100.0)
              aElli(i,j).CovarianceParams(aC[0],aC[1],aC[2],aW(j,0),aW(j,1),aW(j,2),aN);
            else
              aElli(i,j) = 0.0;
            }
          }
        onProgress(-1,"");
        break;
      }
    };

  /** Inverse 3D Elliptic Transform for spectrum ******************************/
  /* @article{Diallo2005GP,
      author = {M. S. Diallo and M. Kulesh and M. Holschneider and F. Scherbaum},
      title = {Instantaneous polarization attributes in the time-frequency domain and wavefield separation},
      journal = {Geophysical Prospecting},
      volume = {53},
      number = {5},
      pages = {723-731},
      year = {2005},
      doi = {10.1111/j.1365-2478.2005.00500.x},
      urllink = {http://www.blackwell-synergy.com/doi/abs/10.1111/j.1365-2478.2005.00500.x},
      eprint = {http://www.blackwell-synergy.com/doi/pdf/10.1111/j.1365-2478.2005.00500.x}
    } */
  void eval3DSpectrParInverse(PPPSpectrContainer<PPPcomplex> &aWg, PPPSpectrContainer<PPPellipse3D> &aElli) {
    strstream str;
    str << PPPTRANSELLI_CALCINV << endl << "  " << WTEMorozov << ", "
        << PPPSIGNALCONTAINER_NAME << "[3]" << ends;
    onMessage(str.str());
    aWg.resize(aElli.voices(),aElli.points(),3);
    aWg.setTime(aElli.getTime());
    aWg.setFreq(aElli.getFreq());
    aWg.setNames(PPPTRANSELLI_INV,PPPSPECTRCONTAINER_TIME,PPPSPECTRCONTAINER_FREQ);
    // Inverse
    PPPVectorContainer<PPPcomplex> p(3);
    unsigned i,j,k;
    for(i=0; i<aElli.voices(); i++) for(j=0; j<aElli.points(); j++)
      {
      aElli(i,j).MorozovParamsInverse(p[0],p[1],p[2]);
      for (k=0; k<3; k++) aWg(i,j,k) = p[k];
      }
    };

  /** Linery Filter (3D method for signal) ************************************/
  void eval3DFilter(PPPSignalContainer<double> &aDest, PPPSignalContainer<double> &aSource,
  PPPSignalContainer<PPPellipse3D> &aElli) {
    if(aSource.channels() != 3)
      onError(PPPTRANSELLI_ERRSOURS+string("eval3DFilter"));
    aDest.prepare(aSource);
    onMessage(PPPTRANSELLI_SIGLINF);
    for(unsigned i=0; i<aElli.points(); i++)
      {
      double rmax = aElli(i).absmax();
      for(unsigned k=0; k<3; k++)
        if(rmax == 0.0)
          aDest(i,k) = 0.0;
        else
          aDest(i,k) = aSource(i,k)*(1.0 - aElli(i).ratio())*fabs(aElli(i).vecmax(k)/rmax);
      }
    };

  /** 3D Notrmal-Ratio Filter *************************************************/
  void eval3DFilter(PPPSpectrContainer<PPPcomplex> &aWgFiltered, PPPSpectrContainer<PPPcomplex> &aWg,
    PPPSpectrContainer<PPPellipse3D> &aElli, PPPEllipse3Dfilter &aFilter, unsigned aChann) {
    register unsigned i,j,k,numb=0;
    for(i=0; i<aElli.voices(); i++)
      for(j=0; j<aElli.points(); j++)
        {
        if(aFilter.getType() == PPPEllipse3Dfilter::NRFNorm)
          {
          aWgFiltered(i,j,3*aChann+0) = aWg(i,j,0)*aElli(i,j).planarity_cosx()/(M_PI/2.0);
          aWgFiltered(i,j,3*aChann+1) = aWg(i,j,1)*aElli(i,j).planarity_cosy()/(M_PI/2.0);
          aWgFiltered(i,j,3*aChann+2) = aWg(i,j,2)*aElli(i,j).planarity_cosz()/(M_PI/2.0);
          }
        else if(aFilter.isContent(aElli(i,j)))
          {
          aWgFiltered(i,j,3*aChann+0) = aWg(i,j,0);
          aWgFiltered(i,j,3*aChann+1) = aWg(i,j,1);
          aWgFiltered(i,j,3*aChann+2) = aWg(i,j,2);
          }
        else
          {
          for (k=0; k<3; k++) aWgFiltered(i,j,3*aChann+k) = 0.0;
          numb++;
          }
        }
    strstream str;
    str << aFilter.getInfo();
    if(aFilter.getType() != PPPEllipse3Dfilter::NRFNorm)
      {
      long cnt = aElli.voices()*aElli.points();
      str << endl << "  " << PPPTRANSELLI_DEL << ": " << numb << "/" << cnt << " (" << 100.0*numb/cnt << "%)";
      }
    str << ends;  
    onMessage(str.str());
    };

  }; // end of object

#endif
