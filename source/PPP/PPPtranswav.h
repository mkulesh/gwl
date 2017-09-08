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

#ifndef _PPPTRANSWAVELET
#define _PPPTRANSWAVELET

#define PPPTRANSWAVELET_NAME      "wavelet transformation"
#define PPPTRANSWAVELET_CALC      "calculation wavelet transform"
#define PPPTRANSWAVELET_NSWT      "time integral"
#define PPPTRANSWAVELET_NFWT      "fft based"
#define PPPTRANSWAVELET_NCWT      "multi-convolution method"
#define PPPTRANSWAVELET_CURR      "current frequency"
#define PPPTRANSWAVELET_CURRP     "current point"
#define PPPTRANSWAVELET_SWTTNAME  "mother wavelet time axis"
#define PPPTRANSWAVELET_SWTFNAME  "mother wavelet frequency axis"
#define PPPTRANSWAVELET_SWTTRACE  "trace information"
#define PPPTRANSWAVELET_IWT       "calculation inverse wavelet transform:"
#define PPPTRANSWAVELET_IFWT      "calculation fast inverse wavelet transform:"
#define PPPTRANSWAVELET_RIDCALC   "calculation Ridge function"
#define PPPTRANSWAVELET_RIDNAME   "Ridge function"
#define PPPTRANSWAVELET_PAUNW2D   "calculation 2-D phase unwrapping"
#define PPPTRANSWAVELET_LIMCALC   "calculation of local intermittency measure"
#define PPPTRANSWAVELET_LIMNAME   "local intermittency measure"
#define PPPTRANSWAVELET_ERRWT     "undefined transform's type in procedure: "
#define PPPTRANSWAVELET_POS       "positive spectrum"
#define PPPTRANSWAVELET_NEG       "negative spectrum"
#define PPPTRANSWAVELET_FULL      "full spectrum"
#define PPPTRANSWAVELET_WAVPAR    "parameter of wavelet"
#define PPPTRANSWAVELET_WAVCUTOFF "wavelet cutoff"
#define PPPTRANSWAVELET_AMPLF     "amplitude factor"

/************************************************************************
 * PPPTransWavelet
 * in this class was implemented procedures to Fourie spectrum
 * calculate and processing
 ***********************************************************************/
template<class AType> class PPPTransWavelet : public PPPTransFour<AType>
  {
  private:

     char _currFreq[100];
     string _tracename;
     PPPConvolutor<PPPcomplex>* _convolutor;

  public:

     typedef enum {ATtimeInt, ATfft, ATconv} CWTAlgType;
     typedef enum {TFTtime, TFTfreq} TimeFreqType;

    PPPTransWavelet(void) {
      _tracename = "NULL";
      PPPBaseTemplate<AType>::setObjectName(PPPTRANSWAVELET_NAME);
      };

     void setTraceName(const string &aName) {
      _tracename = aName;
      };

   /** Wavelet Transformation **************************************************/
    void WT(PPPSpectrContainer<PPPcomplex> &aCWT, PPPSignalContainer<AType> &aVoice, PPPSpectrParams &aT){
      // wavelet
      PPPWavelet *wav = aT.getWavelet(PPPSpectrParams::WSTdirect);
      // log
      strstream str;
      switch(aT.getTransformType())
        {
        case ATtimeInt:
          str << PPPTRANSWAVELET_CALC << " (" << PPPTRANSWAVELET_NSWT <<"): "
              << PPPBaseTemplate<AType>::getTypeName() << "[" << aVoice.channels() << "]" << endl
              << wav->getInfo() << endl
              << "  " << PPPTRANSWAVELET_WAVCUTOFF << " = " << aT.getCutoffPrec();
          if(wav->hasParams())
          str << endl << "  " << PPPTRANSWAVELET_WAVPAR << " = " << aT.getWaveletPar(PPPSpectrParams::WSTdirect);
          str << ends;
          break;
        case ATfft:
          str << PPPTRANSWAVELET_CALC << " (" << PPPTRANSWAVELET_NFWT <<"): "
              << PPPBaseTemplate<AType>::getTypeName() << "[" << aVoice.channels() << "]" << endl
              << wav->getInfo();
          if(wav->hasParams())
          str << endl << "  " << PPPTRANSWAVELET_WAVPAR << " = " << aT.getWaveletPar(PPPSpectrParams::WSTdirect);
          str << ends;
          break;
        case ATconv:
          str << PPPTRANSWAVELET_CALC << " (" << PPPTRANSWAVELET_NCWT <<"): "
              << PPPBaseTemplate<AType>::getTypeName() << "[" << aVoice.channels() << "]" << endl
              << wav->getInfo();
          if(wav->hasParams())
          str << endl << "  " << PPPTRANSWAVELET_WAVPAR << " = " << aT.getWaveletPar(PPPSpectrParams::WSTdirect);
          str << ends;
          break;
        default:
          PPPBaseObject::onError(PPPTRANSWAVELET_ERRWT+string("WT"));
        }
      PPPBaseObject::onMessage(str.str());
      // clock
      clock_t cThen1 = clock();
      // wavelet container
      aCWT.prepare(aT.voices(),aVoice.points(),aVoice.channels(),aVoice.getAxis(),aT.getFreq(),PPPSPECTRCONTAINER_NAME);
      // trace
      bool aTrace = (_tracename != "NULL");
      PPPSignalContainer<double> _trace;
      if(aTrace)
        _trace.prepare(aCWT.voices(), 1, aCWT.getFreq(), PPPTRANSWAVELET_SWTTRACE);
      // Fourier spectrum of source signal
      PPPSignalContainer<PPPcomplex> aFour;
      // wavelet axis
      PPPAxis wt;
      if(aT.getTransformType() == ATtimeInt)
        {
        _evalTimeAxis(wt,aCWT.getTime());
        }
      else if(aT.getTransformType() == ATfft)
        {
        bool aSaveStatus = PPPBaseTemplate<AType>::getShowMessage();
        PPPBaseTemplate<AType>::setShowMessage(false);
        PPPTransFour<AType>::FT(aFour,aVoice);
        PPPBaseTemplate<AType>::setShowMessage(aSaveStatus);
        _evalFrequencyAxis(wt, aFour.getAxis());
        }
      else
        {
        _evalTimeAxis(wt,aCWT.getTime());
        int nsig0 = 0;
        int nsig1 = (int)aCWT.points()-1;
        int nsag0 = -((int)(wt.size()-1)/2);
        int nsag1 = (int)(wt.size()-1)/2;
        int nsug0 = 0;
        int nsug1 = (int)aCWT.points()-1;
        _convolutor = new PPPConvolutor<PPPcomplex>(
            nsig0,  // (in) lower index of first signal (out) actually needed lower index
            nsig1,  // (in) upper index of first signal (out) actually needed upper index
            nsag0,  // (in) lower index of second signal (out) actually needed lower index
            nsag1,  // (in) upper index of second signal (out) actually needed upper index
            nsug0,  // (in) lower index of result signal
            nsug1,  // (in) upper index of result signal
            !PPPConvolutor<PPPcomplex>::UNKNOWN_LIMIT, // take all of sig that is needed
            !PPPConvolutor<PPPcomplex>::UNKNOWN_LIMIT,  // take all of sag in the limits specified
            aCWT.channels(),
            1);
        for(unsigned c=0; c<aCWT.channels(); c++)
          _convolutor->getSig(c).setSignalView(aVoice.getChannel(c), wt.getSamplingPeriod());
        }
      clock_t cThen2 = clock();
      // transformation
      clock_t cThen3 = 0, cThen4 = 0;
      for(unsigned i=0; i<aCWT.voices(); i++)
        {
        sprintf(_currFreq,"%s: %3d",PPPTRANSWAVELET_CURR,i+1);
        PPPBaseTemplate<AType>::onProgress(100*i/(aCWT.voices()-1),_currFreq);
        switch(aT.getTransformType())
          {
          case ATtimeInt:
               cThen3 = _evalVoiceTimeInt(aCWT, i, aVoice, wt, *wav, aT.getCutoffPrec());
               break;
          case ATfft:
               cThen3 = _evalVoiceConvAnalyt(aCWT, i, aFour, wt, *wav);
               break;
          case ATconv:
               cThen3 = _evalVoiceConvFft(aCWT, i, aVoice, wt, *wav);
               break;
          }
        cThen4 += cThen3;
        if(aTrace)
          _trace(i) = cThen3;
        }
      // clock information
      strstream strcl;
      strcl << "preparation time: " << (cThen2-cThen1) << ", calculation time: " << cThen4 << " ticks" << ends;
      PPPBaseObject::onMessage(strcl.str());
      // clear of convolution object
      if(aT.getTransformType() == ATconv)
        {
        delete _convolutor;
        }
      // finalizing the trace
      if(aTrace)
        _trace.write(_tracename);
      PPPBaseTemplate<AType>::onProgress(-1,"");
      }

    /** Inverse Wavelet Transformation ******************************************/
    void IWT(PPPSignalContainer<AType> &aVoice, PPPSpectrContainer<PPPcomplex> &aCWT, PPPSpectrParams &aT, bool aAmpl=true) {
      if(aT.getWavelet(PPPSpectrParams::WSTinverse)->getType() == PPPWavelet::CWdelta)
        {
        IFWT(aVoice, aCWT, aT, aAmpl);
        return;
        }
      strstream str;
      PPPWavelet *invwav = aT.getWavelet(PPPSpectrParams::WSTinverse);
      PPPcomplex Chm(1.0,0.0);
      if(aAmpl) Chm = aT.getInverseConstant();
      str << PPPTRANSWAVELET_IWT << " "
          << PPPBaseTemplate<AType>::getTypeName() << "[" << aCWT.channels() << "]" << endl
          << aT.getWavelet(PPPSpectrParams::WSTinverse)->getInfo() << endl
          << "  " << PPPTRANSWAVELET_WAVPAR << " = " << aT.getWaveletPar(PPPSpectrParams::WSTdirect) << endl
          << "  " << PPPTRANSWAVELET_AMPLF << " = " << Chm << endl
          << "  " << PPPTRANSWAVELET_WAVCUTOFF << " = " << aT.getCutoffPrec() << ends;
      PPPBaseObject::onMessage(str.str());
      aVoice.prepare(aCWT.points(),1,aCWT.getTime(),PPPTRANSFOUR_INV);
      unsigned i,j,k,n,k1,k2;
      double f,df,db;
      PPPcomplex s0,s1,s2,s3;
      db = aCWT.getTime().getSamplingPeriod();
      // transformation
      for(unsigned c=0; c<aVoice.channels(); c++)
        for(j=0; j<aCWT.points(); j++)
          {
          sprintf(_currFreq,"%s: [%d] %3d",PPPTRANSWAVELET_CURRP,c+1,j+1);
          PPPBaseTemplate<AType>::onProgress(100*j/(aCWT.points()-1),_currFreq);
          for(s2=0.0, i=1; i<aCWT.voices(); i++)
            {
            f = aCWT.getFreq(i);
            df = fabs(aCWT.getFreq(i-1) - f);
            invwav->setFrequency(f);
            if(aT.getCutoffPrec() == 0.0)
              n = aCWT.points();
            else  
              n = (unsigned)(invwav->getCutoffTimeNorm(aT.getCutoffPrec())/db);
            k1 = (j<n)? 0 : j-n;
            k2 = (j+n>aCWT.points())? aCWT.points() : j+n;
            for(s1=0.0, k=k1; k<k2; k++)
              {
              s0 = invwav->evalCmplTime(aCWT.getTime(j)-aCWT.getTime(k));
              s1 += aCWT(i,k,c)*s0/f;
              }
            s2 += s1*df;
            }
          s3 = s2*db/Chm;
          if(aT.getFreq().getSign()!=PPPAxis::ASfull) s3 = s3*2.0;
          PPPBaseTemplate<AType>::cmplConvert(aVoice(j,c),s3);
          }
      PPPBaseTemplate<AType>::onProgress(-1,"");
      };

    /** Fast Inverse Wavelet Transformation *************************************/
    void IFWT(PPPSignalContainer<AType> &aVoice, PPPSpectrContainer<PPPcomplex> &aCWT, PPPSpectrParams &aT, bool aAmpl=true) {
      strstream str;
      PPPcomplex Chm(1.0,0.0);
      if(aAmpl) Chm = aT.getInverseConstant();
      str << PPPTRANSWAVELET_IFWT << " "
          << PPPBaseTemplate<AType>::getTypeName() << "[" << aCWT.channels() << "]" << endl
          << "  " << PPPTRANSWAVELET_AMPLF << " = " << Chm << ends;
      PPPBaseObject::onMessage(str.str());
      aVoice.prepare(aCWT.points(),aCWT.channels(),aCWT.getTime(),PPPTRANSFOUR_INV);
      unsigned k,j;
      double f,df;
      AType val;
      PPPcomplex z1;
      for(unsigned c=0; c<aCWT.channels(); c++)
        for(k=0; k<aCWT.points(); k++)
          {
          for(aVoice(k,c)=0, j=1; j<aCWT.voices(); j++)
            {
            f = aCWT.getFreq(j);
            df = fabs(aCWT.getFreq(j-1) - f);
            if(f>0.0) z1 = 1.0*aCWT(j,k,c)*df/(f*Chm);
                 else z1 = -1.0*aCWT(j,k,c)*df/(f*Chm);
            PPPBaseTemplate<AType>::cmplConvert(val, z1);
            aVoice(k,c) += val;
            }
          if(aT.getFreq().getSign()!=PPPAxis::ASfull) aVoice(k,c) = 2.0*aVoice(k,c);
          }
      };


    /** Cating and separating of wavelet spectrum *******************************/
    void WaveletCat(PPPSpectrContainer<AType> &aDest, PPPSpectrContainer<AType> &aPlus, PPPSpectrContainer<AType> &aMinus) {
      if ((aMinus.voices() != aPlus.voices()) || (aMinus.points() != aPlus.points()) || (aMinus.channels() != aPlus.channels()))
        PPPBaseObject::onError(ARG_VALUE+string("WaveletCat"));
      PPPAxis ax;
      ax.assign(aPlus.getFreq());
      ax.fillsign(PPPAxis::ASfull);
      aDest.prepare(2*aPlus.voices(),aPlus.points(),aPlus.channels(),aPlus.getTime(),ax,PPPTRANSWAVELET_FULL);
      unsigned i,j;
      for(unsigned c=0; c<aPlus.channels(); c++)
        {
        for(i=0; i<aMinus.voices(); i++)
          for (j=0; j<aMinus.points(); j++)
            aDest(i,j,c) = aMinus(aMinus.voices()-1-i,j,c);
        for(i=0; i<aPlus.voices(); i++)
          for (j=0; j<aPlus.points(); j++)
            aDest(i+aDest.voices()/2,j,c) = aPlus(i,j,c);
        }
      return;
      };

    void WaveletSeparate(PPPSpectrContainer<AType> &aPlus, PPPSpectrContainer<AType> &aMinus,
      PPPSpectrContainer<AType> &aSource) {
      if(aSource.getFreq(0)*aSource.getFreq(aSource.voices()-1) >= 0)
        PPPBaseObject::onError(ARG_VALUE+string("WaveletSeparate"));
      unsigned i,j,k;
      unsigned count = aSource.voices()/2;
      // minus spectr
      aMinus.resize(count,aSource.points(),aSource.channels());
      aMinus.setNames(PPPTRANSWAVELET_NEG,PPPSPECTRCONTAINER_TIME,PPPSPECTRCONTAINER_FREQ);
      aMinus.setTime(aSource.getTime());
      aMinus.getFreq().assign(aSource.getFreq(),PPPAxis::ASminus);
      for(unsigned c=0; c<aSource.channels(); c++)
        for(i=0,k=count-1; i<count; i++,k--)
          for (j=0; j<aSource.points(); j++)
            aMinus(k,j,c) = aSource(i,j,c);
      // plus spectr
      aPlus.resize(count,aSource.points(),aSource.channels());
      aPlus.setNames(PPPTRANSWAVELET_POS,PPPSPECTRCONTAINER_TIME,PPPSPECTRCONTAINER_FREQ);
      aPlus.setTime(aSource.getTime());
      aPlus.getFreq().assign(aSource.getFreq(),PPPAxis::ASplus);
      for(unsigned c=0; c<aSource.channels(); c++)
        for(i=count,k=0; i<aSource.voices(); i++,k++)
          for (j=0; j<aSource.points(); j++)
            aPlus(k,j,c) = aSource(i,j,c);
      return;
      };

    /** 2D phase unwpaping ******************************************************/
    void PhaseUnwrap2D(PPPMatrixContainer<double> &aDest, PPPMatrixContainer<double> &aSource) {
      PPPBaseObject::onMessage(PPPTRANSWAVELET_PAUNW2D);
      unsigned i,j,k;
      double d;
      PPPMatrixContainer<double> trans;
      trans.assign(aSource);
      aDest.resize(aSource.rows(), aSource.cols());
      for(i=0; i<trans.rows(); i++)
        for(j=1; j<trans.cols(); j++)
        {
        d = trans(i,j) - trans(i,j-1);
        if(fabs(d) >= M_PI)
          {
          for(k=j; k<trans.cols(); k++)
            trans(i,k) = trans(i,k) - 2.0*M_PI*d/fabs(d);
          }
        aDest(i,j) = trans(i,j);
        }
      };

    /** Ridge Calculation *******************************************************/
    void Ridge(PPPSignalContainer<double> &aRidge, PPPSpectrContainer<PPPcomplex> &aCWT,
      TimeFreqType aType, PPPAxis::AxisSign aT) {
      if(aT == PPPAxis::ASfull) PPPBaseObject::onError(ARG_VALUE+string("Ridge"));
      strstream str;
      double dffmax;
      unsigned i,j,imax;
      switch(aType)
        {
        case TFTtime:
          str << PPPTRANSWAVELET_RIDCALC << ": " << PPPSIGNALCONTAINER_AXISNAME << "-" << PPPSPECTRCONTAINER_FREQ << ends;
          PPPBaseObject::onMessage(str.str());
          aRidge.prepare(aCWT.points(),1,aCWT.getTime(),PPPTRANSWAVELET_RIDNAME);
          for(j=0; j<aCWT.points(); j++)
            {
            dffmax = abs(aCWT(0,j));
            imax = 0;
            for(i=1; i<aCWT.voices(); i++)
              {
              if(abs(aCWT(i,j)) > dffmax)
                {
                dffmax = abs(aCWT(i,j));
                imax = i;
                }
              }
            aRidge(j) = aCWT.getFreq(imax);
            }
          break;
        case TFTfreq:
          str << PPPTRANSWAVELET_RIDCALC << ": " << PPPSPECTRCONTAINER_FREQ << "-" << PPPSIGNALCONTAINER_AXISNAME << ends;
          PPPBaseObject::onMessage(str.str());
          aRidge.resize(aCWT.voices());
          aRidge.setNames(PPPTRANSWAVELET_RIDNAME,PPPSPECTRCONTAINER_FREQ);
          for(j=0; j<aCWT.voices(); j++)
            {
            dffmax = abs(aCWT(j,0));
            imax = 0;
            for(i=1; i<aCWT.points(); i++)
              {
              if(abs(aCWT(j,i)) > dffmax)
                {
                dffmax = abs(aCWT(j,i));
                imax = i;
                }
              }
            aRidge.getAxis()[j] = (aT == PPPAxis::ASplus)? aCWT.getFreq(j) : aCWT.getFreq(aCWT.voices()-1-j);
            aRidge(j) = aCWT.getTime(imax);
            }
          break;
        }
      };

    /** Computation of instantaneous frequency from wavelet transform ***********/
    void LIM(PPPSpectrContainer<double> &aDest, PPPSignalContainer<double> &aWavPower,
      PPPSpectrContainer<PPPcomplex> &aCWT, unsigned aCat = 0) {
      PPPBaseObject::onMessage(PPPTRANSWAVELET_LIMCALC);
      aWavPower.prepare(aCWT.voices(),1,aCWT.getFreq(),PPPTRANSFOUR_FOUR);
      unsigned i,j;
      for(j=0; j<aCWT.voices(); j++)
        {
        double p = 0.0;
        for(i=aCat; i<(aCWT.points()-aCat); i++) p += pow(abs(aCWT(j,i)), 2.0);
        aWavPower(j) = p/((double)(aCWT.voices()));
        }
      aDest.prepare(aCWT);
      aDest.setObjectName(PPPTRANSWAVELET_LIMNAME);
      for(j=0; j<aCWT.voices(); j++)
        for(i=0; i<aCWT.points(); i++)
          aDest(j,i) = pow(abs(aCWT(j,i)), 2.0)/aWavPower(j);
      return;
      };


  private:

    void _evalTimeAxis(PPPAxis &aWaxis, PPPAxis &aTime) {
      unsigned slength = aTime.size();
      aWaxis.resize(2*slength-1);
      aWaxis.setObjectName(PPPTRANSWAVELET_SWTTNAME);
      for(unsigned j=0; j<slength; j++)
        {
        aWaxis[slength-1-j] = aTime[0]-aTime[j];
        aWaxis[slength-1+j] = aTime[j]-aTime[0];
        }
      aWaxis.setParams(PPPAxis::ATlin);
      return;
      };

    clock_t _evalVoiceTimeInt(PPPSpectrContainer<PPPcomplex> &aCWT, unsigned const aInd,
      PPPSignalContainer<AType> &aVoice, PPPAxis &aWaxis, PPPWavelet &aWav, double aCutOff=0.0)
      {
      clock_t cThen1 = clock();
      PPPVectorContainer<PPPcomplex> w;
      aWav.setFrequency(aCWT.getFreq(aInd));
      aWav.evalTimeRepresentation(w, aWaxis, aCutOff);
      for(unsigned c=0; c<aCWT.channels(); c++)
        {
        for(unsigned j=0; j<aCWT.points(); j++)
          aCWT(aInd,j,c) = _evalIntegralInTime(aVoice.getChannel(c), j, w, aWaxis.getSamplingPeriod());
        }
      clock_t cThen2 = clock();
      return cThen2-cThen1;
      };

    PPPcomplex _evalIntegralInTime(PPPVectorContainer<AType> &aSig, unsigned const aInd,
      PPPVectorContainer<PPPcomplex> aWav, double aDeltaTime) {
      unsigned slength = aSig.size(); // the length of the source signal
      unsigned nmax = aWav.size(); // the length of the wavelet
      unsigned wcen = nmax/2; // the index of the centrum of the wavelet
      PPPcomplex s1(0.0,0.0);
      int sind;
      register unsigned n,k;
      for(n=0; n<nmax; n++)
        {
        sind = aInd-wcen+n;
        if(sind>=0 && sind<(int)slength)
          {
          k = (unsigned)sind;
          s1 += aSig[k]*conj(aWav[n])*aDeltaTime;
          }
        }
      return s1;
      };

    void _evalFrequencyAxis(PPPAxis &aWaxis, PPPAxis &aFreq) {
      unsigned slength = aFreq.size();
      aWaxis.resize(slength);
      aWaxis.setObjectName(PPPTRANSWAVELET_SWTFNAME);
      for(unsigned j=0; j<=slength/2; j++)
        {
        if(j>0) aWaxis[slength-j] = -aFreq[j];
        aWaxis[j] = aFreq[j];
        }
      aWaxis.setParams(PPPAxis::ATlin);
      return;
      };

    clock_t _evalVoiceConvAnalyt(PPPSpectrContainer<PPPcomplex> &aCWT, unsigned const aInd,
      PPPSignalContainer<PPPcomplex> &aFour, PPPAxis &aWaxis, PPPWavelet &aWav)
      {
      clock_t cThen1 = clock();
      PPPVectorContainer<PPPcomplex> w;
      PPPVectorContainer<PPPcomplex> Hn(aCWT.points()),Sn;
      PPPFft<PPPcomplex> trans;
      aWav.setFrequency(aCWT.getFreq(aInd));
      aWav.evalFreqRepresentation(w, aWaxis);
      for(unsigned c=0; c<aCWT.channels(); c++)
        {
        aCWT.getChannel(c).getRow(Sn,aInd);
        for(unsigned j=0; j<aCWT.points(); j++) Hn[j] = aFour(j,c)*conj(w[j]);
        trans.ifft(Sn,Hn);
        Sn.unlink();
        }
      clock_t cThen2 = clock();
      return cThen2-cThen1;
      };

    clock_t _evalVoiceConvFft(PPPSpectrContainer<PPPcomplex> &aCWT, unsigned const aInd,
      PPPSignalContainer<AType> &aVoice, PPPAxis &aWaxis, PPPWavelet &aWav)
      {
      clock_t cThen1 = clock();
      PPPVectorContainer<PPPcomplex> Sn;
      aWav.setFrequency(aCWT.getFreq(aInd));
      aWav.evalTimeRepresentation(_convolutor->getSagView(), aWaxis, 0.0, false);
      _convolutor->finishSag();
      for(unsigned c=0; c<aCWT.channels(); c++)
        {
        aCWT.getChannel(c).getRow(Sn,aInd);
        _convolutor->doConvolution();
        _convolutor->fillSug(Sn,PPPConvolutor<PPPcomplex>::ZERO_PADDING,c);
        Sn.unlink();
        }
      clock_t cThen2 = clock();
      return cThen2-cThen1;
      };

  }; // end of object


#endif
