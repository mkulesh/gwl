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
#ifndef _PPPTRANSFOUR
#define _PPPTRANSFOUR

#define PPPTRANSFOUR_NAME      "Fourier transformation"
#define PPPTRANSFOUR_FOUR      "Fourier spectrum"
#define PPPTRANSFOUR_CURR      "current frequency"
#define PPPTRANSFOUR_INV       "inverse signal"
#define PPPTRANSFOUR_CROSS     "cross-correlation spectrum"
#define PPPTRANSFOUR_ICROSS    "inverse cross-correlation"
#define PPPTRANSFOUR_FT        "calculation Fourier integral:"
#define PPPTRANSFOUR_IFT       "calculation inverse Fourier integral:"
#define PPPTRANSFOUR_FFT       "calculation fast Fourier transform:"
#define PPPTRANSFOUR_IFFT      "calculation inverse fast Fourier transform:"
#define PPPTRANSFOUR_FFILT     "calculation band pass filter:"
#define PPPTRANSFOUR_FHT       "calculation Hilbert transform:"
#define PPPTRANSFOUR_CCF       "calculation Cross-correlation spectrum:"
#define PPPTRANSFOUR_ICCF      "calculation inverce transform of cross-correlation spectrum:"
#define PPPTRANSFOUR_CCS       "calculation cross-correlation seismogram:"
#define PPPTRANSFOUR_ICCS      "calculation inverce transform of cross-correlation seismogram:"
#define PPPTRANSFOUR_NOTIDENT  "input signals must have equal dimensions in procedure: "

/************************************************************************
 * PPPTransFour
 * in this class was implemented procedures to Fourie spectrum
 * calculate and processing
 ***********************************************************************/
template<class AType> class PPPTransFour : public PPPBaseTemplate<AType>
{
protected:

    unsigned _transType;
    PPPMathFunc _math;

public:
    
    typedef enum
    {
        FTTauto, FTTcft, FTTfft
    } FTType;

    PPPTransFour (void)
    {
        PPPBaseObject::setObjectName(PPPTRANSFOUR_NAME);
        _transType = FTTauto;
    }
    
    PPPTransFour (FTType aType)
    {
        PPPBaseObject::setObjectName(PPPTRANSFOUR_NAME);
        _transType = aType;
    }
    
    /** Fourier transform ********************************************************/
    void FT (PPPSignalContainer<PPPcomplex> &aFour, PPPSignalContainer<AType> &aVoice)
    {
        PPPSignalContainer<PPPcomplex> aVoiceTmp;
        switch (_transType)
        {
        case FTTauto:
            if (_math.isPowerOfTwo(aVoice.points()))
                FFT(aFour, aVoice);
            else
                CFT(aFour, aVoice);
            break;
        case FTTcft:
            CFT(aFour, aVoice);
            break;
        case FTTfft:
            FFT(aFour, aVoice);
            break;
        };
        return;
    }
    
    /** Inverse Fourier transform ************************************************/
    void IFT (PPPSignalContainer<AType> &aVoice, PPPSignalContainer<PPPcomplex> &aFour, double aMin,
              double aMax)
    {
        switch (_transType)
        {
        case FTTauto:
            if (_math.isPowerOfTwo(aFour.points()))
                IFFT(aVoice, aFour, aMin);
            else
                ICFT(aVoice, aFour, aMin);
            break;
        case FTTcft:
            ICFT(aVoice, aFour, aMin);
            break;
        case FTTfft:
            IFFT(aVoice, aFour, aMin);
            break;
        };
        return;
    }
    
    /** Modify Fourier Spectra **************************************************/
    void FTShift (PPPSignalContainer<AType> &aFour, bool aModifyAxis = true)
    {
        if (aModifyAxis)
        {
            double aMax = aFour.getAxis().getMax();
            if (aFour.getAxis().getMin() >= 0)
                aFour.getAxis().fill(-aMax / 2.0, aMax / 2.0, PPPAxis::ATlin);
            else
                aFour.getAxis().fill(0.0, 2.0 * aMax, PPPAxis::ATlin);
        }
        unsigned halb = aFour.points() / 2;
        for (unsigned i = 0; i < halb; i++)
            for (unsigned c = 0; c < aFour.channels(); c++)
                SWAP(aFour(i, c), aFour(i + halb, c));
    }
    
    void FTSeparate (PPPSignalContainer<AType> &aPlus, PPPSignalContainer<AType> &aMinus,
                     PPPSignalContainer<AType> &aFour)
    {
        double aMax = aFour.getAxis().getMax() / 2.0;
        PPPAxis ax(aFour.points() / 2, 0.0, aMax, PPPAxis::ATlin, PPPSIGNALCONTAINER_AXISNAME);
        aPlus.prepare(aFour.points() / 2, aFour.channels(), ax, PPPTRANSFOUR_FOUR);
        aMinus.prepare(aFour.points() / 2, aFour.channels(), ax, PPPTRANSFOUR_FOUR);
        for (unsigned c = 0; c < aFour.channels(); c++)
            for (unsigned i = 0; i < aFour.points() / 2; i++)
            {
                aPlus(i, c) = aFour(i, c);
                if (i == 0)
                    aMinus(i, c) = 0.0;
                else
                    aMinus(i, c) = aFour(aFour.points() - i, c);
            }
        return;
    }
    
    void FTCut (PPPSignalContainer<AType> &aFour, PPPSignalContainer<AType> &aPlus,
                PPPSignalContainer<AType> &aMinus)
    {
        if (aPlus.points() != aMinus.points() || aPlus.channels() != aMinus.channels()) PPPBaseTemplate<
                AType>::onError(PPPTRANSFOUR_NOTIDENT + string("FTCut"));
        unsigned flag = (2 * (unsigned) (aFour.points() / 2) == aFour.points()) ? 0 : 1;
        for (unsigned c = 0; c < aFour.channels(); c++)
            for (unsigned i = 0; i < aFour.points() / 2; i++)
            {
                aFour(i, c) = aPlus(i, c);
                if (i == 0)
                    aFour(aFour.points() / 2 + i, c) = 0.0;
                else
                    aFour(aFour.points() / 2 + i + flag, c) = aMinus(aFour.points() / 2 - i, c);
            }
        return;
    }
    
    /** Hilbert Transforms ******************************************************/
    void HT (PPPSignalContainer<PPPcomplex> &aVoice1, PPPSignalContainer<AType> &aVoice,
             PPPAxis::AxisSign aT)
    {
        strstream str;
        str << PPPTRANSFOUR_FHT << " " << PPPBaseTemplate<AType>::getTypeName() << "["
        << aVoice.channels() << "]" << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        bool aSaveStatus = PPPBaseTemplate<AType>::getShowMessage();
        PPPBaseTemplate<AType>::setShowMessage(false);
        unsigned count = aVoice.points();
        PPPVectorContainer<double> h(count);
        h.assign(0.0);
        switch (aT)
        {
        case PPPAxis::ASplus:
            h[0] = h[count / 2] = 1.0;
            for (unsigned i = 1; i <= count / 2 - 1; i++)
                h[i] = 2.0;
            break;
        case PPPAxis::ASminus:
            h[count - 1] = h[count / 2 + 1] = 1.0;
            for (unsigned i = count / 2 + 2; i <= count - 2; i++)
                h[i] = 2.0;
            break;
        }
        PPPSignalContainer<PPPcomplex> aFour;
        FT(aFour, aVoice);
        for (unsigned c = 0; c < aFour.channels(); c++)
            for (unsigned k = 0; k < aFour.points(); k++)
                aFour(k, c) = aFour(k, c) * h[k];
        PPPTransFour<PPPcomplex> trans_cmpl;
        trans_cmpl.setShowMessage(false);
        trans_cmpl.IFT(aVoice1, aFour, aVoice.getAxis().getMin(), aVoice.getAxis().getMax());
        PPPBaseTemplate<AType>::setShowMessage(aSaveStatus);
        return;
    }
    
    /** Band Pass Filter **********************************************************/
    void BandPassFilter (PPPSignalContainer<AType> &aVoice1, PPPSignalContainer<AType> &aVoice,
                         double aMin, double aMax)
    {
        strstream str;
        str << PPPTRANSFOUR_FFILT << " " << PPPBaseTemplate<AType>::getTypeName() << "[" << aMin
        << "..." << aMax << "]Hz" << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        bool aSaveStatus = PPPBaseTemplate<AType>::getShowMessage();
        PPPBaseTemplate<AType>::setShowMessage(false);
        PPPSignalContainer<PPPcomplex> aFour;
        FT(aFour, aVoice);
        PPPTransFour<PPPcomplex> trans_cmpl;
        trans_cmpl.setShowMessage(false);
        trans_cmpl.FTShift(aFour);
        for (unsigned c = 0; c < aFour.channels(); c++)
            for (unsigned i = 0; i < aFour.points(); i++)
                if (fabs(aFour.getAxis(i)) < aMin || fabs(aFour.getAxis(i)) > aMax) aFour(i, c) =
                        0.0;
        trans_cmpl.FTShift(aFour);
        IFT(aVoice1, aFour, aVoice.getAxis().getMin(), aVoice.getAxis().getMax());
        PPPBaseTemplate<AType>::setShowMessage(aSaveStatus);
        return;
    }
    
    /** Cross-correlation of two signals ****************************************/
    void CCF (PPPSignalContainer<PPPcomplex> &aDest, PPPSignalContainer<AType> &aS1,
              PPPSignalContainer<AType> &aS2)
    {
        strstream str;
        str << PPPTRANSFOUR_CCF << " " << PPPBaseTemplate<AType>::getTypeName() << "["
        << aS1.channels() << "]" << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        if (aS1.points() != aS2.points() || aS1.channels() != aS2.channels()) PPPBaseTemplate<AType>::onError(
                PPPTRANSFOUR_NOTIDENT + string("CCF"));
        bool aSaveStatus = PPPBaseTemplate<AType>::getShowMessage();
        PPPBaseTemplate<AType>::setShowMessage(false);
        // resize of seismoframms
        PPPSignalContainer<AType> aSig1(aS1), aSig2(aS2);
        aSig1.toPowerOfTwo(2 * aS1.points());
        aSig2.toPowerOfTwo(2 * aS1.points());
        // Fourier transforms
        PPPSignalContainer<PPPcomplex> aFour1, aFour2;
        aFour1.setObjectName(PPPTRANSFOUR_CROSS);
        aFour2.setObjectName(PPPTRANSFOUR_CROSS);
        FT(aFour1, aSig1);
        FT(aFour2, aSig2);
        aDest.prepare(aFour1);
        for (unsigned c = 0; c < aFour1.channels(); c++)
            for (unsigned k = 0; k < aFour1.points(); k++)
                aDest(k) = aFour1(k, c) * conj(aFour2(k, c));
        PPPBaseTemplate<AType>::setShowMessage(aSaveStatus);
        return;
    }
    
    /** Inverce Cross-correlation procedure of signal ***************************/
    void ICCF (PPPSignalContainer<AType> &aDest, PPPSignalContainer<PPPcomplex> &aSource,
               double aMin, double aMax)
    {
        strstream str;
        str << PPPTRANSFOUR_ICCF << " " << PPPBaseTemplate<AType>::getTypeName() << "["
        << aSource.channels() << "]" << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        bool aSaveStatus = PPPBaseTemplate<AType>::getShowMessage();
        PPPBaseTemplate<AType>::setShowMessage(false);
        IFT(aDest, aSource, aMin, aMax);
        FTShift(aDest, false);
        aDest.setNames(PPPTRANSFOUR_ICROSS, PPPSIGNALCONTAINER_AXISNAME);
        PPPBaseTemplate<AType>::setShowMessage(aSaveStatus);
        return;
    }
    
    /** Cross-correlation of seismogramm ****************************************/
    void CCS (PPPSignalContainer<PPPcomplex> &aDest, PPPSignalContainer<AType> &aSource,
              PPPMatrixContainer<unsigned> &aNumb)
    {
        strstream str;
        str << PPPTRANSFOUR_CCS << " " << PPPBaseTemplate<AType>::getTypeName() << " ["
        << aNumb.rows() << "]" << endl << "  " << aNumb.matrixToStr() << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        bool aSaveStatus = PPPBaseTemplate<AType>::getShowMessage();
        PPPBaseTemplate<AType>::setShowMessage(false);
        if (aNumb.cols() != 2) PPPBaseTemplate<AType>::onError(ARG_VALUE + string("CCS"));
        PPPSignalContainer<AType> aS1, aS2;
        PPPSignalContainer<PPPcomplex> aSRes;
        aDest.resize(2 * aSource.points(), aNumb.rows());
        for (unsigned i = 0; i < aNumb.rows(); i++)
        {
            aS1.assign(aSource, aNumb(i, 0));
            aS2.assign(aSource, aNumb(i, 1));
            CCF(aSRes, aS1, aS2);
            if (i == 0) aDest.setAxis(aSRes.getAxis());
            aDest.getChannel(i).assign(aSRes.getChannel(0));
        }
        PPPBaseTemplate<AType>::setShowMessage(aSaveStatus);
        return;
    }
    
private:
    
    /** Full Fourie ingegral  ***************************************************/
    void CFT (PPPSignalContainer<PPPcomplex> &aFour, PPPSignalContainer<AType> &aVoice)
    {
        strstream str;
        str << PPPTRANSFOUR_FT << " " << PPPBaseTemplate<AType>::getTypeName() << "["
        << aVoice.channels() << "]" << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        PPPAxis ax(aVoice.points(), 0.0, aVoice.getAxis().getSamplingFreq(), PPPAxis::ATlin,
                   PPPSPECTRCONTAINER_FREQ);
        aFour.prepare(aVoice.points(), aVoice.channels(), ax, PPPTRANSFOUR_FOUR);
        char currFreq[100];
        for (unsigned c = 0; c < aVoice.channels(); c++)
            for (unsigned k = 0; k < aVoice.points(); k++)
            {
                sprintf(currFreq, "%s: [%d] %d", PPPTRANSFOUR_CURR, c + 1, k);
                PPPBaseObject::onProgress(100 * c * k / (aVoice.points() * aVoice.channels()),
                                          currFreq);
                PPPcomplex c2 = 0.0;
                for (unsigned j = 0; j < aVoice.points(); j++)
                    c2 += aVoice(j, c) * exp(
                            PPPcomplex(
                                    0.0,
                                    (-2.0 * M_PI * (double) j * (double) k) / (double) (aVoice
                                            .points())));
                aFour(k, c) = c2;
            }
        PPPBaseTemplate<AType>::onProgress(-1, "");
        return;
    }
    
    /** Inverse full Fourie ingegral ********************************************/
    void ICFT (PPPSignalContainer<AType> &aVoice, PPPSignalContainer<PPPcomplex> &aFour,
               double aMin)
    {
        strstream str;
        str << PPPTRANSFOUR_IFT << " " << PPPBaseTemplate<AType>::getTypeName() << "["
        << aFour.channels() << "]" << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        double aMax = aMin + (double) (aFour.points() - 1) / aFour.getAxis().getMax();
        PPPAxis ax(aFour.points(), aMin, aMax, PPPAxis::ATlin, PPPSIGNALCONTAINER_AXISNAME);
        aVoice.prepare(aFour.points(), aFour.channels(), ax, PPPTRANSFOUR_INV);
        char currFreq[100];
        for (unsigned c = 0; c < aFour.channels(); c++)
            for (unsigned k = 0; k < aFour.points(); k++)
            {
                sprintf(currFreq, "  %s: [%d] %d", PPPTRANSFOUR_CURR, c + 1, k);
                PPPBaseTemplate<AType>::onProgress(
                        100 * c * k / (aFour.points() * aFour.channels()), currFreq);
                PPPcomplex c2 = 0.0;
                for (unsigned j = 0; j < aFour.points(); j++)
                    c2 += aFour(j, c) * exp(
                            PPPcomplex(
                                    0.0,
                                    (2.0 * M_PI * (double) j * (double) k) / (double) aFour.points()))
                          / ((double) aFour.points());
                PPPBaseTemplate<AType>::cmplConvert(aVoice(k, c), c2);
            }
        PPPBaseTemplate<AType>::onProgress(-1, "");
        return;
    }
    
    /** Fast Fourie Transform  **************************************************/
    void FFT (PPPSignalContainer<PPPcomplex> &aFour, PPPSignalContainer<AType> &aVoice)
    {
        PPPFft<AType> trans;
        strstream str;
        str << PPPTRANSFOUR_FFT << " " << trans.getMethodName() << ", "
        << PPPBaseTemplate<AType>::getTypeName() << "[" << aVoice.channels() << "]" << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        PPPAxis ax(aVoice.points(), 0.0, aVoice.getAxis().getSamplingFreq(), PPPAxis::ATlin,
                   PPPSPECTRCONTAINER_FREQ);
        aFour.prepare(aVoice.points(), aVoice.channels(), ax, PPPTRANSFOUR_FOUR);
        for (unsigned c = 0; c < aVoice.channels(); c++)
            trans.fft(aFour.getChannel(c), aVoice.getChannel(c));
        return;
    }
    
    /** Inverse Fast Fourie Transform  ******************************************/
    void IFFT (PPPSignalContainer<AType> &aVoice, PPPSignalContainer<PPPcomplex> &aFour,
               double aMin)
    {
        PPPFft<AType> trans;
        strstream str;
        str << PPPTRANSFOUR_IFFT << " " << trans.getMethodName() << ", "
        << PPPBaseTemplate<AType>::getTypeName() << "[" << aFour.channels() << "]" << ends;
        double aMax = aMin + (double) (aFour.points() - 1) / aFour.getAxis().getMax();
        PPPAxis ax(aFour.points(), aMin, aMax, PPPAxis::ATlin, PPPSIGNALCONTAINER_AXISNAME);
        aVoice.prepare(aFour.points(), aFour.channels(), ax, PPPTRANSFOUR_INV);
        for (unsigned c = 0; c < aFour.channels(); c++)
            trans.ifft(aVoice.getChannel(c), aFour.getChannel(c));
        return;
    }
    
};
// end of object

#endif

