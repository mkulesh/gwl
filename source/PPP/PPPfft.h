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
#ifndef _PPPFTH
#define _PPPFTH

#define PPPFFT_NAME        "Fast fourier transformation"
#define PPPFFT_INVALIDSIZE "invalid size has been detected for input or output vector in procedure: "
#define PPPFFT_ERRINIT     "FFT initialization failure in procedure: "

/************************************************************************
 * PPPFft: Interafce object for fast fourier transform
 * Any implementation of this FFT must implement the following methods
 *   1. initforward(...) or initbackward(...)
 *   2. execute()
 *   3. clear()
 * We can also use fft() or ifft() procedures which include these methods 
 ***********************************************************************/
template<class AType> class PPPFft : public PPPBaseTemplate<AType>
{
private:
    
    PPPMathFunc _math;
    UTFFT _trans;

public:
    
    PPPFft (void)
    {
        PPPBaseTemplate<AType>::setObjectName(PPPFFT_NAME);
    }
    
    void initforward (PPPVectorContainer<PPPcomplex> &aFour, PPPVectorContainer<AType> &aVoice)
    {
        if (aVoice.size() == 0) PPPBaseTemplate<AType>::onError(
                PPPFFT_INVALIDSIZE + string("initforward()"));
        if (aVoice.size() != aFour.size()) PPPBaseTemplate<AType>::onError(
                PPPFFT_INVALIDSIZE + string("initforward()"));
        if (!_math.isPowerOfTwo(aVoice.size())) PPPBaseTemplate<AType>::onError(
                ARG_POW2 + string("initforward()"));
        bool status = false;
        if (aVoice.isReal())
            status = _trans.init((double*) aVoice.begin(), (double*) aFour.begin(), aVoice.size(),
                                 _trans.R2C);
        else
            status = _trans.init((double*) aVoice.begin(), (double*) aFour.begin(), aVoice.size(),
                                 _trans.C2CFORWARD);
        if (!status) PPPBaseTemplate<AType>::onError(PPPFFT_ERRINIT + string("initforward()"));
    }
    
    void initbackward (PPPVectorContainer<AType> &aVoice, PPPVectorContainer<PPPcomplex> &aFour)
    {
        if (aFour.size() == 0) PPPBaseTemplate<AType>::onError(
                PPPFFT_INVALIDSIZE + string("initbackward()"));
        if (aVoice.size() != aFour.size()) PPPBaseTemplate<AType>::onError(
                PPPFFT_INVALIDSIZE + string("initbackward()"));
        if (!_math.isPowerOfTwo(aVoice.size())) PPPBaseTemplate<AType>::onError(
                ARG_POW2 + string("initbackward()"));
        bool status = false;
        if (aVoice.isReal())
            status = _trans.init((double*) aFour.begin(), (double*) aVoice.begin(), aVoice.size(),
                                 _trans.C2R);
        else
            status = _trans.init((double*) aFour.begin(), (double*) aVoice.begin(), aVoice.size(),
                                 _trans.C2CBACKWARD);
        if (!status) PPPBaseTemplate<AType>::onError(PPPFFT_ERRINIT + string("initbackward()"));
    }
    
    bool execute (void)
    {
        return _trans.execute();
    }
    
    bool clear (void)
    {
        return _trans.clear();
    }
    
    inline string &getMethodName (void)
    {
        return _trans.getMethodName();
    }
    
    void fft (PPPVectorContainer<PPPcomplex> &aFour, PPPVectorContainer<AType> &aVoice)
    {
        initforward(aFour, aVoice);
        execute();
        clear();
    }
    
    void ifft (PPPVectorContainer<AType> &aVoice, PPPVectorContainer<PPPcomplex> &aFour)
    {
        initbackward(aVoice, aFour);
        execute();
        clear();
    }
    
};
// end of object

#endif /* _PPPFTH */
