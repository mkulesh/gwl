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
#ifndef _PPPWAVELETS
#define _PPPWAVELETS

#define PPPWAVELETS_ERRCUTOFF1  "epsilon should be between 0 and 1 in procedure: "
#define PPPWAVELETS_ERRCUTOFF2  "epsilon should be positiv in procedure: "
#define PPPWAVELETS_ERRFREQNULL "frequency is 0 in procedure: "
#define PPPWAVELETS_NONE        "wavelet not defined"
#define PPPWAVELETS_MORLET      "complex Morlet wavelet"
#define PPPWAVELETS_MORLETRE    "real Morlet wavelet"
#define PPPWAVELETS_CAUCHY      "complex Cauchy wavelet"
#define PPPWAVELETS_SHANON      "complex Shanon wavelet"
#define PPPWAVELETS_HAAR        "real HAAR wavelet"
#define PPPWAVELETS_DELTA       "delta-function"
#define PPPWAVELETS_ERRWAV      "wavelet undefined in procedure: "
#define PPPWAVELETS_ERRSIG      "size of input signal is invalid in procedure: "
#define PPPWAVELETS_DELTARANGE  0.5

/************************************************************************
 * PPPWavelet
 ***********************************************************************/
class PPPWavelet : public PPPBaseObject
{
public:
    
    typedef enum
    {
        CWnone, 	  // no wavelet
        CWmorlet, 	  // complex Morlet wavelet
        CWmorletre, 	  // real Morlet wavelet
        CWcauchy, 	  // complex Cauchy wavelet
        CWshanon, 	  // complec Shanon wavelet
        CWhaar, 	  // real HAAR wavelet
        CWdelta		  // delta function as indicator for fast inverse transform
    } WaveletType;

protected:
    
    WaveletType _type;   // type of wavelet
    double _f0;     // central frequency of the wavelet
    double _freq;   // frequency of dilation
    double _pos;    // position of translation
    
public:
    
    PPPWavelet ()
    {
        _type = CWnone;
        _freq = 1.0;
        _pos = 0.0;
        _f0 = 1.0;
    }
    
    /**
     *  public virtual methods to be implemeted in a child object
     */
    virtual void setFrequencyResolution (double dfreq)
    {
        return;
    }
    
    // r is the native time variable of the wavelet without transition and dilation
    // we assume that the wavelet has the maximum in the time domain by r=0
    virtual double evalRealTime (double r) const
    {
        return 0.0;
    }
    
    virtual double evalImagTime (double r) const
    {
        return 0.0;
    }
    
    // om is native circular frequency of the wavelet without dilation
    // we assume that the wavelet has the maximum in the frequency domain
    // by the central frequency _f0
    virtual double evalRealFreq (double om) const
    {
        return 0.0;
    }
    
    virtual double evalImagFreq (double om) const
    {
        return 0.0;
    }
    
    // eps is the y-value of the mother wavelet or its Fourier spectrum without dilation
    virtual double getCutoffTime (double eps) const
    {
        return 0.0;
    }
    
    virtual double getCutoffFreq (double eps) const
    {
        return 0.0;
    }
    
    // the wavelet has the maximum in the time domain by r=0 or t=_pos
    virtual double evalMaxAmplitude (void) const
    {
        return abs(evalCmplTime(_pos));
    }
    
    // the wavelet has the maximum in the frequency domain by the
    // central frequency _f0 or by f=_freq*_f0
    virtual double evalMaxFourier (void) const
    {
        return abs(evalCmplFreq(_freq * _f0));
    }
    
    /**
     *  general properties
     */
    inline WaveletType getType (void)
    {
        return _type;
    }
    
    inline double getFrequency (void) const
    {
        return _freq;
    }
    
    inline double getPosition (void) const
    {
        return _pos;
    }
    
    inline double getCentralFrequency (void) const
    {
        return _f0;
    }
    
    void setFrequency (double aFreq)
    {
        _freq = aFreq / _f0;
    }
    
    void setPosition (double aPos)
    {
        _pos = aPos;
    }
    
    inline bool hasParams ()
    {
        return ((_type != CWdelta) && (_type != CWhaar));
    }
    
    /**
     *  evaluation of the wavelet values in the time and frequency domains
     *  Because we defined the wavelet in virtual methods in its native form,
     *  we have to implement here the dilation and translation rules:
     *  1) wavelet translation is given by: Tg(t) = _freq * g(_freq*(t-_pos))
     *  2) wavelet dilation is given by: Dg(f) = \hat g(2.0*M_PI*f/_freq)
     *  At that the central frequency is already considered in the procedure setFrequency()
     */
    // the input variable t is the time
    PPPcomplex evalCmplTime (double t) const
    {
        double r = _freq * (t - _pos);
        PPPcomplex val(evalRealTime(r), evalImagTime(r));
        return _freq * val;
    }
    
    // the input variable f is the physical frequency (Hz)
    PPPcomplex evalCmplFreq (double f) const
    {
        double om = 2.0 * M_PI * f / _freq;
        return PPPcomplex(evalRealFreq(om), evalImagFreq(om));
    }
    
    unsigned evalTimeRepresentation (PPPVectorContainer<PPPcomplex> &aDest, PPPAxis &aTime,
                                     double aCutOff, bool aReSize = true)
    {
        unsigned nmax = 0; // length of the wavelet
        unsigned slength = (aTime.size() + 1) / 2; // lenght og the source signal
        double db = aTime.getSamplingPeriod();
        // for given frequensy, we calculate the wavelet width
        if (aCutOff == 0.0)
            nmax = slength;
        else
            nmax = (unsigned) (getCutoffTimeNorm(aCutOff) / db) + 2;
        if (nmax > slength) nmax = slength;
        double DeltaOverPeriod = fabs(getDeltaTime() - aTime[slength - 1 + nmax - 2]);
        // next, we fill the array aDest with wavelet values
        if (aReSize)
            aDest.realloc(2 * nmax - 1);
        else if (aDest.size() != (2 * nmax - 1)) PPPBaseObject::onError(
                PPPWAVELETS_ERRSIG + string("evalTimeRepresentation()"));
        aDest.setObjectName(getObjectName());
        double sign = (getFrequency() >= 0.0) ? 1.0 : -1.0;
        for (unsigned n = 0; n < aDest.size(); n++)
            aDest[n] = sign * evalCmplTime(aTime[slength - nmax + n]);
        // If we transform using Haar wavelet, then we must correct the
        // numerical error at the borders of wavelet array aDest because the
        // change of the length of the last and the first wavelet intervals
        // could be less then sampling period. We use the variable DeltaOverPeriod
        if (getType() == PPPWavelet::CWhaar && DeltaOverPeriod < db && aDest.size() > 1)
        {
            aDest[0] = aDest[1] * DeltaOverPeriod / db;
            aDest[aDest.size() - 1] = aDest[aDest.size() - 2] * DeltaOverPeriod / db;
        }
        return nmax;
    }
    
    void evalTimeRepresentation (PPPVectorContainer<PPPcomplex> &aDest, PPPAxis &aTime)
    {
        aDest.realloc(aTime.size());
        aDest.setObjectName(getObjectName());
        for (unsigned i = 0; i < aTime.size(); i++)
            aDest[i] = evalCmplTime(aTime[i]);
    }
    
    void evalTimeRepresentation (PPPSignalContainer<PPPcomplex> &aDest, PPPAxis &aTime)
    {
        aDest.prepare(aTime.size(), 1, aTime, getObjectName());
        evalTimeRepresentation(aDest.getChannel(0), aTime);
    }
    
    void evalFreqRepresentation (PPPVectorContainer<PPPcomplex> &aDest, PPPAxis &aFreq)
    {
        aDest.realloc(aFreq.size());
        aDest.setObjectName(getObjectName());
        for (unsigned i = 0; i < aFreq.size(); i++)
            aDest[i] = evalCmplFreq(aFreq[i]);
    }
    
    void evalFreqRepresentation (PPPSignalContainer<PPPcomplex> &aDest, PPPAxis &aFreq)
    {
        aDest.prepare(aFreq.size(), 1, aFreq, getObjectName());
        evalFreqRepresentation(aDest.getChannel(0), aFreq);
    }
    
    /**
     *  width of the wavelet in the time and frequency domains
     */
    // wavelet width in the time domain
    double getCutoffTimeNorm (double eps) const
    {
        if (eps >= 1.0 || eps <= 0.0) onError(PPPWAVELETS_ERRCUTOFF1 + string("getCutoffTimeNorm"));
        return fabs(getCutoffTime(evalMaxAmplitude() * eps / fabs(_freq)) / _freq);
    }
    
    double getDeltaTime () const
    {
        return getCutoffTimeNorm(PPPWAVELETS_DELTARANGE);
    }
    
    // wavelet width in the frequency domain
    double getCutoffFreqNorm (double eps) const
    {
        if (eps >= 1.0 || eps <= 0.0) onError(PPPWAVELETS_ERRCUTOFF1 + string("getCutoffFreqNorm"));
        return fabs(_freq * getCutoffFreq(evalMaxFourier() * eps) / (2.0 * M_PI));
    }
    
    double getDeltaFreq () const
    {
        return getCutoffFreqNorm(PPPWAVELETS_DELTARANGE);
    }
    
    /**
     *  Chm coefficient used in the inverse wavelet transform
     */
    PPPcomplex Chm (PPPWavelet *inv, double aMax)
    {
        return sqrt(2.0 * M_PI) * Chm_IntegralTrapez(0, aMax, 15, inv);
    }
    
    /**
     *  string interface
     */
    WaveletType strToWavelet (const string &aName)
    {
        if (aName == "morlet") return CWmorlet;
        if (aName == "morletre") return CWmorletre;
        if (aName == "cauchy") return CWcauchy;
        if (aName == "shanon") return CWshanon;
        if (aName == "haar") return CWhaar;
        if (aName == "delta") return CWdelta;
        onError(PPPWAVELETS_ERRWAV + string("strToWavelet"));
        return PPPWavelet::CWnone;
    }
    
    const char *getInfo (void)
    {
        strstream aDest;
        aDest << getType() << ":" << endl;
        aDest << "  time-frequency position (" << getPosition() << "," << getFrequency() << ")"
        << endl;
        aDest << "  time maximum = " << evalMaxAmplitude() << ", time width = " << getDeltaTime()
        << endl;
        aDest << "  frequency maximum = " << evalMaxFourier() << ", frequency width = "
        << getDeltaFreq() << ends;
        onNotation(aDest.str());
        return getNotation();
    }
    
    friend ostream& operator << (ostream& aDest, WaveletType aSour)
    {
        switch (aSour)
        {
        case CWnone:
            aDest << PPPWAVELETS_NONE;
            break;
        case CWmorlet:
            aDest << PPPWAVELETS_MORLET;
            break;
        case CWmorletre:
            aDest << PPPWAVELETS_MORLETRE;
            break;
        case CWcauchy:
            aDest << PPPWAVELETS_CAUCHY;
            break;
        case CWshanon:
            aDest << PPPWAVELETS_SHANON;
            break;
        case CWhaar:
            aDest << PPPWAVELETS_HAAR;
            break;
        case CWdelta:
            aDest << PPPWAVELETS_DELTA;
            break;
        }
        return aDest;
    }
    
private:
    
    PPPcomplex Chm_func (double f, PPPWavelet *inv)
    {
        PPPcomplex z1 = conj(evalCmplFreq(f));
        PPPcomplex z2 = inv->evalCmplFreq(f);
        PPPcomplex z3 = conj(evalCmplFreq(-f));
        PPPcomplex z4 = inv->evalCmplFreq(-f);
        return (z1 * z2 + z3 * z4) / f;
    }
    
    /** Trapezoidal rule ***********************************************/
    /** aSize is count of point = power of two, f.e 128 **/
    PPPcomplex Chm_IntegralTrapez (double a, double b, unsigned aSize, PPPWavelet *inv)
    {
        double x, tnm, del;
        PPPcomplex s, sum;
        unsigned it, j;
        if (aSize == 1)
        {
            return (s = 0.5 * (b - a) * (Chm_func(a, inv) + Chm_func(b, inv)));
        }
        else
        {
            for (it = 1, j = 1; j < aSize - 1; j++)
                it <<= 1;
            tnm = it;
            del = (b - a) / tnm;
            x = a + 0.5 * del;
            for (sum = 0.0, j = 1; j <= it; j++, x += del)
                sum += Chm_func(x, inv);
            s = ((b - a) * sum / tnm);
            return s;
        }
    }
    
};
// end of object

#undef PPPWAVELETS_DELTARANGE
#endif

