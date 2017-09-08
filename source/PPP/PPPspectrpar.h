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
#ifndef _PPPSPECTRPARAMS
#define _PPPSPECTRPARAMS

#define PPPSPECTRPARAMS_OBJVER    "WP1.4"
#define PPPSPECTRPARAMS_NAME      "wavelet transform parameters"
#define PPPSPECTRPARAMS_FREQ      "frequency"
#define PPPSPECTRPARAMS_ERRWAV    "wavelet undefined in procedure: "

/************************************************************************/
/** PPPSpectrParams                                                      */
/************************************************************************/
class PPPSpectrParams : public PPPBaseObject
{
private:

    unsigned _transformType;                    // slow or fast wavelet transform
    double _cutoffprec;                         // wavelet cutoff for slow wavelet transform
    // characterizations of frequency axis
    PPPAxis _freq;                              // the value of the frequencies
    PPPAxis _time;                              // the time axis of transform
    // wavelet of direct transform
    PPPWavelet::WaveletType _waveletType;       // the type of the wavelet of direct transform
    double _waveletPar;                         // a single parameter for the wavelet
    PPPWavelet * _wavelet;                      // the wavelt itself
    // wavelet of inverse transform
    PPPWavelet::WaveletType _invwaveletType;    // the type of wavelet for the inverse transform
    double _invwaveletPar;                      // a single paramter for the inverse transform
    PPPWavelet * _invwavelet;                   // the reconstruction wavelet
    
public:
    
    typedef enum
    {
        WSTdirect,        // direct transform
        WSTinverse        // inverse transform
    } TransformDir;

    PPPSpectrParams (void) :
            _wavelet(NULL),
            _invwavelet(NULL),
            _cutoffprec(0.01)
    {
        setObjectVer(PPPSPECTRPARAMS_OBJVER);
        setObjectName(PPPSPECTRPARAMS_NAME);
        PPPAxis aAxis(128, 1, 50.0, PPPAxis::ATlin, PPPSPECTRPARAMS_FREQ);
        initialize(0, aAxis, "morlet", 1.0, "delta", 1.0);
    }
    
    ~PPPSpectrParams ()
    {
        if (_wavelet != NULL) delete _wavelet;
        if (_invwavelet != NULL) delete _invwavelet;
    }
    
    void initialize (unsigned atransformType, PPPAxis &aAxis, string awavName, double awaveletPar,
                     string ainvwavName, double ainvwaveletPar, double acutoffprec = 0.01)
    {
        PPPWavelet tmpwav;
        _transformType = atransformType;
        _waveletType = tmpwav.strToWavelet(awavName);
        _waveletPar = awaveletPar;
        _invwaveletType = tmpwav.strToWavelet(ainvwavName);
        _invwaveletPar = ainvwaveletPar;
        _cutoffprec = acutoffprec;
        _prepareWavelets();
        _freq.assign(aAxis);
    }
    
    void setInverseParams (string ainvwavName, double ainvwaveletPar, double acutoffprec = 0.01)
    {
        PPPWavelet tmpwav;
        _invwaveletType = tmpwav.strToWavelet(ainvwavName);
        _invwaveletPar = ainvwaveletPar;
        _cutoffprec = acutoffprec;
        _prepareWavelets();
    }
    
    inline unsigned getTransformType () const
    {
        return _transformType;
    }
    
    inline unsigned voices () const
    {
        return _freq.size();
    }
    
    inline PPPAxis& getFreq (void)
    {
        return _freq;
    }
    
    inline PPPAxis& getTime (void)
    {
        return _time;
    }
    
    inline int getWriteFreq () const
    {
        return -1;
    }
    
    inline double getCutoffPrec () const
    {
        return _cutoffprec;
    }
    
    inline PPPWavelet::WaveletType getWaveletType (const TransformDir aType) const
    {
        return (aType == WSTdirect) ? _waveletType : _invwaveletType;
    }
    
    inline PPPWavelet * getWavelet (const TransformDir aType)
    {
        return (aType == WSTdirect) ? _wavelet : _invwavelet;
    }
    
    inline double getWaveletPar (const TransformDir aType)
    {
        return (aType == WSTdirect) ? _waveletPar : _invwaveletPar;
    }
    
    PPPcomplex getInverseConstant ()
    {
        _clearWaveletsPosition();
        return getWavelet(WSTdirect)->Chm(getWavelet(WSTinverse), 100);
    }
    
    /**
     *  file stream operations
     */
    void fwrite (FILE *stream)
    {
        unsigned val;
        fwrite_streaminfo(stream, getObjectVer(), sizeof(double));
        std::fwrite((void*) &_transformType, sizeof(_transformType), 1, stream);
        val = _waveletType;
        std::fwrite((void*) &val, sizeof(val), 1, stream);
        std::fwrite((void*) &_waveletPar, sizeof(_waveletPar), 1, stream);
        val = _invwaveletType;
        std::fwrite((void*) &val, sizeof(val), 1, stream);
        std::fwrite((void*) &_invwaveletPar, sizeof(_invwaveletPar), 1, stream);
        _freq.fwrite(stream);
        PPPBaseObject::fwrite(stream);
    }
    
    void fread (FILE *stream)
    {
        unsigned val;
        fread_streaminfo(stream, getObjectVer(), sizeof(double));
        std::fread((void*) &_transformType, sizeof(_transformType), 1, stream);
        std::fread((void*) &val, sizeof(val), 1, stream);
        _waveletType = (PPPWavelet::WaveletType) val;
        std::fread((void*) &_waveletPar, sizeof(_waveletPar), 1, stream);
        std::fread((void*) &val, sizeof(val), 1, stream);
        _invwaveletType = (PPPWavelet::WaveletType) val;
        std::fread((void*) &_invwaveletPar, sizeof(_invwaveletPar), 1, stream);
        _freq.fread(stream);
        PPPBaseObject::fread(stream);
        _prepareWavelets();
    }
    
    PPPWavelet *createWavelet (PPPWavelet::WaveletType aWavelet, double aPar)
    {
        PPPWavelet *wav;
        switch (aWavelet)
        {
        case PPPWavelet::CWmorlet:
            wav = new PPPWaveletMorlet(aPar, 1.0, 0.0);
            break;
        case PPPWavelet::CWmorletre:
            wav = new PPPWaveletMorletRe(aPar, 1.0, 0.0);
            break;
        case PPPWavelet::CWcauchy:
            wav = new PPPWaveletCauchy(aPar, 1.0, 0.0);
            break;
        case PPPWavelet::CWshanon:
            wav = new PPPWaveletShanon(aPar, 1.0, 0.0);
            break;
        case PPPWavelet::CWhaar:
            wav = new PPPWaveletHaar(1.0, 0.0);
            break;
        case PPPWavelet::CWdelta:
            wav = new PPPWaveletDelta;
            break;
        default:
            onError(PPPSPECTRPARAMS_ERRWAV + string("createWavelet"));
        }
        if (wav == NULL) onError(MEM_ERRALLOC + string("createWavelet"));
        return wav;
    }
    
private:
    
    void _prepareWavelets (void)
    {
        if (_wavelet != NULL) delete _wavelet;
        _wavelet = createWavelet(_waveletType, _waveletPar);
        if (_invwavelet != NULL) delete _invwavelet;
        _invwavelet = createWavelet(_invwaveletType, _invwaveletPar);
    }
    
    void _clearWaveletsPosition ()
    {
        getWavelet(WSTdirect)->setPosition(0.0);
        getWavelet(WSTdirect)->setFrequency(1.0);
        getWavelet(WSTinverse)->setPosition(0.0);
        getWavelet(WSTinverse)->setFrequency(1.0);
    }
    
};
// end of object

#endif
