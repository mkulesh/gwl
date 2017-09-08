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
#ifndef _PPPTRANSFK
#define _PPPTRANSFK

#define PPPTRANSFK_NAME    "Frequency-wavenumber transformation"
#define PPPTRANSFK_ERRINT  "type of interpolation is invalid in procedure: "
#define PPPTRANSFK_ERRCOR  "type of component is invalid in procedure: "

/************************************************************************
 * PPPTransFK
 * @article{Kulesh2007PaAG,
 *   author = {M. Kulesh and M. Holschneider and M. Ohrnberger and E. L\"ueck},
 *   title = {Modeling of wave dispersion using continuous wavelet transforms {II}: wavelet based frequency-velocity analysis},
 *   journal = {Pure and Applied Geophysics},
 *   year = {2007},
 *   volume = {162},
 *   pages = {in press},
 * }
 ************************************************************************/
class PPPTransFK : public PPPBaseObject
{
public:
    
    typedef enum
    {
        MMIindex, MMIspline
    } InterpType;

    typedef enum
    {
        MMCabs, MMCarg, MMCcphase
    } CorrType;

private:
    
    InterpType _inttype;
    CorrType _corrtype;
    double _filter;
    bool _norm;
    double _distance;

public:
    
    PPPTransFK (void) :
            _inttype(MMIindex),
            _corrtype(MMCarg),
            _filter(0.0),
            _norm(true),
            _distance(0.0)
    {
        setObjectName(PPPTRANSFK_NAME);
    }
    
    friend ostream& operator << (ostream& aDest, InterpType aSour)
    {
        switch (aSour)
        {
        case MMIindex:
            aDest << "indexed shift";
            break;
        case MMIspline:
            aDest << "spline interpolation";
            break;
        }
        return aDest;
    }
    
    friend ostream& operator << (ostream& aDest, CorrType aSour)
    {
        switch (aSour)
        {
        case MMCabs:
            aDest << "correlation of modulus";
            break;
        case MMCarg:
            aDest << "correlation of phase";
            break;
        case MMCcphase:
            aDest << "correlation of complex phase";
            break;
        }
        return aDest;
    }
    
    void setIntType (InterpType const aType)
    {
        _inttype = aType;
    }
    
    void setCorrType (CorrType const aType)
    {
        _corrtype = aType;
    }
    
    void setFilter (double const aFilter)
    {
        _filter = aFilter;
    }
    
    void setNormVoice (bool const aType)
    {
        _norm = aType;
    }
    
    void setDistance (double adX)
    {
        _distance = adX;
    }
    
    void evalSpectrCorrelation (PPPSpectrContainer<double> &aFK,
                                PPPSpectrContainer<PPPcomplex> &aWg, PPPAxis &aVel, int aTestVoice =
                                        -1)
    {
        PPPSignalContainer<PPPcomplex> aCmplSource, aShiftSource;
        PPPVectorContainer<double> aMaxMod;
        evalMaximumModulus(aMaxMod, aWg);
        aFK.prepare(aWg.voices(), aVel.size(), 1, aVel, aWg.getFreq(), "FK image");
        strstream aDest;
        aDest << "calculation of spectrum correlation:" << endl << "  " << aFK.getTime().getInfo()
        << endl << "  " << aFK.getFreq().getInfo() << endl << "  type of interpolation: "
        << _inttype << endl << "  type of correlation: " << _corrtype << endl
        << "  value of energy filter = " << _filter << endl;
        if (_norm)
            aDest << "  normalization of voices: on" << ends;
        else
            aDest << "  normalization of voices: off" << ends;
        onMessage(aDest.str());
        char currFreq[100];
        double shift, amax;
        unsigned in_start = (aTestVoice >= 0) ? aTestVoice : 0;
        unsigned in_end = (aTestVoice >= 0) ? aTestVoice + 1 : aWg.voices();
        for (unsigned i = in_start; i < in_end; i++)
        {
            sprintf(currFreq, "%s: %3d", PPPTRANSWAVELET_CURR, i + 1);
            onProgress(100 * i / (aWg.voices() - 1), currFreq);
            aWg.getVoices(aCmplSource, i);
            amax = 0.0;
            if (_inttype == MMIindex)
                for (unsigned n = 0; n < aFK.points(); n++)
                {
                    unsigned ind = (unsigned) aFK.getTime(aFK.points() - 1 - n);
                    aFK(i, n) = evalVoicesCorrelation(aCmplSource, ind, aMaxMod);
                    if (aFK(i, n) > amax) amax = aFK(i, n);
                }
            else if (_inttype == MMIspline)
                for (unsigned n = 0; n < aFK.points(); n++)
                {
                    shift = _distance / (aFK.getTime(n) * aWg.getTime().getSamplingPeriod());
                    shiftVoices(aShiftSource, aCmplSource, MMIspline, shift);
                    aFK(i, n) = evalVoicesCorrelation(aShiftSource, 0, aMaxMod);
                    if (aFK(i, n) > amax) amax = aFK(i, n);
                }
            else
                onError(PPPTRANSFK_ERRINT + string("evalSpectrCorrelation()"));
            if (_norm && amax != 0.0) for (unsigned n = 0; n < aFK.points(); n++)
                aFK(i, n) = aFK(i, n) / amax;
        }
    }
    
private:
    
    void evalMaximumModulus (PPPVectorContainer<double> &aMod, PPPSpectrContainer<PPPcomplex> &aWg)
    {
        aMod.realloc(aWg.channels());
        PPPcomplexAbs aWgAbs;
        for (unsigned k = 0; k < aWg.channels(); k++)
            aMod[k] = aWgAbs(aWg.getChannel(k).getMaxValue(aWgAbs));
        return;
    }
    
    void shiftVoices (PPPSignalContainer<PPPcomplex> &aDest,
                      PPPSignalContainer<PPPcomplex> &aSource, InterpType aType, double aInd)
    {
        aDest.prepare(aSource);
        register unsigned j, k;
        if (aType == MMIindex)
        {
            for (k = 0; k < aSource.channels(); k++)
                for (j = 0; j < aSource.points(); j++)
                {
                    unsigned ni = j + (unsigned) aInd * k;
                    if (ni < aSource.points())
                        aDest(j, k) = aSource(ni, k);
                    else
                        aDest(j, k) = 0.0;
                }
        }
        else if (aType == MMIspline)
        {
            for (k = 0; k < aSource.channels(); k++)
            {
                aSource.prepareSpline(k);
                for (j = 0; j < aSource.points(); j++)
                {
                    double nx = aSource.getAxis().getSamplingPeriod()
                            * ((double) j + aInd * (double) k);
                    aDest(j, k) = aSource.evalSplineInterpolation(nx, k);
                }
            }
        }
        else
            onError(PPPTRANSFK_ERRINT + string("shiftVoices()"));
    }
    
    double evalVoicesCorrelation (PPPSignalContainer<PPPcomplex> &aSource, unsigned aInd,
                                  PPPVectorContainer<double> &aMod)
    {
        register unsigned j, k, m, ni;
        double val, summ;
        PPPcomplex cmpl, prod, _i(0.0, 1.0);
        summ = 0.0;
        switch (_corrtype)
        {
        case MMCabs:
        case MMCarg:
            for (j = 0; j < aSource.points(); j++)
            {
                prod = 1.0;
                for (k = 0; k < aSource.channels(); k++)
                {
                    ni = j + aInd * k;
                    if (ni < aSource.points())
                    {
                        cmpl = aSource(ni, k);
                    }
                    else
                    {
                        prod = 0.0;
                        break;
                    }
                    if (_corrtype == MMCabs)
                        val = abs(aSource(ni, k)) / aMod[k];
                    else if (_corrtype == MMCarg
                            && 100.0 * abs(aSource(ni, k)) / aMod[k] >= _filter)
                        val = arg(aSource(ni, k));
                    else
                    {
                        prod = 0.0;
                        break;
                    }
                    prod *= val;
                }
                summ += abs(prod);
            }
            break;
        case MMCcphase:
            for (j = 0; j < aSource.points(); j++)
            {
                prod = 0.0;
                for (k = 0; k < aSource.channels(); k++)
                    for (m = 0; m < aSource.channels(); m++)
                    {
                        PPPcomplex val1 = 0.0, val2 = 0.0;
                        cmpl = aSource(j, k);
                        if (100.0 * abs(cmpl) / aMod[k] >= _filter) val1 = exp(_i * arg(cmpl));
                        cmpl = aSource(j, m);
                        if (100.0 * abs(cmpl) / aMod[m] >= _filter) val2 = exp(_i * arg(cmpl));
                        prod = prod + val1 * conj(val2);
                    }
                summ += abs(prod);
            }
            break;
        default:
            onError(PPPTRANSFK_ERRCOR + string("evalVoicesCorrelation()"));
        }
        return summ;
    }
    
};
// end of object

#endif
