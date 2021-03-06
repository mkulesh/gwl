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
#ifndef _PPPWAVELETMORLET
#define _PPPWAVELETMORLET

/************************************************************************
 * PPPWaveletMorlet
 ***********************************************************************/
class PPPWaveletMorlet : public PPPWavelet
{
private:
    
    double _sigma;  // wavelet parameter
    
public:
    
    PPPWaveletMorlet (double sigma = 1.0, double freq = 1.0, double pos = 0.0)
    {
        _type = CWmorlet;
        setObjectName (PPPWAVELETS_MORLET);
        setFrequency(freq);
        setPosition(pos);
        _sigma = sigma;
        _f0 = 1.0;
    }
    
    void setFrequencyResolution (double dfreq)
    {
        _sigma = 1.0 / dfreq;
    }
    
private:
    
    double evalRealTime (double r) const
    {
        return cos(2.0 * M_PI * r) * exp(-(r * r) / (2.0 * _sigma * _sigma));
    }
    
    double evalImagTime (double r) const
    {
        return sin(2.0 * M_PI * r) * exp(-(r * r) / (2.0 * _sigma * _sigma));
    }
    
    double evalRealFreq (double om) const
    {
        double a = _sigma * (om - 2.0 * M_PI);
        return (om <= 0.0) ? 0.0 : _sigma * sqrt(2.0 * M_PI) * exp(-a * a / 2.0);
    }
    
    double evalImagFreq (double om) const
    {
        return 0.0;
    }
    
    double getCutoffTime (double eps) const
    {
        return _sigma * sqrt(2.0 * log(1.0 / eps));
    }
    
    double getCutoffFreq (double eps) const
    {
        return sqrt(-2.0 * log(eps / (sqrt(2.0 * M_PI) * _sigma))) / _sigma;
    }
    
};
// end of object

#endif

