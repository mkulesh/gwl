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
#ifndef _PPPWAVELETHAAR
#define _PPPWAVELETHAAR

/************************************************************************
 * PPPWaveletHaar
 ***********************************************************************/
class PPPWaveletHaar : public PPPWavelet
{
public:
    
    PPPWaveletHaar (double freq = 1.0, double pos = 0.0)
    {
        _type = CWhaar;
        setObjectName (PPPWAVELETS_HAAR);
        setFrequency(freq);
        setPosition(pos);
        _f0 = 0.75971;
    }
    
private:
    
    double evalRealTime (double r) const
    {
        if (r < -0.5 || r > 0.5 || r == 0.0) return 0.0;
        if (r < 0.0) return -1.0;
        return 1.0;
    }
    
    double evalImagTime (double r) const
    {
        return 0.0;
    }
    
    double evalRealFreq (double om) const
    {
        return 0.0;
    }
    
    double evalImagFreq (double om) const
    {
        return (om == 0.0) ? 0.0 : (2.0 / om) * (cos(om / 2.0) - 1.0);
    }
    
    double getCutoffTime (double eps) const
    {
        return 0.5;
    }
    
    double getCutoffFreq (double eps) const
    {
        return 4.0 / eps;
    }
    
    double evalMaxAmplitude (void) const
    {
        return fabs(_freq);
    }
    
};
// end of object

#endif

