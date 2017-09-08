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
#ifndef _PPPWAVELETDELTA
#define _PPPWAVELETDELTA

/************************************************************************
 * PPPWaveletDelta
 ***********************************************************************/
class PPPWaveletDelta : public PPPWavelet
{
public:
    
    PPPWaveletDelta (void)
    {
        _type = CWdelta;
        setObjectName (PPPWAVELETS_DELTA);
    }
    
private:
    
    double evalRealTime (double r) const
    {
        return (r == 0.0) ? 1.0 / _freq : 0.0;
    }
    
    double evalImagTime (double r) const
    {
        return 0.0;
    }
    
    double evalRealFreq (double om) const
    {
        return 1.0 / sqrt(2.0 * M_PI);
    }
    
    double evalImagFreq (double om) const
    {
        return 0.0;
    }
    
};
// end of object

#endif

