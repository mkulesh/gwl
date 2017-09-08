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

#ifndef _PPPWAVELETMORLETRE
#define _PPPWAVELETMORLETRE

/************************************************************************
 * PPPWaveletMorletRe
 ***********************************************************************/
class PPPWaveletMorletRe : public PPPWavelet
  {
  private:

    double _sigma;  // wavelet parameter

  public:

    PPPWaveletMorletRe(double sigma=1.0, double freq=1.0, double pos=0.0) {
      _type=CWmorletre;
      setObjectName(PPPWAVELETS_MORLETRE);
      setFrequency(freq);
      setPosition(pos);
      _sigma = sigma;
      _f0 = 1.0;
      };

    void setFrequencyResolution(double dfreq) {
      _sigma = 1.0/dfreq;
      };

  private:

    double evalRealTime(double r) const {
      return cos(2.0*M_PI*r)*exp(-(r*r)/(2.0*_sigma*_sigma));
      };

    double evalImagTime(double r) const {
      return 0;
      };

    double evalRealFreq(double om) const {
      double a1 = _sigma*(om - 2.0*M_PI);
      double a2 = _sigma*(-om - 2.0*M_PI);
      return (om==0.0)? 0.0 : _sigma*sqrt(2.0*M_PI)*(exp(-a1*a1/2.0)+exp(-a2*a2/2.0))/2.0;
      };

    double evalImagFreq(double om) const {
      return 0;
      };

    double getCutoffTime(double eps) const {
      return _sigma*sqrt(2.0*log(1.0/eps));
      };

    double getCutoffFreq(double eps) const {
      return sqrt(-2.0*log(2.0*eps/(sqrt(2.0*M_PI)*_sigma)))/_sigma;
      };

  }; // end of object



#endif

