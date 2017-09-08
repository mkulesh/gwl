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

#ifndef _PPPWAVELETSHANON
#define _PPPWAVELETSHANON

/************************************************************************
 * PPPWaveletShanon
 ***********************************************************************/
class PPPWaveletShanon : public PPPWavelet
  {
  private:

    double _sigma;  // wavelet parameter

  public:

    PPPWaveletShanon(double sigma=1.0, double freq=1.0, double pos=0.0) {
      _type = CWshanon;
      setObjectName(PPPWAVELETS_SHANON);
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
      if(r == 0.0) return _sigma;
      return cos(2.0*M_PI*r)*sin(M_PI*_sigma*r)/(M_PI*r);
      };

    double evalImagTime(double r) const {
      if(r == 0.0) return 0.0;
      return sin(2.0*M_PI*r)*sin(M_PI*_sigma*r)/(M_PI*r);
      };

    double evalRealFreq(double om) const {
      return (H(om-2.0*M_PI+_sigma*M_PI) - H(om-2.0*M_PI-_sigma*M_PI));
      };

    double evalImagFreq(double om) const {
      return 0.0;
      };

    double getCutoffTime(double eps) const {
      return 1.0/(M_PI*eps);
      };

    double getCutoffFreq(double eps) const {
      return _sigma*M_PI;
      };

    inline double H(double w) const { return ((w>0.0) ? 1.0 : 0.0); }; /* Heaviside step function */

  }; // end of object


#endif

