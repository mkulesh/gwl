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

#ifndef _PPPWAVELETCAUCHY
#define _PPPWAVELETCAUCHY

/************************************************************************
 * PPPWaveletCauchy
 ***********************************************************************/
class PPPWaveletCauchy: public PPPWavelet
  {
  private:

    double _p;          // wavelet parameter
    double _fourscale;  // coefficient of Fourier transform

  public:

    PPPWaveletCauchy(double p=1.0, double freq=1.0, double pos=0.0) {
      _type = CWcauchy;
      setObjectName(PPPWAVELETS_CAUCHY);
      setFrequency(freq);
      setPosition(pos);
      _p = p;
      _fourscale = getFourScale();
      _f0 = 1.0;
      };

    void setFrequencyResolution (double dfreq) {
      _p = (int)(1.0/(dfreq)) + 1;
      _fourscale = getFourScale();
      };

  private:

    double evalRealTime(double r) const{
      double tt = 2.0*M_PI*r/(_p-1.0);
      return pow(1.0 + tt*tt, -_p/2.0)*cos(_p*atan(tt));
      };

    double evalImagTime(double r) const {
      double tt = 2.0*M_PI*r/(_p-1.0);
      return pow(1.0 + tt*tt, -_p/2.0)*sin(_p*atan(tt));
      };

    double evalRealFreq(double om) const {
      double f1 = om/(2.0*M_PI);
      return (om<=0.0)? 0.0 : _fourscale*pow(f1,_p-1.0)*exp(-f1*(_p-1.0));
      };

    double evalImagFreq(double om) const{
      return 0.0;
      };

    double getFourScale() {
      PPPMathFunc MFunc;
      return pow(_p-1.0,_p)/MFunc.fact((int)_p-1);
      };

    double getCutoffTime(double eps) const {
      return (_p-1.0)*sqrt(pow(eps, -2.0/_p)-1.0)/(2.0*M_PI);
      };

    double getCutoffFreq(double eps) const {
      return 4.0*M_PI*sqrt(-1.0-log(eps/fabs(_fourscale))/(_p-1.0));
      };

  }; // end of object



#endif

