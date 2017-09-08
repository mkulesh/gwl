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

#ifndef _PPPWAVELETDELTA
#define _PPPWAVELETDELTA

/************************************************************************
 * PPPWaveletDelta
 ***********************************************************************/
class PPPWaveletDelta : public PPPWavelet
  {
  public:

    PPPWaveletDelta(void) {
      _type=CWdelta;
      setObjectName(PPPWAVELETS_DELTA);
      };

  private:

    double evalRealTime(double r) const {
      return (r == 0.0) ? 1.0/_freq : 0.0;
      };

    double evalImagTime(double r) const {
      return 0.0;
      };

    double evalRealFreq(double om) const {
      return 1.0/sqrt(2.0*M_PI);
      };

    double evalImagFreq(double om) const {
      return 0.0;
      };

  }; // end of object



#endif

