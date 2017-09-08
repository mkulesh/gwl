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

#ifndef _PPPAPPROXGAUSS
#define _PPPAPPROXGAUSS

/************************************************************************
 * PPPApproximateGauss
 ***********************************************************************/
class PPPApproximateGauss : public PPPApproximate
  {
  public:

    PPPApproximateGauss(unsigned aSize) {
      setObjectName(PPPAPPROXIMATE_GAU);
      _type = APTgauss;
      if(aSize != 3) onError(ARG_VALUE+string("PPPApproximateGauss"));
      _params.resize(aSize);
      };

    const char *getInfo(void) {
      strstream str;
      str << "Func(f) = a0*f + a1*f*exp(-f*f/(2.0*a2*a2))" << endl;
      str << "  a[i] = " << _params.vectorToStr() << ends;
      onNotation(str.str());
      return getNotation();
      };

    double Func(double f) {
      return _params[0]*f + _params[1]*f*exp(-f*f/(2.0*_params[2]*_params[2]));
      };

    double FuncDivX(double f) {
      return _params[0] + _params[1]*(1.0-f*f/(_params[2]*_params[2]))*exp(-f*f/(2.0*_params[2]*_params[2]));
      };

    double FuncDivXPar(double f, unsigned Ai) {
      double res = 0.0;
      double ff  = f*f/(_params[2]*_params[2]);
      if(Ai == 0) res=1.0;
      if(Ai == 1) res=(1.0-ff)*exp(-ff/2.0);
      if(Ai == 2) res=_params[1]*ff*(3.0-ff)*exp(-ff/2.0)/_params[2];
      return res;
      };

    double FuncDivPar(double f, unsigned Ai) {
      double res = 0.0;
      double ff  = f*f/(_params[2]*_params[2]);
      if(Ai == 0) res=f;
      if(Ai == 1) res=f*exp(-ff/2.0);
      if(Ai == 2) res=_params[1]*f*f*f*exp(-ff/2.0)/(2*_params[2]*_params[2]);
      return res;
      };

  }; // end of object


#endif
