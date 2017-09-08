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

#ifndef _PPPAPPROXCOLECOLE
#define _PPPAPPROXCOLECOLE

/************************************************************************
 * PPPApproximateColeCole
 * The Cole-Cole model
 * Y.Wang, J.Guo, Modifed Kolsky model for seismic attenuation and dispersion.
 * J.Geophys.E., No. 1 (2004), 187-196
 ***********************************************************************/
 class PPPApproximateColeCole : public PPPApproximate
  {
  private:
    PPPcomplex j;

  public:

    PPPApproximateColeCole(void) {
      setObjectName(PPPAPPROXIMATE_COL);
      _type = APTcolecole;
      j = PPPcomplex(0.0, 1.0);
      };

    inline PPPcomplex CmplPhase(double om) {
      PPPcomplex M = _params[0]*(1.0+pow(j*om*_params[2], _params[1]))/(1.0+pow(j*om*_params[3], _params[1]));
      return om/sqrt(M);
      };

    inline PPPcomplex CmplPhasePrim(double om) {
      PPPcomplex M = _params[0]*(1.0+pow(j*om*_params[2], _params[1]))/(1.0+pow(j*om*_params[3], _params[1]));
      PPPcomplex Mprim = (_params[0]*_params[1]/om)*
         (pow(j*om*_params[2], _params[1])-pow(j*om*_params[3], _params[1]))/
         ((1.0+pow(j*om*_params[3], _params[1]))*(1.0+pow(j*om*_params[3], _params[1])));
      return (1.0/sqrt(M))*(1.0-om*Mprim/(2.0*M));
      };

  }; // end of object


class PPPApproximateColeColePhi : public PPPApproximateColeCole
  {
  public:

    PPPApproximateColeColePhi(unsigned aSize) {
      if(aSize != 4) onError(ARG_VALUE+string("PPPApproximateColeColePhi"));
      _params.resize(aSize);
      };

    const char *getInfo(void) {
      strstream str;
      str << "Func(f) = real(om/sqrt[a0*(1+(j*om*a2)^a1)/(1+(j*om*a3)^a1)])" << endl;
      str << "  a[i] = " << _params.vectorToStr() << ends;
      onNotation(str.str());
      return getNotation();
      };

    double Func(double f) {
      return real(CmplPhase(2.0*M_PI*f))/(2.0*M_PI);
      };

    double FuncDivX(double f) {
      return real(CmplPhasePrim(2.0*M_PI*f));
      };

    double FuncDivXPar(double f, unsigned Ai) {
      onError(VIRT_NOTDEF+string("FuncDivXPar"));
      return 0.0;
      };

    double FuncDivPar(double f, unsigned Ai) {
      onError(VIRT_NOTDEF+string("FuncDivPar"));
      return 0.0;
      };

  }; // end of object

class PPPApproximateColeColeAtn : public PPPApproximateColeCole
  {
  public:

    PPPApproximateColeColeAtn(unsigned aSize) {
      if(aSize != 4) onError(ARG_VALUE+string("PPPApproximateColeColeAtn"));
      _params.resize(aSize);
      };

    const char *getInfo(void) {
      strstream str;
      str << "Func(f) = imag(om/sqrt[a0*(1+(j*om*a2)^a1)/(1+(j*om*a3)^a1)])" << endl;
      str << "  a[i] = " << _params.vectorToStr() << ends;
      onNotation(str.str());
      return getNotation();
      };

    double Func(double f) {
      return -1.0*imag(CmplPhase(2.0*M_PI*f));
      };

    double FuncDivX(double f) {
      return -2.0*M_PI*imag(CmplPhasePrim(2.0*M_PI*f));
      };

    double FuncDivXPar(double f, unsigned Ai) {
      onError(VIRT_NOTDEF+string("FuncDivXPar"));
      return 0.0;
      };

    double FuncDivPar(double f, unsigned Ai) {
      onError(VIRT_NOTDEF+string("FuncDivPar"));
      return 0.0;
      };

  }; // end of object

#endif
