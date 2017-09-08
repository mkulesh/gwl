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
#ifndef _PPPAPPROXTWOGAUSS
#define _PPPAPPROXTWOGAUSS

/************************************************************************
 * PPPApproximateTwoGauss
 ***********************************************************************/
 class PPPApproximateTwoGauss : public PPPApproximate
  {
  private:
    PPPcomplex j;
    double ATN_HUGE_VALUE;

  public:

    PPPApproximateTwoGauss(void) {
      setObjectName("Two gauss function");
      _type = APTtwogauss;
      j = PPPcomplex(0.0, 1.0);
      ATN_HUGE_VALUE = 1.0E10;
      };

    inline PPPcomplex CmplPhase(double f) {
      double aDist = _params[6];
      if(aDist == 0.0)
        onError("Distance is zero in procedure: CmplPhase");
      double V1 = _params[0] + _params[1]*exp(-f*f/(_params[2]*_params[2]));
      double V2 = _params[3] + _params[4]*exp(-f*f/(_params[5]*_params[5]));
      double Phic = M_PI*aDist*(f/V1 + f/V2);
      double acos = 2.0*cos(M_PI*aDist*(f/V1 - f/V2));
      PPPcomplex Atnc = 0.0;
      if(acos == 0.0) Atnc  = -ATN_HUGE_VALUE;
      else if(acos<0.0) Atnc = log(-acos) + M_PI*j;
      else Atnc = log(acos);
      PPPcomplex aPhi = (Phic + j*Atnc)/aDist;
      return aPhi;
      };

    inline PPPcomplex CmplPhasePrim(double f) {
      double aDist = _params[6];
      if(aDist == 0.0)
        onError("Distance is zero in procedure: CmplPhasePrim");
      double V1 = _params[0] + _params[1]*exp(-f*f/(_params[2]*_params[2]));
      double V2 = _params[3] + _params[4]*exp(-f*f/(_params[5]*_params[5]));
      double V1p = -2.0*_params[1]*f*exp(-f*f/(_params[2]*_params[2]))/(_params[2]*_params[2]);
      double V2p = -2.0*_params[4]*f*exp(-f*f/(_params[5]*_params[5]))/(_params[5]*_params[5]);
      double dV1 = (1.0/V1 + 1.0/V2 - f*V1p/(V1*V1) - f*V2p/(V2*V2));
      double dV2 = (1.0/V1 - 1.0/V2 - f*V1p/(V1*V1) + f*V2p/(V2*V2));
      return M_PI*(dV1 - j*dV2*tan(M_PI*aDist*(f/V1 - f/V2)));
      };

  }; // end of object


class PPPApproximateTwoGaussPhi : public PPPApproximateTwoGauss
  {
  public:

    PPPApproximateTwoGaussPhi(unsigned aSize) {
      if(aSize != 7) onError(ARG_VALUE+string("PPPApproximateTwoGaussPhi"));
      _params.resize(aSize);
      };

    const char *getInfo(void) {
      strstream str;
      str << "Func(f) = real((f/V1 + f/V2)/2.0 - j*log(2.0*cos((f/V1 - f/V2)/2.0)))" << endl;
      str << "  a[i] = " << _params.vectorToStr() << ends;
      onNotation(str.str());
      return getNotation();
      };

    double Func(double f) {
      return real(CmplPhase(f))/(2.0*M_PI);
      };

    double FuncDivX(double f) {
      return real(CmplPhasePrim(f))/(2.0*M_PI);
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

class PPPApproximateTwoGaussAtn : public PPPApproximateTwoGauss
  {
  public:

    PPPApproximateTwoGaussAtn(unsigned aSize) {
      if(aSize != 7) onError(ARG_VALUE+string("PPPApproximateTwoGaussAtn"));
      _params.resize(aSize);
      };

    const char *getInfo(void) {
      strstream str;
      str << "Func(f) = -imag((f/V1 + f/V2)/2.0 - j*log(2.0*cos((f/V1 - f/V2)/2.0)))" << endl;
      str << "  a[i] = " << _params.vectorToStr() << ends;
      onNotation(str.str());
      return getNotation();
      };

    double Func(double f) {
      return -1.0*imag(CmplPhase(f));
      };

    double FuncDivX(double f) {
      return -1.0*imag(CmplPhasePrim(f));
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
