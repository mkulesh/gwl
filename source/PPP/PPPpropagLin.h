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

#ifndef _PPPPROPAGATORLIN
#define _PPPPROPAGATORLIN

#define PPPPROPAGATORLIN_NAME    "linear wavelet diffeomorphism"
#define PPPPROPAGATORLIN_WAV     "calculation of linear diffeomorphism in wavelet space"

/************************************************************************
 * PPPPropagatorLin
 * @TechReport{Xie2003DFG,
 *   author = {Q. Xie and M. Holschneider and M. Kulesh},
 *   title = {Some remarks on linear diffeomorphisms in wavelet space},
 *   institution = {Preprint series of the DFG priority program 1114 "Mathematical methods for time series analysis and digital image processing"},
 *   month = {July},
 *   number = {37},
 *   year = {2003},
 *   urllink = {http://www.math.uni-bremen.de/zetem/DFG-Schwerpunkt/preprints/pdf/37.pdf}
 * }
************************************************************************/
class PPPPropagatorLin : public PPPBaseObject
  {
  private:
    PPPVectorContainer<double>  _par;
    PPPWavelet                  *_wavelet1, *_wavelet2;
    double                      _alpha,_beta,_delta,_sign;

  public:

    PPPPropagatorLin(void) {
      setObjectName(PPPPROPAGATORLIN_NAME);
      _par.resize(6);
      _par[0] = _par[3] = 1.0;
      _par[1] = _par[2] = _par[4] = _par[5] = 0.0;
      _alpha = _beta = _delta = 0.0;
      _sign = 1.0;
      };

    const char *getInfo(void) {
      return getNotation();
      };

    inline PPPVectorContainer<double> & getParams(void) {
      return _par;
      };

    void evalWaveletPropag(PPPSpectrContainer<PPPcomplex> &aDest,
      PPPSpectrContainer<PPPcomplex> &aSour, PPPSpectrParams &aTransPar) {
      onMessage(PPPPROPAGATORLIN_WAV);
      PPPcomplex CmplPhi = evalInverseConst(aTransPar,_par[0],_par[2],_par[3]);
      aDest.assign(aSour);
      double f,b;
      unsigned i,j;
      for(i=0; i<aDest.voices(); i++) for(j=0; j<aDest.points(); j++)
        {
        f = aSour.getFreq(i);
        b = aSour.getTime(j);
        aDest(i,j) = CmplPhi*aSour.Get(fprim(f,b),tprim(f,b));
        }
      // Notation
      strstream str;
      char s1[100], s2[100];
      sprintf(s1,"t` = %f/f %+f t %+f",_par[2],_par[3],_par[5]);
      sprintf(s2,"f` = 1/(%f/f %+f t %+f)",_par[0],_par[1],_par[4]);
      str << "  Wg2(t,f) = CmplPhi*Wg1(t`,f`)" << endl;
      str << "  CmplPhi = " << CmplPhi << endl;
      str << "  " << s1 << endl;
      str << "  " << s2 << ends;
      onNotation(str.str());
      return;
      };

  private:

    double tprim(double f, double b) {
      return _par[2]/f + _par[3]*b + _par[5];
      };

    double fprim(double f, double b) {
      return 1.0/(_par[0]/f + _par[1]*b + _par[4]);
      };

    /** Trapezoidal rule ***********************************************/
    /** aSize is count of point = power of two, f.e 128 **/
    PPPcomplex evalIntegralTrapez(double a, double b, int aSize) {
      double x,tnm,del;
      PPPcomplex s, sum;
      int it,j;
      if (aSize == 1)
        {
        return (s=0.5*(b-a)*(evalKernel(a)+evalKernel(b)));
        }
      else
        {
        for (it=1,j=1;j<aSize-1;j++) it <<= 1;
        tnm=it;
        del=(b-a)/tnm;
        x=a+0.5*del;
        for (sum=0.0,j=1;j<=it;j++,x+=del) sum += evalKernel(x);
        s=((b-a)*sum/tnm);
        return s;
        }
      };

    PPPcomplex evalKernel(double aOmega) {
      PPPcomplex g = _wavelet1->evalCmplFreq(_sign*_delta*aOmega/(2.0*M_PI));
      PPPcomplex h = _wavelet2->evalCmplFreq(_sign*_alpha*aOmega/(2.0*M_PI));
      PPPcomplex z(0.0, _sign*_beta*aOmega);
      return conj(g)*h*exp(z)/aOmega;
      };

    PPPcomplex evalInverseConst(PPPSpectrParams &aT, double aAlpha, double aBeta, double aDelta) {
      _alpha = aAlpha;
      _beta = aBeta;
      _delta = aDelta;
      _wavelet1 = aT.getWavelet(PPPSpectrParams::WSTdirect);
      _wavelet2 = aT.getWavelet(PPPSpectrParams::WSTinverse);
      _sign = 1.0;
      PPPcomplex z1 = evalIntegralTrapez(0,30,12);
      PPPcomplex z2 = 1.0/(2.0*z1);
      return z2;
      };

  }; // end of object


#endif
