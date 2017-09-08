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

#ifndef _PPPAPPROXBISPLINE
#define _PPPAPPROXBISPLINE

/************************************************************************
 * PPPApproximateBspline
 ***********************************************************************/
class PPPApproximateBspline : public PPPApproximate
  {
  public:

    PPPApproximateBspline(unsigned aSize) {
      setObjectName(PPPAPPROXIMATE_SPL);
      _type = APTbspline;
      _params.resize(aSize);
      };

    const char *getInfo(void) {
      strstream str;
      str << "Func(f) = sum((a[i]*Bispline(3.0+(f-_xmin-i*dF)/dF), i=0..." << (size()-1) << endl;
      str << "  a[i] = " << _params.vectorToStr() << ends;
      onNotation(str.str());
      return getNotation();
      };

    double Func(double f) {
      if(f<_xmin || f>_xmax || size()==0) return 0.0;
      double res = 0.0;
      double dF = (_xmax-_xmin)/((double)size()-3.0);
      for(unsigned i=0; i<size(); i++)
        res += _params[i]*bsp4(3.0+(f-_xmin-(double)i*dF)/dF);
      return res;
      };

    double FuncDivX(double f) {
      if(f<_xmin || f>_xmax || size()==0) return 0.0;
      double res = 0.0;
      double dF = (_xmax-_xmin)/((double)size()-3.0);
      for(unsigned i=0; i<size(); i++)
        res += _params[i]*bsp4prim(3.0+(f-_xmin-(double)i*dF)/dF)/dF;
      return res;
      };

    double FuncDivXPar(double f, unsigned Ai) {
      if(f<_xmin || f>_xmax || Ai>=size()) return 0.0;
      double dF = (_xmax-_xmin)/((double)size()-3.0);
      double res = bsp4prim(3.0+(f-_xmin-(double)Ai*dF)/dF)/dF;
      return res;
      };

    double FuncDivPar(double f, unsigned Ai) {
      if(f<_xmin || f>_xmax || Ai>=size()) return 0.0;
      double dF = (_xmax-_xmin)/((double)size()-3.0);
      double res = bsp4(3.0+(f-_xmin-(double)Ai*dF)/dF);
      return res;
      };

  private:

    double bsp1(double b) {
      double res=0.0;
      if(b>=0.0 && b<=1.0) res = b;
      else if(b>=1.0 && b<=2.0) res = 2.0-b;
           else res=0.0;
      return res;
      };

    double bsp1prim(double b) {
      double res=0.0;
      if(b>=0.0 && b<=1.0) res = 1;
      else if(b>=1.0 && b<=2.0) res = -1;
           else res=0.0;
      return res;
      };

    double bsp3(double b) {
      double res=0.0;
      if(b>=0.0 && b<=1.0) res = 0.5*b*b;
      else if(b>=1.0 && b<=2.0)
           res = 0.5*(-2.0*b*b+6.0*b-3.0);
           else if(b>=2.0 && b<=3.0) res = 0.5*(3.0-b)*(3.0-b);
                else res=0.0;
      return res;
      };

    double bsp3prim(double b) {
      double res=0.0;
      if(b>=0.0 && b<=1.0) res = b;
      else if(b>=1.0 && b<=2.0) res = 0.5*(-4.0*b+6.0);
           else if(b>=2.0 && b<=3.0) res = -0.5*(2*b-6);
                else res=0.0;
      return res;
      };

    double bsp4(double b) {
      double res=0.0;
      if(b>=0.0 && b<=1.0) res = b*b*b/6.0;
      else if(b>=1.0 && b<=2.0) res = (-3.0*b*b*b+12.0*b*b-12.0*b+4.0)/6.0;
           else if(b>=2.0 && b<=3.0) res = (3.0*b*b*b-24.0*b*b+60.0*b-44.0)/6.0;
                else if (b>=3.0 && b<=4.0) res=(4.0-b)*(4.0-b)*(4.0-b)/6.0;
                     else res=0.0;
      return 1.5*res;
      };

    double bsp4prim(double b) {
      double res=0.0;
      if(b>=0.0 && b<=1.0) res = b*b/2.0;
      else if(b>=1.0 && b<=2.0) res = (-9.0*b*b+24.0*b-12.0)/6.0;
           else if(b>=2.0 && b<=3.0) res = (9.0*b*b-48.0*b+60.0)/6.0;
                else if(b>=3.0 && b<=4.0) res=-3.0*(4.0-b)*(4.0-b)/6.0;
                     else res=0.0;
      return 1.5*res;
      };

  }; // end of object


#endif
