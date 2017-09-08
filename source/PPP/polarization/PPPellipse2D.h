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
#ifndef _PPPELLIPSE2D
#define _PPPELLIPSE2D

/************************************************************************
 * PPPellipse2D
 ************************************************************************/
class  PPPellipse2D : public PPPEllipse<2>
  {
  public:

    PPPellipse2D(void) {
      set(0.0);
      };

    PPPellipse2D(double const b) {
      set(b);
      };

    void operator= (double const b) {
      set(b);
      };

    friend PPPellipse2D operator+ (PPPellipse2D &a, PPPellipse2D &b) {
      PPPellipse2D r;
      r.sum(a,b);
      return r;
      };

    friend PPPellipse2D operator* (PPPellipse2D &a, double b) {
      PPPellipse2D r;
      r.mult(a,b);
      return r;
      };

    friend PPPellipse2D operator* (double a, PPPellipse2D &b) {
      PPPellipse2D r;
      r.mult(b,a);
      return r;
      };

    friend ostream& operator << (ostream &aDest, PPPellipse2D const &aSour) {
      aDest << showpos << scientific << aSour.absmin();
      aDest << showpos << scientific << " " << aSour.absmax();
      aDest << showpos << scientific << " " << aSour.tilt_rmax();
      aDest << showpos << scientific << " " << aSour.phasediff();
      aDest << showpos << scientific << " " << aSour.ratio();
      aDest << showpos << scientific << " " << aSour.instant(0);
      aDest << showpos << scientific << " " << aSour.instant(1);
      return aDest;
      };


  };  // end of object


/************************************************************************
 * PPPellipse2D Transforms
 ************************************************************************/
class PPPellipse2Drmin : public unary_function<PPPellipse2D, double>
  {
  public:
  inline double operator()(const PPPellipse2D &aVal) const { return aVal.absmin(); };
  inline unsigned getComponentType(void) const { return TErmin; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_RMIN; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERMIN; };
  };

class PPPellipse2Drmax : public unary_function<PPPellipse2D, double>
  {
  public:
  inline double operator()(const PPPellipse2D &aVal) const { return aVal.absmax(); };
  inline unsigned getComponentType(void) const { return TErmax; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_RMAX; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERMAX; };
  };

class PPPellipse2Dtilt : public unary_function<PPPellipse2D, double>
  {
  private:
    bool isDegree;
  public:
  PPPellipse2Dtilt (bool aDegree = false) {
    isDegree = aDegree;
    };
  inline double operator()(const PPPellipse2D &aVal) const {
    if(isDegree) return (180.0*aVal.tilt_rmax()/M_PI);
    else return aVal.tilt_rmax();
    };
  inline unsigned getComponentType(void) const { return TEtilt; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_TILT; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODETILT; };
  };

class PPPellipse2Dphasediff : public unary_function<PPPellipse2D, double>
  {
  private:
    bool isDegree;
  public:
  PPPellipse2Dphasediff (bool aDegree = false) {
    isDegree = aDegree;
    };
  inline double operator()(const PPPellipse2D &aVal) const {
    if(isDegree) return (180.0*aVal.phasediff()/M_PI);
    else return aVal.phasediff();
    };
  inline unsigned getComponentType(void) const { return TEphasediff; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_PHDIFF; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPHDIFF; };
  };

class PPPellipse2Dratio : public unary_function<PPPellipse2D, double>
  {
  public:
  inline double operator()(const PPPellipse2D &aVal) const { return aVal.ratio(); };
  inline unsigned getComponentType(void) const { return TEratio; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_RATIO; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERATIO; };
  };

class PPPellipse2Dwx : public unary_function<PPPellipse2D, double>
  {
  public:
  inline double operator()(const PPPellipse2D &aVal) const { return aVal.instant(0); };
  inline unsigned getComponentType(void) const { return TEwx; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_WX; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEWX; };
  };

class PPPellipse2Dwy : public unary_function<PPPellipse2D, double>
  {
  public:
  inline double operator()(const PPPellipse2D &aVal) const { return aVal.instant(1); };
  inline unsigned getComponentType(void) const { return TEwy; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_WY; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEWY; };
  };

class PPPellipse2Dphase1 : public unary_function<PPPellipse2D, double>
  {
  public:
  inline double operator()(const PPPellipse2D &aVal) const { return aVal.phase(1); };
  inline unsigned getComponentType(void) const { return TEphase1; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_PHASE1; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPHASE1; };
  };

class PPPellipse2Dphase2 : public unary_function<PPPellipse2D, double>
  {
  public:
  inline double operator()(const PPPellipse2D &aVal) const { return aVal.phase(2); };
  inline unsigned getComponentType(void) const { return TEphase2; };
  inline const char *getComponentName(void) const { return PPPELLIPSE_PHASE2; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPHASE2; };
  };

#endif
