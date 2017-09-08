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
#ifndef _PPPELLI3D
#define _PPPELLI3D


/************************************************************************
 * PPPellipse3D
 ************************************************************************/
class  PPPellipse3D : public PPPEllipse<3>
  {
  public:

    PPPellipse3D(void) {
      set(0.0);
      };

    PPPellipse3D(double const b) {
      set(b);
      };

    void operator= (double const b) {
      set(b);
      };

    friend PPPellipse3D operator* (PPPellipse3D &a, double b) {
      PPPellipse3D r;
      r.mult(a,b);
      return r;
      };

    friend PPPellipse3D operator* (double a, PPPellipse3D &b) {
      PPPellipse3D r;
      r.mult(b,a);
      return r;
      };

    friend ostream& operator << (ostream &aDest, PPPellipse3D const &aSour) {
        aDest << showpos << scientific << " " << aSour.planarity_cosx();
        aDest << showpos << scientific << " " << aSour.planarity_cosy();
        aDest << showpos << scientific << " " << aSour.planarity_cosz();
        aDest << showpos << scientific << " " << aSour.absmin();
        aDest << showpos << scientific << " " << aSour.absmax();
        aDest << showpos << scientific << " " << aSour.ratio();
        aDest << showpos << scientific << " " << aSour.absmidd();
        aDest << showpos << scientific << " " << aSour.ratio1();
        aDest << showpos << scientific << " " << aSour.ratio2();
        aDest << showpos << scientific << " " << aSour.instant(0);
        aDest << showpos << scientific << " " << aSour.instant(1);
        aDest << showpos << scientific << " " << aSour.instant(2);
      return aDest;
      };


  };  // end of object


/************************************************************************
 * PPPellipse3D Transforms
 ************************************************************************/
class PPPellipse3Drmin : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.absmin(); };
  inline unsigned getComponentType(void) const { return TErmin; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_RMIN; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERMIN; };
  };

class PPPellipse3Drmax : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.absmax(); };
  inline unsigned getComponentType(void) const { return TErmax; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_RMAX; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERMAX; };
  };

class PPPellipse3Drmidd : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.absmidd(); };
  inline unsigned getComponentType(void) const { return TErmidd; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_RMIDD; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERMIDD; };
  };

class PPPellipse3Dratio : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.ratio(); };
  inline unsigned getComponentType(void) const { return TEratio; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_RATIO; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERATIO; };
  };

class PPPellipse3Dratio1 : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.ratio1(); };
  inline unsigned getComponentType(void) const { return TEratio1; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_RATIO1; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERATIO1; };
  };

class PPPellipse3Dratio2 : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.ratio2(); };
  inline unsigned getComponentType(void) const { return TEratio2; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_RATIO2; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODERATIO2; };
  };

class PPPellipse3Dwx : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.instant(0); };
  inline unsigned getComponentType(void) const { return TEwx; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_WX; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEWX; };
  };

class PPPellipse3Dwy : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.instant(1); };
  inline unsigned getComponentType(void) const { return TEwy; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_WY; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEWY; };
  };

class PPPellipse3Dwz : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.instant(2); };
  inline unsigned getComponentType(void) const { return TEwz; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_WZ; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEWZ; };
  };

class PPPellipse3Dplanarx : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.planarity_x(); };
  inline unsigned getComponentType(void) const { return TEplanarx; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_PLANARX; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPLANARX; };
  };

class PPPellipse3Dplanary : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.planarity_y(); };
  inline unsigned getComponentType(void) const { return TEplanary; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_PLANARY; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPLANARY; };
  };

class PPPellipse3Dplanarz : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.planarity_z(); };
  inline unsigned getComponentType(void) const { return TEplanarz; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_PLANARZ; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPLANARZ; };
  };

class PPPellipse3Dplancosx : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.planarity_cosx(); };
  inline unsigned getComponentType(void) const { return TEplancosx; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_PLANCOSX; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPLANCOSX; };
  };

class PPPellipse3Dplancosy : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.planarity_cosy(); };
  inline unsigned getComponentType(void) const { return TEplancosy; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_PLANCOSY; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPLANCOSY; };
  };

class PPPellipse3Dplancosz : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.planarity_cosz(); };
  inline unsigned getComponentType(void) const { return TEplancosz; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_PLANCOSZ; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEPLANCOSZ; };
  };

class PPPellipse3Dsignratio : public unary_function<PPPellipse3D, double>
  {
  public:
  inline double operator()(const PPPellipse3D &aVal) const { return aVal.sign_ratio(); };
  inline unsigned getComponentType(void) const { return TEsignratio; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_SRATIO; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODESRATIO; };
  };

class PPPellipse3Ddip : public unary_function<PPPellipse3D, double>
  {
  private:
    bool isDegree;
  public:
  PPPellipse3Ddip (bool aDegree = false) {
    isDegree = aDegree;
    };
  inline double operator()(const PPPellipse3D &aVal) const {
    if(isDegree) return (180.0*aVal.dip_rmax()/M_PI);
    else return aVal.dip_rmax();
    };
  inline unsigned getComponentType(void) const { return TEdip; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_DIP; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEDIP; };
  };

class PPPellipse3Dazimuth : public unary_function<PPPellipse3D, double>
  {
  private:
    bool isDegree;
    bool isModPi;
  public:
  PPPellipse3Dazimuth (bool aDegree = false, bool aModPi = false) {
    isDegree = aDegree;
    isModPi = aModPi;
    };
  inline double operator()(const PPPellipse3D &aVal) const {
    double val = aVal.azimuth_rmax();
    if(isModPi)
      {
      if (val > M_PI/2.0) val = val - M_PI;
      if (val < -M_PI/2.0) val = val + M_PI;
      }
    if(isDegree) return (180.0*val/M_PI);
    else return val;
    };
  inline unsigned getComponentType(void) const { return TEazimuth; };
  inline const char *getComponentName(void) const {  return PPPELLIPSE_AZIMUTH; };
  inline const char *getComponentCode(void) const { return PPPELLIPSE_CODEAZIMUTH; };
  };



#endif
