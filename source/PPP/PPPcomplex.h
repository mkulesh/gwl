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

#ifndef _PPPCOMPLEX
#define _PPPCOMPLEX

static const unsigned        TCall  = 0;
static const unsigned        TCre   = 1;  //  real part
static const unsigned        TCim   = 2;  //  imaginary part
static const unsigned        TCabs  = 3;  //  modulus
static const unsigned        TCarg  = 4;  //  argument
static const unsigned        TCargf = 5;  //  argument with cut-off
static const unsigned        TCabsm = 6;  //  normalized modulus
static const unsigned        TCargm = 7;  //  modulated argument

#define CMPL_REAL       "real part"
#define CMPL_IMGR       "imaginary part"
#define CMPL_ABS        "modulus"
#define CMPL_ARG        "argument"
#define CMPL_ARGF       "argument with cut-off"
#define CMPL_ABSM       "normalized modulus"
#define CMPL_ARGM       "modulated argument"

#define CMPL_CODEREAL   "(re)"
#define CMPL_CODEIMGR   "(im)"
#define CMPL_CODEABS    "(mod)"
#define CMPL_CODEARG    "(arg)"
#define CMPL_CODEARGF   "(arg)"
#define CMPL_CODEABSM   "(mmod)"
#define CMPL_CODEARGM   "(marg)"

/************************************************************************
 * PPPcomplex
 ************************************************************************/
typedef complex<double> PPPcomplex;

ostream& operator << (ostream &aDest, PPPcomplex const &aSour) {
  aDest << showpos << scientific << aSour.real() << " " << aSour.imag();
  return aDest;
  };

istringstream& operator >> (istringstream &aDest, PPPcomplex &aSour) {
  double a,b;
  aDest >> a;
  aDest >> b;
  aSour = PPPcomplex(a,b);
  return aDest;
  };

class PPPcomplexRe : public unary_function<PPPcomplex, double>
  {
  public:
  inline double operator()(const PPPcomplex &aVal) const { return real(aVal); };
  inline unsigned getComponentType(void) const { return TCre; };
  inline const char *getComponentName(void) const { return CMPL_REAL; };
  inline const char *getComponentCode(void) const { return CMPL_CODEREAL; };
  };

class PPPcomplexIm : public unary_function<PPPcomplex, double>
  {
  public:
  inline double operator()(const PPPcomplex &aVal) const { return imag(aVal); };
  inline unsigned getComponentType(void) const { return TCim; };
  inline const char *getComponentName(void) const { return CMPL_IMGR; };
  inline const char *getComponentCode(void) const { return CMPL_CODEIMGR; };
  };

class PPPcomplexAbs : public unary_function<PPPcomplex, double>
  {
  public:
  inline double operator()(const PPPcomplex &aVal) const { return abs(aVal); };
  inline unsigned getComponentType(void) const { return TCabs; };
  inline const char *getComponentName(void) const { return CMPL_ABS; };
  inline const char *getComponentCode(void) const { return CMPL_CODEABS; };
  };

class PPPcomplexArg : public unary_function<PPPcomplex, double>
  {
  private:
    bool isDegree;
  public:
  PPPcomplexArg (bool aDegree = false) {
    isDegree = aDegree;
    };
  inline double operator()(const PPPcomplex &aVal) const {
    if(isDegree) return (180.0*arg(aVal)/M_PI);
    else return arg(aVal);
    };
  inline unsigned getComponentType(void) const { return TCarg; };
  inline const char *getComponentName(void) const { return CMPL_ARG; };
  inline const char *getComponentCode(void) const { return CMPL_CODEARG; };
  };

class PPPcomplexArgFilter : public unary_function<PPPcomplex, double>
  {
  private:
    bool isDegree;
    double maxmodulus;
    double percent;
    double defaultval;
  public:
  PPPcomplexArgFilter(const double aMax, const double aPer, const double aDefaultval = -M_PI,bool aDegree = false) {
    maxmodulus = aMax;
    percent = aPer;
    defaultval = aDefaultval;
    isDegree = aDegree;
    };
  inline double operator()(const PPPcomplex &aVal) const {
    double val = arg(aVal);
    if(100.0*abs(aVal)/maxmodulus < percent) val = defaultval;
    if(isDegree) return (180.0*val/M_PI);
    else return val;
    };
  inline unsigned getComponentType(void) const { return TCargf; };
  inline const char *getComponentName(void) const { return CMPL_ARGF; };
  inline const char *getComponentCode(void) const { return CMPL_CODEARGF; };
  };

class PPPcomplexAbsM : public unary_function<PPPcomplex, double>
  {
  private:
    double maxmodulus;

  public:
  PPPcomplexAbsM(const double aMax) {
    maxmodulus = aMax;
    };
  inline double operator()(const PPPcomplex &aVal) const {
    return abs(aVal)/maxmodulus;
    };
  inline unsigned getComponentType(void) const { return TCabsm; };
  inline const char *getComponentName(void) const { return CMPL_ABSM; };
  inline const char *getComponentCode(void) const { return CMPL_CODEABSM; };
  };


class PPPcomplexArgM : public unary_function<PPPcomplex, double>
  {
  private:
    double maxmodulus;

  public:
  PPPcomplexArgM(const double aMax) {
    maxmodulus = aMax;
    };
  inline double operator()(const PPPcomplex &aVal) const {
    return arg(aVal)*abs(aVal)/maxmodulus;
    };
  inline unsigned getComponentType(void) const { return TCargm; };
  inline const char *getComponentName(void) const { return CMPL_ARGM; };
  inline const char *getComponentCode(void) const { return CMPL_CODEARGM; };
  };

#endif
