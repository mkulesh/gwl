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

#ifndef _PPPOPTI
#define _PPPOPTI

#define PPPOPTI_NAME        "virtual optimization object"

/***********************************************************************
 * PPPOpti
 ************************************************************************/
class PPPOpti : public PPPBaseObject
  {
  private:
    bool         _isInitialized;

  protected:
    unsigned     _points;           // Number of function points
    unsigned     _fixed;            // Number of fixed parameters
    unsigned     _calcfuncCount;    // Number of cost function calculation
    double       _chisq;            // Current and minimal value of cost function

  public:
    typedef enum {OTleven, OTsimann} OptiType;
    PPPVectorContainer<double>    optpar;     // Best parameter set
    PPPVectorContainer<int>       fixpar;     // Massiv of 0 and 1 for fixed and nonfixed parameters

    PPPOpti() {
      setObjectName(PPPOPTI_NAME);
      _isInitialized = false;
      };

    void resize(unsigned const aPointCount, unsigned const aParamCount) {
      _isInitialized = true;
      _points = aPointCount;
      optpar.resize(aParamCount);
      fixpar.resize(aParamCount);
      fixpar.assign(1);
      evalFixed();
      };
            
    inline bool         isInitialized() const { return _isInitialized; };
    inline unsigned     points() const { return _points; };
    inline unsigned     params() const { return optpar.size(); };
    inline unsigned     fixed() const { return _fixed; };
    inline unsigned     calcfuncCount() const { return _calcfuncCount; };
    inline double       chisq() const { return _chisq; };
    virtual void        setParams(PPPVectorContainer<unsigned> &aPar) = 0;
    virtual void        setParams(PPPVectorContainer<double> &aPar) = 0;
    virtual void        setParams(PPPMatrixContainer<double> &aPar) = 0;
    virtual void        calcfunc(PPPVectorContainer<double> &) = 0;
    virtual double      getcostfunc(unsigned) = 0;
    virtual double      getder(unsigned,unsigned) = 0;
    virtual OptiType    getOptiType() = 0;
    virtual void        Minimize() = 0;

    unsigned evalFixed(void) {
      _fixed=0;
      for(unsigned j=0;j<params();j++) if(fixpar[j]) _fixed++;
      return _fixed;
      };

    unsigned evalFixed(bool aFixPhi, bool aFixAtn, unsigned aSize) {
      fixpar.assign(1);
      if(aFixAtn) for(unsigned i=aSize;i<params();i++)
        fixpar[i]=0;
      if(aFixPhi) fixpar[0]=0;
      return evalFixed();
      };

  }; // end of object


#endif

 
