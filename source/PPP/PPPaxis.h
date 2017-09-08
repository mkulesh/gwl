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
#ifndef _PPPAXIS
#define _PPPAXIS

#define PPPAXIS_OBJVER  "AX1.4"
#define PPPAXIS_NAME    "axis"
#define PPPAXIS_SUND    "undefined scaling"
#define PPPAXIS_SLIN    "linear scaling"
#define PPPAXIS_SLOG    "logarithmic scaling"
#define PPPAXIS_SLOGE   "exponent logarithmic scaling"
#define PPPAXIS_POS     "positive axis"
#define PPPAXIS_NEG     "negative axis"
#define PPPAXIS_FULL    "positive-negative axis"
#define PPPAXIS_ERRVAL  "minimum or maximum value is zero in procedure: "
#define PPPAXIS_ERRTYPE "type of axis is incorrect in procedure: "

/************************************************************************
 * PPPAxis
 ***********************************************************************/
class PPPAxis : public PPPVectorContainer<double>
  {
  public:
    typedef enum {
      ATnone,      // without defined scaling
      ATlin,       // linear scaling
      ATlog,       // logarithmic scaling
      ATloge       // exponent logarithmic scaling
      } AxisType;

    typedef enum {
      ASplus,      // positive frequencies
      ASminus, 	// negative frequencies
      ASfull	// positive and negative frequencies
      } AxisSign;

  private:
    AxisType    _axistype;              // linear / logarithmic scaling
    AxisSign    _axissign;		// positive / negative / both frequencies
    double      _min;			// minimal value ( included )
    double      _max;			// maximal value ( included )
    double      _delta;			// spacing between samples
    double      _sample;		// sampling frequency

  public:

    PPPAxis(void) {
      _setdefault();
      };

    PPPAxis(unsigned const aSize) {
      _setdefault();
      resize(aSize);
      };

    PPPAxis(const PPPAxis &aSour) {
      _setdefault();
      assign(aSour);
      };

    PPPAxis(unsigned const aSize, const double aMin, const double aMax,
      const AxisType aAxistype, string const &aName) {
      _setdefault();
      setObjectName(aName);
      realloc(aSize);
      fill(aMin, aMax, aAxistype);
      };

    friend ostream& operator << (ostream& aDest, AxisType aSour) {
      switch(aSour)
        {
        case PPPAxis::ATnone: aDest << PPPAXIS_SUND; break;
        case PPPAxis::ATlin:  aDest << PPPAXIS_SLIN; break;
        case PPPAxis::ATlog:  aDest << PPPAXIS_SLOG; break;
        case PPPAxis::ATloge: aDest << PPPAXIS_SLOGE; break;
        }
      return aDest;
      };

    friend ostream& operator << (ostream& aDest, AxisSign aSour) {
      switch(aSour)
        {
        case PPPAxis::ASplus:  aDest << PPPAXIS_POS; break;
        case PPPAxis::ASminus: aDest << PPPAXIS_NEG; break;
        case PPPAxis::ASfull:  aDest << PPPAXIS_FULL; break;
        }
      return aDest;
      };

    void resize(const unsigned aSize) {
      PPPVectorContainer<double> :: resize(aSize);
      _setdefault();
      };

    void assign(const PPPAxis &aSour, unsigned aStart=0, unsigned aSize=0, int aOffset=1) {
      if(aSize != 0)
        realloc(aSize);
      else
        realloc(aSour.size());
      if(size() == 0)
        _setdefault();
      else
        {
        for(unsigned i=0; i<size(); i++) (*this)[i] = aSour[aStart+i*aOffset];
        setParams(aSour.getType());
        }
      _axissign = aSour.getSign();
      setObjectName(aSour.getObjectName());
      return;
      };

    void assign(PPPAxis &aSour, AxisSign aSign) {
      if(aSour.getSign() == ASfull && aSign == ASplus)
        {
        assign(aSour,aSour.size()/2,aSour.size()/2);
        _axissign = aSign;
        }
      else if(aSour.getSign() == ASfull && aSign == ASminus)
        {
        assign(aSour,aSour.size()/2-1,aSour.size()/2,-1);
        _axissign = aSign;
        }
      else assign(aSour);
      return;
      };

    inline AxisType getType(void) const { return _axistype; };
    inline AxisSign getSign() const     { return _axissign; };
    inline double getSamplingPeriod(void) const { return _delta; };
    inline double getSamplingFreq(void) const   { return _sample; };

    inline double getMin(void) const            { return _min; };
    inline double getMax(void) const            { return _max; };

    inline double getAbsMin() const {
      switch ( getSign() ) {
      	case ASplus : return getMin();
     	case ASminus: return -getMax();
      	case ASfull : return (*this)[ size()/2 ];
      	}
      return 0.0;
      };

    inline double getAbsMax() const {
      switch ( getSign() ) {
      	case ASplus : return getMax();
      	case ASminus: return -getMin();
      	case ASfull : return getMax();
	}
      return 0.0;
      };

    void setParams(const AxisType aAxistype) {
      if(size() == 0)
        {
        _setdefault();
        return;
        }
      _min = (*this)[0];
      _max = (*this)[size()-1];
      if(size() == 1)
        {
        _delta = 1;
        return;
        }
      switch(aAxistype)
        {
        case ATnone: case ATlin:  _delta = (*this)[1]-(*this)[0];  break;
        case ATlog:  _delta = log(_max/_min)/(log(2.0)*(double)(size()-1));  break;
        case ATloge: _delta = log(_max/_min)/(M_E-1.0); break;
        default:
          onError(PPPAXIS_ERRTYPE+string("setParams()"));
          break;
        }
      _axistype = aAxistype;
      _sample = (_axistype<2 && _delta!=0.0)? 1.0/_delta : 0.0;
      return;
      };

    void fill(const double aMin, const double aMax, const AxisType aAxistype) {
      unsigned i;
      switch(aAxistype)
        {
        case ATnone: break;
        case ATlin:
          _delta = (size()==1)? 1.0 : (aMax - aMin)/((double)(size()-1));
          for(i=0; i<size(); i++)
            (*this)[i] = aMin + (double)i*_delta;
          break;
        case ATlog:
          if(aMin == 0.0 || aMax == 0.0)
            onError(PPPAXIS_ERRVAL+string("fill()"));
          _delta = log(aMax/aMin)/(log(2.0)*(double)(size()-1));
          for(i=0; i<size(); i++)
            (*this)[i] = aMin*pow(2.0,(double)i*_delta);
          break;
        case ATloge:
          if(aMin == 0.0 || aMax == 0.0)
            onError(PPPAXIS_ERRVAL+string("fill()"));
          _delta = log(aMax/aMin)/(M_E-1.0);
          for(i=0; i<size(); i++)
            (*this)[i] = aMin*exp(_delta*(exp((double)i/(double)(size()-1))-1.0));
          break;
        default:
          onError(PPPAXIS_ERRTYPE+string("fill()"));
          break;
        }
      setParams(aAxistype);
      return;
      };

    void toPowerOfTwo(unsigned aBase=0) {
      unsigned k = size();
      double i,aMax;
      PPPVectorContainer<double> :: toPowerOfTwo(aBase);
      switch(getType())
        {
        case ATnone:
          for(i=1.0; k<size(); k++, i+=1.0)
            (*this)[k] = getMax() + getSamplingPeriod()*i;
          setParams(ATnone);
          break;
        case ATlin:
          aMax = getMin()+(double)(size()-1)*getSamplingPeriod();
          fill(getMin(), aMax, getType());
          break;
        case ATlog:
          aMax = getMin()*pow(2.0,(double)(size()-1)*getSamplingPeriod());
          fill(getMin(), aMax, getType());
          break;
        }
      return;
      };

    /**
     *  fill linearly the axis
     */
    void filllin(const double aSampleFreq, const double aMin = 0.0) {
      fill(aMin,aMin+((double)(size()-1))/aSampleFreq,ATlin);
      return;
      };

    void fillsign(const AxisSign aType) {
      PPPAxis aAxis;
      aAxis.assign((*this));
      switch(aType)
        {
        case ASplus: break;
        case ASminus:
             for(unsigned i=0; i<aAxis.size(); i++) (*this)[i] = -aAxis[aAxis.size()-1-i];
             _min = -aAxis.getMax();
             _max = -aAxis.getMin();
             break;
        case ASfull:
             resize(2*aAxis.size());
             for(unsigned i=0; i<aAxis.size(); i++) (*this)[i] = -aAxis[aAxis.size()-1-i];
             for(unsigned i=0; i<aAxis.size(); i++) (*this)[aAxis.size()+i] = aAxis[i];
             _min = -aAxis.getMax();
             _max = aAxis.getMax();
             _axistype = aAxis.getType();
             _delta = aAxis.getSamplingPeriod();
             _sample = aAxis.getSamplingFreq();
             break;
        }
      _axissign = aType;
      setObjectName(aAxis.getObjectName()); 
      };

    /**
     * Search an Ordered Table
     * Given an double array [0..n-1], and given a value aX, returns a value j
     * such that xx[j]<aX<=xx[j+1]. Array must be monotonic, either increasing
     * or decreasing. j=0 or j=n-1 is returned to indicate that aX is out of range.
     * Numerical Recipe in C: the art of scientific computing, 1992, Cambridge University Press
     */
    unsigned locateFloor(const double  aX) {
      unsigned count = size();
      int jl=0, ju=count, jm, ascnd=((*this)[count-1]>=(*this)[0]);
      while((ju-jl)>1)
        {
        jm=(ju+jl)>>1;
        if((aX>=(*this)[jm])==ascnd) jl=jm;
        else ju=jm;
        }
      return ((aX==(*this)[0])? 0 : ((aX==(*this)[count-1])? count-1 : jl));
      };

    const char *getInfo(void) {
      strstream aDest;
      aDest<<getObjectName()<<": "<<getTypeName()<<"["<<size()<<"]:";
      aDest<<" min="<<getMin()<<", max="<<getMax()<<", sampling period="<<getSamplingPeriod()<<", sampling frequency="<<getSamplingFreq();
      aDest<<" ("<<getType()<<"; "<<getSign()<<")"<<ends;
      onNotation(aDest.str());
      return getNotation();
      };

    /**
     *  file stream operations
     */
    void fwrite(FILE *stream) {
      fwrite_streaminfo(stream, getObjectVer(), sizeof(double));
      unsigned val;
      val=_axistype; std :: fwrite((void*)&val, sizeof(val), 1, stream);
      val=_axissign; std :: fwrite((void*)&val,sizeof(val),1,stream);
      std :: fwrite((void*)&_min, sizeof(double), 1, stream);
      std :: fwrite((void*)&_max, sizeof(double), 1, stream);
      std :: fwrite((void*)&_delta, sizeof(double), 1, stream);
      std :: fwrite((void*)&_sample, sizeof(double), 1, stream);
      setObjectVer(PPPVECTORCONTAINER_OBJVER);
      PPPVectorContainer<double> :: fwrite(stream);
      setObjectVer(PPPAXIS_OBJVER);
      };

    void fread(FILE *stream) {
      fread_streaminfo(stream, getObjectVer(), sizeof(double));
      unsigned val;
      std :: fread((void*)&val, sizeof(unsigned), 1, stream); _axistype = (AxisType)val;
      std :: fread((void*)&val,sizeof(val),1,stream); _axissign = (AxisSign)val;
      std :: fread((void*)&_min, sizeof(double), 1, stream);
      std :: fread((void*)&_max, sizeof(double), 1, stream);
      std :: fread((void*)&_delta, sizeof(double), 1, stream);
      std :: fread((void*)&_sample, sizeof(double), 1, stream);
      setObjectVer(PPPVECTORCONTAINER_OBJVER);
      PPPVectorContainer<double> :: fread(stream);
      setObjectVer(PPPAXIS_OBJVER);
      };


  private:

    void _setdefault(void) {
      setObjectVer(PPPAXIS_OBJVER);
      setObjectName(PPPAXIS_NAME);
      _axistype = ATnone;
      _axissign = ASplus;
      _delta = _sample = _min = _max = 0.0;
      };

  };  // end of object


#endif
