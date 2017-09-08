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

#ifndef _PPPELLIPSE2DFILTER
#define _PPPELLIPSE2DFILTER

#define PPPELLIPSE2DFILTER_NAME     "2D polarization filter"
#define PPPELLIPSE2DFILTER_LINHOR   "linear-horisontal filter"
#define PPPELLIPSE2DFILTER_LINVERT  "linear-vertical filter"
#define PPPELLIPSE2DFILTER_ELLIHOR  "elliptic-horisontal filter"
#define PPPELLIPSE2DFILTER_ELLIVERT "elliptic-vertical filter"
#define PPPELLIPSE2DFILTER_LINHORS  "signed linear-horisontal filter"

/************************************************************************
 * PPPEllipse2Dfilter
 ************************************************************************/
class PPPEllipse2Dfilter : public PPPBaseObject
  {
  public:

    typedef enum {TRFLinHor, TRFLinVert, TRFElliHor, TRFElliVert, TRFLinHorS} FilterType;

  private:
    FilterType          _type;
    double              _tilt, _ratio;
    PPPellipse2Dratio   _compRatio;
    PPPellipse2Dtilt    _compTilt;

  public:
    PPPEllipse2Dfilter(FilterType aType,double aTilt,double aRatio) {
      setObjectName(PPPELLIPSE2DFILTER_NAME);
      _type = aType;
      _tilt = aTilt;
      _ratio = aRatio;
      };

    friend ostream& operator << (ostream& aDest, FilterType aSour) {
      switch(aSour)
        {
        case PPPEllipse2Dfilter::TRFLinHor:    aDest << PPPELLIPSE2DFILTER_LINHOR; break;
        case PPPEllipse2Dfilter::TRFLinVert:   aDest << PPPELLIPSE2DFILTER_LINVERT; break;
        case PPPEllipse2Dfilter::TRFElliHor:   aDest << PPPELLIPSE2DFILTER_ELLIHOR; break;
        case PPPEllipse2Dfilter::TRFElliVert:  aDest << PPPELLIPSE2DFILTER_ELLIVERT; break;
        case PPPEllipse2Dfilter::TRFLinHorS:   aDest << PPPELLIPSE2DFILTER_LINHORS; break;
        }
      return aDest;
      };

    const char *getInfo(void) {
      strstream str;
      str << getType() << ": " << PPPELLIPSE_TILT << "=" << getTilt() << ", " << PPPELLIPSE_RATIO << "=" << getRatio() << ends;
      onNotation(str.str());
      return getNotation();
      };

    inline FilterType getType() { return _type; };
    inline double getTilt() { return _tilt; };
    inline double getRatio() { return _ratio; };

    bool isContent(PPPellipse2D &aVal) {
      double ratio = _compRatio(aVal);
      double tilt = fabs(_compTilt(aVal));
      switch(getType())
        {
        case TRFLinHor:
             if((ratio<(1.0-getRatio()) && ratio>getRatio()) || tilt>getTilt()) return false;
             break;
        case TRFLinVert:
             if((ratio<(1.0-getRatio()) && ratio>getRatio()) || tilt<getTilt()) return false;
             break;
        case TRFElliHor:
             if((ratio<getRatio() || ratio>1.0-getRatio()) || tilt>getTilt()) return false;
             break;
        case TRFElliVert:
             if((ratio<getRatio() || ratio>1.0-getRatio()) || tilt<getTilt()) return false;
             break;
        case TRFLinHorS:
             if(fabs(ratio)>fabs(getRatio()) || fabs(tilt)>fabs(getTilt())) return false;
             break;
        }
      return true;
      };

  };  // end of object

#endif
