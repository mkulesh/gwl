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
#ifndef _PPPELLIPSE3DFILTER
#define _PPPELLIPSE3DFILTER

#define PPELLIPSE3DFILTER_NAME       "3D polarization filter"
#define PPELLIPSE3DFILTER_LINXY      "X-Y plane linear filter"
#define PPELLIPSE3DFILTER_LINXZ      "X-Z plane linear filter"
#define PPELLIPSE3DFILTER_LINYZ      "Y-Z plane linear filter"
#define PPELLIPSE3DFILTER_ELLIXY     "X-Y plane elliptic filter"
#define PPELLIPSE3DFILTER_ELLIXZ     "X-Z plane elliptic filter"
#define PPELLIPSE3DFILTER_ELLIYZ     "Y-Z plane elliptic filter"
#define PPELLIPSE3DFILTER_NORMTA     "normalization of spectrum to tilt-angle"
#define PPELLIPSE3DFILTER_PROGRADE   "prograde filter"
#define PPELLIPSE3DFILTER_RETROGRADE "retrograde filter"
#define PPELLIPSE3DFILTER_LIN        "linear filter"
#define PPELLIPSE3DFILTER_ELLI       "elliptic filter"
#define PPELLIPSE3DFILTER_AZIMUTH    "azimuth filter"
#define PPELLIPSE3DFILTER_DIP        "abs dip angle filter"
#define PPELLIPSE3DFILTER_DEL        "number of erased point"
#define PPELLIPSE3DFILTER_NANGL      "normal angle"
#define PPELLIPSE3DFILTER_AZINT      "angle interval"

/************************************************************************
 * PPPEllipse3Dfilter
 ************************************************************************/
class PPPEllipse3Dfilter : public PPPBaseObject
  {
  public:

    typedef enum {NRFLinXY, NRFLinXZ, NRFLinYZ,
                  NRFElliXY, NRFElliXZ, NRFElliYZ, NRFNorm,
                  NRFprog, NRFretro, NRFlin, NRFelli,
                  NRFazimuth, NRFdip} FilterType;

  private:
  
    FilterType            _type;
    double                _tilt, _ratio, _int1, _int2;
    PPPellipse3Dratio     _compRatio;
    PPPellipse3Dsignratio _signRatio;
    PPPellipse3Dplancosx  _compPlancosx;
    PPPellipse3Dplancosy  _compPlancosy;
    PPPellipse3Dplancosz  _compPlancosz;
    PPPellipse3Dazimuth   _compAzimuth;
    PPPellipse3Ddip       _compDip;

  public:

    PPPEllipse3Dfilter(FilterType aType,double aTilt,double aRatio,double aInt1=0.0, double aInt2=0.0) {
      setObjectName(PPELLIPSE3DFILTER_NAME);
      _type = aType;
      _tilt = aTilt;
      _ratio = aRatio;
      _int1 = aInt1;
      _int2 = aInt2;
      };

    friend ostream& operator << (ostream& aDest, FilterType aSour) {
      switch(aSour)
        {
        case PPPEllipse3Dfilter::NRFLinXY:    aDest << PPELLIPSE3DFILTER_LINXY; break;
        case PPPEllipse3Dfilter::NRFLinXZ:    aDest << PPELLIPSE3DFILTER_LINXZ; break;
        case PPPEllipse3Dfilter::NRFLinYZ:    aDest << PPELLIPSE3DFILTER_LINYZ; break;
        case PPPEllipse3Dfilter::NRFElliXY:   aDest << PPELLIPSE3DFILTER_ELLIXY; break;
        case PPPEllipse3Dfilter::NRFElliXZ:   aDest << PPELLIPSE3DFILTER_ELLIXZ; break;
        case PPPEllipse3Dfilter::NRFElliYZ:   aDest << PPELLIPSE3DFILTER_ELLIYZ; break;
        case PPPEllipse3Dfilter::NRFNorm:     aDest << PPELLIPSE3DFILTER_NORMTA; break;
        case PPPEllipse3Dfilter::NRFprog:     aDest << PPELLIPSE3DFILTER_PROGRADE; break;
        case PPPEllipse3Dfilter::NRFretro:    aDest << PPELLIPSE3DFILTER_RETROGRADE; break;
        case PPPEllipse3Dfilter::NRFlin:      aDest << PPELLIPSE3DFILTER_LIN; break;
        case PPPEllipse3Dfilter::NRFelli:     aDest << PPELLIPSE3DFILTER_ELLI; break;
        case PPPEllipse3Dfilter::NRFazimuth:  aDest << PPELLIPSE3DFILTER_AZIMUTH; break;
        case PPPEllipse3Dfilter::NRFdip:      aDest << PPELLIPSE3DFILTER_DIP; break;
        }
      return aDest;
      };

    const char *getInfo(void) {
      strstream str;
      str << getType();
      if(getType()==NRFlin || getType()==NRFelli)
        str << ": " << PPPELLIPSE_RATIO << "=" << getRatio();
      else if(getType()==NRFLinXY || getType()==NRFLinXZ || getType()==NRFLinYZ || getType()==NRFElliXY || getType()==NRFElliXZ || getType()==NRFElliYZ)
        str << ": " << PPELLIPSE3DFILTER_NANGL << "=" << getTilt() << ", " << PPPELLIPSE_RATIO << "=" << getRatio();
      else if(getType()==NRFazimuth || getType()==NRFdip || getType()==NRFprog || getType()==NRFretro)
        str << ": " << PPELLIPSE3DFILTER_AZINT << "=[" << getInt1() << ", " << getInt2() << "]";
      str << ends;
      onNotation(str.str());
      return getNotation();
      };

    inline FilterType getType() { return _type; };
    inline double getTilt() { return _tilt; };
    inline double getRatio() { return _ratio; };
    inline double getInt1() { return _int1; };
    inline double getInt2() { return _int2; };

    bool isContent(PPPellipse3D &aVal) {
      double ratio = _compRatio(aVal);
      double tilt = 0.0, sratio = 0.0, azimuth = 0.0, dip = 0.0;
      switch(getType())
        {
        case NRFLinXY: case NRFElliXY:  tilt = _compPlancosz(aVal); break;
        case NRFLinXZ: case NRFElliXZ:  tilt = _compPlancosy(aVal); break;
        case NRFLinYZ: case NRFElliYZ:  tilt = _compPlancosx(aVal); break;
        }
      switch(getType())
        {
        case NRFlin: case NRFelli:
             if((getType() == NRFlin && ratio>=getRatio()) || (getType() == NRFelli && ratio<getRatio())) return false;
             break;
        case NRFLinXY: case NRFLinXZ: case NRFLinYZ:
             if((ratio<(1.0-getRatio()) && ratio>getRatio()) || tilt>getTilt()) return false;
             break;
        case NRFElliXY: case NRFElliXZ: case NRFElliYZ:
             if(ratio<getRatio() || ratio>1.0-getRatio() || tilt>getTilt()) return false;
             break;
        case NRFprog: case NRFretro:
             sratio = _signRatio(aVal);
             if(((getType() == NRFprog && sratio<0.0) || (getType() == NRFretro && sratio>0.0)) || (ratio<getInt1() || ratio>getInt2())) return false;
             break;
        case NRFazimuth:
             azimuth = _compAzimuth(aVal);
             if(azimuth<getInt1() || azimuth>getInt2()) return false;
             break;
        case NRFdip:
             dip = fabs(_compDip(aVal));
             if(dip<getInt1() || dip>getInt2()) return false;
             break;
        }
      return true;
      };

  };  // end of object

#endif
