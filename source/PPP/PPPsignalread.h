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
#ifndef _PPPSIGNALREAD
#define _PPPSIGNALREAD

#define PPPSIGNALREAD_NAME     "imported signal"
#define PPPSIGNALREAD_GEONAME  "geophone"
#define PPPSIGNALREAD_CHANNAME "channels"
#define PPPSIGNALREAD_TYPE01   "ASCII file (without time column)"
#define PPPSIGNALREAD_TYPE02   "ASCII file (with time column)"
#define PPPSIGNALREAD_TYPE03   "plain ASCII file"
#define PPPSIGNALREAD_TYPE04   "channels of 1D seismogramm ASCII file (without time column)"
#define PPPSIGNALREAD_TYPE05   "channels of 1D seismogramm ASCII file (with time column)"
#define PPPSIGNALREAD_TYPE06   "channels of 2D seismogramm ASCII files (without time column)"
#define PPPSIGNALREAD_TYPE07   "channels of 2D seismogramm ASCII files (with time column)"
#define PPPSIGNALREAD_TYPE08   "channels of 3D seismogramm ASCII files (without time column)"
#define PPPSIGNALREAD_TYPE09   "channels of 3D seismogramm ASCII files (with time column)"
#define PPPSIGNALREAD_ROT      "rotation of channels: "
#define PPPSIGNALREAD_RESAMPLE "resamplin of seismogram: "
#define PPPSIGNALREAD_T02P     "lengthening to power of two points"
#define PPPSIGNALREAD_MULT     "multiplication with constant="
#define PPPSIGNALREAD_LONG     "input signal is too long in procedure: "

/************************************************************************
 * PPPSignalRead
 ************************************************************************/
template<class AType> class PPPSignalRead : public PPPBaseTemplate<AType>
  {
  private:
    string                        _file, _fileX, _fileY, _fileZ;
    double                        _sample, _tmin, _tmax;
    bool                          _toPowerOfTwo;
    unsigned                      _geophone;
    PPPVectorContainer<unsigned>  _channels;
    PPPVectorContainer<double>    _rotation;
    unsigned                      _resample;
    AType                         _mult;
    unsigned                      _bps, _wchannels;

  public:

    PPPSignalRead(void):  _file(""), _fileX(""), _fileY(""), _fileZ(""),
      _toPowerOfTwo(false), _tmin(0.0), _tmax(0.0), _sample(0.0), _mult(1.0),
      _geophone(0), _resample(0), _bps(0), _wchannels(0)
      {
      PPPBaseObject::setObjectName(PPPSIGNALREAD_NAME);
      _rotation.resize(0);
      _channels.resize(2);
      _channels[0] = 0;
      _channels[1] = 1;
      };

    inline double getSamplingFreq(void) const  { return _sample; };
    inline void setSamplingFreq(double asample) { _sample = asample; };

    inline string & getFileName(void) { return _file; };
    inline void setFileName(const string &afile) { _file=afile; };
    inline void setFileName(const string &afilex, const string &afiley) { _fileX=afilex; _fileY=afiley; };

    inline double getTmin(void) const  { return _tmin; };
    inline void setTmin(double atmin) { _tmin=atmin; };

    inline double getTmax(void) const  { return _tmax; };
    inline void setTmax(double atmax) { _tmax=atmax; };

    inline bool getToPowerOfTwo(void) const { return _toPowerOfTwo; };
    inline void setToPowerOfTwo(bool atoPowerOfTwo) { _toPowerOfTwo=atoPowerOfTwo; };

    inline unsigned getGeophone(void) const { return _geophone; };
    inline void setGeophone(unsigned ageophone) { _geophone=ageophone; };

    inline unsigned getResample(void) const { return _resample; };
    inline void setReample(unsigned aresample) { _resample=aresample; };

    inline PPPVectorContainer<unsigned> & getChannels(void) { return _channels; };
    inline PPPVectorContainer<double> & getRotation(void) { return _rotation; };

    inline AType getMult(void) const  { return _mult; };
    inline void setMult(AType amult) { _mult=amult; };

    inline unsigned getBPS(void) const { return _bps; };
    inline unsigned getWChannels(void) const { return _wchannels; };

    void readFunctionASCII(PPPSignalContainer<AType> &aDest, bool aIsTime=true) {
      PPPSignalContainer<double> func;
      func.read(_file,aIsTime,_sample,_tmin,_tmax);
      aDest.prepare(func.points(), 1, func.getAxis(), PPPBaseObject::getObjectName());
      for(unsigned i=0; i<aDest.points(); i++)
        {
        PPPcomplex z = (PPPBaseTemplate<AType>::isComplex())?
          PPPcomplex(func(i,_channels[0]), func(i,_channels[1])) : PPPcomplex(func(i,_channels[0]), 0.0);
        PPPBaseTemplate<AType>::cmplConvert(aDest(i,0), z);
        }
      // Notation
      strstream str;
      if(!aIsTime) str << PPPSIGNALREAD_TYPE01;
              else str << PPPSIGNALREAD_TYPE02;
      if(PPPBaseTemplate<AType>::isReal() && _channels.size()>0)
        str << ": " << PPPSIGNALREAD_CHANNAME << "=" << _channels[0];
      if(!PPPBaseTemplate<AType>::isReal() && _channels.size()>0)
        str << ": " << PPPSIGNALREAD_CHANNAME << "=" << _channels.vectorToStr();
      str << ends;
      PPPBaseObject::onNotation(str.str());
      };

    void readFunctionGolm(PPPSignalContainer<AType> &aDest) {
      PPPSpectrContainer<double> func;
      _readPlaneASCII(func,_file,_tmin,_tmax);
      aDest.prepare(func.points(), 1, func.getTime(), PPPBaseObject::getObjectName());
      for(unsigned i=0; i<aDest.points(); i++)
        {
        PPPcomplex z = (PPPBaseTemplate<AType>::isComplex())?
          PPPcomplex(func(_geophone,i,_channels[0]), func(_geophone,i,_channels[1])) : PPPcomplex(func(_geophone,i,_channels[0]), 0.0);
        PPPBaseTemplate<AType>::cmplConvert(aDest(i,0), z);
        }
      // Notation
      strstream str;
      str << func.getNotation() << endl << "  " << PPPSIGNALREAD_GEONAME << "=" << _geophone << ", ";
      if(PPPBaseTemplate<AType>::isReal() && _channels.size()>0)
        str << PPPSIGNALREAD_CHANNAME << "=" << _channels[0];
      if(!PPPBaseTemplate<AType>::isReal() && _channels.size()>0)
        str << PPPSIGNALREAD_CHANNAME << "=" << _channels.vectorToStr();
      str << ends;
      PPPBaseObject::onNotation(str.str());
      };

    void readFunctionWav(PPPSignalContainer<AType> &aDest) {
      PPPWavIn infile(_file);
      _bps = infile.getBPS();
      _wchannels = infile.channels();
      double _rate = infile.getRate()*infile.channels();
      aDest.setObjectName(PPPBaseObject::getObjectName());
      long int ind1 = (_tmax==0)? 0 : (long int)(_tmin*_rate);
      long int ind2 = (_tmax==0)? infile.points() : (long int)(_tmax*_rate+infile.channels());
      if(ind2 > infile.points()) ind2 = infile.points(); 
      long int maxind = (long int)powl(2l,(long int)(7*sizeof(unsigned)));
      if((ind2-ind1)>maxind)
        PPPBaseObject::onError(PPPSIGNALREAD_LONG + string("readFunctionWav"));
      unsigned length = (unsigned)(ind2-ind1)/infile.channels();
      aDest.realloc(length, 1);
      long int i=ind1;
      PPPcomplex z;
      for(unsigned k=0; k<length; k++, i++)
        {
        aDest.getAxis()[k] = (double)i/_rate;
        if(infile.channels() == 1) z = PPPcomplex(infile.getData(i),0.0);
        else { z = PPPcomplex(infile.getData(i),infile.getData(i++)); }
        PPPBaseTemplate<AType>::cmplConvert(aDest(k,0), z);
        }
      aDest.getAxis().setParams(PPPAxis::ATlin);  
      aDest.getAxis().setObjectName(PPPSIGNALCONTAINER_AXISNAME);
      // Notation
      PPPBaseObject::onNotation(infile.getInfo());
      };

    void readSeisASCII(PPPSignalContainer<AType> &aDest, bool aIsTime=true) {
      PPPSignalContainer<double> func;
      func.read(_file,aIsTime,_sample,_tmin,_tmax);
      unsigned newcount = (_channels.size() == 0)? func.channels() : _channels.size();
      aDest.prepare(func.points(), newcount, func.getAxis(), PPPBaseObject::getObjectName());
      for(unsigned j=0; j<newcount; j++)
        for(unsigned i=0; i<func.points(); i++)
          if(_channels.size() == 0)
            aDest(i,j) = func(i,j);
          else
            aDest(i,j) = func(i,_channels[j]);
      // Notation
      strstream str;
      if(!aIsTime) str << PPPSIGNALREAD_TYPE04;
              else str << PPPSIGNALREAD_TYPE05;
      if(_channels.size()>0)
        str << ": " << PPPSIGNALREAD_CHANNAME << "=" << _channels.vectorToStr();
      str << ends;
      PPPBaseObject::onNotation(str.str());
      };

    void readSeisGolm(PPPSignalContainer<AType> &aDest) {
      PPPSpectrContainer<double> func;
      _readPlaneASCII(func,_file,_tmin,_tmax);
      unsigned newcount = (_channels.size() == 0)? func.channels() : _channels.size();
      aDest.prepare(func.points(), newcount, func.getTime(), PPPBaseObject::getObjectName());
      for(unsigned j=0; j<newcount; j++)
        for(unsigned i=0; i<func.points(); i++)
          if(_channels.size() == 0)
            aDest(i,j) = func(_geophone,i,j);
          else
            aDest(i,j) = func(_geophone,i,_channels[j]);
      // Notation
      strstream str;
      str << func.getNotation() << endl << "  " << PPPSIGNALREAD_GEONAME << "=" << _geophone << ", ";
      if(_channels.size()>0)
        str << PPPSIGNALREAD_CHANNAME << "=" << _channels.vectorToStr();
      str << ends;
      PPPBaseObject::onNotation(str.str());
      };

    void readSeis2Files(PPPSignalContainer<AType> &aDest, bool aIsTime=true) {
      PPPSignalContainer<double> sx,sy;
      sx.read(_fileX,aIsTime,_sample,_tmin,_tmax);
      sy.read(_fileY,aIsTime,_sample,_tmin,_tmax);
      if(sx.channels() != sy.channels()) PPPBaseObject::onError(FILE_ERRFORM+string("read2DSeis"));
      if(sx.points() != sy.points()) PPPBaseObject::onError(FILE_ERRFORM+string("read2DSeis"));
      unsigned newcount = (_channels.size() == 0)? sx.channels() : _channels.size();
      aDest.prepare(sx.points(), newcount, sx.getAxis(), PPPBaseObject::getObjectName());
      for(unsigned j=0; j<newcount; j++)
        for(unsigned i=0; i<sx.points(); i++)
          if(_channels.size() == 0)
            PPPBaseTemplate<AType>::cmplConvert(aDest(i,j), PPPcomplex(sx(i,j),sy(i,j)));
          else
            PPPBaseTemplate<AType>::cmplConvert(aDest(i,j), PPPcomplex(sx(i,_channels[j]),sy(i,_channels[j])));
      // Notation
      strstream str;
      if(!aIsTime) str << PPPSIGNALREAD_TYPE06;
              else str << PPPSIGNALREAD_TYPE07;
      if(_channels.size()>0)
        str << ": " << PPPSIGNALREAD_CHANNAME << "=" << _channels.vectorToStr();
      str << ends;
      PPPBaseObject::onNotation(str.str());
      };

    void postProcessing(PPPSignalContainer<AType> &aDest) {
      if(_resample>0)
        {
        aDest.resamplig(_resample);
        strstream str3;
        str3 << endl << "  " << PPPSIGNALREAD_RESAMPLE << _resample << ends;
        PPPBaseObject::onNotation(str3.str());
        }
      if(_toPowerOfTwo)
        {
        strstream str1;
        aDest.toPowerOfTwo();
        str1 << endl << "  " << PPPSIGNALREAD_T02P << ends;
        PPPBaseObject::onNotation(str1.str());
        }
      if(_rotation.size() == 3 && aDest.channels() > 1)
        {
        unsigned aChan1 = (unsigned)_rotation[0];
        unsigned aChan2 = (unsigned)_rotation[1];
        aDest.rotation(aChan1, aChan2, _rotation[2]);
        strstream str2;
        str2 << endl << "  " << PPPSIGNALREAD_ROT << " [" << aChan1 << "," << aChan2 << "] -> "
             << "[" << aChan1 << "`," << aChan2 << "`], Phi=" << _rotation[2] << ends;
        PPPBaseObject::onNotation(str2.str());
        }
      if(_mult != 1.0)
        {
        strstream str3;
        aDest.multiplication(_mult);
        str3 << endl << "  " << PPPSIGNALREAD_MULT << _mult << ends;
        PPPBaseObject::onNotation(str3.str());
        }
      };

  private:

    void _readPlaneASCII(PPPSpectrContainer<double> &aDest, const string &aName, double aMin=0.0, double aMax=0.0) {
      strstream stat1;
      stat1 << FILE_READSTA << aName << ends;
      PPPBaseObject::onMessage(stat1.str());
      vector<string> lines;
      PPPBaseObject::readLinesFromFile(aName,lines);
      char tegs[5][30]=
       {
        /*0*/   "#STA_CODE",
        /*1*/   "#STA_CHAN",
        /*2*/   "#SAMP_FREQ",
        /*3*/   "#NDAT",
        /*4*/   "#START_TIME"
       };
      char str[200], fn[200], str1[200], SignalNotation[200];
      fn[0]=0;
      unsigned i=0, j=0, k=0, flag1=1;
      double sfreq=1.0;
      for (unsigned m=0; m<lines.size(); ++m)
        {
        strcpy(str,lines[m].c_str());
        if(str[0] != '#') { flag1=0; k++; }
        if(str[0] == '#' && flag1 == 0) {flag1=1; j++; k=0; }
        if(strncmp(str,tegs[1],strlen(tegs[1])) == 0 && fn[0] == 0) strcpy(fn,str);
        if(strncmp(str,tegs[2],strlen(tegs[2])) == 0)
          {
          strcpy(str1,str+strlen(tegs[2])+1);
          sscanf(str1,"%lf",&sfreq);
          }
        if(strncmp(str,tegs[4],strlen(tegs[4])) == 0)
          {
          strcpy(SignalNotation,str+strlen(tegs[4])+1);
          }
        if(strcmp(str,fn) == 0) i++;
        }
      if(i==0) i=1;
      unsigned s1=i; // number of geophone
      unsigned s2=(j+1)/i; // number of channels
      unsigned s3=k; // number of point
      unsigned ind1=0, ind2=s3;
      if(aMin != 0 || aMax !=0)
        {
        ind1 = (unsigned)(aMin*sfreq);
        ind2 = (unsigned)(aMax*sfreq)+1;
        }
      aDest.resize(s1,ind2-ind1,s2);
      aDest.setNames(PPPSIGNALREAD_TYPE03, PPPSIGNALCONTAINER_AXISNAME, PPPSIGNALREAD_GEONAME);
      aDest.getTime().filllin(sfreq,ind1/sfreq);
      aDest.getFreq().filllin(1,1);
      for(j=0; j<aDest.channels(); j++) aDest.getChannel(j).setObjectName("Unknown");
      i = j = k = 0;  flag1=1;
      double val;
      for (unsigned m=0; m<lines.size(); ++m)
        {
        strcpy(str,lines[m].c_str());
        if(str[0] != '#')
          {
          sscanf(str,"%lf",&val);
          if(k<ind2 && k>=ind1) aDest(i,k-ind1,j)=val;
          flag1=0; k++;
          }
        if(str[0] == '#' && flag1 == 0)
          {
          flag1=1; j++; k=0;
          if(j == s2) { i++; j=0; }
          }
        if(strncmp(str,tegs[1],strlen(tegs[1])) == 0)
          aDest.getChannel(j).setObjectName(str);
        }
      strstream stat2;
      stat2 << aDest.getObjectName() << ": " << aDest.voices() << " x ";
      for(j=0; j<aDest.channels(); j++) stat2 << aDest.getChannel(j).getObjectName() << "; ";
      stat2 << ends;
      aDest.onNotation(stat2.str());
      return;
      };

  };  // end of object


#endif
