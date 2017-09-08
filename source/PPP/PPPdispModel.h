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
#ifndef _PPPDISPERSIONMODEL
#define _PPPDISPERSIONMODEL

#define PPPDISPERSIONMODEL_OBJVER     "DM1.3"
#define PPPDISPERSIONMODEL_NAME       "dispersion model"
#define PPPDISPERSIONMODEL_WN         "wave number"
#define PPPDISPERSIONMODEL_ATN        "attenuation"
#define PPPDISPERSIONMODEL_FREQ       "frequency interval"
#define PPPDISPERSIONMODEL_EXT        "model from external file"
#define PPPDISPERSIONMODEL_ERRTYPE    "type of approximation is invalid in procedure: "
#define PPPDISPERSIONMODEL_ERRDEF     "type of approximation is not defined for this model"

/************************************************************************
 * PPPDispersionModel
************************************************************************/
class PPPDispersionModel : public PPPBaseObject
  {
  public:
  
    static const unsigned iWn         = 0;
    static const unsigned iWndiv      = 1;
    static const unsigned iCp         = 2;
    static const unsigned iCg         = 3;
    static const unsigned iAtn        = 4;
    static const unsigned iAtndiv     = 5;

  private:
    bool                                   _isAnalytical;
    PPPApproximate                         *_wn;
    PPPApproximate                         *_atn;
    bool                                   _isCausal;
    double                                 _freqmin, _freqmax;
    PPPSignalContainer<double>             _data;
    PPPSignalContainer<double>::InterpType _intrp;

  public:

    PPPDispersionModel(void) : _wn(NULL), _atn(NULL), _isCausal(false),
      _freqmin(0), _freqmax(0), _intrp(PPPSignalContainer<double>::ITsplin), _isAnalytical(true) {
      setObjectVer(PPPDISPERSIONMODEL_OBJVER);
      setObjectName(PPPDISPERSIONMODEL_NAME);
      };

    ~PPPDispersionModel(void) {
      if(isAtn()) delete _atn;
      if(isWn()) delete _wn;
      };

    const char *getInfo(void) {
      strstream aDest;
      aDest << getObjectName() << endl
            << "  " << PPPDISPERSIONMODEL_FREQ << ": [" << getFrMin() << ", " << getFrMax() << "]";
      if(isAnalytical())
        {
        aDest << endl << "  " << PPPDISPERSIONMODEL_WN << ": " << getWn().getInfo();
        aDest << endl << "  " << PPPDISPERSIONMODEL_ATN << ": " << getAtn().getInfo();
        }
      else
        aDest << endl << "  " << PPPDISPERSIONMODEL_EXT;
      aDest << ends;
      onNotation(aDest.str());
      return getNotation();
      };

    inline bool isAnalytical(void) const { return (_isAnalytical); };
    void setAnalytical(bool const aAnalytical) { _isAnalytical = aAnalytical; };

    void prepare(PPPApproximate::ApprType aWnType, PPPVectorContainer<double> &aWn,
      PPPApproximate::ApprType aAtnType, PPPVectorContainer<double> &aAtn) {
      if(isAtn()) delete _atn;
      if(isWn()) delete _wn;
      unsigned aWnSize = aWn.size();
      unsigned aAtnSize = aAtn.size();
      _isCausal = false;
      if(aWnType == PPPApproximate::APTcolecole || aWnType == PPPApproximate::APTtwogauss)
        {
        _isCausal = true;
        aAtnType = aWnType;
        aAtnSize = aWnSize;
        }
      switch(aWnType) {
        case PPPApproximate::APTvel:      _wn = _createApprox(1,aWnSize); break;
        case PPPApproximate::APTgauss:    _wn = _createApprox(2,aWnSize); break;
        case PPPApproximate::APTpolin:    _wn = _createApprox(3,aWnSize); break;
        case PPPApproximate::APTbspline:  _wn = _createApprox(4,aWnSize); break;
        case PPPApproximate::APTcolecole: _wn = _createApprox(5,aWnSize); break;
        case PPPApproximate::APTtwogauss: _wn = _createApprox(7,aWnSize); break;
        default:  onError(PPPDISPERSIONMODEL_ERRTYPE+string("prepare()"));
        }
      switch(aAtnType) {
        case PPPApproximate::APTvel:      _atn = _createApprox(1,aAtnSize); break;
        case PPPApproximate::APTgauss:    _atn = _createApprox(2,aAtnSize); break;
        case PPPApproximate::APTpolin:    _atn = _createApprox(3,aAtnSize); break;
        case PPPApproximate::APTbspline:  _atn = _createApprox(4,aAtnSize); break;
        case PPPApproximate::APTcolecole: _atn = _createApprox(6,aAtnSize); break;
        case PPPApproximate::APTtwogauss: _atn = _createApprox(8,aAtnSize); break;
        default:  onError(PPPDISPERSIONMODEL_ERRTYPE+string("prepare()"));
        }
      getWn().getParams().assign(aWn);
      if(!isCausal())
        getAtn().getParams().assign(aAtn);
      else
        getAtn().link(_wn);
      };

    void assign(PPPDispersionModel &aSource) {
      prepare(aSource.getWn().getType(),aSource.getWn().getParams(),
              aSource.getAtn().getType(),aSource.getAtn().getParams());
      setFrequencyRange(aSource.getFrMin(), aSource.getFrMax());
      };     

    // Approximation parameters
    inline unsigned params(void) {
      return getWn().size() + getAtn().size();
      };
    void getParams(PPPVectorContainer<double> &aDest) {
      aDest.realloc(params());
      unsigned aSize = getWn().size();
      for(unsigned i=0; i<aDest.size(); i++)
        aDest[i] = (i<aSize)? getWn()[i] : getAtn()[i-aSize];
      };
    void setParams(PPPVectorContainer<double> &aSource) {
      for(unsigned i=0; i<getWn().size(); i++) getWn()[i] = aSource[i];
      for(unsigned i=0; i<getAtn().size(); i++) getAtn()[i] = aSource[getWn().size()+i];
      };

    // Frequency range
    void setFrequencyRange(double aMin, double aMax) {
      _freqmin = aMin;
      _freqmax = aMax;
      if(isAnalytical())
        {
        getWn().setRange(aMin, aMax);
        getAtn().setRange(aMin, aMax);
        }
      };
    inline double getFrMin(void) const { return _freqmin; };
    inline double getFrMax(void) const { return _freqmax; };
    inline PPPAxis & getFreq(void) { return _data.getAxis(); };

    // wavenumber
    inline PPPApproximate & getWn(void) {
      if(!isWn()) onError(PPPDISPERSIONMODEL_ERRDEF);
      return (*_wn);
      }
    inline double Phi(double f)     { return (isAnalytical())? getWn().Func(f) : _data.Get(f,iWn,_intrp); };
    inline double Phidiv(double f)  { return (isAnalytical())? getWn().FuncDivX(f) : _data.Get(f,iWndiv,_intrp); };

    // attenuation
    inline PPPApproximate & getAtn(void) {
      if(!isAtn()) onError(PPPDISPERSIONMODEL_ERRDEF);
      return (*_atn);
      }
    inline double Atn(double f)     { return (isAnalytical())? getAtn().Func(f) : _data.Get(f,iAtn,_intrp); };
    inline double Atndiv(double f)  { return (isAnalytical())? getAtn().FuncDivX(f) : _data.Get(f,iAtndiv,_intrp); };
    inline bool isCausal(void)      { return _isCausal; }

    // velocities
    inline double VSp(double f)     { return f/fabs(Phi(f)); };
    inline double VSg(double f)     { return 1.0/fabs(Phidiv(f)); };

    // signal container
    void evalContainer(PPPSignalContainer<double> &aData, PPPAxis &aAxis) {
      aData.prepare(aAxis.size(), 6, aAxis, PPPDISPERSIONMODEL_NAME);
      for(unsigned i=0; i<aData.points(); i++)
        {
        aData(i,iWn) = Phi(aAxis[i]);
        aData(i,iWndiv) = Phidiv(aAxis[i]);
        aData(i,iCp) = VSp(aAxis[i]);
        aData(i,iCg) = VSg(aAxis[i]);
        aData(i,iAtn) = Atn(aAxis[i]);
        aData(i,iAtndiv) = Atndiv(aAxis[i]);
        }
      };  

    void evalModel(PPPAxis &aAxis) {
      setFrequencyRange(aAxis.getMin(), aAxis.getMax());
      evalContainer(_data, aAxis);
      };

    void write(const string &aName) {
      _data.write(aName);
      };

    void fwrite(FILE *stream) {
      fwrite_streaminfo(stream, getObjectVer(), sizeof(double));
      unsigned aT;
      aT = isAnalytical();  std :: fwrite((void*)&aT, sizeof(unsigned), 1, stream);
      if(isAnalytical())
        {
        aT = getWn().getType();   std :: fwrite((void*)&aT, sizeof(unsigned), 1, stream);
        aT = getAtn().getType();  std :: fwrite((void*)&aT, sizeof(unsigned), 1, stream);
        getWn().getParams().fwrite(stream);
        getAtn().getParams().fwrite(stream);
        }
      std :: fwrite((void*)&_freqmin, sizeof(double), 1, stream);
      std :: fwrite((void*)&_freqmax, sizeof(double), 1, stream);
      _data.fwrite(stream);
      PPPBaseObject :: fwrite(stream);
      };

    void fread(FILE *stream) {
      fread_streaminfo(stream, getObjectVer(), sizeof(double));
      unsigned aT;
      std :: fread((void*)&aT, sizeof(unsigned), 1, stream);
      setAnalytical(aT);
      if(isAnalytical())
        {
        unsigned aT1,aT2;
        PPPVectorContainer<double> aV1,aV2;
        std :: fread((void*)&aT1, sizeof(unsigned), 1, stream);
        std :: fread((void*)&aT2, sizeof(unsigned), 1, stream);
        aV1.fread(stream);
        aV2.fread(stream);
        prepare((PPPApproximate::ApprType)aT1,aV1,(PPPApproximate::ApprType)aT2,aV2);
        }
      double aF1,aF2;
      std :: fread((void*)&aF1, sizeof(double), 1, stream);
      std :: fread((void*)&aF2, sizeof(double), 1, stream);
      setFrequencyRange(aF1,aF2);
      _data.fread(stream);
      PPPBaseObject :: fread(stream);
      };

  private:

    PPPApproximate *_createApprox(unsigned const aWnType, unsigned const aSize) {
      PPPApproximate *appr;
      switch(aWnType)
        {
        case 1: appr = new PPPApproximateVel(aSize);         break;
        case 2: appr = new PPPApproximateGauss(aSize);       break;
        case 3: appr = new PPPApproximatePolinom(aSize);     break;
        case 4: appr = new PPPApproximateBspline(aSize);     break;
        case 5: appr = new PPPApproximateColeColePhi(aSize); break;
        case 6: appr = new PPPApproximateColeColeAtn(aSize); break;
        case 7: appr = new PPPApproximateTwoGaussPhi(aSize); break;
        case 8: appr = new PPPApproximateTwoGaussAtn(aSize); break;
        default: onError(PPPDISPERSIONMODEL_ERRTYPE+string("createApprox"));
        }
      if(appr == NULL) onError(MEM_ERRALLOC+string("createWavelet"));
      return appr;
      };

    inline bool isWn(void)  { return (_wn != NULL); };
    inline bool isAtn(void) { return (_atn != NULL); };

  }; // end of object


#endif
