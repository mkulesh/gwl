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
#ifndef _PPPOPTIMULT
#define _PPPOPTIMULT

#define PPPOPTIMULT_NAME            "multi-channel optimization"
#define PPPOPTIMULT_ERRCHANN        "number of channel must be more than one in procedure: "
#define PPPOPTIMULT_ERROUT          "out of range in procedure: "

/************************************************************************/
/** PPPOptiMult                                                         */
/************************************************************************/
template<class ATypeOpti> class PPPOptiMult : public ATypeOpti
  {
  private:
    unsigned                            _indi, _indj, _chan, _points;
    ATypeOpti*                          _opt;
    PPPSignalContainer<double>          _sn;
    PPPSpectrContainer<PPPcomplex>      _wgn;

  public:

    PPPOptiMult() : _indi(0), _indj(0), _chan(0), _points(0) {
      PPPBaseObject :: setObjectName(PPPOPTIMULT_NAME);
      };

    void optimization() {
      getModel().getParams(ATypeOpti :: optpar);
      ATypeOpti :: Minimize();
      };

    void prepare(PPPSignalContainer<double> &aSource, PPPDispersionModel &aModel,
      double aDist, PPPPropagatorDiss::PropType aPropType) {
      _sn.assign(aSource);
      if(_chan != 0) delete _opt;
      _chan = aSource.channels()-1;
      if(_chan == 0)
        PPPBaseObject :: onError(PPPOPTIMULT_ERRCHANN+string("prepare"));
      _opt = new ATypeOpti[_chan];
      for(unsigned i=0; i<_chan; i++)
        {
        _opt[i].prepare(aSource,aModel,(double)(i+1)*aDist,aPropType,0,i+1);
        if(i == 0) _points = _opt[i].points();
        }
      if(!ATypeOpti :: isInitialized())
        ATypeOpti :: resize(_chan*_points, aModel.params());
      ATypeOpti :: evalFixed((getModel().getWn().getType() == PPPApproximate::APTpolin), false, getModel().getWn().size());
      };

    void prepare(PPPSpectrContainer<PPPcomplex> &aSource, PPPDispersionModel &aModel, double aDist,
      PPPPropagatorDiss::PropType aPropType, PPPSpectrParams *aWtpar, unsigned aCmplType) {
      _wgn.assign(aSource);
      if(_chan != 0) delete _opt;
      _chan = aSource.channels()-1;
      if(_chan == 0)
        PPPBaseObject :: onError(PPPOPTIMULT_ERRCHANN+string("prepare"));
      _opt = new ATypeOpti[_chan];
      for(unsigned i=0; i<_chan; i++)
        {
        _opt[i].prepare(aSource,aModel,(double)(i+1)*aDist,aPropType,aWtpar,aCmplType,0,i+1);
        if(i == 0) _points = _opt[i].points();
        }
      if(!ATypeOpti :: isInitialized())
        ATypeOpti :: resize(_chan*_points, aModel.params());
      ATypeOpti :: evalFixed((getModel().getWn().getType() == PPPApproximate::APTpolin), (aCmplType == TCarg), getModel().getWn().size());
      };

    double getcostfunc(unsigned Ai) {
      _indi = Ai/_points;
      if(_indi>=_chan)
        PPPBaseObject :: onError(PPPOPTIMULT_ERROUT+string("getcostfunc"));
      _indj = Ai - _points*_indi;
      return _opt[_indi].getcostfunc(_indj);
      };

    double getder(unsigned Ai,unsigned Aj) {
      _indj = Aj/_points;
      if(_indj>=_chan)
        PPPBaseObject :: onError(PPPOPTIMULT_ERROUT+string("getder"));
      _indi = Aj - _points*_indj;
      return _opt[_indj].getder(Ai,_indi);
      };

    void calcfunc(PPPVectorContainer<double> &aPar) {
      for(unsigned i=0; i<_chan; i++)
        _opt[i].calcfunc(aPar);
      return;
      };

    const char *getInfo(void) {
      strstream aDest;
      if(PPPBaseObject :: isNotation()) aDest << PPPBaseObject :: getNotation() << ends;
      PPPBaseObject :: onNotation(aDest.str());
      return PPPBaseObject :: getNotation();
      };

    inline PPPDispersionModel & getModel(void) {
      if(_chan == 0)
        PPPBaseObject :: onError(PPPOPTIMULT_ERRCHANN+string("getModel"));
      return _opt[0].getModel();
      };

    inline PPPSignalContainer<double> & getResultSig(void) {
      for(unsigned i=0; i<_chan; i++)
        _sn.getChannel(i+1).assign(_opt[i].getResult().getChannel(0));
      return _sn;
      };

    inline PPPSpectrContainer<PPPcomplex> & getResultSpec(void) {
      for(unsigned i=0; i<_chan; i++)
        _wgn.getChannel(i+1).assign(_opt[i].getResult().getChannel(0));
      return _wgn;
      };

  }; // end of object

#endif

