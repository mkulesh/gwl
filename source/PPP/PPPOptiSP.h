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

#ifndef _PPPOPTISPECTRUM
#define _PPPOPTISPECTRUM

#define PPPOPTISPECTRUM_NAME       "optimization of spectrum"
#define PPPOPTISPECTRUM_DER          "derivative matrix"
#define PPPOPTISPECTRUM_TYPESP       "optimization for: "
#define PPPOPTISPECTRUM_SPEC1NAME    "spectrum of first signal"
#define PPPOPTISPECTRUM_SPEC2NAME    "spectrum of second signal"

/************************************************************************
 * PPPOptiSpectrum
 * @article{Holschneider2005GJI,
 *   author = {M. Holschneider and M. S. Diallo and M. Kulesh and M. Ohrnberger and E. L\"uck and F. Scherbaum},
 *   title = {Characterization of dispersive surface waves using continuous wavelet transforms},
 *   journal = {Geophysical Journal International},
 *   year = {2005},
 *   volume = {163},
 *   number = {2},
 *   pages = {463-478},
 *   doi = {10.1111/j.1365-246X.2005.02787.x},
 *   urllink = {http://www.blackwell-synergy.com/doi/abs/10.1111/j.1365-246X.2005.02787.x},
 *   eprint = {http://www.blackwell-synergy.com/doi/pdf/10.1111/j.1365-246X.2005.02787.x}
 * }
 ************************************************************************/
template<class ATypeOpti> class PPPOptiSpectrum : public ATypeOpti
  {
  private:
    PPPSpectrContainer<PPPcomplex>  _wg1, _wg2, _wgn;
    PPPPropagatorDiss               _propag;
    PPPSpectrParams*                _wtpar;
    PPPPropagatorDiss::PropType     _propType;
    unsigned                        _cmplType;
    PPPMatrixContainer<double>      _der;

  public:

    PPPOptiSpectrum() : _cmplType(TCabs), _propType(PPPPropagatorDiss::STwav) {
      PPPBaseObject :: setObjectName(PPPOPTISPECTRUM_NAME);
      _der.setObjectName(PPPOPTISPECTRUM_DER);
      };

    void prepare(PPPSpectrContainer<PPPcomplex> &aSource, PPPDispersionModel &aModel, double aDist,
      PPPPropagatorDiss::PropType aPropType, PPPSpectrParams *aWtpar, unsigned aCmplType, unsigned aChan1=0, unsigned aChan2=1) {
      // wavelet spectrum
      _wg1.prepare(aSource.voices(), aSource.points(), 1, aSource.getTime(), aSource.getFreq(), PPPOPTISPECTRUM_SPEC1NAME);
      _wg1.getChannel(0).assign(aSource.getChannel(aChan1));
      _wg2.prepare(aSource.voices(), aSource.points(), 1, aSource.getTime(), aSource.getFreq(), PPPOPTISPECTRUM_SPEC2NAME);
      _wg2.getChannel(0).assign(aSource.getChannel(aChan2));
      // wavelet propagatior
      getModel().assign(aModel);
      _propag.setDistance(aDist);
      _propag.setShowMessage(false);
      _propType = aPropType;
      // parameters
      _wtpar = aWtpar;
      _cmplType = aCmplType;
      if(!ATypeOpti :: isInitialized())
        {
          ATypeOpti::resize(_wg1.voices()*_wg1.points(), getModel().params());
        _der.resize(ATypeOpti :: params(),ATypeOpti :: points());
        }
      ATypeOpti::evalFixed((getWn().getType() == PPPApproximate::APTpolin), (_cmplType == TCarg), getWn().size());
      };

    void optimization() {
      getModel().getParams(ATypeOpti :: optpar);
      ATypeOpti :: Minimize();
      };

    double getcostfunc(unsigned Ai) {
      register unsigned _indi, _indj;
      _indi = Ai/_wg2.points();
      _indj = Ai - _wg2.points()*_indi;
      return getComponent(_wg2(_indi,_indj))-getComponent(_wgn(_indi,_indj));
      };

    double getder(unsigned Ai,unsigned Aj) {
      return _der(Ai,Aj);
      };

    void calcfunc(PPPVectorContainer<double> &aPar) {
      register unsigned _indi, _indj;
      double df,dx,der1,sign,f;
      double dX = _propag.getDistance();
      unsigned aAlg = 1;
      // Selbst Funktion
      getModel().setParams(aPar);
      if(_propType == PPPPropagatorDiss::STwav) // optimization mit Wavelet-propag
        {
        _propag.evalWaveletPropag(PPPPropagatorDiss::STwav,_wgn,_wg1,*_wtpar,aAlg);
        sign = 1.0;
        }
      if(_propType == PPPPropagatorDiss::STwavcross) // optimization mit Cross-correation Wavelet-propag
        {
        _propag.evalWaveletPropag(PPPPropagatorDiss::STwavcross,_wgn,_wg1,*_wtpar,aAlg);
        sign = -1.0;
        }
      // Ableitung von der Funktion
      for(unsigned j=0;j<ATypeOpti :: points();j++)
        {
        _indi = j/_wgn.points();
        _indj = j - _wgn.points()*_indi;
        df = (_indj>0)? (getComponent(_wgn(_indi,_indj)) - getComponent(_wgn(_indi,_indj-1))) : 0.0;
        dx = (_indj>0)? (_wgn.getTime(_indj) - _wgn.getTime(_indj-1)) : 1.0;
        f = _wgn.getFreq(_indi);
        for(unsigned m=0; m<ATypeOpti :: params(); m++)
          {
          der1 = (m<getWn().size())?
            -sign*(df/dx)*dX*getWn().FuncDivXPar(f,m) :
            -dX*getAtn().FuncDivPar(f,m-getWn().size())*getComponent(_wgn(_indi,_indj));
          _der(m,j) = der1;
          }
        }
      return;
      };

    const char *getInfo(void) {
      strstream aDest;
      if(PPPBaseObject :: isNotation()) aDest << PPPBaseObject :: getNotation() << endl;
      if(_cmplType == TCabs)
        {
        PPPcomplexAbs trans;
        aDest << "  "<< PPPOPTISPECTRUM_TYPESP << trans.getComponentName();
        }
      if(_cmplType == TCarg)
        {
        PPPcomplexArg trans;
        aDest << "  "<< PPPOPTISPECTRUM_TYPESP << trans.getComponentName();
        }
      aDest << ends;
      PPPBaseObject :: onNotation(aDest.str());
      return PPPBaseObject :: getNotation();
      };

    inline PPPDispersionModel & getModel(void) {
      return _propag.getModel();
      };

    inline PPPSpectrContainer<PPPcomplex> & getResult(void) {
      return _wgn;
      };

  private:

    inline PPPApproximate & getWn(void) {
      return _propag.getModel().getWn();
      };

    inline PPPApproximate & getAtn(void) {
      return _propag.getModel().getAtn();
      };

    double getComponent(PPPcomplex const &aVal)  {
      switch(_cmplType) {
        case TCabs:  return abs(aVal);
        case TCarg:  return arg(aVal);
        default: PPPBaseObject :: onError(PPPBASETEMPLATE_ERR+string("getComponent"));
        }
      return 0;
      };

  }; // end of object




#endif
