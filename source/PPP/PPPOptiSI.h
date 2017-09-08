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
#ifndef _PPPOPTISIGNAL
#define _PPPOPTISIGNAL

#define PPPOPTISIGNAL_NAME        "optimization of signal"
#define PPPOPTISIGNAL_DER         "derivative matrix"
#define PPPOPTISIGNAL_SIG1NAME    "first signal"
#define PPPOPTISIGNAL_SIG2NAME    "second signal"

/************************************************************************/
/** PPPOptiSignal                                                        */
/************************************************************************/
template<class ATypeOpti> class PPPOptiSignal : public ATypeOpti
  {
  private:
    PPPMatrixContainer<double>      _der;
    PPPPropagatorDiss               _propag;
    PPPTransFour<double>            _trans;
    PPPPropagatorDiss::PropType     _propType;
    PPPSignalContainer<double>      _s1, _s2, _sn;
    PPPSignalContainer<PPPcomplex>  _ys1, _ys2, _ysn;

  public:

    PPPOptiSignal() : _propType(PPPPropagatorDiss::STfour) {
      PPPBaseObject :: setObjectName(PPPOPTISIGNAL_NAME);
      _der.setObjectName(PPPOPTISIGNAL_DER);
      };

    void prepare(PPPSignalContainer<double> &aSource, PPPDispersionModel &aModel,
      double aDist, PPPPropagatorDiss::PropType aPropType, unsigned aChan1=0, unsigned aChan2=1) {
      _s1.prepare(aSource.points(), 1, aSource.getAxis(), PPPOPTISIGNAL_SIG1NAME);
      _s1.getChannel(0).assign(aSource.getChannel(aChan1));
      _s2.prepare(aSource.points(), 1, aSource.getAxis(), PPPOPTISIGNAL_SIG2NAME);
      _s2.getChannel(0).assign(aSource.getChannel(aChan2));
      _trans.setShowMessage(true);
      _trans.FT(_ys1, _s1);
      _trans.FT(_ys2, _s2);
      _trans.setShowMessage(false);
      // propagatior
      getModel().assign(aModel);
      _propag.setDistance(aDist);
      _propag.setShowMessage(false);
      _propType = aPropType;
      // parameters
      if(!ATypeOpti::isInitialized())
        {
         ATypeOpti::resize(_s1.points(), getModel().params());
        _der.resize(ATypeOpti::params(),ATypeOpti::points());
        }
      ATypeOpti::evalFixed((getWn().getType() == PPPApproximate::APTpolin), false, getWn().size());
      };

    void optimization() {
      getModel().getParams(ATypeOpti::optpar);
      ATypeOpti::Minimize();
      };

    double getcostfunc(unsigned Ai) {
      return (_s2(Ai)-_sn(Ai));
      };

    double getder(unsigned Ai,unsigned Aj) {
      return _der(Ai,Aj);
      };

    void calcfunc(PPPVectorContainer<double> &aPar) {
      PPPSignalContainer<double> Snprim;
      PPPTransWavelet<PPPcomplex> TC;
      PPPSignalContainer<PPPcomplex> Splus,Sminus,YSnprim;
      PPPcomplex cmpli(0.0, 2.0*M_PI);
      double f;
      double sign = (_propType == PPPPropagatorDiss::STfour)? 1.0 : -1.0;
      double dX = _propag.getDistance();
      getModel().setParams(aPar);
      _propag.evalFouriePropag(_propType,_ysn,_ys1);
      _trans.IFT(_sn,_ysn,_s1.getAxis().getMin(),_s1.getAxis().getMax());
      if(ATypeOpti::getOptiType() == ATypeOpti::OTleven)
        {
        for(unsigned m=0; m<ATypeOpti::params(); m++)
          {
          TC.FTSeparate(Splus,Sminus,_ysn);
          for(unsigned i=0; i<Splus.points(); i++)
            {
            f = Splus.getAxis()[i];
            if(f == 0.0 || f<getModel().getFrMin() || f>getModel().getFrMax())
              {
              Splus(i) = 0.0;
              Sminus(i) = 0.0;
              continue;
              }
            if(m<getWn().size())
              {
              Splus(i) = -sign*cmpli*dX*getWn().FuncDivPar(f,m)*Splus(i);
              Sminus(i) = sign*cmpli*dX*getWn().FuncDivPar(f,m)*Sminus(i);
              }
            else
              {
              Splus(i) = -dX*getAtn().FuncDivPar(f,m-getWn().size())*Splus(i);
              Sminus(i) = -dX*getAtn().FuncDivPar(f,m-getWn().size())*Sminus(i);
              }
            }
          YSnprim.resize(_ysn.points());
          YSnprim.setAxis(_ysn.getAxis());
          TC.FTCut(YSnprim,Splus,Sminus);
          _trans.IFT(Snprim,YSnprim,_s1.getAxis().getMin(),_s1.getAxis().getMax());
          for(unsigned j=0;j<ATypeOpti::points();j++) _der(m,j) = Snprim(j);
          }
        }
      return;
      };

    const char *getInfo(void) {
      strstream aDest;
      if(PPPBaseObject :: isNotation()) aDest << PPPBaseObject :: getNotation() << ends;
      PPPBaseObject :: onNotation(aDest.str());
      return PPPBaseObject :: getNotation();
      };

    inline PPPDispersionModel & getModel(void) {
      return _propag.getModel();
      };

    inline PPPSignalContainer<double> & getResult(void) {
      return _sn;
      };

  private:

    inline PPPApproximate & getWn(void) {
      return getModel().getWn();
      };

    inline PPPApproximate & getAtn(void) {
      return getModel().getAtn();
      };

  }; // end of object

#endif

