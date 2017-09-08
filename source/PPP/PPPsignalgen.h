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
#ifndef _PPPSIGNALGEN
#define _PPPSIGNALGEN

#define PPPSIGNALGEN_NAME        "synthetic signal"
#define PPPSIGNALGEN_ZERO        "zero-padded function"
#define PPPSIGNALGEN_DELTA       "delta function"
#define PPPSIGNALGEN_HARMON      "harmonic function"
#define PPPSIGNALGEN_HARMONROT   "rotating harmonic function"
#define PPPSIGNALGEN_HARMONPHASE "harmonic function with harmonic phase"
#define PPPSIGNALGEN_MORLET      "morlet wavelet"
#define PPPSIGNALGEN_CAUCHY      "cauchy wavelet"
#define PPPSIGNALGEN_SHANON      "shanon wavelet"
#define PPPSIGNALGEN_MODULATION  "windowed modulation of signal"
#define PPPSIGNALGEN_RANDNOISE   "add ramdom noise to signal"
#define PPPSIGNALGEN_TESTDISSIP  "synthetic dispersive signal"
#define PPPSIGNALGEN_ERRTYPE     "type of signal is incorrect in procedure: "
#define PPPSIGNALGEN_ERRPAR      "number of parametes is incorrect in procedure: "

/************************************************************************
 * PPPSignalGen
 ************************************************************************/
template<class AType> class PPPSignalGen : public PPPBaseTemplate<AType>
  {
  private:
    unsigned                    _points;
    double                      _sample;
    PPPVectorContainer<double>  _par;

  public:

    PPPSignalGen(void): _points(1024), _sample(100.0)
      {
      PPPBaseObject :: setObjectName(PPPSIGNALGEN_NAME);
      _par.resize(6);
      };

    inline unsigned points(void) { return _points; };
    inline void setPoints(unsigned apoints) { _points = apoints; };

    inline double getSamplingFreq(void) const  { return _sample; };
    inline void setSamplingFreq(double asample) { _sample = asample; };

    inline PPPVectorContainer<double> & getParams(void) { return _par; };

    void genZero(PPPSignalContainer<AType> &aDest, unsigned aNum = 0) {
      aDest.getChannel(aNum).assign(aDest.nullValue());
      // Notation
      strstream str;
      str << PPPSIGNALGEN_ZERO << ends;
      PPPBaseObject :: onNotation(str.str());
      return;
      };

    void genDelta(PPPSignalContainer<AType> &aDest, unsigned aNum = 0) {
      if(_par.size() != 2)
        PPPBaseObject :: onError(PPPSIGNALGEN_ERRPAR+string("genDelta()"));
      PPPAxis x = aDest.getAxis();
      unsigned pos = x.locateFloor(_par[0]);
      for(unsigned i=0; i<aDest.points(); i++) aDest(i,aNum) = (AType)0.0;
      aDest(pos,aNum) = (AType)_par[1];
      // Notation
      strstream str;
      str << PPPSIGNALGEN_DELTA << ": " << endl << "  " << "Fre(t) = delta(t-"<<_par[0]<<")";
      if(PPPBaseTemplate<AType> :: isComplex()) str << endl << "  " << "Fim(t) = 0";
      str << ends;
      PPPBaseObject :: onNotation(str.str());
      return;
      };

    void genHarmon(PPPSignalContainer<AType> &aDest, unsigned aNum = 0) {
      if(_par.size() != 7)
        PPPBaseObject :: onError(PPPSIGNALGEN_ERRPAR+string("genHarmon()"));
      PPPAxis x = aDest.getAxis();
      for(unsigned i=0; i<aDest.points(); i++)
        {
        double e = (double)i/(double)(aDest.points()-1);
        PPPcomplex f = PPPcomplex(
          (_par[0]+_par[1]*e)*cos(2.0*M_PI*_par[2]*x[i]),
          (_par[3]+_par[4]*e)*sin(2.0*M_PI*_par[5]*x[i] + _par[6])
          );
        PPPBaseTemplate<AType>::cmplConvert(aDest(i,aNum), f);
        }
      // Notation
      strstream str;
      str << PPPSIGNALGEN_HARMON << ": " << endl
          << "  Sre(t) = (" << _par[0] << showpos << _par[1] << noshowpos << "e(t))*cos(2*Pi*" << _par[2] << "*t)" << endl;
      if(PPPBaseTemplate<AType> :: isComplex())
      str << "  Sim(t) = (" << _par[3] << showpos << _par[4] << noshowpos << "e(t))*sin(2*Pi*" << _par[5] << "*t + " << _par[6] << ")" << endl;
      str << "  e(t) = 0...1" << ends;
      PPPBaseObject :: onNotation(str.str());
      return;
      };

    void genHarmonPhase(PPPSignalContainer<AType> &aDest, unsigned aNum = 0) {
      if(_par.size() != 4)
        PPPBaseObject :: onError(PPPSIGNALGEN_ERRPAR+string("genHarmonPhase()"));
      PPPAxis x = aDest.getAxis();
      for(unsigned i=0; i<aDest.points(); i++)
        {
        PPPcomplex f = PPPcomplex(
          _par[0]*sin(2.0*M_PI*_par[2]*x[i]+sin(_par[3]*x[i])),
          _par[1]*sin(2.0*M_PI*_par[2]*x[i])
          );
        PPPBaseTemplate<AType>::cmplConvert(aDest(i,aNum), f);
        }
      // Notation
      strstream str;
      str << PPPSIGNALGEN_HARMONPHASE << ": " << endl
          << "  Sre(t) = " << _par[0] << "*sin(2*Pi*" << _par[2] << "*t+sin(" << _par[3] << "*t))" << endl;
      if(PPPBaseTemplate<AType> :: isComplex())
      str << "  Sim(t) = " << _par[1] << "*sin(2*Pi*" << _par[2] << "*t)" << endl;
      str << ends;
      PPPBaseObject :: onNotation(str.str());
      return;
      };

    void genHarmonRot(PPPSignalContainer<AType> &aDest, unsigned aNum = 0) {
      if(_par.size() != 7)
        PPPBaseObject :: onError(PPPSIGNALGEN_ERRPAR+string("genHarmonRot()"));
      PPPAxis x = aDest.getAxis();
      for(unsigned i=0; i<aDest.points(); i++)
        {
        double e = (double)i/(double)(aDest.points()-1);
        double xt = (_par[0]+_par[1]*e)*cos(2.0*M_PI*_par[4]*x[i]);
        double yt = (_par[2]+_par[3]*e)*sin(2.0*M_PI*_par[5]*x[i]);
        PPPcomplex f = PPPcomplex(
          xt*cos(2.0*M_PI*_par[6]*x[i]) + yt*sin(2.0*M_PI*_par[6]*x[i]),
          -xt*sin(2.0*M_PI*_par[6]*x[i]) + yt*cos(2.0*M_PI*_par[6]*x[i])
          );
        PPPBaseTemplate<AType>::cmplConvert(aDest(i,aNum), f);
        }
      // Notation
      strstream str;
      str << PPPSIGNALGEN_HARMONROT << ": " << endl
          << "  S1(t) = (" << _par[0] << showpos << _par[1] << noshowpos << "e(t))*cos(2*Pi*" << _par[4] << "*t)" << endl
          << "  S2(t) = (" << _par[2] << showpos << _par[3] << noshowpos << "e(t))*sin(2*Pi*" << _par[5] << "*t)" << endl
          << "  Sre(t) = S1(t)*cos(2*Pi*" << _par[6] << "*x)+S2(t)*sin(2*Pi*" << _par[6] << "*t)" << endl;
      if(PPPBaseTemplate<AType> :: isComplex())
      str << "  Sim(t) =-S1(t)*sin(2*Pi*" << _par[6] << "*x)+S2(t)*cos(2*Pi*" << _par[6] << "*t)" << endl;
      str << "  e(t) = 0...1" << ends;
      PPPBaseObject :: onNotation(str.str());
      return;
      };

    void genRickerDissip(PPPSignalContainer<AType> &aDest, unsigned aNum = 0) {
      if(_par.size() != 6)
        PPPBaseObject :: onError(PPPSIGNALGEN_ERRPAR+string("genRickerDissip()"));
      // Frequency axis
      double aSample = aDest.getAxis().getSamplingFreq();
      PPPAxis aFreq(aDest.points(), 0, aSample, PPPAxis::ATlin, PPPSPECTRPARAMS_FREQ);
      // propagation parameters
      double Vr = _par[0], f0 = _par[1], dist = _par[2];
      PPPApproximateVel aAppr(3);
      aAppr.setRange(aFreq.getMin(), aFreq.getMax());
      aAppr.getParams()[0] = _par[3];
      aAppr.getParams()[1] = _par[4];
      aAppr.getParams()[2] = _par[5];
      // Fourier spectrum
      PPPVectorContainer<double> aFour(aDest.points());
      for(unsigned i=0; i<aFour.size(); i++)
        {
        double omega = 2.0*M_PI*aFreq[i];
        double om0 = 2.0*M_PI*f0;
        aFour[i] = (2.0*omega*omega/(Vr*Vr*sqrt(M_PI)))*exp(-(omega-om0)*(omega-om0)/(Vr*Vr));
        }
      double aFourMax = aFour.getMaxValue()/sqrt((double)aFour.size());
      // Propagated signal
      PPPSignalContainer<PPPcomplex> aRick;
      aRick.prepare(aFour.size(),1,aFreq,"Ricker wavelet");
      PPPcomplex _i(0.0, 1.0);
      for(unsigned i=0; i<aRick.points(); i++)
        {
        double omega = aRick.getAxis(i);
        aRick(i) = exp(-2.0*M_PI*_i*dist*aAppr.Func(omega))*aFour[i]/aFourMax;
        }
      PPPTransFour<AType> TTR;
      TTR.IFT(aDest,aRick,0.0,(aRick.points()-1)/aSample);
      // Notation
      strstream str;
      str << PPPSIGNALGEN_TESTDISSIP << ":" << endl;
      str << "  YR2(f) = exp(-dX*Atn(f))*exp(-2.0*Pi*i*dX*Func(f))*YR1(2*Pi*f)" << endl;
      str << "  YR1(om) = 2*om^2*exp(-om^2/Vr^2)/(Vr^2*sqrt(Pi)), Vr="<<Vr<<", f0="<<f0 << endl;
      str << "  " << aAppr.getInfo() << ends;
      PPPBaseObject :: onNotation(str.str());
      };

    void evalWindowedModulation(PPPSignalContainer<AType> &aDest, unsigned aNum = 0) {
      if(_par.size() != 3)
        PPPBaseObject :: onError(PPPSIGNALGEN_ERRPAR+string("evalWindowedModulation()"));
      for(unsigned i=0; i<aDest.points(); i++)
        {
        if(aDest.getAxis(i)>_par[0] && aDest.getAxis(i)<_par[1])
          aDest(i,aNum) = aDest(i,aNum)*_par[2];
        else
          aDest(i,aNum) = 0.0;
        }
      // Notation
      strstream str;
      str << PPPSIGNALGEN_MODULATION
          << ": time interval = [" << _par[0] << ", " << _par[1] << "]"
          << ", amplitude factor = " << _par[2] << ends;
      PPPBaseObject :: onNotation(str.str());
      };

    void addRandomNoise(PPPSignalContainer<AType> &aDest, unsigned aNum = 0) {
      if(_par.size() != 1)
        PPPBaseObject :: onError(PPPSIGNALGEN_ERRPAR+string("addRandomNoise()"));
      double m = _par[0];
      for(unsigned i=0; i<aDest.points(); i++)
        {
        PPPcomplex r = aDest(i,aNum) + PPPcomplex(m*(double)rand()/RAND_MAX-m/2.0, m*(double)rand()/RAND_MAX-m/2.0);
        PPPBaseTemplate<AType>::cmplConvert(aDest(i,aNum), r);
        }
      // Notation
      strstream str;
      str << PPPSIGNALGEN_RANDNOISE << ": noise value = [" << -m/2.0 << ", " << m/2.0 << "]" << ends;
      PPPBaseObject :: onNotation(str.str());
      };

  };  // end of object


#endif
