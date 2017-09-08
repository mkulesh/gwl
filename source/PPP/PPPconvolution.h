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

#ifndef _PPPCONVOLUTOR
#define _PPPCONVOLUTOR

#define PPPCONVOLUTORSIG_NAME   "signal of convolution"
#define PPPCONVOLUTORSIG_ERR1   "signal of convolution is not initialized in procedure: "
#define PPPCONVOLUTORSIG_ERR2   "size of signals is invalid in procedure: "
#define PPPCONVOLUTOR_NAME      "convolution"
#define PPPCONVOLUTOR_CONV      "convolution plan"
#define PPPCONVOLUTOR_ERR1      "not both limits must be unknwon"
#define PPPCONVOLUTOR_ERR2      "empty convoltuion"


/************************************************************************
 * PPPConvolutorSig
 * This class implement a spetial vector container and its Fourier transform
 * used in the PPPConvolutor object defined below
 ***********************************************************************/
template <class AType> class PPPConvolutorSig : public PPPBaseTemplate<AType>
  {
  private:
    typedef enum {DTnone, DTforward, DTbackward} DirectionType;

    DirectionType                  _type;
    int                            _n0;            // the actual lower index of sig
    int                            _n1;            // the actual upper index of sig
    PPPVectorContainer<AType>      _sig;           // full signal
    PPPVectorContainer<PPPcomplex> _fousig;        // Fourier transform
    PPPVectorContainer<AType>      _sigview;       // a part of the signal
    PPPFft<AType>                  _psigfsig;      // Fourier transformation object

  public:

    PPPConvolutorSig(int an0, int an1): _type(DTnone), _n0(an0), _n1(an1) {
      PPPBaseObject::setObjectName(PPPCONVOLUTORSIG_NAME);
      };

    void clear() {
      if(_type == DTforward) _sigview.unlink();
      _psigfsig.clear();
      };

    inline PPPVectorContainer<AType> & getSignal() { return _sig; };
    inline int getFirstIndex() const { return _n0; };
    inline int getLastIndex() const { return _n1; };

    void initforward(int nn) {
      _type = DTforward;
      _sig.resize(nn);
      _fousig.resize(nn);
      _psigfsig.initforward(_fousig, _sig);
      _sigview.link(_sig, 0, _n1 - _n0 + 1);
      };

    void initbackward(int nn) {
      _type = DTbackward;
      _sig.resize(nn);
      _fousig.resize(nn);
      _psigfsig.initbackward(_sig, _fousig);
      };

    inline PPPVectorContainer<AType> & getSignalView() {
      return _sigview;
      };

    template <class AType1>
    void setSignalView(PPPVectorContainer<AType1> &vec, double aCoeff = 1.0) {
      if(vec.size() != _sigview.size())
        PPPBaseObject::onError(PPPCONVOLUTORSIG_ERR2+string("setSignal()"));
      for(unsigned i=0; i<vec.size(); i++) _sigview[i] = aCoeff*vec[i];
      finish();
      };

    void finish() {
      if(_n0 == 0 && _n1 == 0)
        PPPBaseObject::onError(PPPCONVOLUTORSIG_ERR1+string("finish()"));
      if( _sig.size() > (unsigned)(_n1-_n0+1)) // 0 padding if necessary
        {
        PPPVectorContainer<AType> tmp;
        tmp.link(_sig, _n1-_n0+1);
        tmp.assign( static_cast <AType> (0.) );
        tmp.unlink();
        }
      _psigfsig.execute();
      };

    template <class AType1, class AType2>
    void doConvolution(PPPConvolutorSig<AType1> &s1, PPPConvolutorSig<AType2> &s2) {
      for(unsigned i=0; i<_fousig.size(); ++i )
        _fousig[i] =  s1._fousig[i] * s2._fousig[i];
      _psigfsig.execute();
      };

  };  // end of object


/************************************************************************
 * PPPConvolutor
 * This class allows efficient convolutions
 ***********************************************************************/
template <class AType> class PPPConvolutor : public PPPBaseTemplate<AType>
  {
  private:

    vector< PPPConvolutorSig<AType>* >      _sig;  // the first signal
    vector< PPPConvolutorSig<AType>* >      _sag;  // the second signal
    vector< PPPConvolutorSig<PPPcomplex>* > _sug;  // the result
    PPPMatrixContainer<unsigned>            _convplan;

  public:

    static const bool UNKNOWN_LIMIT = true;
    static const bool KNOWN_LIMIT = false;
    static const bool ZERO_PADDING = true;

    PPPConvolutor() {
      PPPBaseObject::setObjectName(PPPCONVOLUTOR_NAME);
      };

    /**
    *   to get the convolution of
    *   sig with support [nsig0, nsig1] with
    *   sag with support [nsag0, nsag1]
    *   in the window [nsug0,nsug1]
    */
    PPPConvolutor(
        int & nsig0,  // (in) lower index of first signal (out) actually needed lower index
        int & nsig1,  // (in) upper index of first signal (out) actually needed upper index
        int & nsag0,  // (in) lower index of second signal (out) actually needed lower index
        int & nsag1,  // (in) upper index of second signal (out) actually needed upper index
        int   nsug0,  // (in) lower index of result signal (the actual lower limit is nsig0+nsag0)
        int   nsug1,  // (in) upper index of result signal (the actual upper limit is nsig1+nsag1)
        bool unknownLimitSig = !UNKNOWN_LIMIT, // take all of sig that is needed
        bool unknownLimitSag = !UNKNOWN_LIMIT, // take all of sag that is needed
        unsigned aSigCount = 1,
        unsigned aSagCount = 1
      ) {
      PPPBaseObject::setObjectName(PPPCONVOLUTOR_NAME);
      int nn = _evalIndexies(nsig0,nsig1,nsag0,nsag1,nsug0,nsug1,unknownLimitSig,unknownLimitSag);
      // first signal
      _sig.resize(aSigCount);
      for(unsigned i=0; i<getSigSize(); i++)
        {
        if ((_sig[i] = new PPPConvolutorSig<AType>(nsig0,nsig1)) == NULL) 
          PPPBaseObject::onError(MEM_ERRALLOC+string("PPPConvolutor"));
        _sig[i]->initforward(nn);
        }
      // second signal
      _sag.resize(aSagCount);
      for(unsigned i=0; i<getSagSize(); i++)
        {
        if ((_sag[i] = new PPPConvolutorSig<AType>(nsag0,nsag1)) == NULL)
          PPPBaseObject::onError(MEM_ERRALLOC+string("PPPConvolutor"));
        _sag[i]->initforward(nn);
        }
      // convolution plan and result
      _convplan.resize(getSigSize()*getSagSize(),2);
      _convplan.setObjectName(PPPCONVOLUTOR_CONV);
      unsigned k=0;
      for(unsigned i=0; i<getSigSize(); i++)
        for(unsigned j=0; j<getSagSize(); j++)
          {
          _convplan(k,0) = i;
          _convplan(k,1) = j;
          k++;
          }
      _sug.resize(_convplan.rows());
      for(unsigned i=0; i<getSugSize(); i++)
        {
        if ((_sug[i] = new PPPConvolutorSig<AType>(nsug0,nsug1)) == NULL)
          PPPBaseObject::onError(MEM_ERRALLOC+string("PPPConvolutor"));
        _sug[i]->initbackward(nn);
        }
      };

    ~PPPConvolutor() {
      // delete first signal
      for(unsigned i=0; i<getSigSize(); i++)
        {
        _sig[i]->clear();
        delete _sig[i];
        }
      _sig.clear();
      // delete second signal
      for(unsigned i=0; i<getSagSize(); i++)
        {
        _sag[i]->clear();
        delete _sag[i];
        }
      _sag.clear();
      // delete result
      for(unsigned i=0; i<getSugSize(); i++)
        {
        _sug[i]->clear();
        delete _sug[i];
        }
      _sug.clear();
      };

    // this part should be filled with sig (nsig0), sig (nsig0+1)... , sig (nsig1).
    inline unsigned getSigSize() const {
      return _sig.size();
      };

    inline PPPConvolutorSig<AType> &getSig(unsigned aInd=0) {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aInd, getSigSize());
      #endif
      return (*_sig[aInd]);
      };

    inline PPPVectorContainer<AType> &getSigView(unsigned aInd=0) {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aInd, getSigSize());
      #endif
      return _sig[aInd]->getSignalView();
      };

    void finishSig(unsigned aInd=0) {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aInd, getSigSize());
      #endif
      _sig[aInd]->finish();
      };

    // this part should be filled with sag (nsag0), sag (nsag0+1), ... , sag (nsag1).
    inline unsigned getSagSize() const {
      return _sag.size();
      };

    inline PPPConvolutorSig<AType> &getSag(unsigned aInd=0) {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aInd, getSagSize());
      #endif
      return (*_sag[aInd]);
      };

    inline PPPVectorContainer<AType> &getSagView(unsigned aInd=0) {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aInd, getSagSize());
      #endif
      return _sag[aInd]->getSignalView();
      };

    void finishSag(unsigned aInd=0) {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aInd, getSagSize());
      #endif
      _sag[aInd]->finish();
      };

    //  resulting signal
    inline unsigned getSugSize() const {
      return _sug.size();
      };

    void doConvolution () {
      for(unsigned i=0; i<getSugSize(); i++)
        _sug[i]->doConvolution(*_sig[_convplan(i,0)], *_sag[_convplan(i,1)]);
      };

    void fillSug(PPPVectorContainer <PPPcomplex> &aDest, bool aZero = false, unsigned aInd=0) {
      unsigned aInd1 = _convplan(aInd,0);
      unsigned aInd2 = _convplan(aInd,1);
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aInd1, getSigSize());
      _checkindex(aInd2, getSagSize());
      _checkindex(aInd,  getSugSize());
      #endif
      int _s0 = max( _sig[aInd1]->getFirstIndex() + _sag[aInd2]->getFirstIndex(), _sug[aInd]->getFirstIndex());
      int _s1 = min( _sig[aInd1]->getLastIndex() + _sag[aInd2]->getLastIndex(), _sug[aInd]->getLastIndex());
      if(aZero)
        {
        for( int i=0; i<_s0-_sug[aInd]->getFirstIndex(); ++i )
          aDest[i] = 0;
        for( int i=_s1-_sug[aInd]->getFirstIndex()+1; i < _sug[aInd]->getFirstIndex() + _sug[aInd]->getLastIndex(); ++i )
          aDest[i] = 0;
        }
      PPPVectorContainer <PPPcomplex> tmp;
      tmp.link(aDest, _s0-_sug[aInd]->getFirstIndex(), _s1-_s0+1 );
      tmp.assign(_sug[aInd]->getSignal(),_s0-(_sig[aInd1]->getFirstIndex()+_sag[aInd2]->getFirstIndex()),_s1-_s0+1);
      tmp.unlink();
      };

  private:

    int _evalIndexies(
        int & nsig0,  // (in) lower index of first signal (out) actually needed lower index
        int & nsig1,  // (in) upper index of first signal (out) actually needed upper index
        int & nsag0,  // (in) lower index of second signal (out) actually needed lower index
        int & nsag1,  // (in) upper index of second signal (out) actually needed upper index
        int & nsug0,  // (in) lower index of result signal (the actual lower limit is nsig0+nsag0)
        int & nsug1,  // (in) upper index of result signal (the actual upper limit is nsig1+nsag1)
        bool unknownLimitSig, // take all of sig that is needed
        bool unknownLimitSag  // take all of sag that is needed
      ) {
      int nnsig0 = nsig0;
      int nnsig1 = nsig1;
      int nnsag0 = nsag0;
      int nnsag1 = nsag1;
      if ( unknownLimitSig && unknownLimitSag )
        PPPBaseObject::onError(PPPCONVOLUTOR_ERR1);
      else if ( unknownLimitSag ) // first guess for limits
        {
        nsag0 = nsug0 - nsig1;
        nsag1 = nsug1 - nsig0;
        }
      else if ( unknownLimitSig )
        {
        nsig0 = nsug0 - nsag1;
        nsig1 = nsug1 - nsag0;
        }
      if (( nsig0 + nsag0 > nsug1 ) && ( nsig1 + nsag1 < nsug0 ))
        PPPBaseObject::onError(PPPCONVOLUTOR_ERR2);
      if ( nsig0 + nsag1 < nsug0 ) nnsig0 = nsug0 - nsag1; // lower bounds of sig
      if ( nsag0 + nsig1 < nsug0 ) nnsag0 = nsug0-nsig1; // lower bounds of sag
      if ( nsig1 + nsag0 > nsug1 ) nnsig1 = nsug1 - nsag0; // upper bounds of sig
      if ( nsag1 + nsig0 > nsug1 ) nnsag1 = nsug1 - nsig0; // upper bounds of sag
      nsig0 = nnsig0;
      nsig1 = nnsig1;
      nsag0 = nnsag0;
      nsag1 = nnsag1;
      if (( nsig0 + nsag0 > nsug1 ) | ( nsig1 + nsag1 < nsug0 ) | ( nsig0 > nsig1 ) | ( nsag0 > nsag1 ))
        PPPBaseObject::onError(PPPCONVOLUTOR_ERR2);
      int nn = nsig1 - nsig0 + nsag1 - nsag0 + 1;
      nn = (int) pow(2., ceil(log((double)nn)/log(2.)) ) ; // largest power of 2
      return nn;
      };

    void _checkindex(const unsigned aInd, const unsigned aSize) const {
      if(aInd >= aSize)
        {
        strstream str;
        str << MEM_ERRINDEX << ", size = " << aSize << ", index = " << aInd << ends;
        PPPBaseObject::onError(str.str());
        }
      };

  }; // end of object

#endif

