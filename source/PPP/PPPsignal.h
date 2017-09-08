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

#ifndef _PPPSIGNALCONTAINER
#define _PPPSIGNALCONTAINER

#define PPPSIGNALCONTAINER_OBJVER      "SI1.3"
#define PPPSIGNALCONTAINER_NAME        "seismogram"
#define PPPSIGNALCONTAINER_SEC         "seismic section"
#define PPPSIGNALCONTAINER_AXISNAME    "time"
#define PPPSIGNALCONTAINER_CHANNAME    "channel"
#define PPPSIGNALCONTAINER_CHANNOT     "channel's notation:"
#define PPPSIGNALCONTAINER_CHANNUMB    "number of channels"
#define PPPSIGNALCONTAINER_POINT       "number of time points"
#define PPPSIGNALCONTAINER_SAMPL       "sampling frequency"
#define PPPSIGNALCONTAINER_ERRINTERP   "undefined interpolation type in procedure: "

/************************************************************************
 *  PPPSignalContainer
 ***********************************************************************/
template<class AType> class PPPSignalContainer : public PPPBaseTemplate<AType>
  {
  public:

    typedef enum {ITnone, ITlin, ITpolin, ITsplin} InterpType;

  private:
    PPPAxis                             _axis;
    vector< PPPVectorContainer<AType> > _data;
    PPPVectorContainer<AType>           _splineTmp;

  public:

    PPPSignalContainer(void) {
      PPPBaseTemplate<AType>::setObjectVer(PPPSIGNALCONTAINER_OBJVER);
      setNames(PPPSIGNALCONTAINER_NAME,PPPSIGNALCONTAINER_AXISNAME);
      resize(0,0);
      };

    PPPSignalContainer(const unsigned aPoints, const unsigned aChannels=1) {
      PPPBaseTemplate<AType>::setObjectVer(PPPSIGNALCONTAINER_OBJVER);
      setNames(PPPSIGNALCONTAINER_NAME,PPPSIGNALCONTAINER_AXISNAME);
      resize(aPoints,aChannels);
      };

    /**
     *  this is the part, with is depended of type of "_data" property
     */
    inline unsigned points() const {
      return (_data.size() == 0)? 0 : _data[0].size();
      };

    inline unsigned channels() const {
      return _data.size();
      };

    inline AType& operator() (const unsigned aPoint, const unsigned aChannel=0) const {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aChannel);
      #endif
      return (_data[aChannel])[aPoint];
      };

    void realloc(const unsigned aPoints, const unsigned aChannels=1) {
      _axis.realloc(aPoints);
      _data.resize(aChannels);
      for(unsigned i=0; i<channels(); i++)
        {
        _data[i].realloc(aPoints);
        _data[i].setObjectName(PPPSIGNALCONTAINER_CHANNAME);
        }
      };

    void resize(const unsigned aPoints, const unsigned aChannels=1) {
      _axis.resize(aPoints);
      _data.resize(aChannels);
      for(unsigned i=0; i<channels(); i++)
        {
        _data[i].resize(aPoints);
        _data[i].setObjectName(PPPSIGNALCONTAINER_CHANNAME);
        }
      };

    template <class AType1>
    void prepare(PPPSignalContainer<AType1> &aSour) {
      prepare(aSour.points(), aSour.channels(), aSour.getAxis(), aSour.getObjectName());
      }

    void prepare(const unsigned aPoints, const unsigned aChannels, PPPAxis &aAxis, const string &aN1) {
      realloc(aPoints,aChannels);
      setAxis(aAxis);
      setNames(aN1,aAxis.getObjectName());
      };

    void assign(PPPSignalContainer<AType> &aSour) {
      setNames(aSour.getObjectName(),aSour._axis.getObjectName());
      _axis.assign(aSour._axis);
      _data.resize(aSour.channels());
      for(unsigned i=0; i<channels(); i++)
        {
        _data[i].assign(aSour._data[i]);
        _data[i].setObjectName(PPPSIGNALCONTAINER_CHANNAME);
        }
      };

    void assign(PPPSignalContainer<AType> &aSour, const unsigned aChannel) {
      prepare(aSour.points(),1,aSour.getAxis(),aSour.getObjectName());
      getChannel(0).assign(aSour.getChannel(aChannel));
      };

    /**
     *  generates a view
     */
    inline PPPVectorContainer<AType> & getChannel(const unsigned aChannel) {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aChannel);
      #endif
      return _data[aChannel];
      };

    void getChannel (PPPVectorContainer <AType> &aDest, const unsigned aChannel) {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aChannel);
      #endif
      aDest.link(_data[aChannel]);
      };

    /**
     *  Axis property
     */
    inline PPPAxis & getAxis(void) {
      return _axis;
      };

    inline double getAxis(const unsigned aPoints) {
      return _axis[aPoints];
      };

    inline void setAxis(PPPAxis &aAxis) {
      _axis.assign(aAxis);
      };

    /**
     *  this is the part, with is independed from type of "_data" property
     */
    void setNames(const string &aN1, const string &aN2) {
      PPPBaseTemplate<AType>::setObjectName(aN1);
      _axis.setObjectName(aN2);
      };

    const char *getInfo(void) {
      strstream aDest;
      aDest<<PPPBaseTemplate<AType>::getObjectName()<<": "<<PPPBaseTemplate<AType>::getTypeName()<<" "<<PPPSIGNALCONTAINER_SEC<<endl;
      aDest<<"  "<<getAxis().getInfo()<<endl;
      aDest<<"  "<<PPPSIGNALCONTAINER_CHANNUMB<<": "<<channels() << endl;
      aDest<<"  "<<PPPSIGNALCONTAINER_POINT<<": "<<points() << endl;
      aDest<<"  "<<PPPSIGNALCONTAINER_SAMPL<<": "<<getAxis().getSamplingFreq()<<ends;
      PPPBaseTemplate<AType>::onNotation(aDest.str());
      return PPPBaseTemplate<AType>::getNotation();
      };

    template <class DestType, class TransType>
    void compTransform(PPPSignalContainer<DestType> &aDest, TransType &aTrans)
      {
      aDest.prepare(points(),channels(),getAxis(),aDest.getObjectName());
      for(unsigned j=0; j<channels(); j++)
        for(unsigned i=0; i<points(); i++)
          aDest(i,j) = aTrans((*this)(i,j));
      }

    /*
     * File reading and writing procedures
     */
    void read(const string &aName, bool aIsTime=true, double aSampleFreq=0.0, double aMin=0.0, double aMax=0.0) {
      PPPBaseTemplate<AType>::onMessage(FILE_READASC+aName);
      if(!aIsTime && aSampleFreq == 0.0) PPPBaseObject::onError(ARG_VALUE+string("Read"));
      PPPMatrixContainer<double>  Snew;
      Snew.setShowMessage(false);
      Snew.read(aName.c_str());
      PPPAxis xnew(Snew.rows());
      for(unsigned i=0; i<Snew.rows(); i++)
        xnew[i] = (aIsTime)? Snew(i,0) : (double)i/aSampleFreq;
      unsigned ind1=0, ind2=Snew.rows();
      if(aMin != 0 || aMax !=0)
        {
        ind1 = xnew.locateFloor(aMin);
        ind2 = xnew.locateFloor(aMax);
        }
      if(ind2<ind1 || ind1>Snew.rows() || ind2>Snew.rows())
        PPPBaseObject::onError(ARG_VALUE+string("Read"));
      unsigned count1 = (aIsTime)? Snew.cols()-1 : Snew.cols();
      unsigned count2 = ind2-ind1;
      resize(count2, count1);
      unsigned i,j,k;
      for(i=ind1; i<ind2; i++)
        {
        _axis[i-ind1] = xnew[i];
        for(j=k=0; j<Snew.cols(); j++)
          {
          if(aIsTime && j==0) continue;
          else
            {
            (*this)(i-ind1,k) = Snew(i,j);
            k++;
            }
          }
        }
      _axis.setParams((aIsTime)? PPPAxis::ATnone : PPPAxis::ATlin);
      return;
      };

    template <class TransType>
    void write(const string &aName, TransType aTrans, const bool aIsTime=false, const bool aRepName=false) {
      string newName(aName);
      if(aRepName)
        newName.replace(newName.find("."),1,string(aTrans.getComponentCode())+".");
      strstream str;
      str << FILE_WRITEASC << " " << newName <<" ["<<points()<<"]["<<channels()<<"], "<<aTrans.getComponentName()<<ends;
      PPPBaseTemplate<AType>::onMessage(str.str());
      remove(newName.c_str());
      fstream outfile(newName.c_str(),ios_base::out);
      if(!outfile) PPPBaseObject::onError(FILE_ERROPEN+newName);
      for(unsigned i=0; i<points(); i++)
        {
        outfile.precision(PPPBaseTemplate<AType>::getOutPrecision());
        if(aIsTime) outfile << scientific << _axis[i] << " ";
        for(unsigned j=0; j<channels(); j++)
          outfile << showpos << scientific << aTrans((*this)(i,j)) << " ";
        outfile << endl;
        }
      outfile.close();
      return;
      }

    void write(const string &aName) {
      write(aName,PPPNullTransform<AType>(),true,false);
      };

    /*
     * Math procedures
     */
    void resamplig(const unsigned aBase) {
      PPPSignalContainer<AType> aDest;
      unsigned newsize = points()/aBase;
      PPPAxis newaxis(newsize);
      newaxis.setObjectName(getAxis().getObjectName());
      aDest.prepare(newsize,channels(),newaxis,PPPBaseObject::getObjectName());
      for(unsigned i=0; i<newsize; i++)
        {
        aDest.getAxis()[i] = getAxis(i*aBase);
        for(unsigned k=0; k<channels(); k++)
          aDest(i,k) = (*this)(i*aBase,k);
        }
      aDest.getAxis().setParams(getAxis().getType());
      (*this).assign(aDest);
      };

    void toPowerOfTwo(const unsigned aBase=0) {
      _axis.toPowerOfTwo(aBase);
      for(unsigned i=0; i<channels(); i++) _data[i].resize(_axis.size());
      };

    void rotation(const unsigned aChan1, const unsigned aChan2, const double aPhi) {
      if(aChan1 >= channels() || aChan2 >= channels())
        PPPBaseObject::onError(MEM_ERRINDEX+string("Rotation"));
      AType v1,v2;
      for(unsigned i=0; i<points(); i++)
        {
        v1 =      (*this)(i,aChan1)*cos(aPhi) + (*this)(i,aChan2)*sin(aPhi);
        v2 = -1.0*(*this)(i,aChan1)*sin(aPhi) + (*this)(i,aChan2)*cos(aPhi);
        (*this)(i,aChan1) = v1;
        (*this)(i,aChan2) = v2;
        }
      };

    void multiplication(AType aPar) {
      for(unsigned chann=0; chann<channels(); chann++)
        for(unsigned i=0; i<points(); i++)
          (*this)(i,chann) = (*this)(i,chann)*aPar;
      };

    void integration(void) {
      AType n1, n2, n3;
      for(unsigned chann=0; chann<channels(); chann++)
        {
        n1 = (*this)(0,chann);
        n2 = 0.0;
        for(unsigned i=1; i<points()-1; i++)
          {
          n3 = n2 + (n1 + 4.0*(*this)(i,chann) + (*this)(i+1,chann))/(6.0*getAxis().getSamplingFreq());
          n1 = (*this)(i,chann);
          n2 = n3;
          (*this)(i,chann) = n3;
          }
        (*this)(0,chann) = 0.0;
        (*this)(points()-1,chann) = 0.0;
        }
      };

    void differentiation(PPPSignalContainer<AType> &aDest) {
      aDest.prepare((*this));
      for(unsigned chann=0; chann<channels(); chann++)
        {
        aDest(0,chann) = 0.0;
        for(unsigned i=1; i<points()-1; i++)
          aDest(i,chann) = ((*this)(i,chann) - (*this)(i-1,chann))/(getAxis(i)-getAxis(i-1));
        }
      };

    // Differentiation of phase for complex function
    // A.E.Barnes. The calculation of instantaneous frequency and instantaneous
    // bandwidth // Geophysics, Vol. 57, No. 11, P. 1520-1524 (1992)
    void PhaseDifferentiation(PPPSignalContainer<double> &aDest) {
      PPPcomplex f1,f2,f3;
      double dt;
      PPPVectorContainer<AType> ys;
      PPPVectorContainer<double> yd;
      aDest.prepare(points(),channels(),getAxis(),PPPBaseTemplate<AType>::getObjectName());
      for(unsigned chann=0; chann<channels(); chann++)
        {
        ys.link(getChannel(chann));
        yd.link(aDest.getChannel(chann));
        yd[0] = yd[points()-1] = 0.0;
        for(unsigned i=1; i<points()-1; i++)
          {
          f1 = ys[i-1];
          f2 = ys[i];
          f3 = ys[i+1];
          dt = getAxis()[i+1]-getAxis()[i-1];
          if(abs(f2) == 0.0)
            yd[i] = 0.0;
          else
            yd[i] = (f2.real()*(f3.imag()-f1.imag())-f2.imag()*(f3.real()-f1.real()))/(dt*abs(f2)*abs(f2));
          }
        yd[0] = yd[1];
        yd[points()-1] = yd[points()-2];
        }
      return;
      };

    AType Get(double aX, const unsigned aChannel=0, InterpType aC=ITlin) {
      unsigned count = points();
      PPPVectorContainer<AType> y;
      y.link(getChannel(aChannel));
      unsigned i,j;
      double p;
      AType s;
      s = 0.0;
      if(aX<getAxis()[0] || aX>getAxis()[count-1]) return s;
      switch(aC)
        {
        case ITlin:
          s = y[getAxis().locateFloor(aX)];
          break;
        case ITpolin:
          for(i=0; i<count; i++)
            {
            for(p=1,j=0; j<count; j++)
             {
             if(i == j) continue;
             p *= (aX - getAxis()[j])/(getAxis()[i] - getAxis()[j]);
             }
            s += y[i]*p;
            }
          break;
        case ITsplin:
          if(aC == ITsplin) prepareSpline(aChannel);
          s = evalSplineInterpolation(aX,aChannel);
          break;
        default: PPPBaseObject::onError(PPPSIGNALCONTAINER_ERRINTERP+string("Get"));
        }
      return s;
      };

    void prepareSpline(const unsigned aChannel=0) {
      AType aYp1 = 0.0,aYpn = 0.0;
      unsigned count = points();
      PPPVectorContainer<double> x; x.link(getAxis());
      PPPVectorContainer<AType> y;  y.link(getChannel(aChannel));
      PPPVectorContainer<AType> u(count);
      unsigned i,k;
      double sig;
      AType p,un,qn;
      PPPcomplex z1(-0.5,-0.5), z2(0.5,0.5);
      _splineTmp.resize(count);
      PPPBaseTemplate<AType>::cmplConvert(qn,z1);
      _splineTmp[0] = qn;
      u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-aYp1);
      for(i=2; i<count; i++)
        {
        sig=(x[i-1]-x[i-2])/(x[i]-x[i-2]);
        p=sig*_splineTmp[i-2]+2.0;
        _splineTmp[i-1]=(sig-1.0)/p;
        u[i-1]=(y[i]-y[i-1])/(x[i]-x[i-1]) - (y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
        u[i-1]=(6.0*u[i-1]/(x[i]-x[i-2])-sig*u[i-2])/p;
        }
      PPPBaseTemplate<AType>::cmplConvert(qn,z2);
      un = (3.0/(x[count-1]-x[count-2]))*(aYpn-(y[count-1]-y[count-2])/(x[count-1]-x[count-2]));
      _splineTmp[count-1]=(un-qn*u[count-2])/(qn*_splineTmp[count-2]+1.0);
      for (k=count-1;k>=1;k--) _splineTmp[k-1]=_splineTmp[k-1]*_splineTmp[k]+u[k-1];
      return;
      };

    AType evalSplineInterpolation(double aX, const unsigned aChannel=0) {
      if(aX<getAxis().getMin() || aX>getAxis().getMax()) return (AType)0.0;
      unsigned count = points();
      PPPVectorContainer<AType> y;
      y.link(getChannel(aChannel));
      unsigned klo=1;
      unsigned khi=count;
      unsigned k;
      double h,b,a;
      while (khi-klo > 1)
        {
        k=(khi+klo) >> 1;
        if (getAxis()[k-1] > aX) khi=k;
        else klo=k;
        }
      h=getAxis()[khi-1]-getAxis()[klo-1];
      if (h == 0.0) PPPBaseObject::onError(ARG_VALUE+string("Get"));
      a=(getAxis()[khi-1]-aX)/h;
      b=(aX-getAxis()[klo-1])/h;
      return (a*y[klo-1]+b*y[khi-1]+((a*a*a-a)*_splineTmp[klo-1]+(b*b*b-b)*_splineTmp[khi-1])*(h*h)/6.0);
      };

    /**
     *  file stream operations
     */
    void fwrite(FILE *stream) {
      PPPBaseObject :: fwrite_streaminfo(stream, PPPBaseTemplate<AType>::getObjectVer(), sizeof(AType));
      getAxis().fwrite(stream);
      unsigned newsize = channels();
      std :: fwrite((void*)&newsize,sizeof(newsize),1,stream);
      for (unsigned i=0; i<channels(); i++)
        getChannel(i).fwrite(stream);
      PPPBaseObject :: fwrite(stream);
      };

    void fread(FILE *stream) {
      PPPBaseObject :: fread_streaminfo(stream, PPPBaseTemplate<AType>::getObjectVer(), sizeof(AType));
      getAxis().fread(stream);
      unsigned newsize;
      std :: fread((void*)&newsize,sizeof(newsize),1,stream);
      _data.resize(newsize);
      for (unsigned i=0; i<channels(); i++)
        getChannel(i).fread(stream);
      PPPBaseObject :: fread(stream);
      };

  private:

    void _checkindex(const unsigned aChannel) const {
      if(aChannel>=channels()) PPPBaseObject::onError(MEM_ERRINDEX);
      };

  }; // end of object


#endif




