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
#ifndef _PPPNNPRED
#define _PPPNNPRED

#define PPPNNPRED_NAME    "NN prediction"
#define PPPNNPRED_ADNAME  "ANN prediction"
#define PPPNNPRED_TREND   "polynomial trend"

/************************************************************************
 * PPPNearesNeighbors
 * @TechReport{2006-DFG1114-139,
 *   author = {Kurennaya, Kristina and Kulesh, Michail and Holschneider, Matthias},
 *   title = {Adaptive metrics in the nearest neighbours method},
 *   url = {http://www.math.uni-bremen.de/zetem/DFG-Schwerpunkt/preprints/prep139.pdf},
 *   institution = {Preprint series of the DFG priority program 1114 ``Mathematical methods for time series analysis and digital image processing''},
 *   month = {April},
 *   year = {2006},
 *   keywords = {time series, prediction, nearest neighbours method, metrics, optimization},
 *   isbn = {},
 *   dbn = {},
 *   file = {2006-DFG1114-139.pdf},
 * }
 ***********************************************************************/
class PPPNearesNeighbors : public PPPBaseObject
  {
  private:
    PPPVectorContainer<double>  _trendCoeff;            // coeffitients for polinomial trend
    PPPVectorContainer<double>  _data;                  // source time serie
    unsigned                    _dim, _step;            // the embedding dimension end the time delay
    double                      _maxElem, _minElem;


  public:
    PPPNearesNeighbors(void) {
      setObjectName(PPPNNPRED_NAME);
      _dim = _step = 0;
      };

    void prepare(PPPVectorContainer<double> &aData, unsigned aDim, unsigned aStep=1) {
      _data.resize(aData.size());
      for (unsigned i=0; i<aData.size(); i++) _data[i] = aData[i];
      _dim = aDim;
      _step = aStep;
      };

    void evalPrediction(PPPVectorContainer<double> &aData, unsigned aPredcount, unsigned aNNcount,
      unsigned aL=2, double aKoeff=1.0, double aSub=0.0) {
      unsigned leng = _data.size();
      if(aL == 3)
        setObjectName(PPPNNPRED_ADNAME);
      else
        setObjectName(PPPNNPRED_NAME);
      strstream str;
      str << endl
        << "  Size of source time serie = " << leng << endl
        << "  Count of predicted points = " << aPredcount << endl
        << "  Count of nearest neighbors = "  << aNNcount << endl
        << "  Embedding dimension = " << _dim << endl
        << "  Time delay = " << _step << endl
        << "  Norm L" << aL << endl
        << "  Amplitude factor = " << aKoeff;
      if(_data.getMinValue() == 0.0 && aSub != 0.0)
        {
        str << endl << "  Substituation in the case of zero minimal element = " << aSub;
        for (unsigned i=0; i<leng; i++) if(_data[i] == 0) _data[i]=aSub;
        }
      str << ends;  
      onMessage(str.str());
      aData.resize(aPredcount);
      char numpoint[100];
      setShowProgress(true);
      PPPVectorContainer<double> errvals;
      errvals.setOutPrecision(2);
      for (unsigned i=0; i<aPredcount; i++)
        {
        double val = aKoeff*_predictPoint(aNNcount,aL,errvals);
        aData[i] = val;
        _data.push_back(val);
        if(i>0)
          {
          sprintf(numpoint, "Point %3d, norm values=%s (%d%%)", i+1, errvals.vectorToStr(), 100*i/(aPredcount-1));
          onMessage(numpoint);
          }
        }
      _data.resize(leng);
      };

    void extractTrend(PPPVectorContainer<double> &aSig, unsigned aSize, bool aMod=true) {
      if(aSig.size() == 0)
        onError("lenght of input signal is 0 in procedure: extractTrend()");
      if(aSize == 0)
        onError("degree of trend is 0 in procedure: extractTrend()");
      unsigned j,k;
      _trendCoeff.realloc(aSize);
      PPPVectorContainer<double> S(aSig);
      PPPLinAlg<double> lin;
      // A(k,j) = (k)^(j);
      PPPMatrixContainer<double> A, At;
      A.realloc(aSig.size(),_trendCoeff.size());
      At.realloc(_trendCoeff.size(),aSig.size());
      for(k=0; k<aSig.size(); k++)
        for(j=0; j<_trendCoeff.size(); j++)
          {
          A(k,j) = pow((double)k,(double)j);
          At(j,k) = A(k,j);
          }
      // B=A'*A;
      PPPMatrixContainer<double> B;
      lin.multiplication(B,A,At);
      // D = (A'*aSource);
      PPPVectorContainer<double> D;
      lin.multiplication(D,At,S);
      // P = inv(B)*D;
      PPPVectorContainer<double> C(_trendCoeff.size()),P;
      lin.matrixInversion(B,C,C.size());
      lin.multiplication(P,B,D);
      for(k=0; k<_trendCoeff.size(); k++) _trendCoeff[k]=P[k];
      if(aMod) for(k=0; k<aSig.size(); k++)  aSig[k] -= approximateTrend(k);
      strstream str;
      str << "Polynomial trend: " << P.vectorToStr() << ends;
      onMessage(str.str());
      return;
      };

    double approximateTrend(double aVal) {
      double sum = 0;
      for(unsigned k=0; k<_trendCoeff.size(); k++) sum += _trendCoeff[k]*pow(aVal,(double)k);
      return sum;
      };

    void restoreTend(PPPVectorContainer<double> &aSig, int aStart=0) {
      for(unsigned k=0; k<aSig.size(); k++)  aSig[k] += approximateTrend(k+aStart);
      };

  private:

    double _predictPoint(unsigned aNNcount, unsigned aL, PPPVectorContainer<double> &aErrVals) {
      PPPVectorContainer<double> ref,cand;
      PPPVectorContainer<double> errs,nearval,predval;
      double min,max,min1,val,pred;
      unsigned i,j;
      // analysis of all candidates and calculation of errors
      _getCandidate(ref, _data.size());
      _maxElem = _data.getMaxValue();
      _minElem = _data.getMinValue();
      if(_minElem == 0)
        onError("minimal element of serie is null in procedure: predictPoint()");
      for (i=_data.size()-_dim*_step; i>=_dim*_step; i--)
        {
        _getCandidate(cand, i);
        val = _getDistance(ref,cand,aL,_data[i],pred);
        if(i==(_data.size()-_dim*_step)) min=max=val;
        if(val<min) min=val;
        if(val>max) max=val;
        errs.push_back(val);
        predval.push_back(pred);
        }
      // search of aNNcount nearest neighbors
      aErrVals.resize(0);
      for (i=0; i<aNNcount; i++)
        {
        for (j=0; j<errs.size(); j++)
          {
          if(errs[j]==min)
            {
            nearval.push_back(predval[j]);
            aErrVals.push_back(errs[j]);
            errs[j] = max;
            }
          if(j==0) min1=errs[j];
          if(errs[j]<min1) min1=errs[j];
          }
        min = min1;
        }
      return nearval.getSumm()/((double)aNNcount);
      };

    void _getCandidate(PPPVectorContainer<double> &aCand, unsigned aInd) {
      if(aInd<_dim*_step)
        {
        strstream str;
        str << "Invalid index for candidate prepare in procedure getCandidate()" << endl
            << "  the embedding dimension=" << _dim << "; the time delay=" << _step << "; index=" << aInd << endl;
        onError(str.str());
        }
      aCand.resize(_dim);
      for (unsigned i=0; i<_dim; i++) aCand[i] = _data[aInd-_step*i-1];
      };

    double _getDistance(PPPVectorContainer<double> &aVec1, PPPVectorContainer<double> &aVec2, int aL, double &aData, double &aPred) {
      if(aVec1.size() != aVec2.size())
        {
        strstream str;
        str << "Invalid dimensions of vectors in procedure getDistance()" << endl
            << "  length of 1th vector=" << aVec1.size() << "; length of 2th vector=" << aVec2.size() << endl;
        onError(str.str());
        }
      double summ=0;
      unsigned i,k;
      // L1 norm
      if(aL == 1)
        {
        for (i=0; i<_dim; i++) summ += fabs(aVec1[i]-aVec2[i]);
        aPred = aData;
        return summ;
        }
      // L2 norm
      if(aL == 2)
        {
        for (i=0; i<_dim; i++) summ += (aVec1[i]-aVec2[i])*(aVec1[i]-aVec2[i]);
        aPred = aData;
        return sqrt(summ);
        }
      // L3 norm = L1 adaptive norm
      if(aL == 3)
        {
        double aLam=1,aMu=0,SumByLam,Med;
        int Prec = 10;
        unsigned aLamParts = (int)(Prec*fabs(_maxElem/_minElem));
        unsigned MedInd  = (_dim&1)?(_dim/2):(_dim/2-1);
        PPPVectorContainer<double> M(_dim);
        for (k=0; k<aLamParts+1; k++)
          {
           SumByLam = 0;
           for (i=0; i<_dim; i++) M[i]= aVec1[i]-k/Prec*aVec2[i];
           Med = med_kth_smallest(M,MedInd);
           for (i=0; i<_dim; i++) SumByLam += fabs(aVec1[i]-k/Prec*aVec2[i]-Med);
           if (k==0 || SumByLam<summ)
              {
              summ = SumByLam;
              aLam = k/Prec;
              aMu = Med;
              }
          }
          aPred = aLam*aData + aMu;
          return summ;
        }
      onError("Invalid value for norm parameter in procedure getDistance()");
      return 0.0;
      };

    // implements N. Wirth's method.
    // It is actually a tool which finds the kth smallest element from a list
    // of values, the special case of a median is implemented through a macro.
    // http://ndevilla.free.fr/median/median/node20.html
    double med_kth_smallest(PPPVectorContainer<double> &a, int k){
      register int i,j,l,m ;
      int n=a.size();
      double x;
      l=0; m=n-1;
      while (l<m)
       {
       x=a[k];
       i=l;
       j=m;
       do
        {
         while (a[i]<x) i++ ;
         while (x<a[j]) j-- ;
         if (i<=j)
           {
            SWAP(a[i],a[j]) ;
            i++ ; j-- ;
           }
        }
       while (i<=j);
       if (j<k) l=i;
       if (k<i) m=j;
       }
      return a[k] ;
      };

  }; // end of object

#endif
