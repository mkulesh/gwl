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

#ifndef _PPPOPTILM
#define _PPPOPTILM

#define PPPOPTILM_NAME        "Levenberg-Marquardt optimization"
#define PPPOPTILM_POINTS      "number of points"
#define PPPOPTILM_PARAMS      "number of params"
#define PPPOPTILM_PREC        "maximal precision"
#define PPPOPTILM_ITER        "iteration"
#define PPPOPTILM_SUBITER     "subiteration"
#define PPPOPTILM_COST        "cost function"
#define PPPOPTILM_SUMM        "overall results"
#define PPPOPTILM_FIXPAR      "number of fixed params"
#define PPPOPTILM_FAXPAR1     "fixed params"
#define PPPOPTILM_ITER1       "number of iterations"

/***********************************************************************
 * PPPOptiLevenbergMarq
 * optimization using Levenberg-Marquardt method
 ************************************************************************/
class PPPOptiLevenbergMarq : public PPPOpti
  {
  private:
    double MaxEps,alamda,ochisq;
    PPPVectorContainer<double> atry,beta,minda,oneda;
    PPPMatrixContainer<double> covar,alpha;

  public:

    PPPOptiLevenbergMarq() {
      setObjectName(PPPOPTILM_NAME);
      MaxEps = 1E-08;
      };

    OptiType getOptiType() { return OTleven; };

    void setParams(PPPVectorContainer<unsigned> &aPar) {
      onError(VIRT_NOTDEF+string("setParams")); return;
      };

    void setParams(PPPVectorContainer<double> &aPar) {
      MaxEps = aPar[0];
      };

    void setParams(PPPMatrixContainer<double> &aPar) {
      onError(VIRT_NOTDEF+string("setParams")); return;
      };

    void Minimize(void) {
      char ss[200];
      sprintf(ss,"%s\n  %s=%d, %s=%d, %s=%4.3e",PPPOPTILM_NAME,PPPOPTILM_PARAMS,params(),
              PPPOPTILM_POINTS,points(),PPPOPTILM_PREC,MaxEps);
      onMessage(ss);
      double oochisq;
      _calcfuncCount = 0;
      _chisq = alamda = 0.0;
      LevenbergMarq_prepare();
      for(unsigned itst=0,k=1;;k++)
        {
        sprintf(ss,"%s: %d, %s: %d, %s: %8.7e",PPPOPTILM_ITER,k,PPPOPTILM_SUBITER,itst,PPPOPTILM_COST,_chisq);
        if(getShowProgress()) onMessage(ss);
        oochisq = _chisq;
        LevenbergMarq_min();
        if(_chisq <= oochisq && fabs(oochisq-_chisq) < MaxEps) itst++; else itst=0;
        if(itst>params()+1) break;
        }
      alamda=0.0;
      LevenbergMarq_min();
      // Notation
      strstream str;
      str << PPPOPTILM_SUMM << ":" << endl;
      str << "  " << PPPOPTILM_POINTS << " = " << points() << endl;
      str << "  " << PPPOPTILM_PARAMS << " = " << params() << endl;
      str << "  " << PPPOPTILM_FIXPAR << " = " << (params()-fixed()) << endl;
      if((params()-fixed()) > 0)
      str << "  " << PPPOPTILM_FAXPAR1 << " = " << fixpar.vectorToStr() << endl;
      str << "  " << PPPOPTILM_ITER1 << " = " << calcfuncCount() << endl;
      str << "  " << PPPOPTILM_COST << " = " << chisq() << endl;
      str << "  " << PPPOPTILM_PREC << " = " << MaxEps << ends;
      onNotation(str.str());
      return;
      };

  private:

    void LevenbergMarq_prepare() {
      covar.resize(params(),params());
      alpha.resize(params(),params());
      atry.resize(params());
      beta.resize(params());
      minda.resize(params());
      alamda = 0.001;
      evalFixed();
      oneda.resize(fixed());
      _chisq = LevenbergMarq_cof(optpar,alpha,beta);
      ochisq = _chisq;
      atry = optpar;
      return;
      };

    void LevenbergMarq_covsrt(PPPMatrixContainer<double> &aCovar) {
      unsigned i,j,k;
      for(i=fixed()+1;i<=params();i++)  for(j=1;j<=i;j++)
        {
        aCovar(i-1,j-1) = 0.0;
        aCovar(j-1,i-1) = 0.0;
        }
      k=fixed();
      for(j=params();j>=1;j--)
        {
        if(fixpar[j-1])
          {
          for(i=1;i<=params();i++) SWAP(aCovar(i-1,k-1),aCovar(i-1,j-1));
          for(i=1;i<=params();i++) SWAP(aCovar(k-1,i-1),aCovar(j-1,i-1));
          k--;
          }
        }
      return;
      };

    double LevenbergMarq_cof(PPPVectorContainer<double> &aPar,PPPMatrixContainer<double> &aAlpha,PPPVectorContainer<double> &aBeta) {
      unsigned j,k,l,m;
      double sig2i,dy;
      for(j=1;j<=fixed();j++)
        {
        for(k=1;k<=j;k++) aAlpha(j-1,k-1) = 0.0;
        aBeta[j-1]=0.0;
        }
      _calcfuncCount++;
      calcfunc(aPar);
      double ychisq=0.0,der1;
      for(unsigned i=1;i<=points();i++)
        {
        dy = getcostfunc(i-1);
        sig2i = 1.0;
        for(j=0,l=1;l<=params();l++)
          {
          if (fixpar[l-1])
            {
            j++;
            for(k=0,m=1;m<=l;m++) if (fixpar[m-1])
              {
              k++;
              der1 = aAlpha(j-1,k-1) + getder(l-1,i-1)*getder(m-1,i-1)*sig2i;
              aAlpha(j-1,k-1) = der1;
              }
            aBeta[j-1] += dy*getder(l-1,i-1)*sig2i;
            }
          }
        ychisq += dy*dy*sig2i;
        }
      for(j=2;j<=fixed();j++) for(k=1;k<j;k++) aAlpha(k-1,j-1) = aAlpha(j-1,k-1);
      return ychisq;
      };

    void LevenbergMarq_min(void)
      {
      unsigned j,k,l;
      for(j=1;j<=fixed();j++)
        {
        for(k=1;k<=fixed();k++) covar(j-1,k-1) = alpha(j-1,k-1);
        covar(j-1,j-1) = alpha(j-1,j-1)*(1.0+alamda);
        oneda[j-1]=beta[j-1];
        }
      PPPLinAlg<double> linalg;
      linalg.matrixInversion(covar,oneda,fixed());
      for(j=1;j<=fixed();j++) minda[j-1]=oneda[j-1];
      if (alamda == 0.0)
        {
        LevenbergMarq_covsrt(covar);
        LevenbergMarq_covsrt(alpha);
        return;
        }
      for(j=0,l=1;l<=params();l++) if(fixpar[l-1]) atry[l-1]=optpar[l-1]+minda[++j-1];
      _chisq = LevenbergMarq_cof(atry,covar,minda);
      if (_chisq < ochisq)
        {
        alamda *= 0.1;
        ochisq = _chisq;
        for(j=1;j<=fixed();j++) for(k=1;k<=fixed();k++)
          alpha(j-1,k-1) = covar(j-1,k-1);
        beta = minda;
        optpar = atry;
        }
      else
        {
        alamda *= 10.0;
        _chisq = ochisq;
        }
      return;
      };

  }; // end of object

#endif


