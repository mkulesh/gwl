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
#ifndef _PPPLINALG
#define _PPPLINALG

#define PPPLINALG_NAME     "linear algebra"
#define PPPLINALG_ERRQUAD    "matrix is not quadratic in procedure: "
#define PPPLINALG_ERRSING    "matrix is singular in procedure: "
#define PPPLINALG_ERRSYMM    "matrix is not symmetric in procedure: "
#define PPPLINALG_ERRHERM    "matrix is not hermetian in procedure: "
#define PPPLINALG_ERRITER    "too many iterations in routine: "
#define PPPLINALG_ERRSIZE    "invalid size of input matrix in procedure: "
#define PPPLINALG_ERRPAR     "invalid value of input parameter in procedure: "

/************************************************************************
 * PPPLinAlg
 ***********************************************************************/
template<class AType> class PPPLinAlg : public PPPBaseTemplate<AType>
  {
  public:

    PPPLinAlg(void) {
      PPPBaseTemplate<AType>::setObjectName(PPPLINALG_NAME);
      };

    // Gauss-Jordan matrix inversion and linear equation solution
    // Linear equation solution by Gauss-Jordan elimination.
    // aMat[1..n][1..n] is the input matrix.
    // aB is input containing the right-hand side vector.
    // On output, a is replaced by its matrix inverse, and b is
    // replaced by the corresponding set of solution vectors.
    void matrixInversion(PPPMatrixContainer<AType> &aMat, PPPVectorContainer<AType> &aB, unsigned aSize) {
      if(aMat.cols() != aMat.rows())
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRQUAD+string("MatrixInversion"));
      unsigned i,icol,irow,j,k,l,ll;
      AType big,dum,pivinv;
      PPPVectorContainer<unsigned> indxc(aSize),indxr(aSize),ipiv(aSize);
      for(i=1;i<=aSize;i++)
        {
        big=0.0;
        for(j=1;j<=aSize;j++)
          if (ipiv[j-1] != 1)
            for(k=1;k<=aSize;k++)
              {
              if (ipiv[k-1] == 0)
                {
                if (fabs(aMat(j-1,k-1)) >= big)
                  {
                  big=fabs(aMat(j-1,k-1));
                  irow=j;
                  icol=k;
                  }
                }
              }
        ++(ipiv[icol-1]);
        if (irow != icol)
          {
          for(l=1;l<=aSize;l++) SWAP(aMat(irow-1,l-1),aMat(icol-1,l-1));
          SWAP(aB[irow-1],aB[icol-1]);
          }
        indxr[i-1]=irow;
        indxc[i-1]=icol;
        if (aMat(icol-1,icol-1) == 0.0)
          PPPBaseTemplate<AType>::onError(PPPLINALG_ERRSING+string("MatrixInversion"));
        pivinv=1.0/aMat(icol-1,icol-1);
        aMat(icol-1,icol-1)=1.0;
        for(l=1;l<=aSize;l++) aMat(icol-1,l-1) *= pivinv;
        aB[icol-1] *= pivinv;
        for(ll=1;ll<=aSize;ll++) if (ll != icol)
          {
          dum=aMat(ll-1,icol-1);
          aMat(ll-1,icol-1)=0.0;
          for(l=1;l<=aSize;l++) aMat(ll-1,l-1) -= aMat(icol-1,l-1)*dum;
          aB[ll-1] -= aB[icol-1]*dum;
          }
        }
      for(l=aSize;l>=1;l--)
        {
        if (indxr[l-1] != indxc[l-1]) for(k=1;k<=aSize;k++)
          SWAP(aMat(k-1,indxr[l-1]-1),aMat(k-1,indxc[l-1]-1));
        }
      return;
      };

    // multiplication of two matrix
    void multiplication(PPPMatrixContainer<AType> &aDest, PPPMatrixContainer<AType> &aM1, PPPMatrixContainer<AType> &aM2) {
      if(aM1.rows() != aM2.cols() || aM1.cols() != aM2.rows())
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRSIZE+string("Multiplication"));
      aDest.resize(aM1.cols(),aM1.rows());
      AType sum;
      unsigned i,j,k;
      for(k=0; k<aM1.cols(); k++)
        for(j=0; j<aM1.rows(); j++)
          {
          sum = 0.0;
          for(i=0; i<aM1.cols(); i++) sum = sum + aM1(j,i)*aM2(i,k);
          aDest(j,k) = sum;
          }
      };

    // multiplication of matrix and vector
    void multiplication(PPPVectorContainer<AType> &aDest, PPPMatrixContainer<AType> &aM1, PPPVectorContainer<AType> &aM2) {
      if(aM1.cols() != aM2.size())
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRSIZE+string("Multiplication"));
      aDest.resize(aM1.rows());
      double sum;
      unsigned i,j;
      for(j=0; j<aM1.rows(); j++)
        {
        sum = 0;
        for(i=0; i<aM1.cols(); i++) sum += aM1(j,i)*aM2[i];
        aDest[j] = sum;
        }
      };

    // check, if matrix aMat is Square
    inline bool isSquare(const PPPMatrixContainer<AType> &aMat) const {
      return (aMat.cols()==aMat.rows());
      };

    // check, if matrix aMat is Symetric
    bool isSymetric(PPPMatrixContainer<AType> &aMat) const {
      unsigned i,j;
      if(aMat.cols()!=aMat.rows()) return false;
      for (i=0; i<aMat.cols(); i++)
        for (j=i+1; j<aMat.rows(); j++)
          if (aMat(i,j) != aMat(j,i)) return false;
      return true;
      };

    // check, if matrix aMat is Hermetian
    bool isHermetian(PPPMatrixContainer<AType> &aMat) const {
      unsigned i,j;
      if(aMat.cols()!=aMat.rows()) return false;
      for (i=0; i<aMat.cols(); i++)
        {
        if(imag(aMat(i,i)) != 0) return false;
        for (j=i+1; j<aMat.rows(); j++)
          if (real(aMat(i,j)) != real(aMat(j,i)) || imag(aMat(i,j)) != -imag(aMat(j,i))) return false;
        }
      return true;
      };

    // Computes all eigenvalues and eigenvectors of a real symmetric matrix
    // aMat[1..n][1..n]. On output, elements of aMat above the diagonal are destroyed.
    // aEVal[1..n] returns the eigenvalues of aMat.
    // aEVec[1..n][1..n] is a matrix whose columns contain, on output, the normalized
    // eigenvectors of aMat.
    void jacobi(PPPMatrixContainer<double> &aMat, PPPVectorContainer<double> &aEVal, PPPMatrixContainer<double> &aEVec)
      {
      if(!isSquare(aMat))
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRQUAD+string("Jacobi"));
      if(!isSymetric(aMat))
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRSYMM+string("Jacobi"));
      unsigned j,iq,ip,i;
      double tresh,theta,tau,t,sm,s,g,h,c;
      unsigned n = aMat.rows();
      aEVal.resize(n);
      aEVec.resize(n,n);
      PPPVectorContainer<double> b(n),z(n);
      for (ip=1;ip<=n;ip++)
        {
        for (iq=1;iq<=n;iq++) aEVec(ip-1,iq-1)=0.0;
        aEVec(ip-1,ip-1)=1.0;
        }
      for (ip=1;ip<=n;ip++)
        {
        b[ip-1]=aEVal[ip-1]=aMat(ip-1,ip-1);
        z[ip-1]=0.0;
        }
      unsigned nrot = 0;
      for (i=1;i<=50;i++)
        {
        sm=0.0;
        for (ip=1;ip<=n-1;ip++)
          for (iq=ip+1;iq<=n;iq++)
              sm += fabs(aMat(ip-1,iq-1));
        if (sm == 0.0) return;
        if (i < 4) tresh=0.2*sm/(n*n);
        else       tresh=0.0;
        for (ip=1;ip<=n-1;ip++)
          {
          for (iq=ip+1;iq<=n;iq++)
            {
            g=100.0*fabs(aMat(ip-1,iq-1));
            if (i > 4
                 && (double)(fabs(aEVal[ip-1])+g) == (double)fabs(aEVal[ip-1])
                 && (double)(fabs(aEVal[iq-1])+g) == (double)fabs(aEVal[iq-1]))
              aMat(ip-1,iq-1)=0.0;
            else
              if (fabs(aMat(ip-1,iq-1)) > tresh)
                {
                h=aEVal[iq-1]-aEVal[ip-1];
                if ((double)(fabs(h)+g) == (double)fabs(h))
                  t=(aMat(ip-1,iq-1))/h;
                else
                  {
                  theta=0.5*h/(aMat(ip-1,iq-1));
                  t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                  if (theta < 0.0) t = -t;
                  }
                c=1.0/sqrt(1+t*t);
                s=t*c;
                tau=s/(1.0+c);
                h=t*aMat(ip-1,iq-1);
                z[ip-1] -= h;
                z[iq-1] += h;
                aEVal[ip-1] -= h;
                aEVal[iq-1] += h;
                aMat(ip-1,iq-1)=0.0;
                for (j=1;j<=ip-1;j++)    _jacobiRot(aMat,j,ip,j,iq,s,tau);
                for (j=ip+1;j<=iq-1;j++) _jacobiRot(aMat,ip,j,j,iq,s,tau);
                for (j=iq+1;j<=n;j++)    _jacobiRot(aMat,ip,j,iq,j,s,tau);
                for (j=1;j<=n;j++)       _jacobiRot(aEVec,j,ip,j,iq,s,tau);
                nrot++;
                }
            }
          }
        for (ip=1;ip<=n;ip++)
          {
          b[ip-1] += z[ip-1];
          aEVal[ip-1]=b[ip-1];
          z[ip-1]=0.0;
          }
        }
      PPPBaseTemplate<AType>::onError(PPPLINALG_ERRITER+string("Jacobi"));
      }


    void eigComplex(PPPMatrixContainer<PPPcomplex>& aMat, PPPVectorContainer<double>& aEVal, PPPMatrixContainer<PPPcomplex>& aEVec){
      if(!isHermetian(aMat))
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRHERM+string("EigComplex"));
      PPPMatrixContainer<double> M;
      M.resize(2*aMat.rows(),2*aMat.cols());
      unsigned i,j,k;
      for (i=0; i<aMat.rows(); i++)  for (j=0; j<aMat.cols(); j++)
        {
        M(i,j) = real(aMat(i,j));
        M(i+aMat.rows(),j+aMat.cols()) = real(aMat(i,j));
        M(i+aMat.rows(),j) = imag(aMat(i,j));
        M(i,j+aMat.cols()) = -imag(aMat(i,j));
        }
      PPPMatrixContainer<double> V;
      PPPVectorContainer<double> E;

      jacobi(M,E,V);

      aEVal.resize(aMat.rows());
      aEVec.resize(aMat.rows(),aMat.cols());
      PPPVectorContainer<unsigned> ind(aMat.rows());
      for(i=0; i<aMat.rows(); i++) ind[i] = -1;
      bool flag;
      for(i=0; i<E.size(); i++)
        {
        flag = false;
        for(j=0; ind[j]>=0; j++)
          if(fabs(E[i]-aEVal[ind[j]])<1E-10)
            {
            flag = true;
            break;
            }
        if(!flag)
          {
          ind[j] = j;
          aEVal[ind[j]] = E[i];
          for(k=0; k<aMat.rows(); k++)
            {
            real(aEVec(k,ind[j])) = -V(k,i);
            imag(aEVec(k,ind[j])) = -V(k+aMat.cols(),i);
            }
          }
        }
      };


    // Mean values of matrix columnwise (l==1) or rowwise (l==2)
    void mean(unsigned l,PPPVectorContainer<AType>& v, PPPMatrixContainer<AType>& m) {
      unsigned i,j;
      switch(l)
        {
	case 1: v.resize(m.cols());
	   for(i=0;i<m.cols();i++) for(j=0;j<m.rows();j++)
             v[i]=v[i]+m(i,j)/m.rows();
	   break;
	case 2: v.resize(m.rows()); // TODO: please check there was rows() only
	   for(i=0;i<m.rows();i++) for(j=0;j<m.cols();j++)
             v[i]=v[i]+m(j,i)/m.cols();
           break;
        default: PPPBaseTemplate<AType>::onError(PPPLINALG_ERRPAR+string("mean"));
        }
      };  


    unsigned upDiagonal(PPPMatrixContainer<AType>& a) {
      if(!isSquare(a))
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRQUAD+string("upDiagonal"));
      unsigned i,j,k,step=0;
      AType s;
      for(i=0; i<a.rows()-1; i++)
      if( abs(a(i,i)) != 0.0 )
        {
        for(j=i+1; j<a.rows(); j++)
           {
           s=a(j,i)/a(i,i);
           for(k=i; k<a.rows(); k++) a(j,k)-=a(i,k)*s;
           }
        }
       else
        {
        for(j=i; abs(a(j,i)) == 0.0; j++) if(j == a.rows()-1) break;
        if(j == a.rows()-1) continue;
        for(k=0; k<a.rows(); k++) { SWAP(a(i,k), a(j,k)); }
	step++;
        }
      return step;  
      };

    AType evalDeterminant(PPPMatrixContainer<AType> &aMat)
      {
      if(!isSquare(aMat))
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRQUAD+string("evalDeterminant"));
      if(aMat.rows() == 1) return aMat(0,0);
      if(aMat.rows() == 2) return aMat(0,0)*aMat(1,1) - aMat(1,0)*aMat(0,1);
      PPPMatrixContainer<AType> a(aMat);
      unsigned step = upDiagonal(a);
      AType d=1.0, p=pow(-1.0,(double)step);
      for(unsigned i=0; i<a.rows(); i++) d*=a(i,i);
      return d*p;
      };

    AType evalDeterminantCmpl(PPPMatrixContainer<AType> &aMat)
      {
      if(!isSquare(aMat))
        PPPBaseTemplate<AType>::onError(PPPLINALG_ERRQUAD+string("evalDeterminant"));
      if(aMat.rows() == 1) { return aMat(0,0); }
      if(aMat.rows() == 2) { return aMat(0,0)*aMat(1,1) - aMat(1,0)*aMat(0,1); }
      unsigned step = upDiagonal(aMat);
      AType d(1.0,0.0), p(pow(-1.0,(double)step),0.0);
      for(unsigned i=0; i<aMat.rows(); i++) d*=aMat(i,i);
      return d*p;
      };

  private:


    void _jacobiRot(PPPMatrixContainer<double> &a, unsigned i, unsigned j, unsigned k, unsigned l, double s, double tau){
      double g,h;
      g=a(i-1,j-1);
      h=a(k-1,l-1);
      a(i-1,j-1)=g-s*(h+g*tau);
      a(k-1,l-1)=h+s*(g-h*tau);
      };

  }; // end of object



#endif
