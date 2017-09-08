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
#ifndef _PPPELLIPSE
#define _PPPELLIPSE

#define PPPELLIPSE_RMAX           "instantaneous major axis"
#define PPPELLIPSE_RMIN           "instantaneous minor axis"
#define PPPELLIPSE_RMIDD          "instantaneous second minor axis"
#define PPPELLIPSE_RATIO          "ratio of radius"
#define PPPELLIPSE_RATIO1         "ratio of radius: second minor/minor"
#define PPPELLIPSE_RATIO2         "ratio of radius: second minor/major"
#define PPPELLIPSE_WX             "instantaneous phase x"
#define PPPELLIPSE_WY             "instantaneous phase y"
#define PPPELLIPSE_WZ             "instantaneous phase z"
#define PPPELLIPSE_PHDIFF         "phase difference"
#define PPPELLIPSE_TILT           "tilt angle"
#define PPPELLIPSE_PLANARX        "planarity angle x-component"
#define PPPELLIPSE_PLANARY        "planarity angle y-component"
#define PPPELLIPSE_PLANARZ        "planarity angle z-component"
#define PPPELLIPSE_PLANCOSX       "planarity-X angle"
#define PPPELLIPSE_PLANCOSY       "planarity-Y angle"
#define PPPELLIPSE_PLANCOSZ       "planarity-Z angle"
#define PPPELLIPSE_SRATIO         "sign of ratio"
#define PPPELLIPSE_DIP            "dip angle of rmax"
#define PPPELLIPSE_AZIMUTH        "azimuth of rmax"
#define PPPELLIPSE_PHASE1         "phase 1"
#define PPPELLIPSE_PHASE2         "phase 2"

#define PPPELLIPSE_CODERMAX       "(rmax)"
#define PPPELLIPSE_CODERMIN       "(rmin)"
#define PPPELLIPSE_CODERMIDD      "(rmidd)"
#define PPPELLIPSE_CODERATIO      "(ratio)"
#define PPPELLIPSE_CODERATIO1     "(ratio1)"
#define PPPELLIPSE_CODERATIO2     "(ratio2)"
#define PPPELLIPSE_CODEWX         "(wx)"
#define PPPELLIPSE_CODEWY         "(wy)"
#define PPPELLIPSE_CODEWZ         "(wz)"
#define PPPELLIPSE_CODEPHDIFF     "(phdiff)"
#define PPPELLIPSE_CODETILT       "(tilt)"
#define PPPELLIPSE_CODEPLANARX    "(planarx)"
#define PPPELLIPSE_CODEPLANARY    "(planary)"
#define PPPELLIPSE_CODEPLANARZ    "(planarz)"
#define PPPELLIPSE_CODEPLANCOSX   "(plancosx)"
#define PPPELLIPSE_CODEPLANCOSY   "(plancosy)"
#define PPPELLIPSE_CODEPLANCOSZ   "(plancosz)"
#define PPPELLIPSE_CODESRATIO     "(signratio)"
#define PPPELLIPSE_CODEDIP        "(dip)"
#define PPPELLIPSE_CODEAZIMUTH    "(azimuth)"
#define PPPELLIPSE_CODEPHASE1     "(phase1)"
#define PPPELLIPSE_CODEPHASE2     "(phase2)"

#define PPPELLIPSE_RENE       "Rene method"
#define PPPELLIPSE_MOROZOV    "Morozov method"
#define PPPELLIPSE_SUMCOVAR   "amount covariance method"
#define PPPELLIPSE_COVAR      "adaptive covariance method"
#define PPPELLIPSE_COMPLEX    "complex ellipse method"
#define PPPELLIPSE_MLINET     "time maximum line method"
#define PPPELLIPSE_MLINEF     "frequency maximum line method"

/************************************************************************
 * Constants
 ***********************************************************************/
#define PPPELLIPSE_TYNI 1E-10
static const unsigned TEall        = 0;
static const unsigned TErmax       = 1;  // instantaneous major axis
static const unsigned TErmin       = 2;  // instantaneous minor axis
static const unsigned TErmidd      = 3;  // 3D only: instantaneous second minor axis
static const unsigned TEratio      = 4;  // ratio of radius
static const unsigned TEratio1     = 5;  // 3D only: ratio of radius: second minor/minor
static const unsigned TEratio2     = 6;  // 3D only: ratio of radius: second minor/major
static const unsigned TEwx         = 7;  // instantaneous phase x
static const unsigned TEwy         = 8;  // instantaneous phase y
static const unsigned TEwz         = 9;  // 3D only: instantaneous phase z
static const unsigned TEplanarx    = 10; // 3D only: planarity angle x-component
static const unsigned TEplanary    = 11; // 3D only: planarity angle y-component
static const unsigned TEplanarz    = 12; // 3D only: planarity angle z-component
static const unsigned TEphasediff  = 13; // 2D only: phase difference
static const unsigned TEtilt       = 14; // 2D only: tilt angle
static const unsigned TEplancosx   = 15; // 3D only: directional cosine: planarity-X angle
static const unsigned TEplancosy   = 16; // 3D only: directional cosine: planarity-Y angle
static const unsigned TEplancosz   = 17; // 3D only: directional cosine: planarity-Z angle
static const unsigned TEsignratio  = 18; // 3D only: sign of ratio
static const unsigned TEdip        = 19; // 3D only: dip angle of rmax
static const unsigned TEazimuth    = 20; // 3D only: azimuth of rmax
static const unsigned TEphase1     = 21; // 2D only: phase 1
static const unsigned TEphase2     = 22; // 2D only: phase 2

typedef enum {WTERene, WTEMorozov, WTESumCovar, WTECovar, WTEComplex, WTEMlinet, WTEMlinef} PPPConstETmachine;

ostream& operator << (ostream& aDest, PPPConstETmachine aSour) {
  switch(aSour)
    {
    case WTERene:     aDest << PPPELLIPSE_RENE; break;
    case WTEMorozov:  aDest << PPPELLIPSE_MOROZOV; break;
    case WTESumCovar: aDest << PPPELLIPSE_SUMCOVAR; break;
    case WTECovar:    aDest << PPPELLIPSE_COVAR; break;
    case WTEComplex:  aDest << PPPELLIPSE_COMPLEX; break;
    case WTEMlinet:   aDest << PPPELLIPSE_MLINET; break;
    case WTEMlinef:   aDest << PPPELLIPSE_MLINEF; break;
    }
  return aDest;
  };

/************************************************************************
 * PPPEllipse
 ***********************************************************************/
template<unsigned N> class PPPEllipse
  {
  private:
    double rmin[N];   // rmin axis
    double rmax[N];   // rmax axis
    double rmidd[N];  // additional radius for ellipsoide
    double inst[N];   // istantaneuou frequency
    double tilt,phi,phase1,phase2;

  public:
    PPPEllipse(void) {
      set(0.0);
      };

    void set(double const b) {
      for(unsigned i=0; i<N; i++) rmin[i] = rmax[i] = rmidd[i] = inst[i] = b;
      tilt = phi = phase1 = phase2 = b;
      };

    // eigenvectors
    double vecmax(unsigned aInd)  const { return rmax[aInd]; };
    double vecmin(unsigned aInd)  const { return rmin[aInd]; };
    double vecmidd(unsigned aInd) const { return rmidd[aInd]; };

    // eigenvalues
    double absmin(void) const {
      double sum = 0.0;
      for(unsigned i=0; i<N; i++) sum += rmin[i]*rmin[i];
      return sqrt(sum);
      };

    double absmax(void) const {
      double sum = 0.0;
      for(unsigned i=0; i<N; i++) sum += rmax[i]*rmax[i];
      return sqrt(sum);
      };

    double absmidd(void) const {
      double sum = 0.0;
      for(unsigned i=0; i<N; i++) sum += rmidd[i]*rmidd[i];
      return sqrt(sum);
      };


    // ellipticity ratio
    double ratio(void) const {
      double min = absmin();
      double max = absmax();
      if(max < min) return 0.0;
      if(max < PPPELLIPSE_TYNI) return 0.0;
      return fabs(min/max);
      };

    double ratio1(void) const {
      double min = absmidd();
      double max = absmin();
      if(max < min) return 0.0;
      if(max < PPPELLIPSE_TYNI) return 0.0;
      return fabs(min/max);
      };

    double ratio2(void) const {
      double min = absmidd();
      double max = absmax();
      if(max < min) return 0.0;
      if(max < PPPELLIPSE_TYNI) return 0.0;
      return fabs(min/max);
      };

    double sign_ratio() const {
      double px = planarity_x();
      double py = planarity_y();
      double pz = planarity_z();
      double s = px*rmidd[0] + py*rmidd[1] + pz*rmidd[2];
      if(s > 0.0) return 1.0;
             else return -1.0;
      };

    // planarity vector and it directional cosine
    double planarity_x(void) const {
      return rmax[1]*rmin[2]-rmax[2]*rmin[1];
      };

    double planarity_y(void) const {
      return -(rmax[0]*rmin[2]-rmax[2]*rmin[0]);
      };

    double planarity_z(void) const {
      return rmax[0]*rmin[1]-rmax[1]*rmin[0];
      };

    double planarity(void) const {
      double px = planarity_x();
      double py = planarity_y();
      double pz = planarity_z();
      return sqrt(px*px + py*py + pz*pz);
      };

    double planarity_cosx(void) const {
      if(planarity() == 0.0) return 0.0;
      return (acos(fabs(planarity_x())/planarity()));
      };

    double planarity_cosy(void) const {
      if(planarity() == 0.0) return 0.0;
      return (acos(fabs(planarity_y())/planarity()));
      };

    double planarity_cosz(void) const {
      if(planarity() == 0.0) return 0.0;
      return (acos(fabs(planarity_z())/planarity()));
      };

    // Dip and rise angle of rmax
    double tilt_rmax() const {
      return tilt;
      };

    double dip_rmax() const {
      if(rmax[2] == 0.0) return 0.0;
      return atan(sqrt(rmax[0]*rmax[0]+rmax[1]*rmax[1])/rmax[2]);
      };

    // Azimuth of rmax
    double azimuth_rmax() const {
      if(rmax[0] == 0.0 || rmax[1] == 0)
        return 0.0;
      else
        return atan2(rmax[1],rmax[0]);
      };

    // Phase difference
    double phasediff() const {
      return phi;
      };

    double phase(unsigned const aComp) const {
      return (aComp==1)? phase1 : phase2;
      };

    // Instantaneous frequency
    double instant(unsigned const aComp) const {
      return inst[aComp];
      };

    void sum(PPPEllipse a, PPPEllipse b) {
      set(0.0);
      for(unsigned i=0; i<N; i++)
        {
        rmin[i] = a.rmin[i] + b.rmin[i];
        rmax[i] = a.rmax[i] + b.rmax[i];
        rmidd[i] = a.rmidd[i] + b.rmidd[i];
        inst[i] = a.inst[i] + b.inst[i];
        }
      tilt = a.tilt + b.tilt;
      phi = a.phi + b.phi;
      phase1 = a.phase1 + b.phase1;
      phase2 = a.phase2 + b.phase2;
      };

    void mult(PPPEllipse a, double b) {
      set(0.0);
      for(unsigned i=0; i<N; i++)
        {
        rmin[i] = a.rmin[i]*b;
        rmax[i] = a.rmax[i]*b;
        rmidd[i] = a.rmidd[i]*b;
        inst[i] = a.inst[i]*b;
        }
      tilt = a.tilt*b;
      phi = a.phi*b;
      phase1 = a.phase1*b;
      phase2 = a.phase2*b;
      };


    /* @article{Rene1986Geoph,
        author = {R. M. Rene and J. L. Fitter and P. M. Forsyth and K. Y. Kim and D. J. Murray and J. K. Walters and J. D. Westerman},
        title = {Multicomponent seismic studies using complex trace analysis},
        journal = {Geophysics},
        year = {1986},
        volume = {51},
        number = {6},
        pages = {1235-1251},
        urllink = {http://link.aip.org/link/?GPY/51/1235/1}
      } */
    void ComplexTraceParams(PPPcomplex &aC1, PPPcomplex &aC2) {
      set(0.0);
      PPPcomplex z1,z2;
      double At1 = abs(aC1);
      double At2 = abs(aC2);
      z1 = conj(aC1)*aC2;
      phi = arg(z1);
      double s0 = At1*At1 + At2*At2;
      double s1 = At1*At1 - At2*At2;
      double s2 = 2.0*At1*At2*cos(phi);
      rmax[0] = sqrt((s0+sqrt(s1*s1+s2*s2))/2.0);
      if(s0 > sqrt(s1*s1+s2*s2))
        rmin[0] = sqrt((s0-sqrt(s1*s1+s2*s2))/2.0);
      else
        rmin[0] = 0.0;
      z2 = PPPcomplex(s1, s2);
      tilt = arg(z2)/2.0;
      return;
      };

    /* @article{Morozov1996Geoph,
        author = {I. B. Morozov and S. B. Smithson},
        title = {Instantaneous polarization attributes and directional filtering},
        journal = {Geophysics},
        year = {1996},
        volume = {61},
        number = {3},
        pages = {872-881},
        urllink = {http://link.aip.org/link/?GPY/61/872/1}
      } */
    void MorozovParams(PPPcomplex &aC1, PPPcomplex &aC2) {
      set(0.0);
      PPPcomplex z1,z2,ax,ay,bx,by;
      PPPcomplex A = (aC1*aC1 + aC2*aC2)/2.0;
      PPPcomplex B = (aC1 + aC2)*(aC1 + aC2)/2.0;
      // ellipce phase
      double eps = 0.00001;
      z1 = A + eps*B;
      phase1 = arg(z1)/2.0;
      // instantaneous major vector
      z2 = PPPcomplex(0.0, -phase1);
      ax = aC1*exp(z2);
      ay = aC2*exp(z2);
      rmax[0] = sqrt(ax.real()*ax.real() + ay.real()*ay.real());
      // instantaneous minor vector
      z2 = PPPcomplex(0.0, -phase1-M_PI/2.0);
      bx = aC1*exp(z2);
      by = aC2*exp(z2);
      rmin[0] = sqrt(bx.real()*bx.real() + by.real()*by.real());
      // tilt
      tilt = atan(ay.real()/ax.real());
      return;
      };

    void MorozovParams(PPPcomplex &aC1, PPPcomplex &aC2, PPPcomplex &aC3) {
      set(0.0);
      PPPcomplex z1,z2,z3;
      PPPcomplex A = (aC1*aC1 + aC2*aC2 + aC3*aC3)/2.0;
      PPPcomplex B = (aC1 + aC2 + aC3)*(aC1 + aC2 + aC3)/2.0;
      // ellipce phase
      double eps = 0.00001;
      z1 = A + eps*B;
      phi = arg(z1)/2.0;
      // instantaneous major vector
      z2 = PPPcomplex(0.0, -phi);
      z3 = aC1*exp(z2); rmax[0] = z3.real();
      z3 = aC2*exp(z2); rmax[1] = z3.real();
      z3 = aC3*exp(z2); rmax[2] = z3.real();
      // instantaneous minor vector
      z2 = PPPcomplex(0.0, -phi-M_PI/2.0);
      z3 = aC1*exp(z2); rmin[0] = z3.real();
      z3 = aC2*exp(z2); rmin[1] = z3.real();
      z3 = aC3*exp(z2); rmin[2] = z3.real();
      return;
      };

    void MorozovParamsInverse(PPPcomplex &aC1, PPPcomplex &aC2, PPPcomplex &aC3) {
      aC1 = PPPcomplex(rmax[0]*cos(phi)-rmin[0]*sin(phi), rmax[0]*sin(phi)+rmin[0]*cos(phi));
      aC2 = PPPcomplex(rmax[1]*cos(phi)-rmin[1]*sin(phi), rmax[1]*sin(phi)+rmin[1]*cos(phi));
      aC3 = PPPcomplex(rmax[2]*cos(phi)-rmin[2]*sin(phi), rmax[2]*sin(phi)+rmin[2]*cos(phi));
      return;
      };

    /* @BOOK{Kanasewich1981,
        author = {E. R. Kanasewich},
        title = {Time Sequence Analysis in Geophysics},
        year = {1981},
        publisher = {University of Alberta Press, Edmonton, Alberta}
      } */
    void CovarianceSumParams(PPPVectorContainer<double> &aV1, PPPVectorContainer<double> &aV2) {
      set(0.0);
      unsigned j,k,l;
      PPPMatrixContainer<double> cov(2,2);
      for(k=0; k<2; k++) for(l=0; l<2; l++) cov(k,l) = 0.0;
      // mean value of time window
      double sum1=0.0, sum2=0.0;
      for(j=0; j<aV1.size(); j++)
        {
        sum1 += aV1[j];
        sum2 += aV2[j];
        }
      sum1 /= (double)aV1.size();
      sum2 /= (double)aV1.size();
      // Calculation of covariance matrix
      for(j=0; j<aV1.size(); j++)
        {
        cov(0,0) = cov(0,0) + (aV1[j]-sum1)*(aV1[j]-sum1);
        cov(0,1) = cov(0,1) + (aV1[j]-sum1)*(aV2[j]-sum2);
        cov(1,0) = cov(0,1);
        cov(1,1) = cov(1,1) + (aV2[j]-sum2)*(aV2[j]-sum2);
        }
      for(k=0; k<2; k++) for(l=0; l<2; l++) cov(k,l) = 2.0*cov(k,l)/((double)aV1.size());
      // Calculation of eigen values end eigen vectors
      PPPLinAlg<double> lin;
      PPPVectorContainer<double> eig;
      PPPMatrixContainer<double> vec;
      lin.jacobi(cov,eig,vec);
      eig[0] = fabs(eig[0]);
      eig[1] = fabs(eig[1]);
      // Calculation of params
      rmin[0] = sqrt(eig.getMinValue());
      rmax[0] = sqrt(eig.getMaxValue());
      tilt = (eig[0]<eig[1])? atan(vec(0,0)/vec(0,1)) : atan(vec(1,0)/vec(1,1));
      return;
      };

    void CovarianceSumParams(PPPVectorContainer<double> &aV1, PPPVectorContainer<double> &aV2, PPPVectorContainer<double> &aV3) {
      set(0.0);
      unsigned j,k,l;
      PPPMatrixContainer<double> cov(3,3);
      for(k=0; k<3; k++) for(l=0; l<3; l++) cov(k,l) = 0.0;
      // mean value of time window
      double sum1=0.0, sum2=0.0, sum3=0.0;
      for(j=0; j<aV1.size(); j++)
        {
        sum1 += aV1[j];
        sum2 += aV2[j];
        sum3 += aV3[j];
        }
      sum1 /= (double)aV1.size();
      sum2 /= (double)aV1.size();
      sum3 /= (double)aV1.size();
      // Calculation of covariance matrix
      for(j=0; j<aV1.size(); j++)
        {
        cov(0,0) = cov(0,0) + (aV1[j]-sum1)*(aV1[j]-sum1);
        cov(0,1) = cov(0,1) + (aV1[j]-sum1)*(aV2[j]-sum2);
        cov(0,2) = cov(0,2) + (aV1[j]-sum1)*(aV3[j]-sum3);
        cov(1,0) = cov(0,1);
        cov(1,1) = cov(1,1) + (aV2[j]-sum2)*(aV2[j]-sum2);
        cov(1,2) = cov(1,2) + (aV2[j]-sum2)*(aV3[j]-sum3);
        cov(2,0) = cov(0,2);
        cov(2,1) = cov(1,2);
        cov(2,2) = cov(2,2) + (aV3[j]-sum3)*(aV3[j]-sum3);
        }
      for(k=0; k<3; k++) for(l=0; l<3; l++) cov(k,l) = 2.0*cov(k,l)/((double)aV1.size());
      // Calculation of eigen values end eigen vectors
      PPPLinAlg<double> lin;
      PPPVectorContainer<double> eig;
      PPPMatrixContainer<double> vec;
      lin.jacobi(cov,eig,vec);
      // sort of eigen values
      for(unsigned i=0; i<3; i++) eig[i] = fabs(eig[i]);
      double fmax = eig.getMaxValue();
      unsigned imax = (eig[0]==fmax)? 0 : (eig[1]==fmax)? 1:2;
      double fmin = eig.getMinValue();
      unsigned imin = (eig[0]==fmin)? 0 : (eig[1]==fmin)? 1:2;
      unsigned imidd = 0;
      for(unsigned i=0; i<3; i++) if(i!=imin && i!=imax) imidd=i;
      eig[imidd] = sqrt(eig[imidd]);
      eig[imax] = sqrt(eig[imax]);
      eig[imin] = sqrt(eig[imin]);
      // Normirungs of eigen vectors
      double n;
      for(unsigned i=0; i<3; i++)
        {
        n = sqrt(vec(0,i)*vec(0,i) + vec(1,i)*vec(1,i) + vec(2,i)*vec(2,i));
        for(unsigned j=0; j<3; j++) vec(j,i) = eig[i]*vec(j,i)/n;
        }
      for(unsigned j=0; j<3; j++)
        {
        rmax[j]  = vec(j,imax);
        rmin[j]  = vec(j,imidd);
        rmidd[j] = vec(j,imin);
        }
      return;
      };

    /* @article{Diallo2006GeophB,
        author = {M. S. Diallo and M. Kulesh and M. Holschneider and K. Kurennaya and F. Scherbaum},
        title = {Instantaneous polarization attributes based on an adaptive approximate covariance method},
        journal = {Geophysics},
        year = {2006},
        volume = {71},
        number = {5},
        pages = {V99-V104},
        keywords = {seismic waves; seismology; geophysical signal processing; geophysical techniques},
        urllink = {http://link.aip.org/link/?GPY/71/V99/1}
      } */
    double CovarianceInteg(PPPcomplex aZx, PPPcomplex aZy, double aWx, double aWy, double aT) {
      PPPMathFunc mfunc;
      double val = abs(aZx)*abs(aZy)*(
        mfunc.sinc(aT*(aWy-aWx)/2.0)*cos(arg(aZy)-arg(aZx))+
        mfunc.sinc(aT*(aWy+aWx)/2.0)*cos(arg(aZy)+arg(aZx))
        );
      double mean = (mfunc.sinc(aT*aWx/2.0)*aZx.real())*(mfunc.sinc(aT*aWy/2.0)*aZy.real());
      return val - mean;
      };

    void CovarianceParams(PPPcomplex &aC1, PPPcomplex &aC2, double aW1, double aW2, unsigned aN=1) {
      set(0.0);
      // Calculation of covariance matrix
      inst[0] = aW1/(2.0*M_PI);
      inst[1] = aW2/(2.0*M_PI);
      PPPMatrixContainer<double> cov(2,2);
      double p = 2.0*M_PI*((double)aN);
      cov(0,0) = CovarianceInteg(aC1,aC1,aW1,aW1,p/aW1);
      cov(0,1) = CovarianceInteg(aC1,aC2,aW1,aW2,p*2.0/(aW1+aW2));
      cov(1,0) = cov(0,1);
      cov(1,1) = CovarianceInteg(aC2,aC2,aW2,aW2,p/aW2);
      // Calculation of eigen values end eigen vectors
      PPPLinAlg<double> lin;
      PPPVectorContainer<double> eig;
      PPPMatrixContainer<double> vec;
      lin.jacobi(cov,eig,vec);
      eig[0] = fabs(eig[0]);
      eig[1] = fabs(eig[1]);
      // Calculation of params
      rmin[0] = sqrt(eig.getMinValue());
      rmax[0] = sqrt(eig.getMaxValue());
      tilt = (eig[0]<eig[1])? atan(vec(0,0)/vec(0,1)) : atan(vec(1,0)/vec(1,1));
      return;
      };

    void CovarianceParams(PPPcomplex &aC1, PPPcomplex &aC2, PPPcomplex &aC3,
      double aW1, double aW2, double aW3, unsigned aN=1) {
      set(0.0);
      // Calculation of covariance matrix
      PPPMatrixContainer<double> cov(3,3);
      inst[0] = aW1/(2.0*M_PI);
      inst[1] = aW2/(2.0*M_PI);
      inst[2] = aW3/(2.0*M_PI);
      double p = 2.0*M_PI*((double)aN);
      cov(0,0) = CovarianceInteg(aC1,aC1,aW1,aW1,p/aW1);
      cov(0,1) = CovarianceInteg(aC1,aC2,aW1,aW2,p*2.0/(aW1+aW2));
      cov(0,2) = CovarianceInteg(aC1,aC3,aW1,aW3,p*2.0/(aW1+aW3));
      cov(1,0) = cov(0,1);
      cov(1,1) = CovarianceInteg(aC2,aC2,aW2,aW2,p/aW2);
      cov(1,2) = CovarianceInteg(aC2,aC3,aW2,aW3,p*2.0/(aW2+aW3));
      cov(2,0) = cov(0,2);
      cov(2,1) = cov(1,2);
      cov(2,2) = CovarianceInteg(aC3,aC3,aW3,aW3,p/aW3);
      // Calculation of eigen values end eigen vectors
      PPPLinAlg<double> lin;
      PPPVectorContainer<double> eig;
      PPPMatrixContainer<double> vec;
      lin.jacobi(cov,eig,vec);
      // sort of eigen values
      for(unsigned i=0; i<3; i++) eig[i] = fabs(eig[i]);
      double fmax = eig.getMaxValue();
      unsigned imax = (eig[0]==fmax)? 0 : (eig[1]==fmax)? 1:2;
      double fmin = eig.getMinValue();
      unsigned imin = (eig[0]==fmin)? 0 : (eig[1]==fmin)? 1:2;
      unsigned imidd = 0;
      for(unsigned i=0; i<3; i++) if(i!=imin && i!=imax) imidd=i;
      eig[imidd] = sqrt(eig[imidd]);
      eig[imax] = sqrt(eig[imax]);
      eig[imin] = sqrt(eig[imin]);
      // Normirungs of eigen vectors
      double n;
      for(unsigned i=0; i<3; i++)
        {
        n = sqrt(vec(0,i)*vec(0,i) + vec(1,i)*vec(1,i) + vec(2,i)*vec(2,i));
        for(unsigned j=0; j<3; j++) vec(j,i) = eig[i]*vec(j,i)/n;
        }
      for(unsigned j=0; j<3; j++)
        {
        rmax[j]  = vec(j,imax);
        rmin[j]  = vec(j,imidd);
        rmidd[j] = vec(j,imin);
        }
      return;
      };

    /* @article{Diallo2006GeophA,
        author = {M. S. Diallo and M. Kulesh and M. Holschneider and F. Scherbaum and F. Adler},
        title = {Characterization of polarization attributes of seismic waves using continuous wavelet transforms},
        journal = {Geophysics},
        year = {2006},
        volume = {71},
        number = {3},
        pages = {V67-V77},
        keywords = {seismic waves; wavelet transforms; geophysical signal processing},
        urllink = {http://link.aip.org/link/?GPY/71/V67/1}
      } */
    void CWTParams(PPPcomplex &aC1, PPPcomplex &aC2, double aWx=0.0, double aWy=0.0) {
      set(0.0);
      if(abs(aC1) == 0.0 && abs(aC2) == 0.0) return;
      inst[0] = (aWx+aWy)/(4.0*M_PI);
      inst[1] = (aWy-aWx)/(4.0*M_PI);
      rmax[0] = (abs(aC1) + abs(aC2))/2.0;
      rmin[0] = (abs(aC1) - abs(aC2))/2.0;
      phase1 = (arg(aC1)+arg(aC2))/2.0;
      phase2 = (arg(aC1)-arg(aC2))/2.0;
      PPPcomplex z1,z2;
      z2 = aC1*aC2;
      tilt = arg(z2)/2.0;
      z1 = (aC1+conj(aC2))/(aC1-conj(aC2));
      phi = arg(z1)+M_PI/2.0;
      if(phi > M_PI)  phi -= 2.0*M_PI;
      if(phi < -M_PI) phi += 2.0*M_PI;
      return;
      };

    void CWTParamsInverse(PPPcomplex &aC1, PPPcomplex &aC2) {
      PPPcomplex z1(0,phase1+phase2);
      PPPcomplex z2(0,phase1-phase2);
      aC1 = ((rmax[0]+rmin[0])*exp(z1))/2.0;
      aC2 = ((rmax[0]-rmin[0])*exp(z2))/2.0;
      return;
      };

    /* @TechReport{Kulesh2007DFGc,
        author = {M. Kulesh and M. Holschneider and M. S. Diallo},
        title = {Geophysics Wavelet Library: Applications of the Continuous Wavelet Transform to the Polarization and Dispersion Analysis of Signals},
        institution = {Preprint series of the DFG priority program 1114 ``Mathematical methods for time series analysis and digital image processing''},
        month = {March},
        number = {Preprint 156},
        year = {2007},
        urllink = {http://www.math.uni-bremen.de/zetem/DFG-Schwerpunkt/preprints/pdf/156.pdf}
      } */
    void CWTParamsDefomDiss(double aAtn, double aPhase) {
      rmax[0] *= aAtn;
      rmin[0] *= aAtn;
      phase2 -= aPhase;
      return;
      };

    void CWTParamsDefomDiss(double aAtn, PPPEllipse &aElli2) {
      rmax[0] *= aAtn;
      rmin[0] *= aAtn;
      phase1 = aElli2.phase1;
      phase2 = aElli2.phase2;
      return;
      };

  };   // end of object



#endif
