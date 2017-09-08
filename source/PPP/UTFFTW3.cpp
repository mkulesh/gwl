#include "UTFFT.h"
#include "fftw3.h"

/************************************************************************
 * An interface to fftw3 library
 ***********************************************************************/
class UTFFT::FFTDATA
  {
  public:
    fftw_plan p;

    ~FFTDATA(void) {
//      fftw_destroy_plan (p);
      };

  }; // end of object

/************************************************************************
 * UTFFT
 ***********************************************************************/
UTFFT::UTFFT() {
  _method = "fft_w3";
  data = new UTFFT::FFTDATA;
  }

UTFFT::~UTFFT() {
  delete data;
  }

bool UTFFT::setData(double *in, double *out) {
  _in = in;
  _out = out;
  if ( _algotype == CFFT )
    {
    switch (_transtype )
      {
      case C2CFORWARD:
        data -> p = fftw_plan_dft_1d(_n,
                    (double(*)[2]) in,
                    (double(*)[2]) out,
                    FFTW_FORWARD, FFTW_ESTIMATE);
        break;
      case C2CBACKWARD:
        data -> p = fftw_plan_dft_1d(_n,
                    (double(*)[2]) in,
                    (double(*)[2]) out,
                    FFTW_BACKWARD,FFTW_ESTIMATE);
        break;
      case R2C:
        data -> p = fftw_plan_dft_r2c_1d(_n,
                    (double *) in,
                    (double(*)[2]) out,
                    FFTW_ESTIMATE);
        break;
      case C2R:
        data -> p = fftw_plan_dft_c2r_1d(_n,
                    (double(*)[2]) in,
                    (double *) out,
                    FFTW_PRESERVE_INPUT);
        break;
      default:
        return !success;
        break;
      }
    return success;
    }
  else
    {
    return !success; // TODO:
    }
  }

bool UTFFT::execute(void) {
  fftw_execute(data -> p);
  switch ( _transtype )
    {
    case C2CFORWARD: break;
    case C2CBACKWARD: case C2R:
      {
      double fac = 1./_n;
      for ( double * i = outbegin(); i<outend(); ++i) (*i) *= fac;
      break;
      }
    case R2C:
      {
      double * tmp = outbegin()+2;
      for ( double * i = outend() - 2;  i> outbegin()+_n; i -=2 )
        {
        *i = *tmp;
        (*(i+1)) = - (*(tmp+1));
        tmp += 2;
        }
      }
      break;
    default:
      return !success;
      break;
    }
  return success;
  }

bool UTFFT::clear(void) {
   fftw_destroy_plan(data -> p);
  }

