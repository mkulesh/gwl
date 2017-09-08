#ifndef UTFFT_H_
#define UTFFT_H_

#include <string>

#define success true;

/************************************************************************
 * UTFFT: Interafce object for fast fourier transform
 *   computes
 *   forward:    \sum_{k=0}^{n-1} a_k e^{-ikl/n}
 *   backward:   {1/n} \sum_{k=0}^{n-1} a_k e^{ikl/n}
 *  any implementation of this FFT must implement the following methods
 *   1. init(...) or corresponding constructor
 *   2. execute()
 *   3. clear()
 *  The user has to supply the necessary memory. Size of necessary memory can be
 *  checked with realSizeOfIn() and realSizeOfOut()
 *  Make suhre that in != out
 ***********************************************************************/
class UTFFT {

  public:
    typedef enum {C2CFORWARD, C2CBACKWARD, R2C, C2R} transtype;
    typedef enum {CFFT, RFFT } algotype;

  private:
    double      *_in;
    double      *_out;
    unsigned    _n;
    transtype   _transtype;
    algotype  	_algotype;
    std::string      _method;

    class       FFTDATA;
    FFTDATA     *data;

  public:
                   
    /*
     *  Following procedures have to be realized in implementation section
     */

    /*
     *  Default contructor. If necessary, allocates FFTDATA *data variable
     */
    UTFFT();

    /*
     *  Destructor. If necessary, destroys FFTDATA *data variable
     */
    ~UTFFT();

    /*
     *  Initialization with memory allocated by user
     *  in = malloc ( realSizeOfIn() * sizeof ( double ) )
     *  out = malloc ( realSizeOfOut() * sizeof ( double ) )
     *  If necessary, allocates FFTDATA *data variable
     */
    bool setData(double *in, double *out);

    /*
     *  Calculation FFT
     */
    bool execute(void);

    /*
     *  If FFTDATA *data variable has been allocated in setData(), then destroys FFTDATA *data variable
     */
    bool clear(void);

    /*
     *  Pre-defined interface, do not imlement in implementation section
     */
    inline void setType(unsigned n, transtype trans, algotype algo = CFFT) {
      _n = n;
      _transtype = trans;
      _algotype = algo;
      };

    inline bool init(
        double *in,                    // in: [re0,im0,re1,im1,.....
        double *out,                   // out: [re0,im0,re1,im1,....
        unsigned n, 		       // n: number of points (power of 2)
        transtype trans = C2CFORWARD,  // trans: forward, backwart transformation
        algotype algo = CFFT) {
      setType(n,trans,algo);
      return setData(in, out);
      };

    inline unsigned size(void) const { return _n; };
    inline double * inbegin(void) const { return _in; };
    inline double * inend(void) const { return _in + realSizeOfIn(); };
    inline double * outbegin(void) const { return _out; };
    inline double * outend (void) const { return _out + realSizeOfOut(); };
    inline std::string &getMethodName(void) { return _method; };

    unsigned realSizeOfIn (void) const {
      if ( _algotype == CFFT )
        {
        switch (_transtype)
          {
          case C2CFORWARD: case C2CBACKWARD: case C2R: return 2*_n;
          case (R2C): return _n;
          }
        }
      return 0;
      };

    unsigned realSizeOfOut (void) const {
      if ( _algotype == CFFT )
        {
        switch (_transtype)
          {
          case C2CFORWARD: case C2CBACKWARD: case R2C: return 2*_n;
          case C2R: return _n;
          }
        }
      return 0;
      };

  }; // end of object

#endif /*UTFFT_H_*/
