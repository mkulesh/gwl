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
#include "UTFFT.h"
#include "fftw3.h"

/************************************************************************
 * An interface to fftw3 library
 ***********************************************************************/
class UTFFT::FFTDATA
{
public:
    fftw_plan p;

    ~FFTDATA (void)
    {
    //      fftw_destroy_plan (p);
    }
    
};
// end of object

/************************************************************************
 * UTFFT
 ***********************************************************************/
UTFFT::UTFFT ()
{
    _method = "fft_w3";
    data = new UTFFT::FFTDATA;
}

UTFFT::~UTFFT ()
{
    delete data;
}

bool UTFFT::setData (double *in, double *out)
{
    _in = in;
    _out = out;
    if (_algotype == CFFT)
    {
        switch (_transtype)
        {
        case C2CFORWARD:
            data->p = fftw_plan_dft_1d(_n, (double (*)[2]) in, (double (*)[2]) out,
            FFTW_FORWARD,
                                       FFTW_ESTIMATE);
            break;
        case C2CBACKWARD:
            data->p = fftw_plan_dft_1d(_n, (double (*)[2]) in, (double (*)[2]) out,
            FFTW_BACKWARD,
                                       FFTW_ESTIMATE);
            break;
        case R2C:
            data->p = fftw_plan_dft_r2c_1d(_n, (double *) in, (double (*)[2]) out,
            FFTW_ESTIMATE);
            break;
        case C2R:
            data->p = fftw_plan_dft_c2r_1d(_n, (double (*)[2]) in, (double *) out,
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
        return !success;
    }
}

bool UTFFT::execute (void)
{
    fftw_execute(data->p);
    switch (_transtype)
    {
    case C2CFORWARD:
        break;
    case C2CBACKWARD:
    case C2R:
    {
        double fac = 1. / _n;
        for (double * i = outbegin(); i < outend(); ++i)
            (*i) *= fac;
        break;
    }
    case R2C:
    {
        double * tmp = outbegin() + 2;
        for (double * i = outend() - 2; i > outbegin() + _n; i -= 2)
        {
            *i = *tmp;
            (*(i + 1)) = -(*(tmp + 1));
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

bool UTFFT::clear (void)
{
    fftw_destroy_plan(data->p);
}

