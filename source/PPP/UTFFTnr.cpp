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
#ifndef UTFFTNR_CPP
#define UTFFTNR_CPP

#include "UTFFT.h"

/************************************************************************
 * An interface to FFT routines, proposed by Dominic Mazzoni (September 2000)
 * http://www.koders.com/cpp/fid1C48DD8DBB45CA4B87C9A38DAFC532E793443862.aspx
 *
 * The basic algorithm for his code was based on Numerican Recipes
 * in Fortran.  Dominic Mazzoni optimized his code further by reducing
 * array accesses, caching the bit reversal table, and eliminating
 * float-to-double conversions.
 *
 * Some of this code was based on a free implementation of an FFT
 * by Don Cross, available on the web at:
 * http://www.intersrv.com/~dcross/fft.html
 ***********************************************************************/
class UTFFT::FFTDATA
{
private:
    
    double _zeroval;
    int _fastBits;
    int **_bitTable;
    double *_in, *_out, *_tempimag;
    unsigned _n;
    transtype _transtype;

public:
    
    FFTDATA () :
            _zeroval(0.0),
            _fastBits(16)
    {
        _bitTable = new int *[_fastBits];
        int len = 2;
        for (int b = 1; b <= _fastBits; b++)
        {
            _bitTable[b - 1] = new int[len];
            for (int i = 0; i < len; i++)
                _bitTable[b - 1][i] = reverseBits(i, b);
            len <<= 1;
        }
    }

    int numberOfBitsNeeded (int aPow)
    {
        for (int i = 0;; i++)
            if (aPow & (1 << i)) return i;
    }

    int reverseBits (int index, int aBits)
    {
        int i, rev;
        for (i = rev = 0; i < aBits; i++)
        {
            rev = (rev << 1) | (index & 1);
            index >>= 1;
        }
        return rev;
    }

    int fastReverseBits (int index, int aBits)
    {
        if (aBits <= _fastBits)
            return _bitTable[aBits - 1][index];
        else
            return reverseBits(index, aBits);
    }

    void init (double *in, double *out, unsigned n, transtype trans)
    {
        _in = in;
        _out = out;
        _n = n;
        _transtype = trans;
        if (_transtype == C2R) _tempimag = new double[n];
    }

    void clear ()
    {
        if (_transtype == C2R) delete _tempimag;
    }

    inline const double &getInReal (const unsigned aInd)
    {
        return (_transtype == R2C) ? _in[aInd] : _in[2 * aInd];
    }

    inline const double &getInImag (const unsigned aInd)
    {
        return (_transtype == R2C) ? _zeroval : _in[2 * aInd + 1];
    }

    inline const double &getRealOut (const unsigned aInd)
    {
        return (_transtype == C2R) ? _out[aInd] : _out[2 * aInd];
    }

    inline const double &getImagOut (const unsigned aInd)
    {
        return (_transtype == C2R) ? _tempimag[aInd] : _out[2 * aInd + 1];
    }

    inline double &setOutReal (const unsigned aInd)
    {
        return (_transtype == C2R) ? _out[aInd] : _out[2 * aInd];
    }

    inline double &setOutImag (const unsigned aInd)
    {
        return (_transtype == C2R) ? _tempimag[aInd] : _out[2 * aInd + 1];
    }
    
};
// end of object

/************************************************************************
 * UTFFT
 ***********************************************************************/
UTFFT::UTFFT ()
{
    _method = "fft_nr";
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
    return success;
}

bool UTFFT::execute (void)
{
    register unsigned i, j, k, n;
    unsigned BlockSize, BlockEnd;
    int NumBits = data->numberOfBitsNeeded(_n); /* Number of bits needed to store indices */
    double angle_numerator =
            (_transtype == C2CBACKWARD || _transtype == C2R) ? 2.0 * M_PI : -2.0 * M_PI;
    double tr, ti; /* temp real, temp imaginary */
    data->init(_in, _out, _n, _transtype);
    for (i = 0; i < _n; i++)
    {
        j = data->fastReverseBits(i, NumBits);
        double a = data->getInReal(i);
        data->setOutReal(j) = a;
        double b = data->getInImag(i);
        data->setOutImag(j) = b;
    }
    BlockEnd = 1;
    for (BlockSize = 2; BlockSize <= _n; BlockSize <<= 1)
    {
        double delta_angle = angle_numerator / (double) BlockSize;
        double sm2 = sin(-2.0 * delta_angle);
        double sm1 = sin(-delta_angle);
        double cm2 = cos(-2.0 * delta_angle);
        double cm1 = cos(-delta_angle);
        double w = 2.0 * cm1;
        double ar0, ar1, ar2, ai0, ai1, ai2;
        for (i = 0; i < _n; i += BlockSize)
        {
            ar2 = cm2;
            ar1 = cm1;
            ai2 = sm2;
            ai1 = sm1;
            for (j = i, n = 0; n < BlockEnd; j++, n++)
            {
                ar0 = w * ar1 - ar2;
                ar2 = ar1;
                ar1 = ar0;
                ai0 = w * ai1 - ai2;
                ai2 = ai1;
                ai1 = ai0;
                k = j + BlockEnd;
                tr = ar0 * data->getRealOut(k) - ai0 * data->getImagOut(k);
                ti = ar0 * data->getImagOut(k) + ai0 * data->getRealOut(k);
                data->setOutReal(k) = data->getRealOut(j) - tr;
                data->setOutImag(k) = data->getImagOut(j) - ti;
                data->setOutReal(j) += tr;
                data->setOutImag(j) += ti;
            }
        }
        BlockEnd = BlockSize;
    }
    if (_transtype == C2CBACKWARD || _transtype == C2R)
    {
        double denom = (double) _n;
        for (i = 0; i < _n; i++)
        {
            data->setOutReal(i) /= denom;
            data->setOutImag(i) /= denom;
        }
    }
    data->clear();
    return success;
}

bool UTFFT::clear (void)
{
    return success;
}

#endif /*UTFFTNR_CPP*/
