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
#ifndef _PPPAPPROXPOLINOM
#define _PPPAPPROXPOLINOM

/************************************************************************
 * PPPApproximatePolinom
 ***********************************************************************/
class PPPApproximatePolinom : public PPPApproximate
{
public:
    
    PPPApproximatePolinom (unsigned aSize)
    {
        setObjectName (PPPAPPROXIMATE_POL);
        _type = APTpolin;
        _params.resize(aSize);
    }
    
    const char *getInfo (void)
    {
        strstream str;
        str << "Func(f) = sum((a[i]*f^i)/i), i=0..." << (size() - 1) << endl;
        str << "  a[i] = " << _params.vectorToStr() << ends;
        onNotation(str.str());
        return getNotation();
    }
    
    double Func (double f)
    {
        if (f < _xmin || f > _xmax || size() == 0) return 0.0;
        if (f == 0.0) return _params[0];
        double res = _params[0];
        for (unsigned i = 1; i < size(); i++)
            res += _params[i] * pow(f, (double) i) / ((double) i);
        return res;
    }
    
    double FuncDivPar (double f, unsigned Ai)
    {
        if (f < _xmin || f > _xmax || Ai >= size()) return 0.0;
        if (Ai == 0) return 1.0;
        return pow(f, (double) Ai) / ((double) Ai);
    }
    
    double FuncDivX (double f)
    {
        if (f < _xmin || f > _xmax || size() == 0) return 0.0;
        double res = 0.0;
        for (unsigned i = 1; i < size(); i++)
            res += _params[i] * pow(f, (double) (i - 1));
        return res;
    }
    
    double FuncDivXPar (double f, unsigned Ai)
    {
        if (f < _xmin || f > _xmax || Ai == 0 || Ai >= size()) return 0.0;
        return pow(f, (double) (Ai - 1));
    }
    
};
// end of object

#endif
