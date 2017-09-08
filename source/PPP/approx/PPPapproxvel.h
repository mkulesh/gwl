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
#ifndef _PPPAPPROXVEL
#define _PPPAPPROXVEL

/************************************************************************
 * PPPApproximateVel
 ***********************************************************************/
class PPPApproximateVel : public PPPApproximate
{
public:
    
    PPPApproximateVel (unsigned aSize)
    {
        setObjectName (PPPAPPROXIMATE_VEL);
        _type = APTvel;
        if (aSize != 3) onError(ARG_VALUE + string("PPPApproximateVel"));
        _params.resize(aSize);
    }
    
    const char *getInfo (void)
    {
        strstream str;
        str << "Func(f) = f/VSp(f), FuncDiv(f)=1.0/VSg(f), VSp(f) = a1 + a2*exp(-f*f/(a3*a3))"
        << endl;
        str << "  a[i] = " << _params.vectorToStr() << ends;
        onNotation(str.str());
        return getNotation();
    }
    
    double Func (double f)
    {
        double Vp = _params[0] + _params[1] * exp(-f * f / (_params[2] * _params[2]));
        return f / Vp;
    }
    
    double FuncDivX (double f)
    {
        double Vp = _params[0] + _params[1] * exp(-f * f / (_params[2] * _params[2]));
        double Vg =
                Vp * Vp / (Vp
                        + 2.0 * f * f * _params[1] * exp(-f * f / (_params[2] * _params[2])) / (_params[2]
                                * _params[2]));
        return 1.0 / Vg;
    }
    
    double FuncDivXPar (double f, unsigned Ai)
    {
        onError(VIRT_NOTDEF + string("FuncDivXPar"));
        return 0.0;
    }
    
    double FuncDivPar (double f, unsigned Ai)
    {
        onError(VIRT_NOTDEF + string("FuncDivPar"));
        return 0.0;
    }
    
};
// end of object

#endif
