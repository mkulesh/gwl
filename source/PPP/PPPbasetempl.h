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
#ifndef _PPPBASETEMPLATE
#define _PPPBASETEMPLATE

#define PPPBASETEMPLATE_INT        "integer"
#define PPPBASETEMPLATE_REAL       "real"
#define PPPBASETEMPLATE_LREAL      "long real"
#define PPPBASETEMPLATE_CMPL       "complex"
#define PPPBASETEMPLATE_EPAR2D     "2D elliptic"
#define PPPBASETEMPLATE_EPAR3D     "3D elliptic"
#define PPPBASETEMPLATE_MAPM       "huge prec real"
#define PPPBASETEMPLATE_ERR        "invalid component index in procedure: "
#define PPPBASETEMPLATE_FREECODE   ""

#ifdef __BCPLUSPLUS__
#define TEMPLNULLVALUE NULL
#else
#define TEMPLNULLVALUE 0
#endif

/************************************************************************
 * PPPBaseTemplate
 ************************************************************************/
template<class AType> class PPPBaseTemplate : public PPPBaseObject
{
public:

    PPPBaseTemplate (void)
    {
    }
    
    const char *getTypeName (void)
    {
        PPPNullTransform<AType> trans;
        return trans.getComponentName();
    }
    
    AType nullValue (void)
    {
        return (AType) TEMPLNULLVALUE;
    }
    
    bool isInteger (void)
    {
        return (sizeof(AType) == sizeof(int));
    }
    
    bool isReal (void)
    {
        return (sizeof(AType) == sizeof(double));
    }
    
    bool isComplex (void)
    {
        return (sizeof(AType) == sizeof(PPPcomplex));
    }
    
    double getComponent (PPPcomplex const &aVal, unsigned const aC)
    {
        switch (aC)
        {
        case TCre:
            return real(aVal);
        case TCim:
            return imag(aVal);
        case TCabs:
            return abs(aVal);
        case TCarg:
            return arg(aVal);
        default:
            onError(PPPBASETEMPLATE_ERR + string("getComponent"));
        }
        return 0;
    }
    
    void cmplConvert (double &aVal, PPPcomplex const &aCmpl)
    {
        aVal = real(aCmpl);
    }
    
    void cmplConvert (PPPcomplex &aVal, PPPcomplex const &aCmpl)
    {
        aVal = aCmpl;
    }
    
};
// end of object

#undef TEMPLNULLVALUE
#endif

