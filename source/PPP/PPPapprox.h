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
#ifndef _PPPAPPROXIMATE
#define _PPPAPPROXIMATE

#define PPPAPPROXIMATE_VEL         "phase velocity approximation"
#define PPPAPPROXIMATE_GAU         "gauss approximation"
#define PPPAPPROXIMATE_POL         "polynomial approximation"
#define PPPAPPROXIMATE_SPL         "bispline approximation"
#define PPPAPPROXIMATE_COL         "causal cole-cole model"
#define PPPAPPROXIMATE_ERRVAL      "invalid value of approximation type in procedure: "

/************************************************************************
 * PPPApproximate
 ***********************************************************************/
class PPPApproximate : public PPPBaseObject
{
public:
    
    typedef enum
    {
        APTnone, APTvel, APTgauss, APTpolin, APTbspline, APTcolecole, APTtwogauss
    } ApprType;

protected:
    
    ApprType _type;
    PPPVectorContainer<double> _params;
    double _xmin;
    double _xmax;

public:
    
    PPPApproximate (void)
    {
        _type = APTnone;
        _params.resize(0);
        setRange(0.0, 0.0);
    }
    
    inline ApprType getType (void) const
    {
        return _type;
    }
    
    inline unsigned size (void) const
    {
        return _params.size();
    }
    
    inline double getXmin (void) const
    {
        return _xmin;
    }
    
    inline double getXmax (void) const
    {
        return _xmax;
    }
    
    inline PPPVectorContainer<double> & getParams (void)
    {
        return _params;
    }
    
    inline double & operator [] (const unsigned aInd)
    {
        return _params[aInd];
    }
    
    void setRange (double aXmin, double aXmax)
    {
        _xmin = aXmin;
        _xmax = aXmax;
    }
    
    virtual double Func (double) = 0;
    virtual double FuncDivX (double) = 0;
    virtual double FuncDivPar (double, unsigned) = 0;
    virtual double FuncDivXPar (double, unsigned) = 0;
    virtual const char *getInfo (void) = 0;

    void assign (PPPApproximate *aSource)
    {
        if (_type != aSource->_type) onError(PPPAPPROXIMATE_ERRVAL + string("assign"));
        _params.assign(aSource->_params);
        _xmin = aSource->_xmin;
        _xmax = aSource->_xmax;
    }
    
    void link (PPPApproximate * aAppr)
    {
        _params.link(aAppr->_params);
    }
    
};
// end of object

#endif
