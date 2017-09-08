/* 
 * This file is a part of GWL - Geophysical Wavelet Library
 * Copyright (C) 2007 Mikhail Kulesh and Matthias Holschneider
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * For more information please visit: http://users.math.uni-potsdam.de/~gwl
 * Email: mkulesh@math.uni-potsdam.de
 * ICQ: 103-405-403
 */

#ifndef _PPPSIGNALPLOT
#define _PPPSIGNALPLOT

#define PPPSIGNALPLOT_NAME          "plot of a signal"
#define PPPSIGNALPLOT_ERRTYPE       "incorrect signal data type in input file: "


/************************************************************************
 * PPPSignalPlot
 ***********************************************************************/
class PPPSignalPlot : public PPPBaseObject
  {
  private:

    typedef enum {REAL_SIGNAL, COMPLEX_SIGNAL, NONE} SIGTYPE;
    PPPSignalContainer < double > * _signalR;
    PPPSignalContainer < PPPcomplex > * _signalC;

  public:

    PPPSignalPlot(): _signalR (NULL), _signalC (NULL) {
      setObjectName(PPPSIGNALPLOT_NAME);
      };

    void read(char const * aName) {
      PPPObjectIO io;
      PPPBaseObject * tmp = io.read(aName);
      if ( dynamic_cast < PPPSignalContainer < PPPcomplex > * > (tmp) )
        {
        _signalC = dynamic_cast < PPPSignalContainer < PPPcomplex > * > (tmp);
        }
      else if ( dynamic_cast < PPPSignalContainer < double >  * > (tmp) )
        {
        _signalR = dynamic_cast < PPPSignalContainer < double >  * > (tmp);
        }
      else
        {
        onError(PPPSIGNALPLOT_ERRTYPE+string(aName));
        }
      };

    PPPAxis const & getTime() const {
      if ( _signalC != NULL )
        return _signalC -> getAxis();
      else
        return _signalR -> getAxis();
      };

    void getPlotCurves (unsigned trace, vector < QwtPlotCurve * > & curves, double & minvalue, double & maxvalue) {
      if ( _signalR != NULL )
        {
        maxvalue = _signalR -> getChannel(trace).getMaxValue();
        minvalue = _signalR -> getChannel(trace).getMinValue();
        QwtPlotCurve * pc = new QwtPlotCurve ();
        curves.push_back ( pc );
        pc -> setData(getTime().begin(),
          _signalR -> getChannel(trace).begin(),
          getTime().size());
        }
      else if ( _signalC != NULL )
        {
        PPPVectorContainer<double> tmp;
        // real part
        QwtPlotCurve * pcre = new QwtPlotCurve ();
        curves.push_back ( pcre );
        _signalC -> getChannel(trace).compTransform(tmp, PPPcomplexRe());
          pcre -> setData(getTime().begin(),tmp.begin(),getTime().size());
          maxvalue = tmp.getMaxValue();
          minvalue = -maxvalue;
        // imag part
        QwtPlotCurve * pcim = new QwtPlotCurve ();
        curves.push_back ( pcim );
        _signalC -> getChannel(trace).compTransform(tmp, PPPcomplexIm());
          pcim -> setData(getTime().begin(),tmp.begin(),getTime().size());
          if(tmp.getMaxValue() > maxvalue)
            {
            maxvalue = tmp.getMaxValue();
            minvalue = -maxvalue;
            }
        /*
        // modulus
        QwtPlotCurve * pcab = new QwtPlotCurve ();
        curves.push_back ( pcab );
        _signalC -> getChannel(trace).compTransform(tmp, PPPcomplexAbs());
          pcab -> setData(getTime().begin(),tmp.begin(),getTime().size());
        // negative modulus
        QwtPlotCurve * pcma = new QwtPlotCurve ();
        curves.push_back ( pcma );
        for(unsigned i=0; i<tmp.size(); i++) tmp[i] = -tmp[i];
          pcma -> setData(getTime().begin(),tmp.begin(),getTime().size());
        */
        }
      };

  };  // end of object

#endif

