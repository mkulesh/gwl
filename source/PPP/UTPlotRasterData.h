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
#ifndef PLOTRASTERDATA_H_
#define PLOTRASTERDATA_H_

#include <cmath>
#include <qwt_raster_data.h>

class UTPlotRasterData : public QwtRasterData
{
public:
	typedef enum { LIN, LOG } INTERPTYPE;
	
private:
		
	double const * _data;				// data points
	int _m, _n; 					// dimensions of data
	double _xmin, _xmax, _ymin, _ymax;  	// the point _data[0] has coordinates (_xmin,_ymax)
										// the point _data[_n * _m -1] has coordinates (_xmax,_ymin)
	double _zmin, _zmax;
	
	INTERPTYPE _xinterp, _yinterp;
	
public:
	
	UTPlotRasterData (
		double const * data, 		// the data
		int m, int n ,		// the dimensions
		double xmin, double xmax,	// physical range
		double ymin, double ymax	,	// physical range
		INTERPTYPE xinterp = LIN,
		INTERPTYPE yinterp = LIN
	);
	 
	double xmin () const 
	{
		return _xmin;	
	}
	
	double xmax () const 
	{
		return _xmax;	
	}
	
	double ymin () const 
	{
		return _ymin;	
	}
	
	double ymax () const 
	{
		return _ymax;	
	}
	virtual QwtRasterData * 	copy () const;
	
	/**
	 *  this methods computes the image itself
	 */
	virtual double value (double x, double y) const;

	/**
	 *  the range of values
	 */
	virtual QwtDoubleInterval range ( ) const;

	
	virtual ~UTPlotRasterData();
	
private:
	inline double getAtIntegerValue ( int i, int j ) const 
	{
		if ( (0 <= i ) && ( i < (int)_m) && ( 0 <= j ) && ( j < (int)_n ) ) 
		{
			return _data [ j * _m + i ];	
		}
		else
		{
			return 0.;
		}	
	}
};

/**
 * implementation
 */

#include <algorithm>

UTPlotRasterData::UTPlotRasterData(
		double const * data, 				// the data
		int m, int n ,		// the dimensions
		double xmin, double xmax,	// physical range
		double ymin, double ymax,		// physical range
		INTERPTYPE xinterp,
		INTERPTYPE yinterp
		) :
	QwtRasterData( QwtDoubleRect(xmin, ymin, xmax, ymax ) ),
    _data ( data ),  // better copy ????
    _m(m),
    _n(n),
    _xmin ( xmin ),
    _xmax ( xmax ),
    _xinterp ( xinterp ),
    _yinterp ( yinterp )
{
if(ymax>ymin)
  {
  _ymin = ymin;
  _ymax = ymax;
  }
else
  {
  _ymax = ymin;
  _ymin = ymax;
  }
	_zmin = * std::min_element ( _data, _data + n*m);
	_zmax = * std::max_element ( _data, _data + n*m);
}

QwtRasterData * UTPlotRasterData::copy () const {
	return new UTPlotRasterData ( 
		_data, 
		_m, _n, 
		_xmin,  _xmax, 
		_ymin, _ymax, 
		_xinterp, _yinterp 
	);
}

QwtDoubleInterval UTPlotRasterData::range ( ) const 
{
	return	QwtDoubleInterval ( _zmin, _zmax );
}

double UTPlotRasterData::value ( double x, double y ) const
{
	int i, j; // index for x and y value

	i=j=0;
	
	switch ( _xinterp )
	{
		case ( LIN ) :
		{
			i = (int) (0.5 + (((int)_m)-1) * ( x - _xmin ) / (_xmax - _xmin));
			break;	
		}
		case ( LOG ) :
		{
			i = (int) (0.5 + (((int)_m)-1) * log ( x / _xmin ) / log (_xmax/_xmin));
			break;	
		}
	}// end switch
	
	switch ( _yinterp )
	{
		case ( LIN ) :
		{
			j = (int) (0.5 + (((int)_n)-1) * ( y - _ymin ) / (_ymax - _ymin));
			break;
		}
		case ( LOG ) :
		{
			j = (int) (0.5 + (((int)_n)-1) * log ( y / _ymin ) / log (_ymax/_ymin));
			break;	
		}
	}// end switch

	return getAtIntegerValue (  abs(i), abs(j)  );
}

UTPlotRasterData::~UTPlotRasterData()
{
}

#endif /*PLOTRASTERDATA_H_*/
