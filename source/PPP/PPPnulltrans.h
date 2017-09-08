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

#ifndef _PPPNULLTRANSFORM
#define _PPPNULLTRANSFORM

/************************************************************************
 * PPPNullTransform
 ************************************************************************/
template<class T> class PPPNullTransform : public unary_function<T, T>
  {
  public:
    inline T const &operator()(T const &aVal) const { return aVal; };
    const char* getComponentName(int const &aVal) { return PPPBASETEMPLATE_INT; };
    const char* getComponentName(double const &aVal) { return PPPBASETEMPLATE_REAL; };
    const char* getComponentName(long double const &aVal) { return PPPBASETEMPLATE_LREAL; };
    const char* getComponentName(PPPcomplex const &aVal) { return PPPBASETEMPLATE_CMPL; };
    const char* getComponentName(PPPellipse2D const &aVal) { return PPPBASETEMPLATE_EPAR2D; };
    const char* getComponentName(PPPellipse3D const &aVal) { return PPPBASETEMPLATE_EPAR3D; };
    const char *getComponentCode(void) { return PPPBASETEMPLATE_FREECODE; };
    #ifdef PPPCONF_USEMAPM
    const char* getComponentName(MAPM const &aVal) { return PPPBASETEMPLATE_MAPM; };
    #endif
    const char *getComponentName(void) {
      T testval;
      return getComponentName(testval);
      };
  }; // end of object

#endif  
