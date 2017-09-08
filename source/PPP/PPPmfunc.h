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
#ifndef _PPPMATHFUNC
#define _PPPMATHFUNC

#define PPPMATHFUNC_NAME      "math function"

/************************************************************************
 * PPPMathFunc
 ************************************************************************/
class PPPMathFunc : public PPPBaseObject
  {

  public:
    PPPMathFunc(void) {
      setObjectName(PPPMATHFUNC_NAME);
      };

    /** Power of two **********************************************************/
    int isPowerOfTwo(int aSize) {
      if ( aSize < 2 )  return false;
      if ( aSize & (aSize-1) )  return false;
      return true;
      };

    /** Computing SINC function ***********************************************/
    double sinc(double x){
      if(x==0) return 1.0;
      else return sin(x)/x;
      };

    /** Factorial *************************************************************/
    double fact(int aN) {
      if(aN == 0 || aN == 1) return 1.0;
      double s=1.0;
      for(int i=1; i<=aN; i++) s*=i;
      return s;
      };

  }; // end of object



#endif
