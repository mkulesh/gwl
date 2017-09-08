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

#ifndef _PPPOBJECTIO
#define _PPPOBJECTIO

#define PPPOBJECTIO_NAME      "i/o interface"

/************************************************************************
 * PPPObjectIO
 ***********************************************************************/
class PPPObjectIO : public PPPBaseObject
  {
  public:

    typedef enum  
    {
    	UNDEF, 
    	AXIS, 
    	SIGD, 
    	SIGC, 
    	SIGE2D, 
    	SIGE3D, 
    	SPECD, 
    	SPECC, 
    	SPECE2D, 
    	SPECE3D
    } ObjectType;

  private:
    ObjectType _type;

  public:

    PPPObjectIO(void): _type(UNDEF) {
      setObjectName(PPPOBJECTIO_NAME);
      };

    inline ObjectType getType(void) { return _type; };

    PPPBaseObject *read(const char *aName) {
      PPPBaseObject *aDest=NULL;
      _type = read_binheader(aName);
      switch(_type)
        {
        case AXIS:    aDest = new PPPAxis; break;
        case SIGD:    aDest = new PPPSignalContainer<double>; break;
        case SIGC:    aDest = new PPPSignalContainer<PPPcomplex>; break;
        case SIGE2D:  aDest = new PPPSignalContainer<PPPellipse2D>; break;
        case SIGE3D:  aDest = new PPPSignalContainer<PPPellipse3D>; break;
        case SPECD:   aDest = new PPPSpectrContainer<double>; break;
        case SPECC:   aDest = new PPPSpectrContainer<PPPcomplex>; break;
        case SPECE2D: aDest = new PPPSpectrContainer<PPPellipse2D>; break;
        case SPECE3D: aDest = new PPPSpectrContainer<PPPellipse3D>; break;
        }
      if(aDest != NULL)
        {
        FILE *infile = fopen(aName, "rb");
        if(infile == NULL)
          PPPBaseObject :: onError(FILE_ERROPEN+string(aName));
        aDest->fread(infile);
        fclose(infile);
        }
      return aDest;
      };

    ObjectType read_binheader(const char *aName) {
      FILE *infile = fopen(aName, "rb");
      if(infile == NULL)
        PPPBaseObject :: onError(FILE_ERROPEN+string(aName));
      char name[6];
      unsigned size = 0;
      std :: fread((void*)&name,sizeof(char),5,infile);
      std :: fread((void*)&size,sizeof(size),1,infile);
      name[5] = 0;
      fclose(infile);
      if(strcmp(name,PPPAXIS_OBJVER) == 0) return AXIS;
      if(strcmp(name,PPPSIGNALCONTAINER_OBJVER) == 0 && size==sizeof(double)) return SIGD;
      if(strcmp(name,PPPSIGNALCONTAINER_OBJVER) == 0 && size==sizeof(PPPcomplex)) return SIGC;
      if(strcmp(name,PPPSIGNALCONTAINER_OBJVER) == 0 && size==sizeof(PPPellipse2D)) return SIGE2D;
      if(strcmp(name,PPPSIGNALCONTAINER_OBJVER) == 0 && size==sizeof(PPPellipse3D)) return SIGE3D;
      if(strcmp(name,PPPSPECTRCONTAINER_OBJVER) == 0 && size==sizeof(double)) return SPECD;
      if(strcmp(name,PPPSPECTRCONTAINER_OBJVER) == 0 && size==sizeof(PPPcomplex)) return SPECC;
      if(strcmp(name,PPPSPECTRCONTAINER_OBJVER) == 0 && size==sizeof(PPPellipse2D)) return SPECE2D;
      if(strcmp(name,PPPSPECTRCONTAINER_OBJVER) == 0 && size==sizeof(PPPellipse3D)) return SPECE3D;
      return UNDEF;
      };

  }; // end of object

#endif

