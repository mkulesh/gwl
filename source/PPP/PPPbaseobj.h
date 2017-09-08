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

#ifndef _PPPBASEOBJECT
#define _PPPBASEOBJECT

#define PPPBASEOBJECT_OBJVER     "BO1.3"

#include <cstring>

/************************************************************************
 * PPPBaseObject
 ************************************************************************/
class PPPBaseObject
  {
  public:

    typedef enum {NMnone,NMappend,NMinsert} NotatType;

  private:
    string      _objectver;
    string      _objectname;
    string      _notation;
    bool        _showMessage;
    bool        _showError;
    bool        _showProgress;
    NotatType   _showNotation;
    int         _outPrecision;

  public:

    PPPBaseObject(void) {
      _objectver = PPPBASEOBJECT_OBJVER;
      _objectname = "";
      clearNotation();
      _showMessage = true;
      _showError = true;
      _showProgress = true;
      _showNotation = NMinsert;
      _outPrecision = 6;
      };

    void setObjectVer(const string &aName) {
      _objectver = aName;
      };

    const char *getObjectVer(void) {
      return _objectver.c_str();
      };

    void setObjectName(const string &aName) {
      _objectname = aName;
      };

    const char *getObjectName(void) const {
      return _objectname.c_str();
      };

    int getOutPrecision(void) const {
      return _outPrecision;
      };

    void setOutPrecision(int aPrec) {
      _outPrecision = aPrec;
      };

    bool getShowMessage(void) const {
      return _showMessage;
      };

    void setShowMessage(bool aMess) {
      _showMessage = aMess;
      };

    void onMessage(const string &aNotation) const {
      if(!_showMessage) return;
      string str;
      if(_objectname.empty())
        str = aNotation;
      else
        str = _objectname + ": " + aNotation;
      ConApplication.onMessage(str);
      };

    bool getShowError(void) const {
      return _showError;
      };

    void setShowError(bool aMess) {
      _showError = aMess;
      };

    void onError(const string &aNotation) const {
      if(!_showError) return;
      string str;
      if(_objectname.empty())
        str = aNotation;
      else
        str = _objectname + ": " + aNotation;
      ConApplication.onError(str);
      };

    NotatType getShowNotation(void) const {
      return _showNotation;
      };

    void setShowNotation(NotatType aMess) {
      _showNotation = aMess;
      };

    void onNotation(const string &aNotation) {
      if(_showNotation == NMnone) return;
      if(_showNotation == NMinsert)
        {
        clearNotation();
        _notation.assign(aNotation);
        }
      if(_showNotation == NMappend)
        {
        _notation.append(aNotation);
        }
      };

    void clearNotation() {
      _notation = "";
      };

    bool isNotation() {
      return (_notation.size()>0);
      };

    const char *getNotation(void) const {
      return _notation.c_str();
      };

    bool getShowProgress(void) const {
      return _showProgress;
      };

    void setShowProgress(bool aMess) {
      _showProgress = aMess;
      };

    void onProgress(const int aPercent, const string &aNotation) const {
      if(!_showProgress) return;
      ConApplication.onProgress(aPercent, aNotation);
      };

    void readLinesFromFile(const string &aName, vector<string> &aLines) {
      fstream infile(aName.c_str(),ios_base::in);
      if (!infile) onError(FILE_ERROPEN+aName);
      string line;
      aLines.clear();
      int pos;
      while (getline(infile,line))
        {
        for(pos=line.size()-1; pos>=0; pos--)
          if(line[pos]!='\0' && line[pos]!=' ' && line[pos]!=0x0D) break;
        if(pos<0) continue;
        line.resize(pos+1);
        aLines.push_back(line);
        }
      infile.close();
      return;
      };

    /**
     *  file stream operations
     */
    void fwrite(FILE *stream) {
      unsigned newsize = _objectname.size();
      std :: fwrite((void*)&newsize,sizeof(newsize),1,stream);
      std :: fwrite((void*)_objectname.c_str(),sizeof(char),newsize,stream);
      };

    virtual void fread(FILE *stream) {
      unsigned newsize;
      std :: fread((void*)&newsize,sizeof(newsize),1,stream);
      char *name = new char[newsize+1];
      std :: fread((void*)name,sizeof(char),newsize,stream);
      name[newsize] = 0;
      setObjectName(name);
      delete name;
      };

    void fwrite_streaminfo(FILE *stream, const char *aName, unsigned aSize) {
      std :: fwrite((void*)aName, sizeof(char), 5, stream);
      std :: fwrite((void*)&aSize,sizeof(aSize),1, stream);
      };

    void fread_streaminfo(FILE *stream, const char *aName, unsigned aSize) {
      char name[6];
      unsigned size;
      std :: fread((void*)&name,sizeof(char),5,stream);
      std :: fread((void*)&size,sizeof(size),1,stream);
      name[5] = 0;
      if(strcmp(name,aName) != 0 || size != aSize)
        {
        stringstream str;
        str << FILE_ERRFORM << "fread()\n  " << FILE_ERREXP << aName << " (" << aSize << "b)" << ends;
        onError(str.str());
        }
      };

  }; // end of object

#endif
