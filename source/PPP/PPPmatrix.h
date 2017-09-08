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
#ifndef _PPPMATRIXCONTAINER
#define _PPPMATRIXCONTAINER

#define PPPMATRIXCONTAINER_OBJVER     "MT1.3"
#define PPPMATRIXCONTAINER_NAME       "matrix"
#define PPPMATRIXCONTAINER_ERRPAR     "unable to destroy parent matrix before all childs vectors are not destroyed"
#define PPPMATRIXCONTAINER_HOR        "horisontal direction"
#define PPPMATRIXCONTAINER_VERT       "vertical direction"

/************************************************************************
 *  PPPMatrixContainer
 ***********************************************************************/
template<class AType> class PPPMatrixContainer : public PPPBaseTemplate<AType>
  {
  private:

    unsigned    _childscount;
    unsigned    _rowscount;
    unsigned    _colscount;
    AType       *_data;

  public:

    typedef enum {MDhor, MDvert} MatrixDir;

    PPPMatrixContainer(void) {
      _setdefault();
      };

    PPPMatrixContainer(const unsigned aRowscount, const unsigned aColscount) {
      _setdefault();
      resize(aRowscount, aColscount);
      };

    PPPMatrixContainer(const PPPMatrixContainer &aSour) {
      _setdefault();
      assign(aSour);
      };

    ~PPPMatrixContainer(void) {
      realloc(0,0);
      };

    /**
     *  this is the part, with is depended of type of "_data" property
     */
    inline unsigned cols() const {
      return _colscount;
      };

    inline unsigned rows() const {
      return _rowscount;
      };

    /**
    *  acess to elemetns of matrix
    */
    AType const * begin () {
      return &_data[0];
      };

    AType const * end () const {
      return &_data[rows()*cols()];
      };

    inline AType & operator() (const unsigned aRow, const unsigned aCol) const {
      #ifdef PPPCONF_CHECKINDEX
      _checkindex(aRow,aCol);
      #endif
      return _data[aRow*cols()+aCol];
      };

    inline void getRow (PPPVectorContainer <AType> &aDest,  unsigned aRow ) {
      aDest._link(_data + aRow*cols(), &_childscount, cols(), 1);
      };

    inline void getCol (PPPVectorContainer <AType> &aDest,  unsigned aCol ) {
      aDest._link(_data + aCol, &_childscount, rows(), cols());
      };


    /**
    *  resize and realloc
    */
    void realloc(const unsigned aRowscount, const unsigned aColscount) {
      if(rows() == aRowscount && cols() == aColscount) return;
      if(_childscount > 0) PPPBaseTemplate<AType>::onError(PPPMATRIXCONTAINER_ERRPAR);
      unsigned aSize = aRowscount*aColscount;
      unsigned count = rows()*cols();
      if(count != 0)
        {
        delete _data;
        _rowscount = 0;
        _colscount = 0;
        }
      if(aSize != 0)
        {
        _data = new AType[aSize];
        if(_data == NULL) PPPBaseTemplate<AType>::onError(MEM_ERRALLOC+string("realloc"));
        _rowscount = aRowscount;
        _colscount = aColscount;
        }
      return;
      };

    void resize(const unsigned aRowscount, const unsigned aColscount) {
      if(rows() == aRowscount && cols() == aColscount) return;
      if(_childscount > 0) PPPBaseTemplate<AType>::onError(PPPMATRIXCONTAINER_ERRPAR);
      unsigned aSize = aRowscount*aColscount;
      unsigned count = rows()*cols();
      unsigned i,j;
      if(aSize == 0 && count != 0)
        {
        delete _data;
        _rowscount = 0;
        _colscount = 0;
        return;
        }
      else if(aSize != 0 && count == 0)
        {
        if((_data = new AType[aSize]) == NULL) PPPBaseTemplate<AType>::onError(MEM_ERRALLOC+string("resize"));
        for(i=0; i<aSize; i++) _data[i] = PPPBaseTemplate<AType>::nullValue();
        _rowscount = aRowscount;
        _colscount = aColscount;
        return;
        }
      else if(aSize != 0 && count != 0)
        {
        AType *values1;
        if((values1 = new AType[aSize]) == NULL) PPPBaseTemplate<AType>::onError(MEM_ERRALLOC+string("resize"));
        for(i=0; i<aRowscount; i++)  for(j=0; j<aColscount; j++)
          {
          if(i<rows() && j<cols()) values1[i*aColscount+j] = _data[i*cols()+j];
          else values1[i*aColscount+j] = PPPBaseTemplate<AType>::nullValue();
          }
        delete _data;
        _data = values1;
        _rowscount = aRowscount;
        _colscount = aColscount;
        return;
        }
      };


    /**
     *  this is the part, with is independed from type of "_data" property
     */
    void assign(const AType aVal) {
      unsigned i,j;
      for(i=0; i<rows(); i++) for(j=0; j<cols(); j++) (*this)(i,j) = aVal;
      };

    void assign(const PPPMatrixContainer &aSour,
      unsigned aRowStart=0, unsigned aRows=0, int aRowOffset=1,
      unsigned aColStart=0, unsigned aCols=0, int aColOffset=1) {
      PPPBaseObject::setObjectName(aSour.getObjectName());
      if(aRows != 0 && aCols != 0)
        realloc(aRows,aCols);
      else
        realloc(aSour.rows(),aSour.cols());
      unsigned i,j;
      for(i=0; i<rows(); i++)
        for(j=0; j<cols(); j++)
          (*this)(i,j) = aSour(aRowStart+i*aRowOffset,aColStart+j*aColOffset);
      };


    template <class TransType>
    void compTransform(PPPMatrixContainer<double> &aDest, TransType aTrans) {
      aDest.realloc(rows(), cols());
      aDest.setObjectName(aTrans.getComponentName());
      unsigned i,j;
      for(i=0; i<rows(); i++) for(j=0; j<cols(); j++)
        {
        aDest(i,j) = aTrans((*this)(i,j));
        }
      return;
      };


    /*
     * File reading and writing procedures
     */
    void read(string const &aName) {
      PPPBaseTemplate<AType>::onMessage(FILE_READASC+aName);
      vector<string> lines;
      PPPBaseTemplate<AType>::readLinesFromFile(aName,lines);
      unsigned j,aColscount=0;
      AType v;
      for (unsigned i=0; i<lines.size(); ++i)
        {
        istringstream stringline(lines[i]);
        j = 0;
        while (! stringline.eof())
          {
          stringline >> v;
          ++j;
          }
        aColscount = max(aColscount, j);
        }
      realloc(lines.size(),aColscount);
      for (unsigned i=0; i<lines.size(); ++i)
        {
        istringstream stringline(lines[i]);
        for(j=0; j<aColscount; j++)
          if( !(stringline >> (*this)(i,j)) ) break;
        }
      return;
      };

    template <class TransType>
    void write(string const &aName, const unsigned aRowscount, const unsigned aColscount, TransType &aTrans,
      MatrixDir const aDirect=MDvert, int indj1=-1, int indj2=-1, int indi1=-1, int indi2=-1) {
      strstream str;
      remove(aName.c_str());
      fstream outfile(aName.c_str(),ios_base::out);
      if(!outfile) PPPBaseTemplate<AType>::onError(FILE_ERROPEN+aName);
      int i,j,k1,k2,n1=0,n2=0;
      if(indj1<0) indj1=0;
      if(indj2<0) indj2=rows();
      if(indi1<0) indi1=0;
      if(indi2<0) indi2=cols();
      if(aRowscount != 0) k1 = (indj2-indj1)/aRowscount;
      if(aColscount != 0) k2 = (indi2-indi1)/aColscount;
      if(k1==0) k1=1;
      if(k2==0) k2=1;
      if(aDirect == MDhor)
        {
        for(j=indi1,n2=0; j<indi2; j+=k2, n2++)
          {
          for(i=indj1,n1=0; i<indj2; i+=k1, n1++)
            {
            outfile.precision(PPPBaseTemplate<AType>::getOutPrecision());
            if(!PPPBaseTemplate<AType>::isInteger()) outfile << showpos << scientific;
            outfile << aTrans((*this)(i,j)) << " ";
            }
          outfile << endl;
          }
        }
      if(aDirect == MDvert)
        {
        for(i=indj1,n1=0; i<indj2; i+=k1, n1++)
          {
          for(j=indi1,n2=0; j<indi2; j+=k2, n2++)
            {
            outfile.precision(PPPBaseTemplate<AType>::getOutPrecision());
            if(!PPPBaseTemplate<AType>::isInteger()) outfile << showpos << scientific;
            outfile << aTrans((*this)(i,j)) << " ";
            }
          outfile << endl;
          }
        }
      str << FILE_WRITEASC << " " << aName <<" ["<<n1<<"]["<<n2<<"], "<<aTrans.getComponentName()<<ends;
      outfile.close();
      PPPBaseTemplate<AType>::onMessage(str.str());
      return;
      };

    void write(const string &aName, MatrixDir const aDirect=MDvert) {
      PPPNullTransform<AType> trans;
      write(aName,rows(),cols(),trans,aDirect);
      };

    /*
     * Conversion of string to PPPMatrixContainer
     */
    void strToMatrix(const string &aSour) {
      char *uk0 = (char *)aSour.c_str();
      char *uk1,*uk2,*str,*strval;
      unsigned i,j,flag1,flag2,vallen;
      if(strlen(uk0) == 0) return;
      uk1 = (strchr(uk0,'{')==NULL)? uk1=&uk0[0] : strchr(uk0,'{')+1;
      uk2 = (strrchr(uk0,'}')==NULL)? uk2=&uk0[strlen(uk0)] : strrchr(uk0,'}');
      vallen = (int)(uk2-uk1);
      str = new char[vallen+1];
      strval = new char[vallen+1];
      strncpy(str,uk1,vallen);
      str[vallen] = 0;
      for(i=j=0, uk1=str; *uk1 != 0; uk1++)
        {
        if(*uk1 == '{') i++;
        if(*uk1 == ',') j++;
        }
      j = (j-(i-1))/i + 1;
      resize(i,j);
      uk1 = uk2 = &str[0];
      i = j = 0;
      PPPVectorContainer<AType> vecval;
      while(1)
        {
        for(flag1=flag2=0; ;uk2++)
          {
          if(*uk2 == 0) break;
          if(*uk2 != ' ' && *uk2 != ',' && *uk2 != '{'&& *uk2 != '}') flag1 = 1;
          if(*uk2 == '}') flag2 = 1;
          if(uk2>uk1 && *uk2 == ',' && flag2 == 1) break;
          }
        if(flag1 == 1)
          {
          vallen = (int)(uk2-uk1);
          strncpy(strval,uk1,vallen);
          strval[vallen] = 0;
          vecval.strToVector(strval);
          for(j=0; j<cols(); j++) (*this)(i,j) = vecval[j];
          if(*uk2 == ',') uk2++;
          }
        uk1 = uk2;
        if(*uk1 == 0) break;
        i++;
        }
      delete str;
      delete strval;
      return;
      };

    /*
     * Conversion of PPPMatrixContainer to string
    */
    template <class TransType>
    const char *matrixToStr(TransType &aTrans) {
      string aDest = "{";
      PPPVectorContainer<AType> vec(cols());
      unsigned i,j;
      for(j=0; j<rows(); j++)
        {
        for(i=0; i<cols(); i++) vec[i] = (*this)(j,i);
        aDest.append(vec.vectorToStr(aTrans));
        if(j!=rows()-1) aDest.append(",");
        }
      aDest.append("}");
      PPPBaseTemplate<AType>::onNotation(aDest);
      return PPPBaseTemplate<AType>::getNotation();
      };

    const char *matrixToStr(void) {
      PPPNullTransform<AType> trans;
      return matrixToStr(trans);
      };

    /*
     * Math procedures
     */
    template <class TransType>
    AType getMinValue (TransType &aTrans) {
      return *min_element ( begin(), end(), PPPCompareLess<AType,TransType> ( aTrans )  );
      };

    AType getMinValue (void) {
      return *min_element ( begin(), end() );
      };

    template <class TransType>
    AType getMaxValue (TransType const &aTrans) {
      return *max_element ( begin(), end(), PPPCompareLess<AType,TransType> ( aTrans ) ); //getExtremeValue ( aTrans, compare_greater < double >() );
      };

    AType getMaxValue (void) {
      return *max_element ( begin(), end() );
      };

    /**
     *  file stream operations
     */
    void fwrite(FILE *stream) {
      PPPBaseObject :: fwrite_streaminfo(stream, PPPBaseTemplate<AType>::getObjectVer(), sizeof(AType));
      unsigned newsize1 = rows(), newsize2 = cols();
      std :: fwrite((void*)&newsize1,sizeof(newsize1),1,stream);
      std :: fwrite((void*)&newsize2,sizeof(newsize2),1,stream);
      for (unsigned i=0; i<rows(); i++)
        for (unsigned j=0; j<cols(); j++)
          {
          AType tmp = (*this)(i,j);
          std :: fwrite((void*)&tmp,sizeof(tmp),1,stream);
          }
      PPPBaseObject :: fwrite(stream);
      };

    void fread(FILE *stream) {
      PPPBaseObject :: fread_streaminfo(stream, PPPBaseTemplate<AType>::getObjectVer(), sizeof(AType));
      unsigned newsize1,newsize2;
      std :: fread((void*)&newsize1,sizeof(newsize1),1,stream);
      std :: fread((void*)&newsize2,sizeof(newsize2),1,stream);
      realloc(newsize1,newsize2);
      AType tmp;
      for (unsigned i=0; i<rows(); i++)
        for (unsigned j=0; j<cols(); j++)
          {
          std :: fread((void*)&tmp,sizeof(AType),1,stream);
          (*this)(i,j) = tmp;
          }
      PPPBaseObject :: fread(stream);
      };

  private:

    void _checkindex(const unsigned aRow,const unsigned aCol) const {
      if(aRow>=rows() || aCol>=cols()) PPPBaseTemplate<AType>::onError(MEM_ERRINDEX);
      };

    void _setdefault(void) {
      PPPBaseTemplate<AType>::setObjectVer(PPPMATRIXCONTAINER_OBJVER);
      PPPBaseTemplate<AType>::setObjectName(PPPMATRIXCONTAINER_NAME);
      _rowscount = 0;
      _colscount = 0;
      _childscount = 0;
      return;
      };

  };  // end of object


#endif
