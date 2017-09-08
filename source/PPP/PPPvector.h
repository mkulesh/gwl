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
#ifndef _PPPVECTORCONTAINER
#define _PPPVECTORCONTAINER

#define PPPVECTORCONTAINER_OBJVER  "VC1.3"
#define PPPVECTORCONTAINER_NAME    "vector"
#define PPPVECTORCONTAINER_ERRRES  "unable to resize reference object"
#define PPPVECTORCONTAINER_ERRPAR  "unable to destroy parent vector before all childs are not destroyed"
#define PPPVECTORCONTAINER_ERRUNL  "unable to unlink direct object"
#define PPPVECTORCONTAINER_ERRUNLR "unable to unlink reference object before all childs are not destroyed"

/************************************************************************
 *  PPPVectorContainer
 ***********************************************************************/
template<class AType>
class PPPVectorContainer : public PPPBaseTemplate<AType>
{
    
    friend class PPPMatrixContainer<AType> ;

private:
    
    unsigned _childscount;   // the number of links to *this
    unsigned _size;		// the data is _data[0], _data[_offset], ..., _data[(_size-1)*_offset].
    int _offset;	// the offset
    unsigned *_parent;	// the pointer to the linked object
    AType *_data;		// the actual data
    
public:
    
    PPPVectorContainer ()
    {
        _setdefault();
    }
    
    PPPVectorContainer (const unsigned aSize)
    {
        _setdefault();
        resize(aSize);
    }
    
    PPPVectorContainer (const PPPVectorContainer & aSour)
    {
        _setdefault();
        assign(aSour);
    }
    
    ~PPPVectorContainer (void)
    {
        if (_parent == NULL)
            realloc(0);
        else
            unlink();
    }
    
    /**
     *  this is the part, with is depended of type of "_data" property
     */
    inline unsigned size (void) const
    {
        return _size;
    }
    
    inline AType & operator [] (const unsigned aInd) const
    {
#ifdef PPPCONF_CHECKINDEX
        _checkindex(aInd);
#endif
        return *(_data + _offset * aInd);
    }
    
    void realloc (const unsigned aSize)
    { // reallocate to size aSize
        if (size() == aSize) return;  // same size do nothing
        if (_parent != NULL) PPPBaseTemplate<AType>::onError(PPPVECTORCONTAINER_ERRRES);
        if (_childscount > 0) PPPBaseTemplate<AType>::onError(PPPVECTORCONTAINER_ERRPAR);
        if (size() != 0)
        {
            delete _data;
            _data = NULL;
            _size = 0;
        }
        if (aSize != 0)
        {
            if ((_data = new AType[aSize]) == NULL) // if allocation fails
            PPPBaseTemplate<AType>::onError(MEM_ERRALLOC + string("realloc"));
            _size = aSize;
        }
        return;
    }
    
    void resize (const unsigned aSize)
    {
        if (size() == aSize) return;
        if (_parent != NULL) PPPBaseTemplate<AType>::onError(PPPVECTORCONTAINER_ERRRES);
        if (_childscount > 0) PPPBaseTemplate<AType>::onError(PPPVECTORCONTAINER_ERRPAR);
        if (aSize == 0 && size() != 0)
        {
            delete _data;
            _data = NULL;
            _size = 0;
            return;
        }
        else if (aSize != 0 && size() == 0)
        {
            _data = new AType[aSize];
            if (_data == NULL) PPPBaseTemplate<AType>::onError(MEM_ERRALLOC + string("resize"));
            for (unsigned i = 0; i < aSize; i++)
                _data[i] = PPPBaseTemplate<AType>::nullValue();
            _size = aSize;
            return;
        }
        else if (aSize != 0 && size() != 0)
        {
            AType *_data1;
            if ((_data1 = new AType[aSize]) == NULL) PPPBaseTemplate<AType>::onError(
                    MEM_ERRALLOC + string("resize"));
            for (unsigned i = 0; i < aSize; i++)
            {
                if (i < size())
                    _data1[i] = _data[i];
                else
                    _data1[i] = PPPBaseTemplate<AType>::nullValue();
            }
            delete _data;
            _data = _data1;
            _size = aSize;
            return;
        }
    }
    
    /**
     *  iterators support
     */
    AType * begin () const
    {
        return _data;
    }
    
    /**
     *  generates a view. *this is a view of aSour
     *  after linking of a to b we have
     *  a[0] <--> b[aStart]
     *  a[1] <--> b[aStart + aOffset]
     *  a[2] <--> b[aStart + 2*aOffset]
     *  ...
     *  a[aSize-1] <--> b[aStart + (aSize-1)*aOffset]
     *  
     *  if aStart + (aSize-1)*aOffset is out of range, this generates an error
     */
    void link (PPPVectorContainer & aSour, 		 // the source to which we link
            unsigned aStart = 0, 			 // the first index
            unsigned aSize = 0, 			 // how many elements. if 0 to end
            int aOffset = 1)
    {				 // which offset
        unsigned size = (aSize == 0) ? (aSour.size() - aStart) : aSize;
        int index_beg = aStart;
        int index_end = index_beg + aOffset * (int) (size - 1);
        aSour._checkindex(index_beg);
        aSour._checkindex(index_end);
        _link(aSour._data + aSour._offset * index_beg, &(aSour._childscount), size,
              aOffset * aSour._offset);
        return;
    }
    
    void unlink (void)
    {
        if (_parent == NULL) PPPBaseTemplate<AType>::onError(PPPVECTORCONTAINER_ERRUNL);
        if (_childscount > 0) PPPBaseTemplate<AType>::onError(PPPVECTORCONTAINER_ERRUNLR);
        (*_parent)--;
        _parent = NULL;
        _size = 0;
        _data = NULL;
        _offset = 1;
        return;
    }
    
    void operator = (const PPPVectorContainer<AType> & aSour)
    {
        assign(aSour);
    }
    
    void assign (const AType aVal)
    {
        for (unsigned i = 0; i < size(); i++)
            (*this)[i] = aVal;
        return;
    }
    
    void assign (const PPPVectorContainer & aSour, unsigned aStart = 0, unsigned aSize = 0,
                 int aOffset = 1)
    {
        PPPBaseObject::setObjectName(aSour.getObjectName());
        unsigned newsize = aSize;
        if (newsize == 0 && aOffset == 1) newsize = aSour.size() - aStart;
        if (newsize == 0 && aOffset != 1) newsize = (unsigned) ceil(
                (double) (aSour.size() - aStart) / ((double) aOffset));
        realloc(newsize);
        for (unsigned i = 0; i < size(); i++)
            (*this)[i] = aSour[aStart + i * aOffset];
    }
    
    void push_back (AType aVal)
    {
        resize(size() + 1);
        (*this)[size() - 1] = aVal;
    }
    
    template<class DestType, class TransType>
    void compTransform (PPPVectorContainer<DestType> &aDest, TransType aTrans)
    {
        aDest.realloc(size());
        aDest.setObjectName(aTrans.getComponentName());
        for (unsigned i = 0; i < size(); i++)
            aDest[i] = aTrans((*this)[i]);
        return;
    }
    
    const char *getInfo (void)
    {
        strstream str;
        str << PPPBaseTemplate<AType>::getObjectName() << ": "
        << PPPBaseTemplate<AType>::getTypeName() << "[" << size() << "]";
        if (_childscount > 0) str << ", Childs=" << _childscount;
        if (_offset != 1) str << ", offset=" << _offset;
        if (_parent != NULL) str << ", Reference object";
        str << ends;
        PPPBaseTemplate<AType>::onNotation(str.str());
        return PPPBaseTemplate<AType>::getNotation();
    }
    
    /*
     * File reading and writing procedures
     */
    void read (const string &aName)
    {
        PPPBaseTemplate<AType>::onMessage(FILE_READASC + aName);
        vector < string > lines;
        PPPBaseTemplate<AType>::readLinesFromFile(aName, lines);
        realloc(lines.size());
        for (unsigned i = 0; i < lines.size(); i++)
        {
            istringstream stringline(lines[i]);
            stringline >> (*this)[i];
        }
        return;
    }
    
    template<class TransType>
    void write (const string &aName, TransType &aTrans)
    {
        strstream str;
        str << FILE_WRITEASC << aName << "[" << size() << "], " << aTrans.getComponentName()
        << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        remove(aName.c_str());
        fstream outfile(aName.c_str(), ios_base::out);
        if (!outfile) PPPBaseTemplate<AType>::onError(FILE_ERROPEN + aName);
        for (unsigned i = 0; i < size(); i++)
        {
            outfile.precision(PPPBaseTemplate<AType>::getOutPrecision());
            if (!PPPBaseTemplate<AType>::isInteger()) outfile << showpos << scientific;
            outfile << aTrans((*this)[i]) << endl;
        }
        outfile.close();
        return;
    }
    
    void write (const string &aName)
    {
        PPPNullTransform<AType> trans;
        write(aName, trans);
    }
    
    /*
     * Conversion of string to PPPVectorContainer
     * string must be contain the _data with format: "{val1, val2, val3, ...}"
     */
    void strToVector (const string &aSour)
    {
        char *uk0 = (char *) aSour.c_str();
        char *uk1, *uk2, *str, *strval;
        int size, flag, vallen;
        AType val;
        if (strlen(uk0) == 0) return;
        uk1 = (strchr(uk0, '{') == NULL) ? uk1 = &uk0[0] : strchr(uk0, '{') + 1;
        uk2 = (strrchr(uk0, '}') == NULL) ? uk2 = &uk0[strlen(uk0)] : strrchr(uk0, '}');
        vallen = (int) (uk2 - uk1);
        if (vallen <= 0)
        {
            realloc(0);
            return;
        }
        str = new char[vallen + 2];
        strval = new char[vallen + 2];
        strncpy(str, uk1, vallen);
        str[vallen] = 0;
        for (size = 0, uk1 = str; *uk1 != 0; uk1++)
            if (*uk1 == ',') size++;
        size++;
        resize(size);
        uk1 = uk2 = &str[0];
        size = 0;
        while (1)
        {
            for (flag = 0;; uk2++)
            {
                if (*uk2 == 0) break;
                if (*uk2 != ' ' && *uk2 != ',') flag = 1;
                if (uk2 > uk1 && *uk2 == ',' && *(uk2 - 1) != ',') break;
            }
            if (flag == 1)
            {
                vallen = (int) (uk2 - uk1);
                strncpy(strval, uk1, vallen);
                strval[vallen] = 0;
                if (PPPBaseTemplate<AType>::isInteger()) sscanf(strval, "%d", &val); // TODO check format
                if (PPPBaseTemplate<AType>::isReal()) sscanf(strval, "%lf", &val); // TODO check format
                (*this)[size] = val;
                size++;
                if (*uk2 == ',') uk2++;
            }
            uk1 = uk2;
            if (*uk1 == 0) break;
        }
        delete str;
        delete strval;
    }
    
    /*
     * Conversion of PPPVectorContainer to string
     * string will be contain the _data with format: "{val1, val2, val3, ...}"
     */
    template<class TransType>
    const char *vectorToStr (TransType &aTrans)
    {
        strstream str;
        str << "{";
        for (unsigned i = 0; i < size(); i++)
        {
            str.precision(PPPBaseTemplate<AType>::getOutPrecision());
            if (!PPPBaseTemplate<AType>::isInteger()) str << showpos << scientific;
            str << aTrans((*this)[i]);
            if (i != size() - 1) str << ", ";
        }
        str << "}" << ends;
        PPPBaseTemplate<AType>::onNotation(str.str());
        return PPPBaseTemplate<AType>::getNotation();
    }
    
    const char *vectorToStr (void)
    {
        PPPNullTransform<AType> trans;
        return vectorToStr(trans);
    }
    
    /*
     * Math procedures
     */
    AType getSumm (void)
    {
        AType val = PPPBaseTemplate<AType>::nullValue();
        for (unsigned i = 0; i < size(); i++)
            val = val + (*this)[i];
        return val;
    }
    
    template<class TransType>
    AType getMinValue (TransType &aTrans)
    {
        double res = 0.0, nres;
        AType val;
        for (unsigned i = 0; i < size(); i++)
        {
            nres = aTrans((*this)[i]);
            if (i == 0 || nres < res)
            {
                res = nres;
                val = (*this)[i];
            }
        }
        return val;
    }
    
    AType getMinValue (void)
    {
        PPPNullTransform<AType> trans;
        return getMinValue(trans);
    }
    
    template<class TransType>
    AType getMaxValue (TransType &aTrans)
    {
        double res = 0.0, nres;
        AType val;
        for (unsigned i = 0; i < size(); i++)
        {
            nres = aTrans((*this)[i]);
            if (i == 0 || nres > res)
            {
                res = nres;
                val = (*this)[i];
            }
        }
        return val;
    }
    
    AType getMaxValue (void)
    {
        PPPNullTransform<AType> trans;
        return getMaxValue(trans);
    }
    
    void toPowerOfTwo (unsigned aBase = 0)
    {
        PPPMathFunc math;
        if (aBase == 0)
        {
            if (math.isPowerOfTwo(size())) return;
            aBase = 1;
            for (;;)
            {
                aBase *= 2;
                if (aBase >= size()) break;
            }
        }
        else
        {
            if (!math.isPowerOfTwo(aBase)) PPPBaseTemplate<AType>::onError(
                    ARG_POW2 + string("ToPowerOfTwo"));
        }
        resize(aBase);
        return;
    }
    
    /**
     *  file stream operations
     */
    void fwrite (FILE *stream)
    {
        PPPBaseObject::fwrite_streaminfo(stream, PPPBaseTemplate<AType>::getObjectVer(),
                                         sizeof(AType));
        unsigned newsize = size();
        std::fwrite((void*) &newsize, sizeof(newsize), 1, stream);
        for (unsigned i = 0; i < size(); i++)
        {
            AType tmp = (*this)[i];
            std::fwrite((void*) &tmp, sizeof(tmp), 1, stream);
        }
        PPPBaseObject::fwrite(stream);
    }
    
    void fread (FILE *stream)
    {
        PPPBaseObject::fread_streaminfo(stream, PPPBaseTemplate<AType>::getObjectVer(),
                                        sizeof(AType));
        unsigned newsize;
        std::fread((void*) &newsize, sizeof(newsize), 1, stream);
        realloc(newsize);
        AType tmp;
        for (unsigned i = 0; i < size(); i++)
        {
            std::fread((void*) &tmp, sizeof(AType), 1, stream);
            (*this)[i] = tmp;
        }
        PPPBaseObject::fread(stream);
    }
    
private:
    
    void _checkindex (const unsigned aInd) const
    {
        if (aInd >= size())
        {
            strstream str;
            str << MEM_ERRINDEX << ", size = " << size() << ", index = " << aInd << ends;
            PPPBaseTemplate<AType>::onError(str.str());
        }
    }
    
    void _setdefault (void)
    {
        PPPBaseTemplate<AType>::setObjectVer(PPPVECTORCONTAINER_OBJVER);
        PPPBaseTemplate<AType>::setObjectName(PPPVECTORCONTAINER_NAME);
        _childscount = 0;
        _parent = NULL;
        _data = NULL;
        _size = 0;
        _offset = 1;
        return;
    }
    
    /**
     *  primitive linking method
     */
    void _link (AType *aData, unsigned *aParent, unsigned aSize, int aOffset)
    {
        if (_parent == NULL)
            realloc(0);
        else
            unlink();
        _parent = aParent;
        _size = aSize;
        _data = aData;
        _offset = aOffset;
        (*_parent)++;
        return;
    }
    
};

#endif

