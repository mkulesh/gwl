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
#ifndef _PPPSPECTRCONTAINER
#define _PPPSPECTRCONTAINER

#define PPPSPECTRCONTAINER_OBJVER     "SP1.3"
#define PPPSPECTRCONTAINER_NAME       "spectrum"
#define PPPSPECTRCONTAINER_TIME       "time"
#define PPPSPECTRCONTAINER_FREQ       "frequency"
#define PPPSPECTRCONTAINER_CHANNAME   "spectrum channel"
#define PPPSPECTRCONTAINER_ERRINTERP  "undefined interpolation type in procedure: "

/************************************************************************
 *  PPPSpectrContainer
 ***********************************************************************/
template<class AType>
class PPPSpectrContainer : public PPPBaseTemplate<AType>
{
public:
    
    typedef enum
    {
        ITnone, ITlin
    } InterpType;

private:
    
    vector<PPPMatrixContainer<AType> > _data;
    PPPAxis _time;
    PPPAxis _freq;

public:
    
    PPPSpectrContainer (void)
    {
        PPPBaseObject::setObjectVer(PPPSPECTRCONTAINER_OBJVER);
        setNames(PPPSPECTRCONTAINER_NAME, PPPSPECTRCONTAINER_TIME, PPPSPECTRCONTAINER_FREQ);
        resize(0, 0, 0);
    }
    
    PPPSpectrContainer (const unsigned aVoices, const unsigned aPoints,
                        const unsigned aChannels = 1)
    {
        PPPBaseObject::setObjectVer(PPPSPECTRCONTAINER_OBJVER);
        setNames(PPPSPECTRCONTAINER_NAME, PPPSPECTRCONTAINER_TIME, PPPSPECTRCONTAINER_FREQ);
        resize(aVoices, aPoints, aChannels);
    }
    
    PPPSpectrContainer (PPPSpectrContainer &aSour)
    {
        PPPBaseObject::setObjectVer(PPPSPECTRCONTAINER_OBJVER);
        assign(aSour);
    }
    
    inline unsigned voices () const
    {
        return (_data.size() == 0) ? 0 : _data[0].rows();
    }
    
    inline unsigned points () const
    {
        return (_data.size() == 0) ? 0 : _data[0].cols();
    }
    
    inline unsigned channels () const
    {
        return _data.size();
    }
    
    inline PPPMatrixContainer<AType> & getChannel (const unsigned aChannel = 0)
    {
#ifdef PPPCONF_CHECKINDEX
        _checkindex(aChannel);
#endif
        return _data[aChannel];
    }
    
    inline AType& operator() (const unsigned aVoice, const unsigned aPoint,
                              const unsigned aChannel = 0) const
    {
#ifdef PPPCONF_CHECKINDEX
        _checkindex(aChannel);
#endif
        return (_data[aChannel])(aVoice, aPoint);
    }
    
    void realloc (const unsigned aVoices, const unsigned aPoints, const unsigned aChannels = 1)
    {
        _time.realloc(aPoints);
        _freq.realloc(aVoices);
        _data.resize(aChannels);
        for (unsigned i = 0; i < channels(); i++)
        {
            _data[i].realloc(aVoices, aPoints);
            _data[i].setObjectName(PPPSPECTRCONTAINER_CHANNAME);
        }
    }
    
    void resize (const unsigned aVoices, const unsigned aPoints, const unsigned aChannels = 1)
    {
        _time.resize(aPoints);
        _freq.resize(aVoices);
        _data.resize(aChannels);
        for (unsigned i = 0; i < channels(); i++)
        {
            _data[i].resize(aVoices, aPoints);
            _data[i].setObjectName(PPPSPECTRCONTAINER_CHANNAME);
        }
    }
    
    void setNames (const string &aN1, const string &aN2, const string &aN3)
    {
        PPPBaseObject::setObjectName(aN1);
        _time.setObjectName(aN2);
        _freq.setObjectName(aN3);
    }
    
    // time axis
    inline PPPAxis & getTime (void)
    {
        return _time;
    }
    
    inline double getTime (const unsigned aPoint)
    {
        return _time[aPoint];
    }
    
    inline void setTime (PPPAxis &aAxis)
    {
        _time.assign(aAxis);
    }
    
    // frequency axis
    inline PPPAxis & getFreq (void)
    {
        return _freq;
    }
    
    inline double getFreq (const unsigned aVoice)
    {
        return _freq[aVoice];
    }
    
    inline void setFreq (PPPAxis &aAxis)
    {
        _freq.assign(aAxis);
    }
    
    void getVoices (PPPSignalContainer<AType> &aDest, const unsigned aVoice)
    {
        aDest.prepare(points(), channels(), getTime(), "Spectr voices");
        for (unsigned j = 0; j < channels(); j++)
            for (unsigned k = 0; k < points(); k++)
                aDest(k, j) = (*this)(aVoice, k, j);
    }
    
    template<class AType1>
    void prepare (PPPSpectrContainer<AType1> & aSour)
    {
        setNames(aSour.getObjectName(), aSour.getTime().getObjectName(),
                 aSour.getFreq().getObjectName());
        realloc(aSour.voices(), aSour.points(), aSour.channels());
        setTime(aSour.getTime());
        setFreq(aSour.getFreq());
    }
    
    void prepare (const unsigned aVoices, const unsigned aPoints, const unsigned aChannels,
                  PPPAxis &aAxisTime, PPPAxis &aAxisFreq, const string &aN1)
    {
        setNames(aN1, getTime().getObjectName(), getFreq().getObjectName());
        realloc(aVoices, aPoints, aChannels);
        setTime(aAxisTime);
        setFreq(aAxisFreq);
    }
    
    void assign (PPPSpectrContainer & aSour, unsigned aVoiStart = 0, unsigned aVoices = 0,
                 int aVoiOffset = 1, unsigned aPoiStart = 0, unsigned aPoints = 0, int aPoiOffset =
                         1)
    {
        setNames(aSour.getObjectName(), aSour.getTime().getObjectName(),
                 aSour.getFreq().getObjectName());
        if (aVoices != 0 && aPoints != 0)
            realloc(aVoices, aPoints, aSour.channels());
        else
            realloc(aSour.voices(), aSour.points(), aSour.channels());
        _time.assign(aSour.getTime(), aPoiStart, aPoints, aPoiOffset);
        _freq.assign(aSour.getFreq(), aVoiStart, aVoices, aVoiOffset);
        for (unsigned i = 0; i < channels(); i++)
            getChannel(i).assign(aSour.getChannel(i), aVoiStart, aVoices, aVoiOffset, aPoiStart,
                                 aPoints, aPoiOffset);
        return;
    }
    
    const char *getInfo (void)
    {
        strstream aDest;
        if (channels() == 1)
            aDest << PPPBaseObject::getObjectName() << ": " << PPPBaseTemplate<AType>::getTypeName()
            << "[" << voices() << "][" << points() << "]" << endl;
        else
            aDest << PPPBaseObject::getObjectName() << ": " << PPPBaseTemplate<AType>::getTypeName()
            << "[" << channels() << "][" << voices() << "][" << points() << "]" << endl;
        aDest << "  " << getTime().getInfo() << endl;
        aDest << "  " << getFreq().getInfo() << ends;
        PPPBaseObject::onNotation(aDest.str());
        return PPPBaseObject::getNotation();
    }
    
    AType Get (double aY, double aX, const unsigned aChannel = 0, InterpType aC = ITlin)
    {
        if (aX < getTime().getMin() || aX > getTime().getMax() || aY < getFreq().getMin()
            || aY > getFreq().getMax()) return 0.0;
        unsigned i1 = getFreq().locateFloor(aY);
        unsigned j1 = getTime().locateFloor(aX);
        if (i1 >= voices() - 1 || j1 >= points() - 1) return 0.0;
        if (aC == ITnone)
            return (*this)(i1, j1, aChannel);
        else if (aC == ITlin)
        {
            double T = (aY - getFreq(i1)) / (getFreq(i1 + 1) - getFreq(i1));
            double U = (aX - getTime(j1)) / (getTime(j1 + 1) - getTime(j1));
            return ((1 - T) * (1 - U) * (*this)(i1, j1, aChannel) + T * (1 - U)
                                                                    * (*this)(i1 + 1, j1, aChannel)
                    + T * U * (*this)(i1 + 1, j1 + 1, aChannel)
                    + (1 - T) * U * (*this)(i1, j1 + 1, aChannel));
        }
        else
            PPPBaseObject::onError(PPPSPECTRCONTAINER_ERRINTERP + string("Get"));
        return 0.0;
    }
    
    template<class TransType>
    void write (string const &aName, TransType aTrans, const unsigned aChann = 0,
                const int aVoices = -1, const int aPoints = -1, const bool aRepName = false)
    {
        unsigned int vcount = (aVoices < 0) ? voices() : aVoices;
        unsigned int pcount = (aPoints < 0) ? points() : aPoints;
        string newName(aName);
        if (aRepName) newName.replace(newName.find("."), 1,
                                      string(aTrans.getComponentCode()) + ".");
        getChannel(aChann).write(newName, vcount, pcount, aTrans, PPPMatrixContainer<AType>::MDhor);
    }
    
    void write (string const &aName, const int aVoices = -1, const int aPoints = -1)
    {
        PPPNullTransform<AType> nulltr;
        for (unsigned i = 0; i < channels(); i++)
        {
            string newName(aName);
            if (channels() > 1)
            {
                char num[10];
                sprintf(num, "(%d)", i + 1);
                newName.replace(newName.find("."), 1, string(num) + ".");
            }
            write(newName, nulltr, i, aVoices, aPoints);
        }
    }
    
    void writegnuplot (string const &aName)
    {
        strstream str;
        str << FILE_WRITEASC << " " << aName << " xyz[" << voices() * points() << "]" << ends;
        PPPBaseTemplate<AType>::onMessage(str.str());
        remove(aName.c_str());
        fstream outfile(aName.c_str(), ios_base::out);
        if (!outfile) PPPBaseTemplate<AType>::onError(FILE_ERROPEN + aName);
        unsigned i, j, k;
        for (i = 0; i < voices(); i++)
        {
            for (j = 0; j < points(); j++)
            {
                outfile.precision(PPPBaseTemplate<AType>::getOutPrecision());
                if (!PPPBaseTemplate<AType>::isInteger()) outfile << showpos << scientific;
                outfile << getTime(j) << " " << getFreq(i);
                for (k = 0; k < channels(); k++)
                    outfile << " " << (*this)(i, j, k);
                outfile << endl;
            }
            outfile << endl;
        }
    }
    
    /**
     *  file stream operations
     */
    void fwrite (FILE *stream)
    {
        PPPBaseObject::fwrite_streaminfo(stream, PPPBaseObject::getObjectVer(), sizeof(AType));
        getTime().fwrite(stream);
        getFreq().fwrite(stream);
        unsigned newsize = channels();
        std::fwrite((void*) &newsize, sizeof(newsize), 1, stream);
        for (unsigned i = 0; i < channels(); i++)
            getChannel(i).fwrite(stream);
        PPPBaseObject::fwrite(stream);
    }
    
    void fread (FILE *stream)
    {
        PPPBaseObject::fread_streaminfo(stream, PPPBaseObject::getObjectVer(), sizeof(AType));
        getTime().fread(stream);
        getFreq().fread(stream);
        unsigned newsize;
        std::fread((void*) &newsize, sizeof(newsize), 1, stream);
        _data.resize(newsize);
        for (unsigned i = 0; i < channels(); i++)
            getChannel(i).fread(stream);
        PPPBaseObject::fread(stream);
    }
    
private:
    
    void _checkindex (const unsigned aChannel) const
    {
        if (aChannel >= channels()) PPPBaseObject::onError(MEM_ERRINDEX);
    }
    
};
// end of object

#endif

