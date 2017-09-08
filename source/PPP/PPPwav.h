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
#ifndef _PPPWAVIN
#define _PPPWAVIN

#define PPPWAVIN_NAME      "wav stream"
#define PPPWAVIN_ERRFOR    "incorrect format in wav file: "
#define PPPWAVIN_SAMPLE    "sample rate"
#define PPPWAVIN_LENGTH    "number of samples"
#define PPPWAVIN_BPS       "bits per sample"
#define PPPWAVIN_CHAN      "number of channels"

// header of wav file
typedef struct
{
    char rID[4];            // 'RIFF'
    long int rLen;
    char wID[4];            // 'WAVE'
    char fId[4];            // 'fmt '
    long int pcm_header_len;   // varies...
    short int wFormatTag;
    short int nChannels;      // 1,2 for stereo data is (l,r) pairs
    long int nSamplesPerSec;
    long int nAvgBytesPerSec;
    short int nBlockAlign;
    short int nBitsPerSample;
} WAV_HDR;

// header of wav file
typedef struct
{
    char dId[4];            // 'data' or 'fact'
    long int dLen;
} CHUNK_HDR;

/************************************************************************
 * PPPWavIn
 ************************************************************************/
class PPPWavIn : public PPPBaseObject
{
protected:
    
    int *_data;              // Voice data
    unsigned _bps;           // BitsPerSample
    unsigned _channels;      // Mono (_channels=1) or stereo (_channels=2)
    double _rate;            // SamplesPerSec [Hz]
    long int _points;        // length of _data
    
public:
    
    PPPWavIn (const string &aName) :
            _data(NULL),
            _points(0),
            _bps(0),
            _channels(0),
            _rate(0)
    {
        int i, wbuff_len, sflag;
        FILE *fw;
        unsigned int wstat;
        char obuff[80];
        WAV_HDR *wav;
        CHUNK_HDR *chk;
        short int *uptr;
        unsigned char *cptr;
        long int rmore;
        char *wbuff;
        
        setObjectName(PPPWAVIN_NAME);
        onMessage(FILE_READWAV + aName);
        
        // allocate wav header
        if ((wav = new WAV_HDR) == NULL) onError(MEM_ERRALLOC + string("PPPWavIn"));
        if ((chk = new CHUNK_HDR) == NULL) onError(MEM_ERRALLOC + string("PPPWavIn"));
        
        /* open wav file */
        fw = fopen(aName.c_str(), "rb");
        if (fw == NULL) onError(FILE_ERROPEN + aName);
        
        /* read riff/wav header */
        wstat = std::fread((void *) wav, sizeof(WAV_HDR), (size_t) 1, fw);
        if (wstat != 1) onError(PPPWAVIN_ERRFOR + string("HDR"));
        
        // check format of header 
        for (i = 0; i < 4; i++)
            obuff[i] = wav->rID[i];
        obuff[4] = 0;
        if (strcmp(obuff, "RIFF") != 0) onError(PPPWAVIN_ERRFOR + string("RIFF"));
        
        for (i = 0; i < 4; i++)
            obuff[i] = wav->wID[i];
        obuff[4] = 0;
        if (strcmp(obuff, "WAVE") != 0) onError(PPPWAVIN_ERRFOR + string("WAVE"));
        
        for (i = 0; i < 3; i++)
            obuff[i] = wav->fId[i];
        obuff[3] = 0;
        if (strcmp(obuff, "fmt") != 0) onError(PPPWAVIN_ERRFOR + string("FMT"));
        
        if (wav->wFormatTag != 1) onError(PPPWAVIN_ERRFOR + string("wFormatTag"));
        
        if ((wav->nBitsPerSample != 16) && (wav->nBitsPerSample != 8)) onError(
                PPPWAVIN_ERRFOR + string("BPS"));
        
        // save demographics
        _rate = (double) (wav->nSamplesPerSec);
        _bps = wav->nBitsPerSample;
        _channels = wav->nChannels;
        
        // skip over any remaining portion of wav header
        rmore = wav->pcm_header_len - (sizeof(WAV_HDR) - 20);
        wstat = fseek(fw, rmore, SEEK_CUR);
        if (wstat != 0) onError(PPPWAVIN_ERRFOR + string("SEEK"));
        
        // read chunks until a 'data' chunk is found
        sflag = 1;
        while (sflag != 0)
        {
            // check attempts
            if (sflag > 10) onError(PPPWAVIN_ERRFOR + string("CHUNKS"));
            
            // read chunk header
            wstat = std::fread((void *) chk, sizeof(CHUNK_HDR), (size_t) 1, fw);
            if (wstat != 1) onError(PPPWAVIN_ERRFOR + string("CHUNK"));
            
            // check chunk type
            for (i = 0; i < 4; i++)
                obuff[i] = chk->dId[i];
            obuff[4] = 0;
            if (strcmp(obuff, "data") == 0) break;
            
            // skip over chunk
            sflag++;
            wstat = fseek(fw, chk->dLen, SEEK_CUR);
            if (wstat != 0) onError(PPPWAVIN_ERRFOR + string("SEEK"));
        }
        
        /* find length of remaining data */
        wbuff_len = chk->dLen;
        
        // find number of samples 
        _points = chk->dLen;
        _points /= wav->nBitsPerSample / 8;
        
        /* allocate new buffers */
        wbuff = new char[wbuff_len];
        if (wbuff == NULL) onError(MEM_ERRALLOC + string("PPPWavIn"));
        
        _data = new int[points()];
        if (_data == NULL) onError(MEM_ERRALLOC + string("PPPWavIn"));
        
        /* read signal data */
        wstat = std::fread((void *) wbuff, wbuff_len, (size_t) 1, fw);
        if (wstat != 1) onError(PPPWAVIN_ERRFOR + string("WBUF"));
        
        // convert data
        if (wav->nBitsPerSample == 16)
        {
            uptr = (short *) wbuff;
            for (i = 0; i < points(); i++)
                _data[i] = (int) (uptr[i]);
        }
        else
        {
            cptr = (unsigned char *) wbuff;
            for (i = 0; i < points(); i++)
                _data[i] = (int) (cptr[i]) - 0x80;
        }
        
        // be polite - clean up
        delete wbuff;
        delete wav;
        delete chk;
        fclose(fw);
        return;
    }
    
    ~PPPWavIn ()
    {
        if (_data != NULL) delete _data;
    }
    
    // routine for reading one sample from a (previously loaded) wave file
    // returns current sample as a double
    inline int getData (long int aInd) const
    {
        _checkindex(aInd);
        return (_data[aInd]);
    }
    
    // returns number of samples in file
    inline long int points () const
    {
        return _points;
    }
    
    // reports number of channels (1==mono, 2==stereo)
    inline unsigned channels () const
    {
        return _channels;
    }
    
    // reports the number of bits in each sample
    inline unsigned getBPS () const
    {
        return _bps;
    }
    
    // reports sample rate in Hz
    inline double getRate () const
    {
        return _rate;
    }
    
    const char *getInfo (void)
    {
        char s1[80], s2[80], s3[80], s4[80];
        sprintf(s1, "%s = %1.0lf (Hz)", PPPWAVIN_SAMPLE, getRate());
        sprintf(s2, "%s = %ld", PPPWAVIN_LENGTH, points());
        sprintf(s3, "%s = %d", PPPWAVIN_BPS, getBPS());
        sprintf(s4, "%s = %d", PPPWAVIN_CHAN, channels());
        strstream str;
        str << getObjectName() << ":" << endl << "  " << s1 << endl << "  " << s2 << endl << "  "
        << s3 << endl << "  " << s4 << ends;
        onNotation(str.str());
        return getNotation();
    }
    
private:
    
    void _checkindex (long int aInd) const
    {
        if (aInd >= points())
        {
            strstream str;
            str << MEM_ERRINDEX << ", size = " << points() << ", index = " << aInd << ends;
            PPPBaseObject::onError(str.str());
        }
    }
    
};
// end of object

#endif 

