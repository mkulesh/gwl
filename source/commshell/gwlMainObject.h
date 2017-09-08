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
#ifndef _GWLMAINOBJECT
#define _GWLMAINOBJECT

#define PPPCONF_USEPARSING true      // using command line parsing
//#define PPPCONF_CHECKINDEX true      // index control in all container objects
#ifndef __WIN32__
#define PPPCONF_USEFFTW3 true      // using FFTW3 external library for FFT calculation for linux platform
#endif

#include "PPPtypes.h"
#include <map>

#define gwlMain_ERR01 "Type of output is incorrect in procedure: write()"
#define gwlMain_ERR02 "Format of input file is unknown"

/************************************************************************
 * gwlNameOptions
 ************************************************************************/
class gwlNameOptions
{
private:
    string _name, _default, _fullname;
public:
    gwlNameOptions (const char *aname, const char *adefault) :
            _name(aname),
            _default(adefault)
    {
        _fullname = _name + string(" (by default '") + _default + string("')");
    }
    
    const char * getFullName ()
    {
        return _fullname.c_str();
    }
    
};

/************************************************************************
 * gwlMainObject
 ************************************************************************/
template<class ATypeSource, class ATypeDest> class gwlMainObject : public PPPBaseObject
{
protected:

    UTParsing o_parser;
    UTOption_lit o_nomess;
    UTOption_lit o_progr;
    UTOption_lit o_iscmpl;
    gwlNameOptions o_infile_name;
    UTOption_file o_infile;
    gwlNameOptions o_outfile_name;
    UTOption_file o_outfile;
    UTOption_int o_outtype;
    UTOption_str o_wavelet;
    UTOption_dbl o_wavpar;

    PPPSpectrParams aSpectrumPar;
    vector<string> aFileNames;

public:

    ATypeSource aSource;
    ATypeDest aDest;

    gwlMainObject (const char *aAppName, const char *aModName, const char *aDefIn,
                   const char *aDefOut) :
            o_parser(aModName, aAppName),
            o_infile_name("input file name", aDefIn),
            o_outfile_name("output file name", aDefOut),
            o_nomess("m", "nomess", "if set, no messave will be printed", false),
            o_progr("r", "prog", "if set, show progress indicator during the calculation", false),
            o_iscmpl("c", "iscmpl", "if set, the input/output object will be complex", false),
            o_infile("i", "infile", "<file>", o_infile_name.getFullName(), aDefIn),
            o_outfile("o", "outfile", "<file>", o_outfile_name.getFullName(), aDefOut),
            o_outtype("t", "outtype", "<unsigned>",
                      "type of output: 1-ASCII, 2-binary (by default 2)", 2),
            o_wavelet("w", "wavelet", "<str>",
                      "wavelet name: haar, morlet, morletre, cauchy, shanon (by default 'morlet')",
                      "morlet"),
            o_wavpar("p", "wavpar", "<real>", "parameter of the wavelet (by default 1)", 1.0)
    {
        ConApplication.setAppName(
                string("\n") + string(aAppName) + string(", ") + string(PPPCONF_VERSION));
    }
    
    void parse (int argc, char **argv)
    {
        o_parser.parse(argc, argv);
        ConApplication.setMessageMode(!o_nomess.getValue());
        aFileNames.clear();
        if (o_infile.isParsed()) aFileNames.push_back(o_infile.getValue());
    }
    
    void parseFileNames (void)
    {
        aFileNames.clear();
        string aPar = o_infile.getValue();
        while (1)
        {
            unsigned pos = aPar.find(",");
            if (pos >= aPar.size())
            {
                aFileNames.push_back(aPar);
                break;
            }
            aFileNames.push_back(aPar.substr(0, pos));
            aPar.erase(0, pos + 1);
        }
        if (aFileNames.size() < 2) PPPBaseObject::onError("too few input signals are given");
    }
    
    bool isComplex ()
    {
        return o_iscmpl.isOptionGiven();
    }
    
    PPPObjectIO::ObjectType read_binheader (void)
    {
        if (aFileNames.size() == 0) return PPPObjectIO::UNDEF;
        PPPObjectIO ioobj;
        return ioobj.read_binheader(aFileNames[0].c_str());
    }
    
    template<class AType>
    void read_bindata (AType &aData, const char *aName, bool aCheckType = false)
    {
        FILE *infile = fopen(aName, "rb");
        if (infile == NULL) PPPBaseObject::onError(FILE_ERROPEN + string(aName));
        aData.fread(infile);
        if (aCheckType && strcmp(aData.getObjectVer(), PPPSPECTRCONTAINER_OBJVER) == 0) aSpectrumPar
                .fread(infile);
        fclose(infile);
    }
    
    virtual void read (void)
    {
        read_bindata(aSource, o_infile.getValue(), true);
        strstream str;
        str << aSource.getInfo() << ends;
        PPPBaseObject::onMessage(str.str());
    }
    
    template<class AType>
    void write_bindata (AType &aData, const char *aName)
    {
        onMessage(FILE_WRITEBIN + string(aName));
        FILE *outfile = fopen(aName, "wb");
        aData.fwrite(outfile);
        if (strcmp(aData.getObjectVer(), PPPSPECTRCONTAINER_OBJVER) == 0) aSpectrumPar.fwrite(
                outfile);
        fclose(outfile);
    }
    
    virtual void write ()
    {
        if (o_outtype.getValue() == 1)
        {
            aDest.write(o_outfile.getValue());
        }
        else if (o_outtype.getValue() == 2)
        {
            write_bindata(aDest, o_outfile.getValue());
        }
        else
            PPPBaseObject::onError(gwlMain_ERR01);
    }
    
    virtual void calc (void) = 0;

    void evaluate (void)
    {
        if (o_infile.isParsed()) read();
        calc();
        if (o_outfile.isParsed() && strlen(o_outfile.getValue()) > 0) write();
    }
    
};
// end of object

#endif

