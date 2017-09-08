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
#ifndef _UTOPTION
#define _UTOPTION

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <argtable2.h>
using namespace std;
#define UTOPTION_ERRPARS  "error parsing in procedure: "
#define MEM_ERRALLOC      "allocation memory problem in procedure: "

/************************************************************************
 * UTOption
 ***********************************************************************/
class UTOption
{
private:
    
    bool _parsed;

public:

    /**
     * default constructor: at beginning we are not parsed
     */
    UTOption ()
    {
        _parsed = false;
    }
    
    virtual UTOption *
    getClone (const char * shortoption, const char * longoption) = 0;

    virtual const char *
    getShortOption () const =0;

    virtual const char *
    getLongOption () const = 0;

    virtual const char *
    getDataType () const =0;

    virtual const char *
    getGlossary () const = 0;

    void setParsed ()
    {
        _parsed = true;
    }
    
    bool isParsed () const
    {
        return _parsed;
    }
    
    /**
     * self-documentation
     */
    virtual void print () = 0;

    virtual void printUsage () = 0;

    /**
     * argtable2 compatible pointer to struct
     */
    virtual void * getData () const = 0;

    virtual bool isOptionGiven () const = 0;

    virtual
    void copyValueFrom (UTOption * option) = 0;

    virtual ~UTOption ()
    {
    }
    
};
// end of object

class UTMultiOption
{
    
private:
    vector<UTOption *> _option;		// the options
    
public:
    
    UTMultiOption (UTOption & option, unsigned n)
    {
        
        for (unsigned i = 0; i < n; ++i)
        {
            ostringstream s;
            char * shortoption;
            if (option.getShortOption() == NULL)
            {
                shortoption = NULL;
            }
            else
            {
                s << option.getShortOption() << i;
                shortoption = new char[s.str().length() + 1];
                s.str().copy(shortoption, string::npos);
                shortoption[s.str().length()] = 0;
            }
            
            s.clear();
            
            char * longoption;
            if (option.getLongOption() == NULL)
            {
                longoption = NULL;
            }
            else
            {
                s << option.getLongOption() << i;
                longoption = new char[s.str().length() + 1];
                s.str().copy(longoption, string::npos);
                longoption[s.str().length()] = 0;
            }
            
            _option.push_back(option.getClone(shortoption, longoption));
        }
        
        _option.push_back(&option);
        
    }
    
    UTOption *
    getOption (unsigned n)
    {
        return _option[n];
    }
    
    unsigned getN () const
    {
        return (unsigned) ((int) _option.size() - 1);
    }
};
// end of class UTMultiOption

/************************************************************************
 * wrapper for argtable options
 * UTOption_xxx
 ***********************************************************************/
template<typename D, typename argtabletype>
class UTOption_xxx : public UTOption
{
    
protected:
    argtabletype * _data;

protected:
    D _value;
    bool _set;

public:
    
    typedef UTOption_xxx<D, argtabletype> option_type;

    UTOption_xxx (argtabletype * opt, D const & defaultValue) :
            _data(opt),
            _value(defaultValue),
            _set(false)
    {
    }
    
    UTOption_xxx (argtabletype * opt) :
            _data(opt),
            _set(false)
    {
    }
    
    void print ()
    {
        cout << _data->hdr.longopts << "=" << getValue() << endl;
    }
    
    void printUsage ()
    {
        cout << "wrong usage of option " << _data->hdr.longopts << endl;
        cout << "--" << _data->hdr.longopts << " " << _data->hdr.datatype << " "
        << _data->hdr.glossary << endl;
    }
    
    void setDefault (D const & def)
    {
        _value = def;
    }
    
    bool isOptionGiven () const
    {
        return _data->count > 0;
    }
    
    bool isSet () const
    {
        return _set;
    }
    
    void * getData () const
    {
        return _data;
    }
    
    virtual D getValue ()
    {
        if (!isParsed())
        { // default behavior
            return _value;
        }
        else if (isSet() || !isOptionGiven())
        {
            return _value;
        }
        else
        {
            return getParsed();
        }
    }
    
    virtual
    void setValue (D const & val)
    {
        _value = val;
        _set = true;
    }
    
    const char *
    getShortOption () const
    {
        return (_data->hdr).shortopts;
    }
    
    const char *
    getLongOption () const
    {
        return (_data->hdr).longopts;
    }
    
    const char *
    getDataType () const
    {
        return (_data->hdr).datatype;
    }
    
    const char *
    getGlossary () const
    {
        return (_data->hdr).glossary;
    }
    
    virtual void copyValueFrom (UTOption * option)
    {
        setValue(dynamic_cast<option_type *>(option)->getValue());
    }
    
    virtual D getParsed () = 0;
    
};
// end of object

/************************************************************************
 * double options
 * UTOption_dbl
 ***********************************************************************/
class UTOption_dbl : public UTOption_xxx<double, arg_dbl>
{
private:
    double _unit; // unit 
    
public:
    
    UTOption_dbl (const char* shortOption, const char* longOption, const char* typeInfo,
                  const char* info, double defaultValue) :
            UTOption_xxx<double, arg_dbl>(arg_dbl0(shortOption, longOption, typeInfo, info),
                                          defaultValue),
            _unit(1)
    {
    }
    
    UTOption_dbl (const char* shortOption, const char* longOption, const char* typeInfo,
                  const char* info) :
            UTOption_xxx<double, arg_dbl>(arg_dbl0(shortOption, longOption, typeInfo, info)),
            _unit(1)
    {
    }
    
    UTOption_dbl (const char* shortOption, const char* longOption, const char* typeInfo,
                  double unit, const char* info) :
            UTOption_xxx<double, arg_dbl>(arg_dbl0(shortOption, longOption, typeInfo, info)),
            _unit(unit)
    {
    }
    
    UTOption_dbl (const char* shortOption, const char* longOption, const char* typeInfo,
                  double unit, const char* info, double defaultValue) :
            UTOption_xxx<double, arg_dbl>(arg_dbl0(shortOption, longOption, typeInfo, info),
                                          defaultValue),
            _unit(unit)
    {
    }
    
    double getParsed ()
    {
        if (!isParsed())
        {
            cout << UTOPTION_ERRPARS + string("getParsed") << endl;
            exit(-1);
        }
        return _data->dval[0];
    }
    
    double getValue ()
    {
        return UTOption_xxx<double, arg_dbl>::getValue() * _unit;
    }
    
    void setValue (double const & val)
    {
        _value = val / _unit;
        _set = true;
    }
    
    UTOption *
    getClone (const char * shortoption, const char * longoption)
    {
        return new UTOption_dbl(shortoption, longoption, getDataType(), getGlossary());
    }
    
    virtual ~UTOption_dbl ()
    {
    }
    
};
// end of object

/************************************************************************
 * string options
 * UTOption_str
 ***********************************************************************/
class UTOption_str : public UTOption_xxx<char const*, arg_str>
{
    
public:
    
    UTOption_str (const char* shortOption, const char* longOption, const char* typeInfo,
                  const char* info, const char* defaultValue) :
            UTOption_xxx<char const*, arg_str>(arg_str0(shortOption, longOption, typeInfo, info),
                                               defaultValue)
    {
    }
    
    UTOption_str (const char* shortOption, const char* longOption, const char* typeInfo,
                  const char* info) :
            UTOption_xxx<char const*, arg_str>(arg_str0(shortOption, longOption, typeInfo, info))
    {
    }
    
    char const* getParsed ()
    {
        if (!isParsed())
        {
            cout << UTOPTION_ERRPARS + string("getParsed") << endl;
            exit(-1);
        }
        return _data->sval[0];
    }
    
    UTOption *
    getClone (const char * shortoption, const char * longoption)
    {
        return new UTOption_str(shortoption, longoption, getDataType(), getGlossary());
    }
    
    virtual ~UTOption_str ()
    {
    }
    
};
// end of object

/************************************************************************
 * integer options
 * UTOption_int
 ***********************************************************************/
class UTOption_int : public UTOption_xxx<int, arg_int>
{
    
public:
    
    UTOption_int (const char* shortOption, const char* longOption, const char* typeInfo,
                  const char* info, int defaultValue) :
            UTOption_xxx<int, arg_int>(arg_int0(shortOption, longOption, typeInfo, info),
                                       defaultValue)
    {
    }
    
    UTOption_int (const char* shortOption, const char* longOption, const char* typeInfo,
                  const char* info) :
            UTOption_xxx<int, arg_int>(arg_int0(shortOption, longOption, typeInfo, info))
    {
    }
    
    int getParsed ()
    {
        if (!isParsed())
        {
            cout << UTOPTION_ERRPARS + string("getParsed") << endl;
            exit(-1);
        }
        return _data->ival[0];
    }
    
    UTOption *
    getClone (const char * shortoption, const char * longoption)
    {
        return new UTOption_int(shortoption, longoption, getDataType(), getGlossary());
    }
    
    virtual ~UTOption_int ()
    {
    }
    
};
// end of object

/************************************************************************
 * file option
 * UTOption_file
 ***********************************************************************/
class UTOption_file : public UTOption_xxx<char const *, arg_file>
{
    
public:
    
    UTOption_file (const char* shortOption, const char* longOption, const char* typeInfo,
                   const char* info, const char * defaultValue) :
            UTOption_xxx<char const *, arg_file>(arg_file0(shortOption, longOption, typeInfo, info),
                                                 defaultValue)
    {
    }
    
    UTOption_file (const char* shortOption, const char* longOption, const char* typeInfo,
                   const char* info) :
            UTOption_xxx<char const *, arg_file>(arg_file0(shortOption, longOption, typeInfo, info))
    {
    }
    
    UTOption_file (const char* shortOption, const char* longOption, const char* typeInfo,
                   const char* info, bool obligatory) :
            UTOption_xxx<char const *, arg_file>(arg_file1(shortOption, longOption, typeInfo, info))
    {
    }
    
    char const * getParsed ()
    {
        if (!isParsed())
        {
            cout << UTOPTION_ERRPARS + string("getParsed") << endl;
            exit(-1);
        }
        return _data->filename[0];
    }
    
    UTOption *
    getClone (const char * shortoption, const char * longoption)
    {
        return new UTOption_file(shortoption, longoption, getDataType(), getGlossary());
    }
    
    virtual ~UTOption_file ()
    {
    }
    
};
// end of object

/************************************************************************
 * literal options
 * UTOption_lit
 ***********************************************************************/
class UTOption_lit : public UTOption_xxx<bool, arg_lit>
{
    
public:
    
    UTOption_lit (const char* shortOption, const char* longOption, const char* info,
                  bool defaultValue) :
            UTOption_xxx<bool, arg_lit>(arg_lit0(shortOption, longOption, info), defaultValue)
    {
    }
    
    UTOption_lit (const char* shortOption, const char* longOption, const char* info) :
            UTOption_xxx<bool, arg_lit>(arg_lit0(shortOption, longOption, info))
    {
    }
    
    bool getParsed ()
    {
        if (!isParsed())
        {
            cout << UTOPTION_ERRPARS + string("getParsed") << endl;
            exit(-1);
        }
        return (_data->count > 0);
    }
    
    UTOption *
    getClone (const char * shortoption, const char * longoption)
    {
        return new UTOption_lit(shortoption, longoption, getGlossary());
    }
    
    virtual ~UTOption_lit ()
    {
    }
    
};
// end of object

#endif

