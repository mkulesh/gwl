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
#ifndef _PPPAPPLICATIONCON
#define _PPPAPPLICATIONCON

/************************************************************************
 * PPPApplicationCon
 ***********************************************************************/
class PPPApplicationCon
{
private:

    string _appName;
    int _messageMode;
    string _messageLog, _messageFile;
    int _summaryMode;
    string _summaryLog, _summaryFile;
    time_t _startTime;

public:
    
    PPPApplicationCon (void)
    {
        _appName = "";
        _messageLog = "";
        _summaryLog = "";
        time(&_startTime);
        _messageMode = 1;
        _summaryMode = 0;
        _messageFile = "";
        _summaryFile = "";
    }
    
    ~PPPApplicationCon (void)
    {
        if (_summaryMode > 0 && _summaryFile != "" && _summaryLog != "")
        {
            stringstream str;
            str << "Calculation time: " << getExecTime() << " sec" << ends;
            _summaryLog.append(str.str());
            write(_summaryFile, _summaryLog, _summaryMode);
        }
    }
    
    void write (const string &aFile, const string &aSource, int aFormat)
    {
        fstream FileStream;
        remove(aFile.c_str());
        FileStream.open(aFile.c_str(), ios_base::out);
        if (aFormat == 1)
        {
            FileStream << getAppName() << endl;
            FileStream << aSource << endl;
        }
        if (aFormat == 2 || aFormat == 3)
        {
            stringstream streams;
            streams << getAppName() << endl;
            streams << aSource << endl;
            int buffsize = streams.str().length();
            char *str, *str1;
            if ((str = new char[buffsize]) == NULL) onError(MEM_ERRALLOC + string("write"));
            if ((str1 = new char[buffsize]) == NULL) onError(MEM_ERRALLOC + string("writewrite"));
            istream *tmp;
            if (aFormat == 2) FileStream << "function readme()" << endl;
            while (1)
            {
                str[0] = 0;
                tmp = &streams.getline(str, buffsize);
                if (tmp->fail()) break;
                if (aFormat == 2 || (aFormat == 3 && str[0] != '!'))
                    FileStream << "% " << str << endl;
                else
                    FileStream << (str + 1) << endl;
            }
            delete str;
            delete str1;
        }
        FileStream.close();
        return;
    }
    
    /*
     *  Application parameters
     */
    void setAppName (const string &aName)
    {
        _appName = aName;
        _appName.append(" ");
        _appName.append(__DATE__);
    }
    
    const char *getAppName (void)
    {
        return _appName.c_str();
    }
    
    double getExecTime (void)
    {
        time_t EndTime;
        time(&EndTime);
        return difftime(EndTime, _startTime);
    }
    
    /*
     *  Exceptions
     */
    void onError (const string &aNotation)
    {
        string error;
        error.append("\nFATAL ERROR\n  ");
        error.append(aNotation);
        error.append("\n  Program terminated\n");
        _summaryLog.append(error);
        printf("%s", error.c_str());
        exit(-1);
    }
    
    /*
     *  Messages and progresses
     */
    void setMessageMode (int aMode)
    {
        _messageMode = aMode;
    }
    
    void onMessage (const string &aNotation)
    {
        if (_messageMode == 0) return;
        _messageLog.append(aNotation);
        printf("%s\n", aNotation.c_str());
    }
    
    void onProgress (const int aPercent, const string &aNotation)
    {
        if (_messageMode == 0) return;
        if (aPercent < 0)
        {
            printf("\n");
            return;
        }
        if (aPercent == 0)
        {
            printf("  %s\r", aNotation.c_str());
            return;
        }
        if (aPercent > 0)
        {
            if (aNotation.size() > 0)
                printf("  %s (%d%%)\r", aNotation.c_str(), aPercent);
            else
                printf("  %d%%\r", aPercent);
        }
        return;
    }
    
    /*
     *  Summary
     */
    void setSummaryMode (int aMode, const string &aFile)
    {
        _summaryMode = aMode;
        _summaryFile = aFile;
    }
    
    const char *getSummary (void)
    {
        return _summaryLog.c_str();
    }
    
    void onSummary (const string &aNotation)
    {
        _summaryLog.append(aNotation);
        _summaryLog.append("\n");
    }
    
};
// end of object

#endif

